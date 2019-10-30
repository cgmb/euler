#include <assert.h>
#include <float.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <errno.h>
#include <unistd.h>
#include <math.h>
#include <time.h>

#include <algorithm>
#include <vector>

#include "math/vec2f.h"
#include "math/vec2zu.h"
#include "math/vec3f.h"
#include "misc/debug.h"
#include "misc/terminal.h"
#include "misc/file.h"
#include "misc/int.h"
#include "misc/rng.h"
#include "misc/color.h"
#include "misc/dynarray2d.h"
#include "misc/split_index2f.h"
#include "misc/bool_grid_adapters.h"
#include "misc/util.h"
#include "world_grid.h"

// simulation size
size_t X;
size_t Y;

// display size
size_t g_wx;
size_t g_wy;

struct args_t {
  const char* scenario_file;
  bool rainbow;
  bool pressure;
};

// simulation constants
const float k_s = 1; // side length
const float k_inv_s = 1/k_s;
const float k_d = 1.0; // density
const float k_inv_d = 1/k_d;
const float k_g = -10.f; // gravity

dynarray2d<float> g_u;
dynarray2d<float> g_v;

// grid cell properties from the scenario file
WorldGrid  g_grid;

// lookups for if a sample is valid
PBoolGridAdapter g_is_water_p = g_grid.fluid_grid();
UBoolGridAdapter g_is_water_u(g_is_water_p);
VBoolGridAdapter g_is_water_v(g_is_water_p);

PBoolGridAdapter g_was_water_p = g_grid.old_fluid_grid();
UBoolGridAdapter g_was_water_u(g_was_water_p);
VBoolGridAdapter g_was_water_v(g_was_water_p);

PBoolGridAdapter g_is_solid_p = g_grid.solid_grid();
UBoolGridAdapter g_is_solid_u(g_is_solid_p);
VBoolGridAdapter g_is_solid_v(g_is_solid_p);

// https://en.wikipedia.org/wiki/HSL_and_HSV
float hsv_basis(float t) {
  // periodic function, repeats after t=[0,6]
  // returns value in range of [0,1]
  t -= 6.f*floorf(1.f/6*t);
  if (t < 0.f) {
    t += 6.f;
  }

  if (t < 1.f) {
    return t;
  } else if (t < 3.f) {
    return 1.f;
  } else if (t < 4.f) {
    return 4.f - t;
  } else {
    return 0.f;
  }
}

// color data
bool g_rainbow_enabled;
dynarray2d<float> g_r;
dynarray2d<float> g_g;
dynarray2d<float> g_b;
const float k_source_color_period = 10.f; // seconds
const float k_initial_color_period = 60.f; // grid cells

bool g_view_pressure;
dynarray2d<double> g_p;

// simulation progress
bool g_pause;
uint32_t g_simulate_steps;
uint16_t g_frame_count;

// pressure matrix
struct sparse_entry_t {
  int8_t a_diag;
};

dynarray2d<sparse_entry_t> g_a;

// marker particle data
const size_t k_max_marker_count = INT16_MAX;
std::vector<vec2f> g_markers;
bool g_source_exhausted;
dynarray2d<u8> g_marker_count;

bool is_water(size_t y, size_t x) {
  return g_is_water_p[{x,y}];
}

bool is_water(const vec2zu& index) {
  return g_is_water_p[index];
}

bool was_water(const vec2zu& index) {
  return g_was_water_p[index];
}

bool was_water(size_t y, size_t x) {
  return was_water({x,y});
}

void refresh_marker_counts() {
  g_marker_count.fill(u8(0));
  for (size_t i = 0; i < g_markers.size(); ++i) {
    int x = (int)floorf(g_markers[i].x*k_inv_s);
    int y = (int)floorf(g_markers[i].y*k_inv_s);
    bool in_bounds = x > 0 && x < (int)X && y > 0 && y < (int)Y;
    assert(in_bounds);
    if (in_bounds) {
      vec2zu index(x,y);
      if (g_grid.is_sink(index)) { // remove markers in sinks
        g_markers[i--] = g_markers.back();
        g_markers.pop_back();
      } else {
        assert(!g_grid.is_solid(index));
        if (g_marker_count[index] < 16) {
          g_marker_count[index]++;
        } else { // remove excess markers
          g_markers[i--] = g_markers.back();
          g_markers.pop_back();
        }
      }
    } else { // remove out-of-bounds markers
      g_markers[i--] = g_markers.back();
      g_markers.pop_back();
    }
  }

  g_grid.mark_fluid_as_old();
  for (size_t y = 0; y < Y; ++y) {
    for (size_t x = 0; x < X; ++x) {
      g_grid.set_fluid({x,y}, g_marker_count[{x,y}] > 0);
    }
  }
}

// Find the mean of all valid neighbours.
template <class BoolGrid>
float neighbor_average(const dynarray2d<float>& q,
  const BoolGrid& valid, const vec2zu& i) {
  size_t min_y = (i.y == 0) ? 0 : i.y-1;
  size_t max_y = (i.y == q.height()-1) ? q.height()-1 : i.y+1;
  size_t min_x = (i.x == 0) ? 0 : i.x-1;
  size_t max_x = (i.x == q.width()-1) ? q.width()-1 : i.x+1;
  float total = 0.f;
  size_t count = 0;
  for (size_t y = min_y; y <= max_y; ++y) {
    for (size_t x = min_x; x <= max_x; ++x) {
      vec2zu neighbor(x,y);
      if (valid.at(neighbor)) {
        total += q[neighbor];
        ++count;
      }
    }
  }
  assert(count > 0);
  return total / count;
}

// Add values for cells that became fluid cells this frame due to advection.
template <class BoolGrid>
void extrapolate_new(dynarray2d<float>& q,
  const BoolGrid& was_valid, const BoolGrid& is_valid) {
  for (size_t y = 1; y < q.height()-1; ++y) {
    for (size_t x = 1; x < q.width()-1; ++x) {
      vec2zu index(x,y);
      if (!was_valid[index] && is_valid[index]) {
        q[index] = neighbor_average(q, was_valid, index);
      }
    }
  }
}

void extrapolate_new_p(dynarray2d<float>& q) {
  assert(q.dimensions() == g_was_water_p.dimensions());
  assert(q.dimensions() == g_is_water_p.dimensions());
  extrapolate_new(q, g_was_water_p, g_is_water_p);
}

void extrapolate_new_u(dynarray2d<float>& q) {
  assert(q.dimensions() == g_was_water_u.dimensions());
  assert(q.dimensions() == g_is_water_u.dimensions());
  extrapolate_new(q, g_was_water_u, g_is_water_u);
}

void extrapolate_new_v(dynarray2d<float>& q) {
  assert(q.dimensions() == g_was_water_v.dimensions());
  assert(q.dimensions() == g_is_water_v.dimensions());
  extrapolate_new(q, g_was_water_v, g_is_water_v);
}

template <class BoolGrid>
bool has_valid_neighbor(const BoolGrid& valid, const vec2zu& i) {
  // Unlike most functions, we might actually call this on an edge cell,
  // so we need to handle that case.
  size_t min_y = (i.y == 0) ? 0 : i.y-1;
  size_t max_y = (i.y == valid.height()-1) ? valid.height()-1 : i.y+1;
  size_t min_x = (i.x == 0) ? 0 : i.x-1;
  size_t max_x = (i.x == valid.width()-1) ? valid.width()-1 : i.x+1;
  for (size_t y = min_y; y <= max_y; ++y) {
    for (size_t x = min_x; x <= max_x; ++x) {
      vec2zu index(x,y);
      if (valid.at(index)) {
        return true;
      }
    }
  }
  return false;
}

// Add values for cells adjacent to fluid cells, as that staggered grid means
// that // each U and V cell is overlaps two P cells. At the fluid boundary,
// only one P cell will be filled, and the other will be empty. e.g.
//
//   0   1   2
//   ╎   ╎   ╎
//   ├─x─┼─X─┼╌ 1
//   │  .│   │
//   └─X─┴─X─┴╌ 0
//
// In the diagram above, the marker (.) fills the bottom left P Cell.
// The vertical velocity at the marker's position would be calculated
// using the values at the four Xs, however, the right two Xs are outside
// the fluid region and thus do not have a value.
//
// This is also a problem for advection, which may query a property at up to
// half a cell-length away from the current position of an valid sample point
// or marker particle.
//
// This extrapolation would assign values to the two right Xs equal to the
// average of their valid neighbours. Though, it should probably do linear
// extrapolation instead.
template <class BoolGrid>
void extrapolate_air(dynarray2d<float>& q, const BoolGrid& is_valid) {
  for (size_t y = 0; y < q.height(); ++y) {
    for (size_t x = 0; x < q.width(); ++x) {
      vec2zu index(x,y);
      if (!is_valid.at(index)) {
        if (has_valid_neighbor(is_valid, index)) {
          q[index] = neighbor_average(q, is_valid, index);
        } else {
          q[index] = 0;
        }
      }
    }
  }
}

void extrapolate_air_p(dynarray2d<float>& q) {
  extrapolate_air(q, g_is_water_p);
}

void extrapolate_air_u(dynarray2d<float>& q) {
  extrapolate_air(q, g_is_water_u);
}

void extrapolate_air_v(dynarray2d<float>& q) {
  extrapolate_air(q, g_is_water_v);
}

void colorize() {
  for (size_t y = 0; y < Y; ++y) {
    for (size_t x = 0; x < X; ++x) {
      if (is_water(y,x)) {
        float t = 0.f;
        if (!g_grid.is_source({x,y})) {
          t = (x + y) * 6.f / k_initial_color_period;
        }
        g_r[{x,y}] = hsv_basis(t + 2.f);
        g_g[{x,y}] = hsv_basis(t);
        g_b[{x,y}] = hsv_basis(t - 2.f);
      }
    }
  }
}

// Returns a float in the range [0,1)
float randf() {
  static uint64_t rng_state = 0x9bd185c449534b91;
  uint32_t x = xorshift64_32star(&rng_state);
  return float(x / double(UINT32_MAX));
}

// Returns a pair of floats in the range [0,s)
vec2f rand2f(float s) {
  return s*vec2f(randf(), randf());
}

vec2f index_to_world(const vec2f& i) {
  return k_s * i;
}

void sim_init(args_t in) {
  g_grid = WorldGrid::from_file(in.scenario_file);
  X = g_grid.width();
  Y = g_grid.height();

  g_u.resize(X+1,Y);
  g_v.resize(X,Y+1);

  g_marker_count.resize(X,Y);
  g_a.resize(X,Y);
  g_r.resize(X,Y);
  g_g.resize(X,Y);
  g_b.resize(X,Y);

  g_p.resize(X,Y);

  // setup fluid markers, 4 per cell, jittered
  for (size_t y = 0; y < g_grid.height(); ++y) {
    for (size_t x = 0; x < g_grid.width(); ++x) {
      if (g_grid.is_fluid(vec2zu(x,y))) {
        g_markers.push_back(index_to_world(vec2f(x,      y)      + rand2f(0.5)));
        g_markers.push_back(index_to_world(vec2f(x+0.5f, y)      + rand2f(0.5)));
        g_markers.push_back(index_to_world(vec2f(x,      y+0.5f) + rand2f(0.5)));
        g_markers.push_back(index_to_world(vec2f(x+0.5f, y+0.5f) + rand2f(0.5)));
      }
    }
  }
  refresh_marker_counts();

  // setup color
  if (g_rainbow_enabled) {
    colorize();
  }
}

void update_fluid_sources() {
  // If we ever hit the max number of markers, that's it for adding fluid.
  // The current extrapolation implementation assumes that the fluid is never
  // more than one cell from a cell where fluid was in the previous step. If
  // we stop and then later restart generating fluid, that may not hold true.
  g_source_exhausted |= (g_markers.size() == k_max_marker_count);

  float t = 0.6f / k_source_color_period * g_frame_count;
  for (size_t y = 0; y < Y; ++y) {
    for (size_t x = 0; x < X; ++x) {
      vec2zu index(x,y);
      if (g_grid.is_source(index)) {
        if (!g_source_exhausted && g_marker_count[index] < 4) {
          g_markers.push_back(index_to_world(vec2f(x+randf(), y+randf())));
          g_marker_count[index]++;
          g_source_exhausted |= (g_markers.size() == k_max_marker_count);
        }
        if (g_rainbow_enabled) {
          g_r[index] = hsv_basis(t + 2.f);
          g_g[index] = hsv_basis(t);
          g_b[index] = hsv_basis(t - 2.f);
        }
      }
    }
  }
}

float bilinear(const dynarray2d<float>& q, size_t bottom_left, const vec2f& fraction) {
  const float x = fraction.x;
  const float y = fraction.y;
  const float rx = 1 - x;
  const float ry = 1 - y;

  const size_t bottom_right = bottom_left + 1;
  const size_t top_left = bottom_left + q.width();
  const size_t top_right = top_left + 1;

  return q[bottom_left]*(rx*ry) + q[bottom_right]*(x*ry)
       + q[top_left]*(rx*y)     + q[top_right]*(x*y);
}

float bilinear(const dynarray2d<float>& q, const SplitIndex2f& index) {
  return bilinear(q, q.flatten(index.whole), index.fractional);
}

vec2f to_vec2f(const vec2zu& v) {
  return vec2f(float(v.x), float(v.y));
}

// u
vec2f world_position_to_uindex(const vec2f& position) {
  return k_inv_s * position - vec2f(0,0.5);
}

vec2f uindex_to_world_position(const vec2f& index) {
  return k_s * (index + vec2f(0,0.5));
}

float interpolate_u(const dynarray2d<float>& u, const vec2f& position) {
  vec2f index = world_position_to_uindex(position);
  return bilinear(u, SplitIndex2f(index));
}

// v
vec2f world_position_to_vindex(const vec2f& position) {
  return k_inv_s * position - vec2f(0.5,0);
}

vec2f vindex_to_world_position(const vec2f& index) {
  return k_s * (index + vec2f(0.5,0));
}

float interpolate_v(const dynarray2d<float>& v, const vec2f& position) {
  vec2f index = world_position_to_vindex(position);
  return bilinear(v, SplitIndex2f(index));
}

// p
vec2f world_position_to_pindex(const vec2f& position) {
  return k_inv_s * position - vec2f(0.5,0.5);
}

vec2f pindex_to_world_position(const vec2f& index) {
  return k_s * (index + vec2f(0.5,0.5));
}

float interpolate_p(const dynarray2d<float>& p, const vec2f& position) {
  vec2f index = world_position_to_pindex(position);
  return bilinear(p, SplitIndex2f(index));
}

vec2f velocity_at(const vec2f& pos, const dynarray2d<float>& grid_u, const dynarray2d<float>& grid_v) {
  float u = interpolate_u(grid_u, pos);
  float v = interpolate_v(grid_v, pos);
  return vec2f(u,v);
}

// I've ensured there's no fluid in the p-cells around the boundary, and fluid
// can only travel one cell per step, thus we do not need boundary checks as
// the grid contains extrapolated values for one cell beyond the fluid edge.
// Those extrapolated values are used for advection, but are not advected
// themselves, so we can also ignore them when iterating.
dynarray2d<float> advect_u(const dynarray2d<float>& u, const dynarray2d<float>& v, float dt) {
  dynarray2d<float> out(u.dimensions());
  for (size_t y = 1; y < u.height()-1; ++y) {
    for (size_t x = 1; x < u.width()-1; ++x) {
      vec2zu index(x,y);
      if (g_is_water_u[index]) {
        vec2f pos = uindex_to_world_position(to_vec2f(index));
        vec2f prev_pos = pos - velocity_at(pos, u, v) * dt; // could be optimized
        out[index] = interpolate_u(u, prev_pos);
      }
    }
  }
  return out;
}

dynarray2d<float> advect_v(const dynarray2d<float>& u, const dynarray2d<float>& v, float dt) {
  dynarray2d<float> out(v.dimensions());
  for (size_t y = 1; y < v.height()-1; ++y) {
    for (size_t x = 1; x < v.width()-1; ++x) {
      vec2zu index(x,y);
      if (g_is_water_v[index]) {
        vec2f pos = vindex_to_world_position(to_vec2f(index));
        vec2f prev_pos = pos - velocity_at(pos, u, v) * dt; // could be optimized
        out[index] = interpolate_v(v, prev_pos);
      }
    }
  }
  return out;
}

dynarray2d<float> advect_p(const dynarray2d<float>& q, const dynarray2d<float>& u, const dynarray2d<float>& v, float dt) {
  dynarray2d<float> out(q.dimensions());
  for (size_t y = 1; y < Y-1; ++y) {
    for (size_t x = 1; x < X-1; ++x) {
      vec2zu index(x,y);
      if (g_is_water_p[index]) {
        vec2f pos = pindex_to_world_position(to_vec2f(index));
        vec2f prev_pos = pos - velocity_at(pos, u, v) * dt; // could be optimized
        out[index] = interpolate_p(q, prev_pos);
      }
    }
  }
  return out;
}

float time_to(float p0, float p1, float v) {
  if (fabsf(v) > 0.f) {
    return (p1 - p0) / v;
  } else {
    return FLT_MAX;
  }
}

// We could do something based on:
// "A Fast Voxel Traversal Algorithm for Ray Tracing" (1987)
// However, I know I'm only travelling into a neighbouring cell,
// so we can simply handle all the cases.
//
// This can be extremely painful due to floating-point imprecision if you
// rearrange any equations. Always try to directly test whatever position
// you're moving to. Testing indirections like the intersection time with a
// boundary might not always agree with `is_solid(world_to_index(marker_pos))`
void advect_markers(float dt, const dynarray2d<float>& grid_u, const dynarray2d<float>& grid_v) {
  for (size_t i = 0; i < g_markers.size(); ++i) {
    vec2f p = g_markers[i];
    vec2f v = velocity_at(p, grid_u, grid_v);
    vec2f np = p + v*dt;

    int x_idx = (int)floorf(p.x*k_inv_s);
    int y_idx = (int)floorf(p.y*k_inv_s);
    assert(x_idx > 0 && x_idx < int(X));
    assert(y_idx > 0 && y_idx < int(Y));
    vec2zu index(x_idx, y_idx);

    int nx_idx = (int)floorf(np.x*k_inv_s);
    int ny_idx = (int)floorf(np.y*k_inv_s);
    assert(nx_idx > 0 && nx_idx < int(X));
    assert(ny_idx > 0 && ny_idx < int(Y));
    vec2zu next_index(nx_idx, ny_idx);

    // No possible collisions if we don't leave the cell
    if (index == next_index) {
      g_markers[i] = np;
      continue;
    }

    // Don't move at all if the destination would be solid
    if (g_grid.is_solid(next_index)) {
      continue;
    }

    // If we're going entirely horizontally or vertically,
    // we will go directly from our current cell to our destination.
    // And, since we know the destination isn't solid...
    if (x_idx == nx_idx || y_idx == ny_idx) {
      g_markers[i] = np;
      continue;
    }

    // We're going diagonally to an empty cell.
    // Does it matter which path we take?
    bool nx_y_solid = g_grid.is_solid({size_t(nx_idx), size_t(y_idx)});
    bool x_ny_solid = g_grid.is_solid({size_t(x_idx), size_t(ny_idx)});
    if (!nx_y_solid && !x_ny_solid) {
      // Not if the way is totally clear!
      g_markers[i] = np;
      continue;
    } else if (nx_y_solid && x_ny_solid) {
      // Nor if the way is totally blocked.
      continue;
    }

    // The diagonal is partly blocked.
    // Find out if we cross the horizontal or vertical axis first.
    float np_x = k_s*(x_idx + (v.x > 0 ? 1 : 0));
    float t_x = time_to(p.x, np_x, v.x);

    float np_y = k_s*(y_idx + (v.y > 0 ? 1 : 0));
    float t_y = time_to(p.y, np_y, v.y);

    if (t_x < t_y) {
      if (nx_y_solid) {
        v.x = 0;
      }
    } else {
      if (x_ny_solid) {
        v.y = 0;
      }
    }

    g_markers[i] = p + v*dt;
    continue;
  }
}

void apply_body_forces(dynarray2d<float>& v, float dt) {
  for (size_t y = 0; y < Y-1; ++y) {
    for (size_t x = 0; x < X; ++x) {
      v[{x,y}] += k_g * dt;
    }
  }
}

int8_t nonsolid_neighbor_count(size_t y, size_t x) {
  return 4 - g_grid.is_solid({x-1, y}) - g_grid.is_solid({x+1, y})
           - g_grid.is_solid({x, y-1}) - g_grid.is_solid({x, y+1});
}

int8_t get_a_plus_i(size_t y, size_t x) {
  assert(is_water(y, x));
  return is_water(y, x+1) ? -1 : 0;
}

int8_t get_a_plus_j(size_t y, size_t x) {
  assert(is_water(y, x));
  return is_water(y+1, x) ? -1 : 0;
}

int8_t get_a_minus_i(size_t y, size_t x) {
  assert(is_water(y, x));
  return is_water(y, x-1) ? -1 : 0;
}

int8_t get_a_minus_j(size_t y, size_t x) {
  assert(is_water(y, x));
  return is_water(y-1, x) ? -1 : 0;
}

void apply_preconditioner(const dynarray2d<double>& r, dynarray2d<double>& z) {
  // Incomplete Cholesky
  // A ~= LLᵀ
  // L = F * E_inv + E

  // calculate E_inv (precon)
  dynarray2d<double> precon(X,Y);
  for (size_t y = 0; y < Y; ++y) {
    for (size_t x = 0; x < X; ++x) {
      if (is_water(y, x)) {
        double a = g_a[{x,y}].a_diag;
        double b = sq(get_a_minus_i(y,x) * precon[{x-1,y}]);
        double c = sq(get_a_minus_j(y,x) * precon[{x,y-1}]);

        double e = a - b - c;
        if (e < 0.25*a) {
          e = a != 0 ? a : 1;
        }
        precon[{x,y}] = 1 / sqrt(e);
      }
    }
  }

  // solve Lq = r
  dynarray2d<double> q(X,Y);
  for (size_t y = 0; y < Y; ++y) {
    for (size_t x = 0; x < X; ++x) {
      if (is_water(y, x)) {
        double t = r[{x,y}];
        if (is_water(y,x-1)) {
          t -= get_a_plus_i(y,x-1) * precon[{x-1,y}] * q[{x-1,y}];
        }
        if (is_water(y-1,x)) {
          t -= get_a_plus_j(y-1,x) * precon[{x,y-1}] * q[{x,y-1}];
        }
        q[{x,y}] = t * precon[{x,y}];
      }
    }
  }

  // solve Lᵀz = q
  z.fill(0);
  for (size_t y = Y; y--;) {
    for (size_t x = X; x--;) {
      if (is_water(y, x)) {
        double t = q[{x,y}];
        if (is_water(y,x+1)) {
          t -= get_a_plus_i(y,x) * precon[{x,y}] * z[{x+1,y}];
        }
        if (is_water(y+1,x)) {
          t -= get_a_plus_j(y,x) * precon[{x,y}] * z[{x,y+1}];
        }
        z[{x,y}] = t * precon[{x,y}];
      }
    }
  }
}

double dot(const dynarray2d<double>& a, const dynarray2d<double>& b) {
  double total = 0.f;
  for (size_t y = 0; y < Y; ++y) {
    for (size_t x = 0; x < X; ++x) {
      if (is_water(y, x)) {
        total += a[{x,y}] * b[{x,y}];
      }
    }
  }
  return total;
}

bool all_zero(const dynarray2d<double>& r) {
  for (size_t y = 0; y < Y; ++y) {
    for (size_t x = 0; x < X; ++x) {
      if (is_water(y, x)) {
        if (r[{x,y}] != 0.f) {
          return false;
        }
      }
    }
  }
  return true;
}

double inf_norm(const dynarray2d<double>& r) {
  double max = 0.f;
  for (size_t y = 0; y < Y; ++y) {
    for (size_t x = 0; x < X; ++x) {
      if (is_water(y, x)) {
        double a = fabs(r[{x,y}]);
        if (max < a) {
          max = a;
        }
      }
    }
  }
  return max;
}

void update_search(dynarray2d<double>& s, const dynarray2d<double>& z, double beta) {
  for (size_t y = 0; y < Y; ++y) {
    for (size_t x = 0; x < X; ++x) {
      if (is_water(y, x)) {
        s[{x,y}] = z[{x,y}] + beta * s[{x,y}];
      }
    }
  }
}

void apply_a(const dynarray2d<double>& s, dynarray2d<double>& out) {
  for (size_t y = 0; y < Y; ++y) {
    for (size_t x = 0; x < X; ++x) {
      if (is_water(y, x)) {
        out[{x,y}] = g_a[{x,y}].a_diag * s[{x,y}]
          + (is_water(y, x+1) ? -s[{x+1, y}] : 0)
          + (is_water(y+1, x) ? -s[{x, y+1}] : 0)
          + (is_water(y, x-1) ? -s[{x-1, y}] : 0)
          + (is_water(y-1, x) ? -s[{x, y-1}] : 0);
      }
    }
  }
}

// c += a*b
void fmadd(const dynarray2d<double>& a, double b, dynarray2d<double>& c) {
  for (size_t y = 0; y < Y; ++y) {
    for (size_t x = 0; x < X; ++x) {
      if (is_water(y, x)) {
        c[{x,y}] += a[{x,y}] * b;
      }
    }
  }
}

void project(float dt, const dynarray2d<float>& u, const dynarray2d<float>& v, dynarray2d<float>& uout, dynarray2d<float>& vout) {
  const double c = -k_d*k_s*k_s / dt; // -density * dt^2 / dt
  dynarray2d<double> d0(X,Y); // divergence * c

  // calculate d0
  for (size_t y = 0; y < Y; ++y) {
    for (size_t x = 0; x < X; ++x) {
      if (is_water(y,x)) {
        d0[{x,y}] = c * k_inv_s * (u[{x+1,y}] - u[{x,y}] + v[{x,y+1}] - v[{x,y}]);
      }
    }
  }

  // calculate A
  for (size_t y = 0; y < Y; ++y) {
    for (size_t x = 0; x < X; ++x) {
      if (is_water(y, x)) {
        g_a[{x,y}].a_diag = nonsolid_neighbor_count(y, x);
      }
    }
  }

  const size_t max_iterations = 200;
  const double tol = 1e-6f;

  // conjugate gradient
  dynarray2d<double> p(X,Y);    // pressure guess
  dynarray2d<double> r = d0;    // residual
  if (!all_zero(r)) {
    dynarray2d<double> z(X,Y);  // auxiliary vector
    apply_preconditioner(r, z);
    dynarray2d<double> s = z;   // search vector

    double sigma = dot(z,r);
    for (size_t i = 0; i < max_iterations; ++i) {
      apply_a(s, z);

      double alpha = sigma / dot(z,s);
      fmadd(s, alpha, p);  // p += alpha*s
      fmadd(z, -alpha, r); // r -= alpha*z

      if (inf_norm(r) <= tol) {
        break;
      }

      apply_preconditioner(r, z);

      double sigma_new = dot(z,r);
      double beta = sigma_new / sigma;
      update_search(s,z,beta);
      sigma = sigma_new;
    }
  }

  // clamp pressure to zero or greater
  // This step is not in Bridson's notes, but it fixes some weird artifacts.
  // Without it, the pressure solve sometimes introduced negative pressures
  // to eliminate divergence at solid boundaries, causing extreme stickyness
  for (size_t y = 0; y < Y; ++y) {
    for (size_t x = 0; x < X; ++x) {
      vec2zu index(x,y);
      if (g_is_water_p[index] && p[index] < 0.f) {
        p[index] = 0.f;
      }
    }
  }

  // update horizontal velocities
  for (size_t y = 0; y < u.height(); ++y) {
    for (size_t x = 0; x < u.width(); ++x) {
      vec2zu index(x,y);
      if (g_is_solid_u.at(index)) {
        uout[index] = 0.f;
      } else if (g_is_water_u.at(index)) {
        uout[index] = u[index] - k_inv_d * k_inv_s * dt * (p[{x,y}] - p[{x-1,y}]);
      } else {
        uout[index] = 0.f;
      }
    }
  }

  // update vertical velocities
  for (size_t y = 0; y < v.height(); ++y) {
    for (size_t x = 0; x < v.width(); ++x) {
      vec2zu index(x,y);
      if (g_is_solid_v.at(index)) {
        vout[index] = 0.f;
      } else if (g_is_water_v.at(index)) {
        vout[index] = v[index] - k_inv_d * k_inv_s * dt * (p[{x,y}] - p[{x,y-1}]);
      } else {
        vout[index] = 0.f;
      }
    }
  }

  if (g_view_pressure) {
    g_p = std::move(p);
  }
}

float sq_max(const dynarray2d<float>& q) {
  float max = 0;
  for (size_t i = 0; i < q.size(); ++i) {
    float value = sqf(q[i]);
    if (max < value) {
      max = value;
    }
  }
  return max;
}

template <class BoolGrid>
void zero_bounds(dynarray2d<float>& q,
  const BoolGrid& /*is_water*/, const BoolGrid& is_solid) {
  for (size_t y = 0; y < q.height(); ++y) {
    for (size_t x = 0; x < q.width(); ++x) {
      vec2zu index(x,y);
      if (/*!is_water.at(index) || */is_solid.at(index)) {
        q[index] = 0.f;
      }
    }
  }
}

void zero_bounds_u(dynarray2d<float>& u) {
  zero_bounds(u, g_is_water_u, g_is_solid_u);
}

void zero_bounds_v(dynarray2d<float>& v) {
  zero_bounds(v, g_is_water_v, g_is_solid_v);
}

float calculate_timestep(float frame_time) {
  // Bridson suggests a limit of five for stability, but my implementation of
  // advection and extrapolation assume that new fluid cells are within one
  // grid cell of old fluid cells
  const float m = 63/128.f; // maximum number of cells to traverse in one step

  float dt;
  float max_velocity = sqrtf(sq_max(g_u) + sq_max(g_v));
  if (max_velocity*frame_time < m*k_s) {
    dt = frame_time;
  } else {
    dt = m*k_s / max_velocity;
  }
  return dt;
}

void sim_step() {
  if (g_pause && g_simulate_steps == 0) {
    return;
  }

  const float total_frame_time = 0.1f;
  float frame_time = total_frame_time;
  do {
    float dt = calculate_timestep(frame_time);
    assert((total_frame_time/16) <= dt || frame_time == dt);
    frame_time -= dt;

    extrapolate_air_u(g_u);
    extrapolate_air_v(g_v);

    zero_bounds_u(g_u);
    zero_bounds_v(g_v);

    advect_markers(dt, g_u, g_v);
    refresh_marker_counts();
    update_fluid_sources();
    if (g_rainbow_enabled) {
      extrapolate_new_p(g_r);
      extrapolate_new_p(g_g);
      extrapolate_new_p(g_b);
    }
    extrapolate_new_u(g_u);
    extrapolate_new_v(g_v);

    if (g_rainbow_enabled) {
      extrapolate_air_p(g_r);
      extrapolate_air_p(g_g);
      extrapolate_air_p(g_b);
    }
    extrapolate_air_u(g_u);
    extrapolate_air_v(g_v);

    zero_bounds_u(g_u);
    zero_bounds_v(g_v);

    dynarray2d<float> utmp = advect_u(g_u, g_v, dt);
    dynarray2d<float> vtmp = advect_v(g_u, g_v, dt);

    if (g_rainbow_enabled) {
      g_r = advect_p(g_r, g_u, g_v, dt);
      g_g = advect_p(g_g, g_u, g_v, dt);
      g_b = advect_p(g_b, g_u, g_v, dt);
    }
    apply_body_forces(vtmp, dt);

    zero_bounds_u(utmp);
    zero_bounds_v(vtmp);

    project(dt, utmp, vtmp, g_u, g_v);
  } while (frame_time > 0.f);

  if (g_simulate_steps) {
    g_simulate_steps--;
  }
  g_frame_count++;
}

int float_to_byte_color(float x) {
  x = clampf(0.f, x, nextafterf(1.f, 0.f));
  return (int)256.f*x;
}

float linear_to_sRGB(float x) {
  return powf(x, 1/2.2f); // approximation
}

enum ColorMode {
  CM_FOREGROUND = '3',
  CM_BACKGROUND = '4'
};

void buffer_append_color(buffer* buf, const vec3f& color, ColorMode mode = CM_FOREGROUND) {
  char tmp[20];
  int r_out = float_to_byte_color(linear_to_sRGB(color.x));
  int g_out = float_to_byte_color(linear_to_sRGB(color.y));
  int b_out = float_to_byte_color(linear_to_sRGB(color.z));
  int length = sprintf(tmp, "\x1B[%c8;2;%d;%d;%dm", mode, r_out, g_out, b_out);
  if (length < 0) {
    die("sprintf");
  }
  buffer_append(buf, tmp, length);
}

void buffer_appendz(buffer* buf, const char* s) {
  buffer_append(buf, s, strlen(s));
}

struct run_t {
  int length;
  char type;
};

char type_at(const vec2zu& index) {
  if (g_grid.is_solid(index)) {
    return 'X';
  } else if (g_grid.is_sink(index)) {
    return '=';
  } else {
    int lower = g_grid.is_fluid({index.x, index.y-1}) ? 0x1 : 0;
    int upper = g_grid.is_fluid(index) ? 0x2 : 0;
    return char(lower + upper);
  }
}

vec3f rgb(const vec2zu& index) {
  return vec3f(g_r[index], g_g[index], g_b[index]);
}

vec3f pressure_color(const vec2zu& index) {
  return colormap(clampf(0, g_p[index]/128.0, 1));
}

void draw_fluid(const vec2zu& end, char type, int count, buffer* buf) {
  const char* symbol[] = {" ","▄","▀","█"};

  if (!g_rainbow_enabled && !g_view_pressure) {
    buffer_appendz(buf, T_BLUE);
  }
  const size_t y = end.y;
  for (size_t x = end.x - count; x < end.x; x++) {
    vec2zu index(x,y);
    int symbol_idx = type;
    bool reset_background = false;
    if (g_rainbow_enabled) { // use the background color with half-blocks for double-vertical resolution
      if (type == 3) {
        buffer_append_color(buf, rgb(index), CM_BACKGROUND);
        buffer_append_color(buf, rgb({index.x, index.y-1}), CM_FOREGROUND);
        symbol_idx = 1;
        reset_background = true;
      } else if (type == 2) {
        buffer_append_color(buf, rgb(index), CM_FOREGROUND);
      } else if (type == 1) {
        buffer_append_color(buf, rgb({index.x, index.y-1}), CM_FOREGROUND);
      }
    } else if (g_view_pressure) {
      if (type == 3) {
        buffer_append_color(buf, pressure_color(index), CM_BACKGROUND);
        buffer_append_color(buf, pressure_color({index.x, index.y-1}), CM_FOREGROUND);
        symbol_idx = 1;
        reset_background = true;
      } else if (type == 2) {
        buffer_append_color(buf, pressure_color(index), CM_FOREGROUND);
      } else if (type == 1) {
        buffer_append_color(buf, pressure_color({index.x, index.y-1}), CM_FOREGROUND);
      }
    }
    buffer_appendz(buf, symbol[int(symbol_idx)]);
    if (reset_background) {
      buffer_appendz(buf, T_BG_DEFAULT);
    }
  }
  buffer_appendz(buf, T_RESET);
}

void draw_nonfluid(char type, int count, buffer* buf) {
  const char* symbol;
  if (type == 'X') {
    symbol = "█";
  } else if (type == '=') {
    symbol = "░";
  } else {
    symbol = " ";
  }
  for (int i = 0; i < count; ++i) {
    buffer_appendz(buf, symbol);
  }
}

void draw_run(const run_t& run, const vec2zu& index, buffer* buf) {
  if (run.type < 4) {
    draw_fluid(index, run.type, run.length, buf);
  } else {
    draw_nonfluid(run.type, run.length, buf);
  }
}

void draw_rows(buffer* buf) {
  ssize_t y_cutoff = (Y-2 < 2*g_wy) ? 1 : Y-2 - 2*g_wy;
  for (ssize_t ys = Y-2; ys > y_cutoff; ys -= 2) {
    size_t y = size_t(ys);
    run_t run = { 1, type_at({1,y}) };
    size_t x = 2;
    for (; x < X-1 && x < g_wx+1; x++) {
      char type = type_at({x,y});
      if (type == run.type) {
        run.length++;
      } else {
        draw_run(run, {x,y}, buf);
        run.type = type;
        run.length = 1;
      }
    }
    draw_run(run, {x,y}, buf);
    buffer_appendz(buf, "\x1b[K"); // clear remainer of line
    if (ys > y_cutoff) {
      buffer_appendz(buf, "\r\n");
    }
  }
}

bool g_show_status = false;
long g_frametime;

void draw_status_bar(buffer* buf) {
  if (g_show_status) {
    char tmp[120];
    snprintf(tmp, sizeof(tmp),
      "frame (ms):%4ld | "
      "grid: [%3zu,%3zu] | "
      "markers:%7zu"
      "\n", g_frametime, X, Y, g_markers.size());
    buffer_appendz(buf, tmp);
  }
}

void draw(buffer* buf) {
  buffer_clear(buf);
  reposition_cursor(buf);
  draw_rows(buf);
  draw_status_bar(buf);
  hide_cursor(buf);
  buffer_write(buf);
}

bool process_keypress() {
  char c = '\0';
  if (read(STDIN_FILENO, &c, 1) == -1 && errno != EAGAIN && errno != EINTR) {
    die("read");
  }

  if (c == 'p') {
    g_pause = !g_pause;
  } else if (c == 'f') {
    g_simulate_steps++;
  } else if (c == 'r') {
    if (g_rainbow_enabled) {
      colorize();
    }
  } else if (c == 's') {
    g_show_status = !g_show_status;
  } else if (c == 'q') {
    u_clear_screen();
    return false;
  }
  return true;
}

args_t parse_args(int argc, char** argv) {
  args_t in;
  in.rainbow = false;
  in.pressure = false;

  if (argc < 2) {
    fprintf(stderr, "usage: %s [--rainbow|--pressure] <scenario>\n", argv[0]);
    exit(1);
  }

  for (int i = 1; i < argc - 1; ++i) {
    if (!strcmp(argv[i], "--rainbow")) {
      in.rainbow = true;
    } else if (!strcmp(argv[i], "--pressure")) {
      in.pressure = true;
    } else {
      fprintf(stderr, "Unrecognized input: %s\n", argv[i]);
      exit(1);
    }
  }

  if (in.rainbow && in.pressure) {
    fprintf(stderr, "Cannot display both rainbow and pressure colors "
      "at the same time\n");
    exit(1);
  }

  in.scenario_file = argv[argc-1];
  return in;
}

void update_window_size() {
  int y, x;
  if (get_window_size(&y, &x) == -1) {
    die("get_window_size");
  }
  g_wy = size_t(max(y,0));
  g_wx = size_t(max(x,0));
}

void handle_window_size_changed(int) {
  update_window_size();
  u_clear_screen();
}

timespec subtract(const timespec& lhs, const timespec& rhs) {
  timespec diff;
  diff.tv_sec = lhs.tv_sec - rhs.tv_sec;
  diff.tv_nsec = lhs.tv_nsec - rhs.tv_nsec;
  if (lhs.tv_nsec < rhs.tv_nsec) {
    diff.tv_sec -= 1;
    diff.tv_nsec += 1e9;
  }
  return diff;
}

// wait for up to one second from the given start time
// returns the current time when it exits
timespec wait(long desired_interval_nsec, timespec start) {
  timespec now;
  clock_gettime(CLOCK_MONOTONIC, &now);

  timespec diff = subtract(now, start);
  if (diff.tv_sec == 0) {
    long wait_for = (desired_interval_nsec - diff.tv_nsec) / 1000;
    if (wait_for > 0) {
      usleep(wait_for);
      clock_gettime(CLOCK_MONOTONIC, &now);
    }
  }
  g_frametime = diff.tv_nsec/1000000 + diff.tv_sec*1000;
  return now;
}

int main(int argc, char** argv) {
  enable_fpmath_asserts();

  args_t in = parse_args(argc, argv);
  g_rainbow_enabled = in.rainbow;
  g_view_pressure = in.pressure;

  update_window_size();
  set_window_size_handler(&handle_window_size_changed);

  sim_init(in);

  enable_raw_mode();
  u_clear_screen();
  buffer buf = { 0, 0 };

  draw(&buf);

  timespec interval_start;
  clock_gettime(CLOCK_MONOTONIC, &interval_start);
  while (process_keypress()) {
    sim_step();
    interval_start = wait(1e8, interval_start);
    draw(&buf);
  }

  buffer_free(&buf);
  return 0;
}
