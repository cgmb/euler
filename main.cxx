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
#include "misc/terminal.h"
#include "misc/file.h"
#include "misc/rng.h"
#include "world_grid.h"

// simulation size
const size_t X = 100;
const size_t Y = 40;

// display size
int g_wx;
int g_wy;

struct args_t {
  const char* scenario_file;
  bool rainbow;
};

// simulation constants
const float k_s = 1; // side length
const float k_inv_s = 1/k_s;
const float k_d = 1.0; // density
const float k_inv_d = 1/k_d;
const float k_g = -10.f; // gravity

// All arrays are the same size so functions like bilinear interpolation can
// work on any array. The real size is smaller and is listed in the comments.
float g_u[Y][X]; // [Y][X-1]
float g_v[Y][X]; // [Y-1][X]
float g_utmp[Y][X]; // [Y][X-1]
float g_vtmp[Y][X]; // [Y-1][X]

// grid cell properties from the scenario file
WorldGrid  g_grid;

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
float g_r[Y][X];
float g_g[Y][X];
float g_b[Y][X];
float g_rtmp[Y][X];
float g_gtmp[Y][X];
float g_btmp[Y][X];
const float k_source_color_period = 10.f; // seconds
const float k_initial_color_period = 60.f; // grid cells

// simulation progress
bool g_pause;
uint32_t g_simulate_steps;
uint16_t g_frame_count;

// pressure matrix
struct sparse_entry_t {
  int8_t a_diag;
  int8_t a_plus_i;
  int8_t a_plus_j;
} g_a[Y][X];

// marker particle data
const size_t k_max_marker_count = UINT16_MAX;
std::vector<vec2f> g_markers;
bool g_source_exhausted;
uint16_t g_marker_count[Y][X];

float clampf(float min, float x, float max) {
  if (x < min) {
    return min;
  } else if (x > max) {
    return max;
  } else {
    return x;
  }
}

void refresh_marker_counts() {
  memset(g_marker_count, 0, sizeof(g_marker_count));
  for (size_t i = 0; i < g_markers.size(); ++i) {
    int x = (int)floorf(g_markers[i].x*k_inv_s);
    int y = (int)floorf(g_markers[i].y*k_inv_s);
    bool in_bounds = x > 0 && x < (int)X && y > 0 && y < (int)Y;
    if (in_bounds) {
      if (g_grid.is_sink({size_t(x),size_t(y)})) {
        g_markers[i] = g_markers.back();
        g_markers.pop_back();
      } else {
        g_marker_count[y][x]++;
      }
    } else {
      g_markers[i] = g_markers.back();
      g_markers.pop_back();
    }
  }

  g_grid.mark_fluid_as_old();
  for (size_t y = 0; y < Y; ++y) {
    for (size_t x = 0; x < X; ++x) {
      g_grid.set_fluid({x,y}, g_marker_count[y][x] > 0);
    }
  }
}

bool was_water(size_t y, size_t x) {
  return g_grid.was_fluid({x,y});
}

bool is_water(size_t y, size_t x) {
  return g_grid.is_fluid({x,y});
}

void extrapolate_u(float u[Y][X]) {
  for (size_t y = 0; y < Y; ++y) {
    for (size_t x = 0; x < X-1; ++x) {
      bool was_invalid = !was_water(y,x) && !was_water(y,x+1);
      if (was_invalid && (is_water(y,x) || is_water(y,x+1))) {
        size_t y_lower = y>0 ? y-1 : 0;
        size_t x_lower = x>0 ? x-1 : 0;
        size_t y_upper = y+1<Y ? y+1 : y;
        size_t x_upper = x+2<X ? x+1 : x;
        float total = 0.f;
        size_t count = 0;
        for (size_t y_u = y_lower; y_u <= y_upper; ++y_u) {
          for (size_t x_u = x_lower; x_u <= x_upper; ++x_u) {
            if (was_water(y_u,x_u) || was_water(y_u,x_u+1)) {
              total += u[y_u][x_u];
              count++;
            }
          }
        }
        assert(count > 0);
        u[y][x] = total / count;
      }
    }
  }
}

void extrapolate_v(float v[Y][X]) {
  for (size_t y = 0; y < Y-1; ++y) {
    for (size_t x = 0; x < X; ++x) {
      bool was_invalid = !was_water(y,x) && !was_water(y+1,x);
      if (was_invalid && (is_water(y,x) || is_water(y+1,x))) {
        size_t y_lower = y>0 ? y-1 : 0;
        size_t x_lower = x>0 ? x-1 : 0;
        size_t y_upper = y+2<Y ? y+1 : y;
        size_t x_upper = x+1<X ? x+1 : x;
        float total = 0.f;
        size_t count = 0;
        for (size_t y_v = y_lower; y_v <= y_upper; ++y_v) {
          for (size_t x_v = x_lower; x_v <= x_upper; ++x_v) {
            if (was_water(y_v,x_v) || was_water(y_v+1,x_v)) {
              total += v[y_v][x_v];
              count++;
            }
          }
        }
        assert(count > 0);
        v[y][x] = total / count;
      }
    }
  }
}

void extrapolate_p(float q[Y][X]) {
  // extrapolate a quantity sampled in the pressure cell
  for (size_t y = 0; y < Y; ++y) {
    for (size_t x = 0; x < X; ++x) {
      if (!was_water(y,x) && is_water(y,x)) {
        size_t y_lower = y>0 ? y-1 : 0;
        size_t x_lower = x>0 ? x-1 : 0;
        size_t y_upper = y<Y-1 ? y+1 : y;
        size_t x_upper = x<X-1 ? x+1 : x;
        float total = 0.f;
        size_t count = 0;
        for (size_t y_u = y_lower; y_u <= y_upper; ++y_u) {
          for (size_t x_u = x_lower; x_u <= x_upper; ++x_u) {
            if (was_water(y_u,x_u)) {
              total += q[y_u][x_u];
              count++;
            }
          }
        }
        assert(count > 0);
        q[y][x] = total / count;
      }
    }
  }
}

void colorize() {
  for (size_t y = 0; y < Y; ++y) {
    for (size_t x = 0; x < X; ++x) {
      if (is_water(y,x)) {
        float t = 0.f;
        if (!g_grid.is_source({x,y})) {
          t = (x + y) * 6.f / k_initial_color_period;
        }
        g_r[y][x] = hsv_basis(t + 2.f);
        g_g[y][x] = hsv_basis(t);
        g_b[y][x] = hsv_basis(t - 2.f);
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
  g_grid = WorldGrid::from_file(in.scenario_file, X, Y);

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
      if (g_grid.is_source({x,y})) {
        if (!g_source_exhausted && g_marker_count[y][x] < 4) {
          g_markers.push_back(index_to_world(vec2f(x+randf(), y+randf())));
          g_marker_count[y][x]++;
          g_source_exhausted |= (g_markers.size() == k_max_marker_count);
        }
        g_r[y][x] = hsv_basis(t + 2.f);
        g_g[y][x] = hsv_basis(t);
        g_b[y][x] = hsv_basis(t - 2.f);
      }
    }
  }
}

float bilinear(bool w[2][2], float q[Y][X],
    float x_frac, int x_floori, int x_ceili,
    float y_frac, int y_floori, int y_ceili) {
  // uses bilinear interpolation to combine up to four samples from q
  // at least one sample must be valid, as specified by w
  assert(w[0][0] || w[0][1] || w[1][0] || w[1][1]);

  // note that y1_mid or y2_mid may be wrong if both the points they're comprised of
  // are out of water, but the x interpolation will take care of that.
  float y1_frac;
  if (!w[0][0]) {
    y1_frac = 1.f;
  } else if (!w[1][0]) {
    y1_frac = 0.f;
  } else {
    y1_frac = y_frac;
  }
  float y1_mid = (1-y1_frac)*q[y_floori][x_floori] + y1_frac*q[y_ceili][x_floori];

  float y2_frac;
  if (!w[0][1]) {
    y2_frac = 1.f;
  } else if (!w[1][1]) {
    y2_frac = 0.f;
  } else {
    y2_frac = y_frac;
  }
  float y2_mid = (1-y2_frac)*q[y_floori][x_ceili] + y2_frac*q[y_ceili][x_ceili];

  float x0_frac;
  if (!w[0][0] && !w[1][0]) {
    x0_frac = 1.f;
  } else if (!w[0][1] && !w[1][1]) {
    x0_frac = 0.f;
  } else {
    x0_frac = x_frac;
  }
  return (1-x0_frac)*y1_mid + x0_frac*y2_mid;
}

float interpolate_u(float u[Y][X], float y, float x) {
  // bilinear interpolation
  x = clampf(0, x, X-2);
  y = clampf(0, y, Y-1);

  float x_floor;
  float x_frac = modff(x, &x_floor);
  int x_floori = (int)x_floor;
  int x_ceili = std::min(x_floori+1, int(X)-2);

  float y_floor;
  float y_frac = modff(y, &y_floor);
  int y_floori = (int)y_floor;
  int y_ceili = std::min(y_floori+1, int(Y)-1);

  // bilinearly interpolate between the 4 surrounding points, excluding
  // any points that do not have water
  bool w[2][2];
  w[0][0] = is_water(y_floori,x_floori) || is_water(y_floori,x_floori+1);
  w[0][1] = is_water(y_floori,x_ceili) || is_water(y_floori,x_ceili+1);
  w[1][0] = is_water(y_ceili,x_floori) || is_water(y_ceili,x_floori+1);
  w[1][1] = is_water(y_ceili,x_ceili) || is_water(y_ceili,x_ceili+1);

  return bilinear(w, u, x_frac, x_floori, x_ceili, y_frac, y_floori, y_ceili);
}

float interpolate_v(float v[Y][X], float y, float x) {
  // bilinear interpolation
  x = clampf(0, x, X-1);
  y = clampf(0, y, Y-2);

  float x_floor;
  float x_frac = modff(x, &x_floor);
  int x_floori = (int)x_floor;
  int x_ceili = std::min(x_floori+1, int(X)-1);

  float y_floor;
  float y_frac = modff(y, &y_floor);
  int y_floori = (int)y_floor;
  int y_ceili = std::min(y_floori+1, int(Y)-2);

  // bilinearly interpolate between the 4 surrounding points, excluding
  // any points that do not have water
  bool w[2][2];
  w[0][0] = is_water(y_floori,x_floori) || is_water(y_floori+1,x_floori);
  w[0][1] = is_water(y_floori,x_ceili) || is_water(y_floori+1,x_ceili);
  w[1][0] = is_water(y_ceili,x_floori) || is_water(y_ceili+1,x_floori);
  w[1][1] = is_water(y_ceili,x_ceili) || is_water(y_ceili+1,x_ceili);

  return bilinear(w, v, x_frac, x_floori, x_ceili, y_frac, y_floori, y_ceili);
}

float interpolate_p(float q[Y][X], float y, float x) {
  // bilinear interpolation
  x = clampf(0, x, X-1);
  y = clampf(0, y, Y-1);

  float x_floor;
  float x_frac = modff(x, &x_floor);
  int x_floori = (int)x_floor;
  int x_ceili = std::min(x_floori+1, int(X)-1);

  float y_floor;
  float y_frac = modff(y, &y_floor);
  int y_floori = (int)y_floor;
  int y_ceili = std::min(y_floori+1, int(Y)-1);

  // bilinearly interpolate between the 4 surrounding points, excluding
  // any points that do not have water
  bool w[2][2];
  w[0][0] = is_water(y_floori,x_floori);
  w[0][1] = is_water(y_floori,x_ceili);
  w[1][0] = is_water(y_ceili,x_floori);
  w[1][1] = is_water(y_ceili,x_ceili);

  return bilinear(w, q, x_frac, x_floori, x_ceili, y_frac, y_floori, y_ceili);
}

void advect_u(float u[Y][X], float v[Y][X], float dt, float out[Y][X]) {
  for (size_t y = 0; y < Y; ++y) {
    for (size_t x = 0; x < X-1; ++x) {
      if (is_water(y,x) || is_water(y,x+1)) {
        float dx = u[y][x];
        // we could check if either of the two squares around each
        // vertical component is water, but that would mean checking
        // four squares instead of two, and doing bounds checks.
        // Let's just consider two because I'm not sure that allowing
        // fluid to interact diagonally like that actually helps.
        float div = 0.f;
        float vv[4] = {};
        if (is_water(y,x)) {
          vv[0] = v[y][x];
          div += 1.f;
          if (y > 0) {
            vv[1] = v[y-1][x];
            div += 1.f;
          }
        }
        if (is_water(y,x+1)) {
          vv[2] = v[y][x+1];
          div += 1.f;
          if (y > 0) {
            vv[3] = v[y-1][x+1];
            div += 1.f;
          }
        }
        float dy = (vv[0] + vv[1] + vv[2] + vv[3]) / div;
        float prev_x = x - dx*dt*k_inv_s;
        float prev_y = y - dy*dt*k_inv_s;

        out[y][x] = interpolate_u(g_u, prev_y, prev_x);
      }
    }
  }
}

void advect_v(float u[Y][X], float v[Y][X], float dt, float out[Y][X]) {
  for (size_t y = 0; y < Y-1; ++y) {
    for (size_t x = 0; x < X; ++x) {
      if (is_water(y,x) || is_water(y+1,x)) {
        float dy = v[y][x];
        // we could check if either of the two squares around each
        // horizontal component is water, but that would mean checking
        // four squares instead of two, and doing bounds checks.
        // Let's just consider two because I'm not sure that allowing
        // fluid to interact diagonally like that actually helps.
        float div = 0.f;
        float uu[4] = {};
        if (is_water(y,x)) {
          uu[0] = u[y][x];
          div += 1.f;
          if (x > 0) {
            uu[1] = u[y][x-1];
            div += 1.f;
          }
        }
        if (is_water(y+1,x)) {
          uu[2] = u[y+1][x];
          div += 1.f;
          if (x > 0) {
            uu[3] = u[y+1][x-1];
            div += 1.f;
          }
        }
        float dx = (uu[0] + uu[1] + uu[2] + uu[3]) / div;
        float prev_x = x - dx*dt*k_inv_s;
        float prev_y = y - dy*dt*k_inv_s;

        out[y][x] = interpolate_v(g_v, prev_y, prev_x);
      }
    }
  }
}

void advect_p(float q[Y][X], float u[Y][X], float v[Y][X], float dt, float out[Y][X]) {
  for (size_t y = 0; y < Y; ++y) {
    for (size_t x = 0; x < X; ++x) {
      if (is_water(y,x)) {
        // let's assume we're never advecting something on the grid boundary,
        // as my pressure solve doesn't check those bounds either
        float dy = (v[y][x] + v[y-1][x]) / 2;
        float dx = (u[y][x] + u[y][x-1]) / 2;
        float prev_x = x - dx*dt*k_inv_s;
        float prev_y = y - dy*dt*k_inv_s;

        out[y][x] = interpolate_p(q, prev_y, prev_x);
      }
    }
  }
}

vec2f velocity_at(vec2f pos) {
  float u_x = pos.x*k_inv_s - 1.f;
  float u_y = pos.y*k_inv_s - 0.5f;
  float v_x = pos.x*k_inv_s - 0.5f;
  float v_y = pos.y*k_inv_s - 1.f;

  // out-of-bounds is handled in interpolate
  float x = interpolate_u(g_u, u_y, u_x);
  float y = interpolate_v(g_v, v_y, v_x);
  return vec2f{x,y};
}

float time_to(vec2f p0, vec2f p1, vec2f v, int axis) {
  // p1 = p0 + v*t
  // t = (p1 - p0) / v
  float t;
  if (fabsf(v[axis]) > 0.f) {
    t = (p1[axis] - p0[axis]) / v[axis];
  } else {
    t = FLT_MAX;
  }
  return t;
}

// Collision detection was developed ad-hoc, but I should probably read
// A Fast Voxel Traversal Algorithm for Ray Tracing (1987)
void advect_markers(float dt) {
  for (size_t i = 0; i < g_markers.size(); ++i) {
    vec2f p = g_markers[i];
    vec2f v = velocity_at(p);
    vec2f np;

    int x_idx = (int)floorf(p.x*k_inv_s);
    int y_idx = (int)floorf(p.y*k_inv_s);

    // next horizontal intersect
    int x_dir = v.x > 0 ? 1 : -1;
    int nx_idx = x_idx + (v.x > 0 ? 1 : 0);
    np[0] = nx_idx*k_s;
    float t_x = time_to(p, np, v, 0);
    // at idx = x, we're on the boundary between x-1 and x
    // if we're going left, the pressure cell we care about is x-1
    // if we're going right, the pressure cell we care about is x
    int x_idx_offset = v.x < 0 ? -1 : 0;

    // next vertical intersect
    int y_dir = v.y > 0 ? 1 : -1;
    int ny_idx = y_idx + (v.y > 0 ? 1 : 0);
    np[1] = ny_idx*k_s;
    float t_y = time_to(p, np, v, 1);
    // at idx = y, we're on the boundary between y-1 and y
    // if we're going down, the pressure cell we care about is y-1
    // if we're going up, the pressure cell we care about is y
    int y_idx_offset = v.y < 0 ? -1 : 0;

    float t_prev = 0.f;
    float t_near = std::min(t_x, t_y);
    while (t_near < dt) {
      if (t_x < t_y) {
        // entered new horizontal cell
        if (g_grid.is_solid({nx_idx + x_idx_offset, y_idx})) {
          // hit! we're done going horizontal
          p += v * t_prev;
          dt -= t_prev;
          t_near = 0;
          v[0] = 0.f;
          t_x = FLT_MAX;
          t_y = time_to(p, np, v, 1);
        } else {
          // calculate next intersection
          x_idx = nx_idx;
          nx_idx = x_idx + x_dir;
          np[0] = nx_idx*k_s;
          t_x = time_to(p, np, v, 0);
        }
      } else {
        // entered new vertical cell
        if (g_grid.is_solid({x_idx, y_idx + y_idx_offset})) {
          // hit! we're done going vertical
          p += v * t_prev;
          dt -= t_prev;
          t_near = 0;
          v[1] = 0.f;
          t_y = FLT_MAX;
          t_x = time_to(p, np, v, 0);
        } else {
          // calculate next intersection
          y_idx = ny_idx;
          ny_idx = y_idx + y_dir;
          np[1] = ny_idx*k_s;
          t_y = time_to(p, np, v, 1);
        }
      }
      t_prev = t_near;
      t_near = std::min(t_x, t_y);
    }
    if (t_near < FLT_MAX) {
      g_markers[i] = p + v*dt;
    } else {
      g_markers[i] = p + v*t_prev;
    }
  }
}

void apply_body_forces(float v[Y][X], float dt) {
  for (size_t y = 0; y < Y-1; ++y) {
    for (size_t x = 0; x < X; ++x) {
      v[y][x] += k_g * dt;
    }
  }
}

int8_t nonsolid_neighbor_count(size_t y, size_t x) {
  // this function is only used on fluid cells, and the edge cells
  // should never be fluid, so no bounds checks required
  return 4 - g_grid.is_solid({x-1, y}) - g_grid.is_solid({x+1, y})
           - g_grid.is_solid({x, y-1}) - g_grid.is_solid({x, y+1});
}

int8_t get_a_minus_i(size_t y, size_t x) {
  return x>0 ? g_a[y][x-1].a_plus_i : 0;
}

int8_t get_a_minus_j(size_t y, size_t x) {
  return y>0 ? g_a[y-1][x].a_plus_j : 0;
}

double sq(double x) {
  return x*x;
}

double g_precon[Y][X];
double g_q[Y][X];

void apply_preconditioner(double r[Y][X], double z[Y][X]) {
  // Incomplete Cholesky
  // A ~= LLᵀ
  // L = F * E_inv + E

  // calculate E_inv (precon)
  for (size_t y = 0; y < Y; ++y) {
    for (size_t x = 0; x < X; ++x) {
      if (is_water(y, x)) {
        double a = g_a[y][x].a_diag;
        double b = sq(get_a_minus_i(y,x) * g_precon[y][x-1]);
        double c = sq(get_a_minus_j(y,x) * g_precon[y-1][x]);

        double e = a - b - c;
        if (e < 0.25*a) {
          e = a != 0 ? a : 1;
        }
        g_precon[y][x] = 1 / sqrt(e);
      }
    }
  }

  // solve Lq = r
  memset(g_q, '\0', sizeof(double)*Y*X);
  for (size_t y = 0; y < Y; ++y) {
    for (size_t x = 0; x < X; ++x) {
      if (is_water(y, x)) {
        double t = r[y][x]
          - g_a[y][x-1].a_plus_i * g_precon[y][x-1] * g_q[y][x-1]
          - g_a[y-1][x].a_plus_j * g_precon[y-1][x] * g_q[y-1][x];
        g_q[y][x] = t * g_precon[y][x];
      }
    }
  }

  // solve Lᵀz = q
  memset(z, '\0', sizeof(double)*Y*X);
  for (size_t y = Y; y--;) {
    for (size_t x = X; x--;) {
      if (is_water(y, x)) {
        double t = g_q[y][x]
          - g_a[y][x].a_plus_i * g_precon[y][x] * z[y][x+1]
          - g_a[y][x].a_plus_j * g_precon[y][x] * z[y+1][x];
        z[y][x] = t * g_precon[y][x];
      }
    }
  }
}

double dot(double a[Y][X], double b[Y][X]) {
  double total = 0.f;
  for (size_t y = 0; y < Y; ++y) {
    for (size_t x = 0; x < X; ++x) {
      if (is_water(y, x)) {
        total += a[y][x] * b[y][x];
      }
    }
  }
  return total;
}

bool all_zero(double r[Y][X]) {
  for (size_t y = 0; y < Y; ++y) {
    for (size_t x = 0; x < X; ++x) {
      if (is_water(y, x)) {
        if (r[y][x] != 0.f) {
          return false;
        }
      }
    }
  }
  return true;
}

double inf_norm(double r[Y][X]) {
  double maximum = 0.f;
  for (size_t y = 0; y < Y; ++y) {
    for (size_t x = 0; x < X; ++x) {
      if (is_water(y, x)) {
        double a = fabs(r[y][x]);
        if (a > maximum) {
          maximum = a;
        }
      }
    }
  }
  return maximum;
}

void update_search(double s[Y][X], double z[Y][X], double beta) {
  for (size_t y = 0; y < Y; ++y) {
    for (size_t x = 0; x < X; ++x) {
      if (is_water(y, x)) {
        s[y][x] = z[y][x] + beta * s[y][x];
      }
    }
  }
}

void apply_a(double s[Y][X], double out[Y][X]) {
  for (size_t y = 0; y < Y; ++y) {
    for (size_t x = 0; x < X; ++x) {
      if (is_water(y, x)) {
        out[y][x] = g_a[y][x].a_diag * s[y][x]
                  + g_a[y][x].a_plus_i * (x+1 < X && is_water(y,x+1) ? s[y][x+1] : 0)
                  + g_a[y][x].a_plus_j * (y+1 < Y && is_water(y+1,x) ? s[y+1][x] : 0)
                  + get_a_minus_i(y,x) * (x>0 && is_water(y,x-1) ? s[y][x-1] : 0)
                  + get_a_minus_j(y,x) * (y>0 && is_water(y-1,x) ? s[y-1][x] : 0);
      }
    }
  }
}

// c += a*b
void fmadd(double a[Y][X], double b, double c[Y][X]) {
  for (size_t y = 0; y < Y; ++y) {
    for (size_t x = 0; x < X; ++x) {
      if (is_water(y, x)) {
        c[y][x] += a[y][x] * b;
      }
    }
  }
}

void project(float dt, float u[Y][X], float v[Y][X], float uout[Y][X], float vout[Y][X]) {
  const double c = -k_d*k_s*k_s / dt; // -density * dt^2 / dt
  double d0[Y][X] = {}; // divergence * c

  // calculate d0
  for (size_t y = 0; y < Y; ++y) {
    for (size_t x = 0; x < X; ++x) {
      if (is_water(y,x)) {
        float up = x>0 ? u[y][x-1] : 0;
        float vp = y>0 ? v[y-1][x] : 0;
        d0[y][x] = c * k_inv_s * (u[y][x] - up + v[y][x] - vp);
      }
    }
  }

  // calculate A
  for (size_t y = 0; y < Y; ++y) {
    for (size_t x = 0; x < X; ++x) {
      if (is_water(y, x)) {
        g_a[y][x].a_diag = nonsolid_neighbor_count(y, x);
        g_a[y][x].a_plus_i = x<X-1 && is_water(y,x+1) ? -1 : 0;
        g_a[y][x].a_plus_j = y<Y-1 && is_water(y+1,x) ? -1 : 0;
      }
    }
  }

  const size_t max_iterations = 100;
  const double tol = 1e-6f;

  // conjugate gradient
  double p[Y][X] = {}; // pressure guess
  double r[Y][X];      // residual
  memcpy(r, d0, sizeof(r));
  if (!all_zero(r)) {
    double z[Y][X];      // auxiliary vector
    apply_preconditioner(r, z);
    double s[Y][X];      // search vector
    memcpy(s, z, sizeof(s));

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
      if (is_water(y,x) && p[y][x] < 0.f) {
        p[y][x] = 0.f;
      }
    }
  }

  // update horizontal velocities
  for (size_t y = 0; y < Y; ++y) {
    for (size_t x = 0; x < X-1; ++x) {
      if (g_grid.is_solid({x,y}) || g_grid.is_solid({x+1,y})) {
        uout[y][x] = 0.f;
      } else if (is_water(y,x) || is_water(y,x+1)) {
        uout[y][x] = u[y][x] - k_inv_d * k_inv_s * dt * (p[y][x+1] - p[y][x]);
      } else { // both cells are air
        uout[y][x] = 0.f;
      }
    }
  }

  // update vertical velocities
  for (size_t y = 0; y < Y-1; ++y) {
    for (size_t x = 0; x < X; ++x) {
      if (g_grid.is_solid({x,y}) || g_grid.is_solid({x,y+1})) {
        vout[y][x] = 0.f;
      } else if (is_water(y,x) || is_water(y+1,x)) {
        vout[y][x] = v[y][x] - k_inv_d * k_inv_s * dt * (p[y+1][x] - p[y][x]);
      } else { // both cells are air
        vout[y][x] = 0.f;
      }
    }
  }
}

float maxabs_u(float q[Y][X]) {
  float max = 0;
  for (size_t y = 0; y < Y; ++y) {
    for (size_t x = 0; x < X-1; ++x) {
      float value = fabsf(q[y][x]);
      if (value > max) {
        max = value;
      }
    }
  }
  return max;
}

float maxabs_v(float q[Y][X]) {
  float max = 0;
  for (size_t y = 0; y < Y-1; ++y) {
    for (size_t x = 0; x < X; ++x) {
      float value = fabsf(q[y][x]);
      if (value > max) {
        max = value;
      }
    }
  }
  return max;
}

void zero_bounds_u(float u[Y][X]) {
  for (size_t y = 0; y < Y; ++y) {
    for (size_t x = 0; x < X-1; ++x) {
      // not really necessary to zero air cells, but makes debugging easier
      bool is_air = !is_water(y,x) && !is_water(y,x+1);
      if (is_air || g_grid.is_solid({x,y}) || g_grid.is_solid({x+1,y})) {
        u[y][x] = 0.f;
      }
    }
  }
}

void zero_bounds_v(float v[Y][X]) {
  for (size_t y = 0; y < Y-1; ++y) {
    for (size_t x = 0; x < X; ++x) {
      // not really necessary to zero air cells, but makes debugging easier
      bool is_air = !is_water(y,x) && !is_water(y+1,x);
      if (is_air || g_grid.is_solid({x,y}) || g_grid.is_solid({x,y+1})) {
        v[y][x] = 0.f;
      }
    }
  }
}

float calculate_timestep(float frame_time) {
  // Bridson suggests a limit of five for stability, but my implementation of
  // advection and extrapolation assume that new fluid cells are within one
  // grid cell of old fluid cells
  const float m = 0.75f; // maximum number of cells to traverse in one step

  float dt;
  float max_velocity = sqrtf(sq(maxabs_u(g_u)) + sq(maxabs_v(g_v)));
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
  const float minimum_step = total_frame_time / 8;
  float frame_time = total_frame_time;
  do {
    float dt = calculate_timestep(frame_time);
    assert(minimum_step <= dt || frame_time == dt);
    frame_time -= dt;

    advect_markers(dt);
    refresh_marker_counts();
    if (g_rainbow_enabled) {
      extrapolate_p(g_r);
      extrapolate_p(g_g);
      extrapolate_p(g_b);
    }
    update_fluid_sources();
    extrapolate_u(g_u);
    extrapolate_v(g_v);
    zero_bounds_u(g_u);
    zero_bounds_v(g_v);

    advect_u(g_u, g_v, dt, g_utmp);
    advect_v(g_u, g_v, dt, g_vtmp);
    if (g_rainbow_enabled) {
      advect_p(g_r, g_u, g_v, dt, g_rtmp);
      memcpy(g_r, g_rtmp, sizeof(g_r));

      advect_p(g_g, g_u, g_v, dt, g_gtmp);
      memcpy(g_g, g_gtmp, sizeof(g_g));

      advect_p(g_b, g_u, g_v, dt, g_btmp);
      memcpy(g_b, g_btmp, sizeof(g_b));
    }
    apply_body_forces(g_vtmp, dt);

    zero_bounds_u(g_utmp);
    zero_bounds_v(g_vtmp);

    project(dt, g_utmp, g_vtmp, g_u, g_v);
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

int sprint_color_code(char* buf, float r, float g, float b) {
  // buf must be at least 20 characters
  int r_out = float_to_byte_color(linear_to_sRGB(r));
  int g_out = float_to_byte_color(linear_to_sRGB(g));
  int b_out = float_to_byte_color(linear_to_sRGB(b));
  return sprintf(buf, "\x1B[38;2;%d;%d;%dm", r_out, g_out, b_out);
}

void buffer_appendz(buffer* buf, const char* s) {
  buffer_append(buf, s, strlen(s));
}

void draw_rows(buffer* buf) {
  const char* symbol[4] = {" ","o","O","0"};
  const uint16_t max_symbol_idx = 3;
  const int y_cutoff = std::max((int)Y-1 - g_wy, 1);
  for (int y = Y-1; y-- > y_cutoff;) {
    bool prev_water = false;
    for (int x = 1; x < (int)X-1 && x < g_wx+1; x++) {
      if (g_grid.is_solid({x,y})) {
        if (prev_water) {
          buffer_appendz(buf, T_RESET);
        }
        buffer_appendz(buf, "X");
        prev_water = false;
      } else if (g_grid.is_sink({x,y})) {
        if (prev_water) {
          buffer_appendz(buf, T_RESET);
        }
        buffer_appendz(buf, "=");
      } else {
        uint16_t i = std::min(g_marker_count[y][x], max_symbol_idx);
        bool has_water = i > 0;
        if (!prev_water && has_water && !g_rainbow_enabled) {
          buffer_appendz(buf, T_BLUE);
        } else if (has_water && g_rainbow_enabled) {
          char tmp[20];
          int length = sprint_color_code(tmp, g_r[y][x], g_g[y][x], g_b[y][x]);
          if (length < 0) {
            die("sprintf");
          }
          buffer_append(buf, tmp, length);
        } else if (prev_water && !has_water) {
          buffer_appendz(buf, T_RESET);
        }
        buffer_appendz(buf, symbol[i]);
        prev_water = has_water;
      }
    }
    buffer_appendz(buf, T_RESET);
    buffer_appendz(buf, "\x1b[K"); // clear remainer of line
    if (y > y_cutoff) {
      buffer_appendz(buf, "\r\n");
    }
  }
}

void draw(buffer* buf) {
  buffer_clear(buf);
  reposition_cursor(buf);
  draw_rows(buf);
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
  } else if (c == 'q') {
    u_clear_screen();
    return false;
  }
  return true;
}

args_t parse_args(int argc, char** argv) {
  args_t in;
  in.rainbow = false;

  if (argc < 2) {
    fprintf(stderr, "usage: %s [--rainbow] <scenario>\n", argv[0]);
    exit(1);
  }

  for (int i = 1; i < argc - 1; ++i) {
    if (!strcmp(argv[i], "--rainbow")) {
      in.rainbow = true;
    } else {
      fprintf(stderr, "Unrecognized input: %s\n", argv[i]);
      exit(1);
    }
  }

  in.scenario_file = argv[argc-1];
  return in;
}

void update_window_size() {
  if (get_window_size(&g_wy, &g_wx) == -1) {
    die("get_window_size");
  }
}

void handle_window_size_changed(int signal) {
  (void)signal;
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

  return now;
}

int main(int argc, char** argv) {
  args_t in = parse_args(argc, argv);
  g_rainbow_enabled = in.rainbow;

  update_window_size();
  set_window_size_handler(&handle_window_size_changed);

  enable_raw_mode();
  u_clear_screen();
  buffer buf = { 0, 0 };

  sim_init(in);
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
