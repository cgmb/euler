#include <assert.h>
#include <stdio.h>
#include <string.h>
#include <errno.h>
#include <unistd.h>
#include <math.h>

#include <algorithm>
#include <random>

#include "math/vec2f.h"
#include "misc/terminal.h"

const size_t X = 20;
const size_t Y = 20;

const float k_s = 1; // side length
const float k_invs = 1/k_s;
const float k_d = 1.0; // density
const float k_invd = 1/k_d;

float g_u[Y][X]; // [Y][X-1]
float g_v[Y][X]; // [Y-1][X]
float g_utmp[Y][X]; // [Y][X-1]
float g_vtmp[Y][X]; // [Y-1][X]

bool g_solid[Y][X];

bool g_pause;
unsigned g_simulate_steps;

struct sparse_entry_t {
  int8_t a_diag;
  int8_t a_plus_i;
  int8_t a_plus_j;
} g_a[Y][X];

const size_t N = 8*8*4;
vec2f g_markers[N];
uint8_t g_marker_count[Y][X];
uint8_t g_old_marker_count[Y][X];

const float k_g = -10.f;

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
  memcpy(g_old_marker_count, g_marker_count, sizeof(g_old_marker_count));
  memset(g_marker_count, '\0', sizeof(g_marker_count));
  for (size_t i = 0; i < N; ++i) {
    int x = (int)floorf(g_markers[i].x()*k_invs);
    int y = (int)floorf(g_markers[i].y()*k_invs);
    g_marker_count[y][x]++;
  }
}

bool was_water(size_t y, size_t x) {
  return g_old_marker_count[y][x] > 0;
}

bool is_water(size_t y, size_t x) {
  return g_marker_count[y][x] > 0;
}

void extrapolate_velocity_field() {
  // extrapolate u
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
              total += g_u[y_u][x_u];
              count++;
            }
          }
        }
        assert(count > 0);
        g_u[y][x] = total / count;
      }
    }
  }

  // extrapolate v
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
              total += g_v[y_v][x_v];
              count++;
            }
          }
        }
        assert(count > 0);
        g_v[y][x] = total / count;
      }
    }
  }
}

void sim_init() {
  // setup walls
  for (size_t i = 0; i < X; ++i) {
    g_solid[0][i] = true;
  }
  for (size_t i = 0; i < Y; ++i) {
    g_solid[i][0] = true;
    g_solid[i][X-1] = true;
  }

  std::mt19937 rng_engine(123456789u);
  std::uniform_real_distribution<float> distribution(0.f, 0.5f);
  auto rng = [&](){ return distribution(rng_engine); };

  // setup fluid markers, 4 per cell, jittered
  size_t idx = 0;
  for (size_t i = 1; i < 9; ++i) {
    for (size_t j = 11; j < 19; ++j) {
      for (size_t k = 0; k < 4; ++k) {
        float x = i + (k < 2 ? 0 : 0.5f) + rng();
        float y = j + (k % 2 ? 0 : 0.5f) + rng();
        assert(idx < N);
        g_markers[idx++] = k_s*vec2f{x,y};
      }
    }
  }
  refresh_marker_counts();
}

void advectu(float u[Y][X], float v[Y][X], float dt, float out[Y][X]) {
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
        float prev_x = x - dx*dt*k_invs;
        float prev_y = y - dy*dt*k_invs;


        // bilinear interpolation
        prev_x = clampf(0, prev_x, X-2);
        prev_y = clampf(0, prev_y, Y-1);

        float x_floor;
        float x_frac = modff(prev_x, &x_floor);
        int x_floori = (int)x_floor;
        int x_ceili = std::min(x_floori+1, int(X)-2);

        float y_floor;
        float y_frac = modff(prev_y, &y_floor);
        int y_floori = (int)y_floor;
        int y_ceili = std::min(y_floori+1, int(Y)-1);

        // bilinearly interpolate between the 4 surrounding points, excluding
        // any points that do not have water
        bool w[2][2];
        w[0][0] = is_water(y_floori,x_floori) || is_water(y_floori,x_floori+1);
        w[0][1] = is_water(y_floori,x_ceili) || is_water(y_floori,x_ceili+1);
        w[1][0] = is_water(y_ceili,x_floori) || is_water(y_ceili,x_floori+1);
        w[1][1] = is_water(y_ceili,x_ceili) || is_water(y_ceili,x_ceili+1);

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
        float y1_mid = (1-y1_frac)*u[y_floori][x_floori] + y1_frac*u[y_ceili][x_floori];

        float y2_frac;
        if (!w[0][1]) {
          y2_frac = 1.f;
        } else if (!w[1][1]) {
          y2_frac = 0.f;
        } else {
          y2_frac = y_frac;
        }
        float y2_mid = (1-y2_frac)*u[y_floori][x_ceili] + y2_frac*u[y_ceili][x_ceili];

        float x0_frac;
        if (!w[0][0] && !w[1][0]) {
          x0_frac = 1.f;
        } else if (!w[0][1] && !w[1][1]) {
          x0_frac = 0.f;
        } else {
          x0_frac = x_frac;
        }
        out[y][x] = (1-x0_frac)*y1_mid + x0_frac*y2_mid;
      }
    }
  }
}

void advectv(float u[Y][X], float v[Y][X], float dt,
  float out[Y][X]) {
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
        float prev_x = x - dx*dt*k_invs;
        float prev_y = y - dy*dt*k_invs;


        // bilinear interpolation
        prev_x = clampf(0, prev_x, X-1);
        prev_y = clampf(0, prev_y, Y-2);

        float x_floor;
        float x_frac = modff(prev_x, &x_floor);
        int x_floori = (int)x_floor;
        int x_ceili = std::min(x_floori+1, int(X)-1);

        float y_floor;
        float y_frac = modff(prev_y, &y_floor);
        int y_floori = (int)y_floor;
        int y_ceili = std::min(y_floori+1, int(Y)-2);

        // bilinearly interpolate between the 4 surrounding points, excluding
        // any points that do not have water
        bool w[2][2];
        w[0][0] = is_water(y_floori,x_floori) || is_water(y_floori+1,x_floori);
        w[0][1] = is_water(y_floori,x_ceili) || is_water(y_floori+1,x_ceili);
        w[1][0] = is_water(y_ceili,x_floori) || is_water(y_ceili+1,x_floori);
        w[1][1] = is_water(y_ceili,x_ceili) || is_water(y_ceili+1,x_ceili);

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
        float y1_mid = (1-y1_frac)*v[y_floori][x_floori] + y1_frac*v[y_ceili][x_floori];

        float y2_frac;
        if (!w[0][1]) {
          y2_frac = 1.f;
        } else if (!w[1][1]) {
          y2_frac = 0.f;
        } else {
          y2_frac = y_frac;
        }
        float y2_mid = (1-y2_frac)*v[y_floori][x_ceili] + y2_frac*v[y_ceili][x_ceili];

        float x0_frac;
        if (!w[0][0] && !w[1][0]) {
          x0_frac = 1.f;
        } else if (!w[0][1] && !w[1][1]) {
          x0_frac = 0.f;
        } else {
          x0_frac = x_frac;
        }
        out[y][x] = (1-x0_frac)*y1_mid + x0_frac*y2_mid;
      }
    }
  }
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
  float y1_mid = (1-y1_frac)*u[y_floori][x_floori] + y1_frac*u[y_ceili][x_floori];

  float y2_frac;
  if (!w[0][1]) {
    y2_frac = 1.f;
  } else if (!w[1][1]) {
    y2_frac = 0.f;
  } else {
    y2_frac = y_frac;
  }
  float y2_mid = (1-y2_frac)*u[y_floori][x_ceili] + y2_frac*u[y_ceili][x_ceili];

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
  float y1_mid = (1-y1_frac)*v[y_floori][x_floori] + y1_frac*v[y_ceili][x_floori];

  float y2_frac;
  if (!w[0][1]) {
    y2_frac = 1.f;
  } else if (!w[1][1]) {
    y2_frac = 0.f;
  } else {
    y2_frac = y_frac;
  }
  float y2_mid = (1-y2_frac)*v[y_floori][x_ceili] + y2_frac*v[y_ceili][x_ceili];

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

vec2f velocity_at(vec2f pos) {
  float u_x = pos.x()*k_invs - 1.f;
  float u_y = pos.y()*k_invs - 0.5f;
  float v_x = pos.x()*k_invs - 0.5f;
  float v_y = pos.y()*k_invs - 1.f;

  // out-of-bounds is handled in interpolate
  float x = interpolate_u(g_u, u_y, u_x);
  float y = interpolate_v(g_v, v_y, v_x);
  return vec2f{x,y};
}

// digital differential analyzer collisions
// https://en.wikipedia.org/wiki/Digital_differential_analyzer_(graphics_algorithm)
void advect_markers(float dt) {
  for (size_t i = 0; i < N; ++i) {
    vec2f p = g_markers[i];
    vec2f v = velocity_at(p);
    if (v.x() == 0.f && v.y() == 0.f) {
      continue;
    }
    vec2f e = p + v*dt;
    // step along the line by 1 square at a time, so find the longer direction
    // and scale the vector so that's 1
    const bool shallow_slope = fabsf(v.y()) < fabsf(v.x());
    const vec2f step = shallow_slope ?
      copysignf(1.f,v.x())*vec2f{1,v.y()/v.x()} : copysignf(1.f,v.y())*vec2f{v.x()/v.y(),1};
    const int major_axis = shallow_slope ? 0 : 1; // 0 == X axis, 1 == Y axis
    bool hit = false;
    // the minor velocity direction may have a magnitude of zero,
    // but the major velocity direction cannot be zero
    while ((step[major_axis] > 0.f && p[major_axis] < e[major_axis]) ||
           (step[major_axis] < 0.f && p[major_axis] > e[major_axis])) {
      p += step;
      int x_idx = (int)floorf(p.x()*k_invs);
      int y_idx = (int)floorf(p.y()*k_invs);
      if (g_solid[y_idx][x_idx]) {
        // hit!
        e = p - step;
        hit = true;
        break;
      }
    }
    if (!hit) {
      int x_idx = (int)floorf(e.x()*k_invs);
      int y_idx = (int)floorf(e.y()*k_invs);
      if (g_solid[y_idx][x_idx]) {
        // hit!
        e = p;
      }
    }
    g_markers[i] = e;
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
  return 4 - g_solid[y][x-1] - g_solid[y][x+1]
           - g_solid[y-1][x] - g_solid[y+1][x];
}

int8_t get_a_minus_i(size_t y, size_t x) {
  return x>0 ? g_a[y][x-1].a_plus_i : 0;
}

int8_t get_a_minus_j(size_t y, size_t x) {
  return y>0 ? g_a[y-1][x].a_plus_j : 0;
}

void apply_preconditioner(double r[Y][X], double z[Y][X]) {
  memcpy(z, r, sizeof(double)*Y*X); // skip!
}

double mat_dot(double a[Y][X], double b[Y][X]) {
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

double all_zero(double r[Y][X]) {
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

void get_dense_a(int q[Y*X][Y*X]) {
  for (size_t y = 0; y < Y; ++y) {
    for (size_t x = 0; x < X; ++x) {
      if (is_water(y, x)) {
        q[X*y+x][X*y+x] = g_a[y][x].a_diag;
        if (x > 0)
          q[X*y+x][X*y+x-1] = get_a_minus_i(y,x);
        if (x+1 < X)
          q[X*y+x][X*y+x+1] = g_a[y][x].a_plus_i;
        if (y > 0)
          q[X*y+x][X*(y-1)+x] = get_a_minus_j(y,x);
        if (y+1 < X)
          q[X*y+x][X*(y+1)+x] = g_a[y][x].a_plus_j;
      }
    }
  }
}

void print_fluid_matrix(FILE* f, const char* name, int q[Y*X][Y*X]) {
  fprintf(f, "%s = [", name);
  for (size_t y0 = 0; y0 < Y; ++y0) {
    for (size_t x0 = 0; x0 < X; ++x0) {
      for (size_t y1 = 0; y1 < Y; ++y1) {
        for (size_t x1 = 0; x1 < X; ++x1) {
          if (is_water(y0,x0) && is_water(y1,x1)) {
            fprintf(f, "% 2d ", q[y0*X+x0][y1*X+x1]);
          }
        }
      }
      if (is_water(y0,x0)) {
        fprintf(f, "; ");
      }
    }
  }
  fprintf(f, "]\n");
}

void print_fluid_vector(FILE* f, const char* name, float q[Y][X]) {
  fprintf(f, "%s = [", name);
  for (size_t y = 0; y < Y; ++y) {
    for (size_t x = 0; x < X; ++x) {
      if (is_water(y,x)) {
        fprintf(f, "%f ", q[y][x]);
      }
    }
  }
  fprintf(f, "].'\n");
}

void print_matrix(FILE* f, const char* name, float q[Y][X]) {
  fprintf(f, "%s = [", name);
  for (size_t y = Y; y--;) {
    for (size_t x = 0; x < X; ++x) {
      fprintf(f, "%f ", q[y][x]);
    }
    fprintf(f, ";\n");
  }
  fprintf(f, "]\n");
}

void project(float dt, float u[Y][X], float v[Y][X], float uout[Y][X], float vout[Y][X]) {
  const double c = -k_d*k_s*k_s / dt; // -density * dt^2 / dt
  double d0[Y][X] = {}; // divergence * c
  for (size_t y = 0; y < Y; ++y) {
    for (size_t x = 0; x < X; ++x) {
      if (is_water(y,x)) {
        float up = x>0 ? u[y][x-1] : 0;
        float vp = y>0 ? v[y-1][x] : 0;
        d0[y][x] = c * k_invs * (u[y][x] - up + v[y][x] - vp);
      } else {
        // not really neccessary, but prevents uninitialized reads with memcpy
        d0[y][x] = 0.f;
      }
    }
  }

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
    double z[Y][X];      // auxilliary vector
    apply_preconditioner(r, z);
    double s[Y][X];      // search vector
    memcpy(s, z, sizeof(s));

    double sigma = mat_dot(z,r);
    for (size_t i = 0; i < max_iterations; ++i) {
      apply_a(s, z);

      double alpha = sigma / mat_dot(z,s);
      fmadd(s, alpha, p);  // p += alpha*s
      fmadd(z, -alpha, r); // r -= alpha*z

      if (inf_norm(r) <= tol) {
        break;
      }

      apply_preconditioner(r, z);

      double sigma_new = mat_dot(z,r);
      double beta = sigma_new / sigma;
      update_search(s,z,beta);
      sigma = sigma_new;
    }
  }

  for (size_t y = 0; y < Y; ++y) {
    for (size_t x = 0; x < X-1; ++x) {
      if (g_solid[y][x] || g_solid[y][x+1]) { // todo: unnecessary? probably same as !is_water
        uout[y][x] = 0.f;
      } else if (is_water(y,x) || is_water(y,x+1)) {
        uout[y][x] = u[y][x] - k_invd * k_invs * dt * (p[y][x+1] - p[y][x]);
      } else { // both cells are air
        uout[y][x] = 0.f;
      }
    }
  }

  for (size_t y = 0; y < Y-1; ++y) {
    for (size_t x = 0; x < X; ++x) {
      if (g_solid[y][x] || g_solid[y+1][x]) {
        vout[y][x] = 0.f;
      } else if (is_water(y,x) || is_water(y+1,x)) {
        vout[y][x] = v[y][x] - k_invd * k_invs * dt * (p[y+1][x] - p[y][x]);
      } else { // both cells are air
        vout[y][x] = 0.f;
      }
    }
  }
}

float sq(float x) {
  return x*x;
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

void zero_horizontal_velocity_bounds(float u[Y][X]) {
  for (size_t y = 0; y < Y; ++y) {
    for (size_t x = 0; x < X-1; ++x) {
      // not really neccessary to zero air cells, but makes debugging easier
      bool is_air = !is_water(y,x) && !is_water(y,x+1);
      if (is_air || g_solid[y][x] || g_solid[y][x+1]) {
        u[y][x] = 0.f;
      }
    }
  }
}

void zero_vertical_velocity_bounds(float v[Y][X]) {
  for (size_t y = 0; y < Y-1; ++y) {
    for (size_t x = 0; x < X; ++x) {
      // not really neccessary to zero air cells, but makes debugging easier
      bool is_air = !is_water(y,x) && !is_water(y+1,x);
      if (is_air || g_solid[y][x] || g_solid[y+1][x]) {
        v[y][x] = 0.f;
      }
    }
  }
}

float calculate_timestep(float frame_time, float minimum_step) {
  // Brison suggests a limit of five for stability, but my implementaiton of
  // advection and extrapolation assume that new fluid cells are within one
  // grid cell of old fluid cells
  const float m = 0.75f; // maximum number of cells to traverse in one step

  float dt;
  float max_velocity = sqrtf(sq(maxabs_u(g_u)) + sq(maxabs_v(g_v)));
  if (max_velocity < (m*k_s / frame_time)) {
    dt = frame_time;
  } else {
    dt = std::max(m*k_s / max_velocity, minimum_step);
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
    float dt = calculate_timestep(frame_time, minimum_step);
    frame_time -= dt;

    advect_markers(dt);
    refresh_marker_counts();
    extrapolate_velocity_field();
    zero_horizontal_velocity_bounds(g_utmp);
    zero_vertical_velocity_bounds(g_vtmp);

    advectu(g_u, g_v, dt, g_utmp);
    advectv(g_u, g_v, dt, g_vtmp);
    apply_body_forces(g_vtmp, dt);

    zero_horizontal_velocity_bounds(g_utmp);
    zero_vertical_velocity_bounds(g_vtmp);

    project(dt, g_utmp, g_vtmp, g_u, g_v);
  } while (frame_time > 0.f);

  if (g_simulate_steps) {
    g_simulate_steps--;
  }
}

void draw_rows(struct buffer* buf) {
  const char* symbol[4] = {" ","o","O","0"};
  const uint8_t max_symbol = 3;
  for (int y = Y; y-- > 0;) {
    for (size_t x = 0; x < X; x++) {
      if (g_solid[y][x]) {
        buffer_append(buf, "X", 1);
      } else {
        uint8_t i = std::min(g_marker_count[y][x], max_symbol);
        buffer_append(buf, symbol[i], 1);
      }
    }
    if (y > 0) {
      buffer_append(buf, "\r\n", 2);
    }
  }
}

void refresh_screen(buffer* buf) {
  buffer_clear(buf);

  clear_screen(buf);

  draw_rows(buf);

  buffer_write(buf);
}

void process_keypress() {
  char c = '\0';
  if (read(STDIN_FILENO, &c, 1) == -1 && errno != EAGAIN) {
    die("read");
  }

  if (c == 'p') {
    g_pause = !g_pause;
  } else if (c == 'f') {
    g_simulate_steps++;
  } else if (c == 'q') {
    u_clear_screen();
    exit(0);
  }
}

int main(int argc, char** argv) {
  (void)argc; (void)argv;

  enable_raw_mode();
  buffer buf = { 0, 0 };

  sim_init();

  while (1) {
    refresh_screen(&buf);
    process_keypress();
    sim_step();
  }

  buffer_free(&buf);

  return 0;
}
