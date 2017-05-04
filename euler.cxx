#include <stdio.h>
#include <string.h>
#include <errno.h>
#include <unistd.h>
#include <math.h>

#include <algorithm>
#include <random>

#include "math/vec2f.h"
#include "misc/terminal.h"

const size_t X = 10;
const size_t Y = 10;

const float k_s = 1; // side length
const float k_invs = 1/k_s;
const float k_d = 1.0; // density
const float k_invd = 1/k_d;

float g_p[Y][X];
float g_u[Y][X];
float g_v[Y][X];
float g_utmp[Y][X];
float g_vtmp[Y][X];

bool g_solid[Y][X];

bool g_pause;

struct sparse_entry_t {
  int8_t a_diag;
  int8_t a_plus_i;
  int8_t a_plus_j;
} g_a[Y][X];

const size_t N = 4*4*4;
vec2f g_markers[N];
uint8_t g_marker_count[Y][X];

const float k_g = -9.81;

void sim_init() {
  // setup walls
  for (size_t i = 0; i < X; ++i) {
    g_solid[0][i] = true;
    g_solid[Y-1][i] = true;
  }
  for (size_t i = 0; i < Y; ++i) {
    g_solid[i][0] = true;
    g_solid[i][X-1] = true;
  }

  std::mt19937 rng_engine(123456789u);
  std::uniform_real_distribution<float> distribution(0.f, 0.5f);
  auto rng = [&](){ return distribution(rng_engine); };

  // setup fluid markers
  size_t idx = 0;
  for (size_t i = 1; i < 5; ++i) {
    for (size_t j = Y-2; j > Y-6; --j) {
      for (size_t k = 0; k < 4; ++k) {
        float x = i + (k < 2 ? 0 : 0.5f) + rng();
        float y = j + (k % 2 ? 0 : 0.5f) + rng();
        g_markers[idx++] = k_s*vec2f{x,y};
      }
    }
  }
}

float interpolate(float q[Y][X], float x, float y) {
  float x_floor;
  float x_frac = modff(x, &x_floor);
  int x_floord = (int)x_floor;

  float y_floor;
  float y_frac = modff(y, &y_floor);
  int y_floord = (int)y_floor;

  float y1_mid = (1-y_frac)*q[y_floord][x_floord] + y_frac*q[y_floord+1][x_floord];
  float y2_mid = (1-y_frac)*q[y_floord][x_floord+1] + y_frac*q[y_floord+1][x_floord+1]; // todo: if x or y == Y < 0 || [x,y] > [X,Y]
  return (1-x_frac)*y1_mid + x_frac*y2_mid;
}

void advectu(float u[Y][X], float v[Y][X], float dt,
  float q[Y][X], float out[Y][X]) {
  for (size_t y = 0; y < Y; ++y) {
    for (size_t x = 0; x < X; ++x) {
      float dx = u[y][x];
      float dy = (v[y][x-1] + v[y][x] + v[y+1][x-1] + v[y+1][x]) / 4; // todo: if x == 0, y == Y
      float prev_x = x - dx*dt*k_invs;
      float prev_y = y - dy*dt*k_invs;
      out[y][x] = interpolate(q, prev_x, prev_y);
    }
  }
}

void advectv(float u[Y][X], float v[Y][X], float dt,
  float q[Y][X], float out[Y][X]) {
  for (size_t y = 0; y < Y; ++y) {
    for (size_t x = 0; x < X; ++x) {
      float dx = (u[y-1][x] + u[y][x] + u[y-1][x+1] + u[y][x+1]) / 4; // todo: if x == 0, y == Y
      float dy = v[y][x];
      float prev_x = x - dx*dt*k_invs;
      float prev_y = y - dy*dt*k_invs;
      out[y][x] = interpolate(q, prev_x, prev_y);
    }
  }
}

vec2f velocity_at(vec2f pos) {
  float u_x = pos.x()*k_invs - 0.5f; // todo: x < 0.5
  float u_y = pos.y()*k_invs;        // todo: y < 0
  float v_x = pos.x()*k_invs;        // todo: y < 0
  float v_y = pos.y()*k_invs - 0.5f; // todo: y < 0.5

  float x = interpolate(g_u, u_x, u_y);
  float y = interpolate(g_v, v_x, v_y);
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
  for (size_t y = 0; y < Y; ++y) {
    for (size_t x = 0; x < X; ++x) {
      v[y][x] += k_g * dt;
    }
  }
}

int8_t nonsolid_neighbor_count(size_t y, size_t x) { // todo: X,Y bounds
  return 4 - g_solid[y][x-1] - g_solid[y][x+1]
           - g_solid[y-1][x] - g_solid[y+1][x];
}

bool is_water(size_t y, size_t x) {
  return g_marker_count[y][x] > 0; // todo: X,Y bounds
}

int8_t get_a_minus_i(size_t y, size_t x) {
  return g_a[y][x-1].a_plus_i; // todo: x == 0
}

int8_t get_a_minus_j(size_t y, size_t x) {
  return g_a[y-1][x].a_plus_j; // todo: y == 0
}

void apply_preconditioner(float r[Y][X], float z[Y][X]) {
  memcpy(z, r, sizeof(float)*Y*X); // skip!
}

float mat_dot(float a[Y][X], float b[Y][X]) {
  float total = 0.f;
  for (size_t y = 0; y < Y; ++y) {
    for (size_t x = 0; x < X; ++x) {
      if (is_water(y, x)) {
        total += a[y][x] * b[y][x];
      }
    }
  }
  return total;
}

float all_zero(float r[Y][X]) {
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

float inf_norm(float r[Y][X]) {
  float maximum = r[0][0];
  for (size_t y = 0; y < Y; ++y) {
    for (size_t x = 0; x < X; ++x) {
      if (is_water(y, x)) {
        float a = fabsf(r[y][x]);
        if (a > maximum) {
          maximum = a;
        }
      }
    }
  }
  return maximum;
}

void update_search(float s[Y][X], float z[Y][X], float beta) {
  for (size_t y = 0; y < Y; ++y) {
    for (size_t x = 0; x < X; ++x) {
      if (is_water(y, x)) {
        s[y][x] = z[y][x] + beta * s[y][x];
      }
    }
  }
}

void apply_a(float s[Y][X], float out[Y][X]) {
  for (size_t y = 0; y < Y; ++y) {
    for (size_t x = 0; x < X; ++x) {
      if (is_water(y, x)) {
        out[y][x] = g_a[y][x].a_diag * s[y][x]
                  + g_a[y][x].a_plus_i * s[y][x+1]
                  + g_a[y][x].a_plus_j * s[y+1][x]
                  + get_a_minus_i(y,x) * s[y][x-1]
                  + get_a_minus_j(y,x) * s[y-1][x];
      }
    }
  }
}

// c += a*b
void fmadd(float a[Y][X], float b, float c[Y][X]) {
  for (size_t y = 0; y < Y; ++y) {
    for (size_t x = 0; x < X; ++x) {
      if (is_water(y, x)) {
        c[y][x] += a[y][x] * b;
      }
    }
  }
}

void project(float dt, float u[Y][X], float v[Y][X], float uout[Y][X], float vout[Y][X]) {
  const float c = -k_invd*k_invs*k_invs * dt;
  const float invc = 1/c;
  float d[Y][X] = {}; // divergence / c
  for (size_t y = 0; y < Y; ++y) {
    for (size_t x = 0; x < X; ++x) {
      d[y][x] = -k_invs * invc * (u[y][x] - u[y][x-1] + v[y][x] - v[y-1][x]); // todo: handle x, y < 0
    }
  }
  for (size_t y = 0; y < Y; ++y) {
    for (size_t x = 0; x < X; ++x) {
      if (is_water(y, x)) {
        g_a[y][x].a_diag = nonsolid_neighbor_count(y, x);
        g_a[y][x].a_plus_i = is_water(y,x+1) ? -1 : 0; // todo: x == X
        g_a[y][x].a_plus_j = is_water(y+1,x) ? -1 : 0; // todo: y == Y
      }
    }
  }

  const size_t max_iterations = 100;
  const float tol = 1e-6f;

  // conjugate gradient
  float p[Y][X] = {}; // pressure guess
  float r[Y][X];      // residual
  memcpy(r, d, sizeof(r));
  if (!all_zero(r)) {
    float z[Y][X];      // auxilliary vector
    apply_preconditioner(r, z);
    float s[Y][X];      // search vector
    memcpy(s, z, sizeof(s));

    float sigma = mat_dot(z,r);
    for (size_t i = 0; i < max_iterations; ++i) {
      apply_a(s, z);
      float alpha = sigma / mat_dot(z,s);
      fmadd(s, alpha, p);  // p += alpha*s
      fmadd(z, -alpha, r); // r -= alpha*r
      if (inf_norm(r) <= tol) {
        break;
      }

      apply_preconditioner(r, z);
      memcpy(s, z, sizeof(s));

      float sigma_new = mat_dot(z,r);
      float beta = sigma_new / sigma;
      update_search(s,z,beta);
      sigma = sigma_new;
    }
  }

  for (size_t y = 0; y < Y-1; ++y) { // range right?
    for (size_t x = 0; x < X-1; ++x) { // range right?
      if (is_water(y, x)) {
        uout[y][x] = u[y][x] - k_invd * k_invs * dt * (p[y][x+1] - p[y][x]); // todo: if x == X
        vout[y][x] = v[y][x] - k_invd * k_invs * dt * (p[y+1][x] - p[y][x]); // todo: if y == Y
      }
    }
  }
}

void refresh_marker_counts() {
  memset(g_marker_count, '\0', sizeof(g_marker_count));
  for (size_t i = 0; i < N; ++i) {
    int x = (int)floorf(g_markers[i].x()*k_invs);
    int y = (int)floorf(g_markers[i].y()*k_invs);
    g_marker_count[y][x]++;
  }
}

#include <fenv.h>

void sim_step() {
  if (g_pause) {
    return;
  }

  float dt = 1.f/8;
  advect_markers(dt);
  refresh_marker_counts();

  advectu(g_u, g_v, dt, g_u, g_utmp);
  advectv(g_u, g_v, dt, g_v, g_vtmp);
  apply_body_forces(g_vtmp, dt);
  project(dt, g_utmp, g_vtmp, g_u, g_v);
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
