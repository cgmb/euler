#include <stdio.h>
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

const size_t N = 4*4*4;
vec2f g_markers[N];

const float k_g = -9.81;

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
  std::uniform_real_distribution<float> distribution(0.f, 0.25f);
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
    while (step[major_axis] > 0.f && p[major_axis] < e[major_axis] ||
           step[major_axis] < 0.f && p[major_axis] > e[major_axis]) {
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

void sim_step() {
  float dt = 1.f/8;
  advect_markers(dt);
  apply_body_forces(g_v, dt);
}

void draw_rows(struct buffer* buf) {
  unsigned char marker_count[Y][X] = {};
  for (size_t i = 0; i < N; ++i) {
    int x = (int)floorf(g_markers[i].x()*k_invs);
    int y = (int)floorf(g_markers[i].y()*k_invs);
    marker_count[y][x]++;
  }
  const char* symbol[4] = {" ","o","O","0"};
  const unsigned char max_symbol = 3;
  for (int y = Y; y-- > 0;) {
    for (size_t x = 0; x < X; x++) {
      if (g_solid[y][x]) {
        buffer_append(buf, "X", 1);
      } else {
        unsigned char i = std::min(marker_count[y][x], max_symbol);
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
  hide_cursor(buf);

  draw_rows(buf);

  show_cursor(buf);
  buffer_write(buf);
}

void process_keypress() {
  char c = '\0';
  if (read(STDIN_FILENO, &c, 1) == -1 && errno != EAGAIN) {
    die("read");
  }
  
  if (c == 'q') {
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
