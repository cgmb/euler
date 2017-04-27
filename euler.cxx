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
    for (size_t j = 1; j < 5; ++j) {
      for (size_t k = 0; k < 4; ++k) {
        float x = i + (k < 2 ? 0 : 0.5f) + rng();
        float y = j + (k % 2 ? 0 : 0.5f) + rng();
        g_markers[idx++] = k_s*vec2f{x,y};
      }
    }
  }
}

void sim_step() {
}

void draw_rows(struct buffer* buf) {
  unsigned char marker_count[Y][X] = {};
  for (size_t i = 0; i < N; ++i) {
    size_t x = (size_t)floor(g_markers[i].x());
    size_t y = Y-1-(size_t)floor(g_markers[i].y());
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
