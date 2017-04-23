#include <stdio.h>

#include "math/vec2f.h"
#include "misc/terminal.h"

const size_t X = 10;
const size_t Y = 10;

float g_p[X][Y];
float g_u[X][Y];

float g_q[X][Y];
bool g_solid[X][Y];

const vec2f k_g{0,-9.81};

void sim_init() {
  // setup walls
  for (size_t i = 0; i < X; ++i) {
    g_solid[i][0] = true;
  }
  for (size_t i = 0; i < Y; ++i) {
    g_solid[0][i] = true;
    g_solid[X-1][i] = true;
  }
}

void sim_step() {
}

void draw_rows(struct buffer* buf) {
  for (size_t y = 0; y < Y; y++) {
    buffer_append(buf, "~", 1);
    if (y < Y - 1) {
      buffer_append(buf, "\r\n", 2);
    }
  }
}

#include <errno.h>
#include <unistd.h>

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
  clear_screen(&buf);

  sim_init();

  while (1) {
    refresh_screen(&buf);
    process_keypress();
    sim_step();
  }

  buffer_free(&buf);

  return 0;
}
