#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <signal.h>
#include <termios.h>
#include <unistd.h>
#include <sys/ioctl.h>

#include "terminal.h"

// http://viewsourcecode.org/snaptoken/kilo/02.enteringRawMode.html

static struct termios g_orig_termios;

void u_clear_screen() {
  write(STDIN_FILENO, "\x1b[2J", 4); // clear
  write(STDIN_FILENO, "\x1b[H", 3);  // reposition cursor
}

void clear_screen(buffer* buf) {
  buffer_append(buf, "\x1b[2J", 4); // clear
  buffer_append(buf, "\x1b[H", 3);  // reposition cursor
}

void reposition_cursor(buffer* buf) {
  buffer_append(buf, "\x1b[H", 3);  // reposition cursor
}

void hide_cursor(buffer* buf) {
  buffer_append(buf, "\x1b[?25l", 6);
}

void show_cursor(buffer* buf) {
  buffer_append(buf, "\x1b[?25h", 6);
}

void u_show_cursor() {
  write(STDIN_FILENO, "\x1b[?25h", 6);
}

void die(const char* msg) {
  u_clear_screen();
  perror(msg);
  exit(1);
}

void disable_raw_mode() {
  if (tcsetattr(STDIN_FILENO, TCSAFLUSH, &g_orig_termios) == -1) {
    die("failed to disable raw mode");
  }
}

void enable_raw_mode() {
  atexit(disable_raw_mode);
  atexit(u_show_cursor);
  if (tcgetattr(STDIN_FILENO, &g_orig_termios) == -1) {
    die("failed to enable raw mode");
  }
  struct termios raw = g_orig_termios;
/*
  raw.c_lflag &= ~(ECHO | ICANON);
  raw.c_oflag &= ~(OPOST);
*/
  raw.c_iflag &= ~(BRKINT | ICRNL | INPCK | ISTRIP | IXON);
  raw.c_oflag &= ~(OPOST);
  raw.c_cflag |= (CS8);
  raw.c_lflag &= ~(ECHO | ICANON | IEXTEN | ISIG);
  raw.c_cc[VMIN] = 0;
  raw.c_cc[VTIME] = 0; // tenths of a second
  if (tcsetattr(STDIN_FILENO, TCSAFLUSH, &raw) == -1) {
    die("failed to enable raw mode");
  }
}

void buffer_append(buffer* buf, const char* s, int len) {
  char* data = (char*)realloc(buf->data, buf->len + len);
  if (!data) {
    die("failed to reallocate buffer");
  }
  memcpy(&data[buf->len], s, len);
  buf->data = data;
  buf->len += len;
}

void buffer_append_nchars(buffer* buf, char c, int count) {
  char* data = (char*)realloc(buf->data, buf->len + count);
  if (!data) {
    die("failed to reallocate buffer");
  }
  memset(&data[buf->len], c, count);
  buf->data = data;
  buf->len += count;
}

void buffer_append_char(buffer* buf, char c) {
  buffer_append_nchars(buf, c, 1);
}

void buffer_write(buffer* buf) {
  write(STDIN_FILENO, buf->data, buf->len);
}

void buffer_clear(buffer* buf) {
  buf->len = 0;
}

void buffer_free(buffer* buf) {
  free(buf->data);
}

int get_window_size(int* rows, int* cols) {
  struct winsize ws;
  if (ioctl(STDOUT_FILENO, TIOCGWINSZ, &ws) == -1) {
    return -1;
  } else {
    *cols = ws.ws_col;
    *rows = ws.ws_row;
    return 0;
  }
}

int set_window_size_handler(wshandler_t fn) {
  struct sigaction sa;
  sigemptyset(&sa.sa_mask);
  sa.sa_flags = 0;
  sa.sa_handler = fn;
  return sigaction(SIGWINCH, &sa, 0);
}
