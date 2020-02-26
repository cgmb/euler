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

static void write_stdin(const char* buf, ssize_t count) {
  int result = write(STDIN_FILENO, buf, count);
  (void)result; // todo: handle EINTR
}

static void write_stdinz(const char* buf) {
  write_stdin(buf, strlen(buf));
}

void clear_screen_now() {
  write_stdinz(T_CLEAR T_REPOSITION_CURSOR);
}

void clear_screen(buffer_t* buf) {
  buffer_appendz(buf, T_CLEAR T_REPOSITION_CURSOR);
}

void reposition_cursor(buffer_t* buf) {
  buffer_appendz(buf, T_REPOSITION_CURSOR);
}

void hide_cursor(buffer_t* buf) {
  buffer_appendz(buf, T_HIDE_CURSOR);
}

void show_cursor(buffer_t* buf) {
  buffer_appendz(buf, T_SHOW_CURSOR);
}

void show_cursor_now() {
  write_stdinz(T_SHOW_CURSOR);
}

void die(const char* msg) {
  clear_screen_now();
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
  atexit(show_cursor_now);
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

void buffer_append(buffer_t* buf, const char* s, int len) {
  char* data = (char*)realloc(buf->data, buf->len + len);
  if (!data) {
    die("failed to reallocate buffer");
  }
  memcpy(&data[buf->len], s, len);
  buf->data = data;
  buf->len += len;
}

void buffer_appendz(buffer_t* buf, const char* s) {
  buffer_append(buf, s, strlen(s));
}

void buffer_write(buffer_t* buf) {
  write_stdin(buf->data, buf->len);
}

void buffer_clear(buffer_t* buf) {
  buf->len = 0;
}

void buffer_free(buffer_t* buf) {
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
