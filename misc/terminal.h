#pragma once

typedef struct buffer {
  char* data;
  int len;
} buffer;

void buffer_append(buffer* buf, const char* s, int len);
void buffer_free(buffer* buf);
void buffer_write(buffer* buf);
void buffer_clear(buffer* buf);

void die(const char* msg);

void u_clear_screen();

void clear_screen(buffer* buf);
void reposition_cursor(buffer* buf);
void hide_cursor(buffer* buf);
void show_cursor(buffer* buf);
void enable_raw_mode();
void disable_raw_mode();

typedef void(*wshandler_t)(int);
int get_window_size(int* rows, int* cols);
int set_window_size_handler(wshandler_t fn);

// ANSI color codes
// https://en.wikipedia.org/wiki/ANSI_escape_code
#define T_BLACK   "\x1B[30m"
#define T_RED     "\x1B[31m"
#define T_GREEN   "\x1B[32m"
#define T_YELLOW  "\x1B[33m"
#define T_BLUE    "\x1B[34m"
#define T_MAGENTA "\x1B[35m"
#define T_CYAN    "\x1B[36m"
#define T_WHITE   "\x1B[37m"

#define T_BG_BLACK   "\x1B[40m"
#define T_BG_RED     "\x1B[41m"
#define T_BG_GREEN   "\x1B[42m"
#define T_BG_YELLOW  "\x1B[43m"
#define T_BG_BLUE    "\x1B[44m"
#define T_BG_MAGENTA "\x1B[45m"
#define T_BG_CYAN    "\x1B[46m"
#define T_BG_WHITE   "\x1B[47m"

#define T_RESET   "\x1B[0m"
