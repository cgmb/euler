#ifndef TERMINAL_H
#define TERMINAL_H

struct buffer {
  char* data;
  int len;
};

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

#endif
