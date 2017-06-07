#include <limits.h>
#include <stdio.h>
#include <stdlib.h>

char* load_file(const char* filename, int* length) {
  FILE* f = fopen(filename, "rb");
  char* contents = NULL;
  long fsize = 0; 
  *length = 0;
  if (!f) {
    goto exit;
  }
  if (fseek(f, 0, SEEK_END)) {
    goto close_exit;
  }
  fsize = ftell(f);
  if (fsize > INT_MAX) {
    goto close_exit;
  }
  *length = fsize;
  if (fseek(f, 0, SEEK_SET)) {
    goto close_exit;
  }
  contents = (char*)malloc(fsize + 1);
  if (!contents) {
    goto close_exit;
  }
  if (fread(contents, fsize, 1, f)) {
    goto close_exit;
  }
  contents[fsize] = '\0';
close_exit:
  fclose(f);
exit:
  return contents;
}

void release_file(char* contents) {
  free(contents);
}
