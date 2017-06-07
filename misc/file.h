#ifndef FILE_H
#define FILE_H

char* load_file(const char* filename, int* length);
void release_file(char* contents);

#endif
