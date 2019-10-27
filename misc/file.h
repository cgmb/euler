#pragma once

char* load_file(const char* filename, int* length);
void release_file(char* contents);
