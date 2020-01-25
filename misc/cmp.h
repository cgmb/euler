#pragma once

#include <stdint.h>

inline float min_f(float a, float b) {
  return a < b ? a : b;
}

inline float max_f(float a, float b) {
  return a < b ? b : a;
}

inline int min_i(int a, int b) {
  return a < b ? a : b;
}

inline int max_i(int a, int b) {
  return a < b ? b : a;
}

inline uint8_t min_u8(uint8_t a, uint8_t b) {
  return a < b ? a : b;
}
