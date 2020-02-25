#pragma once

#include <stdint.h>

static inline int min_i(int a, int b) {
  return a < b ? a : b;
}

static inline int max_i(int a, int b) {
  return a < b ? b : a;
}

static inline uint8_t min_u8(uint8_t a, uint8_t b) {
  return a < b ? a : b;
}

static inline float invf(float x) {
  return 1.f / x;
}

static inline float clampf(float min, float x, float max) {
  if (x < min) {
    return min;
  } else if (x > max) {
    return max;
  } else {
    return x;
  }
}
