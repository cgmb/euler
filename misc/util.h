#pragma once

inline float sqf(float x) {
  return x*x;
}

inline double sq(double x) {
  return x*x;
}

inline float clampf(float min, float x, float max) {
  if (x < min) {
    return min;
  } else if (x > max) {
    return max;
  } else {
    return x;
  }
}

inline int min(int a, int b) {
  return a < b ? a : b;
}

inline int max(int a, int b) {
  return a < b ? b : a;
}
