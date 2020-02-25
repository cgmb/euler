#pragma once

#include <math.h>
#include "misc/util.h"

static inline int float_to_byte_color(float x) {
  const float end = nextafterf(256.f, 0.f);
  return (int)clampf(0.f, end*x, end);
}

static inline float linear_to_sRGB(float x) {
  return powf(x, 1/2.2f); // approximation
}

// https://en.wikipedia.org/wiki/HSL_and_HSV
static inline float hsv_basis(float t) {
  // periodic function, repeats after t=[0,6]
  // returns value in range of [0,1]
  t -= 6.f*floorf(1.f/6*t);
  if (t < 0.f) {
    t += 6.f;
  }

  if (t < 1.f) {
    return t;
  } else if (t < 3.f) {
    return 1.f;
  } else if (t < 4.f) {
    return 4.f - t;
  } else {
    return 0.f;
  }
}
