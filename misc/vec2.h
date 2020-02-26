#pragma once

#include <math.h>

// float
typedef struct vec2f {
  float x;
  float y;
} vec2f;

static inline vec2f v2f_mulf(float s, vec2f v) {
  v.x *= s;
  v.y *= s;
  return v;
}

static inline vec2f v2f_add(vec2f a, vec2f b) {
  a.x += b.x;
  a.y += b.y;
  return a;
}

static inline vec2f v2f_sub(vec2f a, vec2f b) {
  a.x -= b.x;
  a.y -= b.y;
  return a;
}

static inline vec2f v2f(float x, float y) {
  return (vec2f){x, y};
}

static inline vec2f modf2f(vec2f arg, vec2f* iptr) {
  vec2f result;
  result.x = modff(arg.x, &iptr->x);
  result.y = modff(arg.y, &iptr->y);
  return result;
}

// integer
typedef struct vec2i {
  int x;
  int y;
} vec2i;

static inline vec2i v2i(float x, float y) {
  return (vec2i){x, y};
}

static inline vec2i up(vec2i i) {
  return (vec2i){ i.x, i.y + 1 };
}

static inline vec2i right(vec2i i) {
  return (vec2i){ i.x + 1, i.y };
}

static inline vec2i to_vec2i(vec2f v) {
  return (vec2i){ (int)v.x, (int)v.y };
}
