#pragma once

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
