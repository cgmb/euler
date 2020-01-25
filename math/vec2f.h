#pragma once

typedef struct vec2f {
  float x;
  float y;
} vec2f;

static inline vec2f mulf_v2f(float s, vec2f v) {
  v.x *= s;
  v.y *= s;
  return v;
}

static inline vec2f add_v2f(vec2f a, vec2f b) {
  a.x += b.x;
  a.y += b.y;
  return a;
}

static inline vec2f sub_v2f(vec2f a, vec2f b) {
  a.x -= b.x;
  a.y -= b.y;
  return a;
}

static inline vec2f make_v2f(float x, float y) {
  vec2f v = {x, y};
  return v;
}
