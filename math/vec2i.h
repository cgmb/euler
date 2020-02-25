#pragma once

typedef struct vec2i {
  int x;
  int y;
} vec2i;

static inline vec2i v2i(float x, float y) {
  return (vec2i){x, y};
}
