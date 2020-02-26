#pragma once

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
