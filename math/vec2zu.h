#pragma once

#include <assert.h>
#include <math.h>
#include <stddef.h>

struct vec2zu {
  vec2zu() = default;
  vec2zu(size_t x, size_t y);

  const size_t& operator[](size_t i) const;
  size_t& operator[](size_t i);

  vec2zu& operator+=(const vec2zu& rhs);
  vec2zu& operator-=(const vec2zu& rhs);

public:
  size_t x;
  size_t y;
};

/* Non-Member Functions */

inline vec2zu operator+(vec2zu lhs, const vec2zu& rhs) {
  return lhs += rhs;
}

inline vec2zu operator-(vec2zu lhs, const vec2zu& rhs) {
  return lhs -= rhs;
}

inline bool operator==(const vec2zu& lhs, const vec2zu& rhs) {
  return lhs.x == rhs.x
    && lhs.y == rhs.y;
}

inline bool operator!=(const vec2zu& lhs, const vec2zu& rhs) {
  return !(lhs == rhs);
}

/* Member Functions */

inline vec2zu::vec2zu(size_t x, size_t y)
  : x(x)
  , y(y)
{}

inline const size_t& vec2zu::operator[](size_t i) const {
  switch (i) {
    case 0: return x;
    case 1: return y;
    default: assert(i < 2);
  }
  __builtin_unreachable();
}

inline size_t& vec2zu::operator[](size_t i) {
  switch (i) {
    case 0: return x;
    case 1: return y;
    default: assert(i < 2);
  }
  __builtin_unreachable();
}

inline vec2zu& vec2zu::operator+=(const vec2zu& rhs) {
  vec2zu& lhs = *this;
  lhs.x += rhs.x;
  lhs.y += rhs.y;
  return lhs;
}

inline vec2zu& vec2zu::operator-=(const vec2zu& rhs) {
  vec2zu& lhs = *this;
  lhs.x -= rhs.x;
  lhs.y -= rhs.y;
  return lhs;
}
