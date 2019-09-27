#pragma once

#include <assert.h>
#include <math.h>
#include <stddef.h>

struct vec2f {
  vec2f() = default;
  vec2f(float x, float y);

  const float& operator[](size_t i) const;
  float& operator[](size_t i);

  vec2f& operator*=(float rhs);
  vec2f& operator+=(const vec2f& rhs);
  vec2f& operator-=(const vec2f& rhs);
  vec2f& operator/=(float rhs);

public:
  float x;
  float y;
};

/* Non-Member Functions */

inline float dot(const vec2f& lhs, const vec2f& rhs) {
  return lhs.x*rhs.x + lhs.y*rhs.y;
}

inline float magnitude_sq(const vec2f& v) {
  return v.x*v.x + v.y*v.y;
}

inline float magnitude(const vec2f& v) {
  return sqrtf(magnitude_sq(v));
}

inline vec2f operator/(vec2f lhs, float rhs) {
  return lhs /= rhs;
}

inline vec2f operator*(vec2f lhs, float rhs) {
  return lhs *= rhs;
}

inline vec2f operator*(float lhs, vec2f rhs) {
  return rhs *= lhs;
}

inline vec2f operator+(vec2f lhs, const vec2f& rhs) {
  return lhs += rhs;
}

inline vec2f operator-(vec2f lhs, const vec2f& rhs) {
  return lhs -= rhs;
}

inline vec2f operator-(vec2f rhs) {
  rhs.x = -rhs.x;
  rhs.y = -rhs.y;
  return rhs;
}

/* Member Functions */

inline vec2f::vec2f(float x, float y)
  : x(x)
  , y(y)
{}

inline const float& vec2f::operator[](size_t i) const {
  switch (i) {
    case 0: return x;
    case 1: return y;
    default: assert(i < 2);
  }
  __builtin_unreachable();
}

inline float& vec2f::operator[](size_t i) {
  switch (i) {
    case 0: return x;
    case 1: return y;
    default: assert(i < 2);
  }
  __builtin_unreachable();
}

inline vec2f& vec2f::operator*=(float rhs) {
  vec2f& lhs = *this;
  lhs.x *= rhs;
  lhs.y *= rhs;
  return lhs;
}

inline vec2f& vec2f::operator+=(const vec2f& rhs) {
  vec2f& lhs = *this;
  lhs.x += rhs.x;
  lhs.y += rhs.y;
  return lhs;
}

inline vec2f& vec2f::operator-=(const vec2f& rhs) {
  vec2f& lhs = *this;
  lhs.x -= rhs.x;
  lhs.y -= rhs.y;
  return lhs;
}

inline vec2f& vec2f::operator/=(float rhs) {
  vec2f& lhs = *this;
  lhs.x /= rhs;
  lhs.y /= rhs;
  return lhs;
}
