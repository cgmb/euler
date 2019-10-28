#pragma once

#include <assert.h>
#include <math.h>
#include <stddef.h>

struct vec3f {
  vec3f() = default;
  vec3f(float x, float y, float z);

  const float& operator[](size_t i) const;
  float& operator[](size_t i);

  vec3f& operator*=(float rhs);
  vec3f& operator+=(const vec3f& rhs);
  vec3f& operator-=(const vec3f& rhs);
  vec3f& operator/=(float rhs);

public:
  float x;
  float y;
  float z;
};

/* Non-Member Functions */

inline float dot(const vec3f& lhs, const vec3f& rhs) {
  return lhs.x*rhs.x + lhs.y*rhs.y + lhs.z+rhs.z;
}

inline float magnitude_sq(const vec3f& v) {
  return v.x*v.x + v.y*v.y + v.z*v.z;
}

inline float magnitude(const vec3f& v) {
  return sqrtf(magnitude_sq(v));
}

inline vec3f operator/(vec3f lhs, float rhs) {
  return lhs /= rhs;
}

inline vec3f operator*(vec3f lhs, float rhs) {
  return lhs *= rhs;
}

inline vec3f operator*(float lhs, vec3f rhs) {
  return rhs *= lhs;
}

inline vec3f operator+(vec3f lhs, const vec3f& rhs) {
  return lhs += rhs;
}

inline vec3f operator-(vec3f lhs, const vec3f& rhs) {
  return lhs -= rhs;
}

inline vec3f operator-(vec3f rhs) {
  rhs.x = -rhs.x;
  rhs.y = -rhs.y;
  rhs.z = -rhs.z;
  return rhs;
}

/* Member Functions */

inline vec3f::vec3f(float x, float y, float z)
  : x(x)
  , y(y)
  , z(z)
{}

inline const float& vec3f::operator[](size_t i) const {
  switch (i) {
    case 0: return x;
    case 1: return y;
    case 2: return z;
    default: assert(i < 3);
  }
  __builtin_unreachable();
}

inline float& vec3f::operator[](size_t i) {
  switch (i) {
    case 0: return x;
    case 1: return y;
    case 2: return z;
    default: assert(i < 3);
  }
  __builtin_unreachable();
}

inline vec3f& vec3f::operator*=(float rhs) {
  vec3f& lhs = *this;
  lhs.x *= rhs;
  lhs.y *= rhs;
  lhs.z *= rhs;
  return lhs;
}

inline vec3f& vec3f::operator+=(const vec3f& rhs) {
  vec3f& lhs = *this;
  lhs.x += rhs.x;
  lhs.y += rhs.y;
  lhs.z += rhs.z;
  return lhs;
}

inline vec3f& vec3f::operator-=(const vec3f& rhs) {
  vec3f& lhs = *this;
  lhs.x -= rhs.x;
  lhs.y -= rhs.y;
  lhs.z -= rhs.z;
  return lhs;
}

inline vec3f& vec3f::operator/=(float rhs) {
  vec3f& lhs = *this;
  lhs.x /= rhs;
  lhs.y /= rhs;
  lhs.z /= rhs;
  return lhs;
}
