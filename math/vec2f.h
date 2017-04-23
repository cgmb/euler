#ifndef VEC2F_H
#define VEC2F_H

#include <array>
#include <cmath>

struct vec2f {
  vec2f() = default;
  vec2f(float x, float y);

  float x() const;
  float y() const;

  const float* data() const;
  float* data();

  const float& operator[](size_t i) const;
  float& operator[](size_t i);

  vec2f& operator*=(float rhs);
  vec2f& operator+=(const vec2f& rhs);
  vec2f& operator-=(const vec2f& rhs);
  vec2f& operator/=(float rhs);

private:
  std::array<float, 2> data_;
};

inline vec2f::vec2f(float x, float y) {
    (*this)[0] = x;
    (*this)[1] = y;
}

inline float vec2f::x() const {
  return (*this)[0];
}

inline float vec2f::y() const {
  return (*this)[1];
}

inline const float* vec2f::data() const {
  return data_.data();
}

inline float* vec2f::data() {
  return data_.data();
}

inline const float& vec2f::operator[](size_t i) const {
  return data_[i];
}

inline float& vec2f::operator[](size_t i) {
  return data_[i];
}

inline vec2f& vec2f::operator*=(float rhs) {
  vec2f& lhs = *this;
  lhs[0] *= rhs;
  lhs[1] *= rhs;
  return lhs;
}

inline vec2f operator*(vec2f lhs, float rhs) {
  return lhs *= rhs;
}

inline vec2f operator*(float lhs, vec2f rhs) {
  return rhs *= lhs;
}

inline vec2f& vec2f::operator+=(const vec2f& rhs) {
  vec2f& lhs = *this;
  lhs[0] += rhs[0];
  lhs[1] += rhs[1];
  return lhs;
}

inline vec2f operator+(vec2f lhs, const vec2f& rhs) {
  return lhs += rhs;
}

inline vec2f& vec2f::operator-=(const vec2f& rhs) {
  vec2f& lhs = *this;
  lhs[0] -= rhs[0];
  lhs[1] -= rhs[1];
  return lhs;
}

inline vec2f operator-(vec2f lhs, const vec2f& rhs) {
  return lhs -= rhs;
}

inline vec2f operator-(vec2f rhs) {
  rhs[0] = -rhs[0];
  rhs[1] = -rhs[1];
  return rhs;
}

inline vec2f& vec2f::operator/=(float rhs) {
  vec2f& lhs = *this;
  lhs[0] /= rhs;
  lhs[1] /= rhs;
  return lhs;
}

inline vec2f operator/(vec2f lhs, float rhs) {
  lhs /= rhs;
  return lhs;
}

inline float dot(const vec2f& lhs, const vec2f& rhs) {
  return lhs[0]*rhs[0] + lhs[1]*rhs[1];
}

inline float magnitude_sq(const vec2f& v) {
  return v[0]*v[0] + v[1]*v[1];
}

inline float magnitude(const vec2f& v) {
  return std::sqrt(magnitude_sq(v));
}

#endif
