#pragma once

#include <assert.h>
#include <stddef.h>

#include "math/vec2zu.h"

class IndexCalculator {
public:
  IndexCalculator(size_t x, size_t y);

  size_t width() const;
  size_t height() const;
  size_t size() const;
  vec2zu dimensions() const;

  vec2zu to2d(size_t i) const;
  size_t to1d(vec2zu i) const;

  size_t up(size_t i) const;
  size_t down(size_t i) const;
  size_t right(size_t i) const;
  size_t left(size_t i) const;

  static bool check_overflow(size_t x, size_t y);

private:
  size_t x_;
  size_t y_;
};

inline IndexCalculator::IndexCalculator(size_t x, size_t y)
  : x_(x)
  , y_(y)
{
  assert(check_overflow(x,y));
}

inline size_t IndexCalculator::width() const {
  return x_;
}

inline size_t IndexCalculator::height() const {
  return y_;
}

inline size_t IndexCalculator::size() const {
  return x_ * y_;
}

inline vec2zu IndexCalculator::dimensions() const {
  return vec2zu(x_, y_);
}

inline vec2zu IndexCalculator::to2d(size_t i) const {
  return vec2zu(i % x_, i / x_);
}

inline size_t IndexCalculator::to1d(vec2zu i) const {
  return i.y*x_ + i.x;
}

inline size_t IndexCalculator::up(size_t i) const {
  assert(i + x_ > x_); // overflow
  assert(i + x_ < size());
  return i + x_;
}

inline size_t IndexCalculator::down(size_t i) const {
  assert(i >= x_);
  return i - x_;
}

inline size_t IndexCalculator::right(size_t i) const {
  assert((i % x_) + 1 < x_);
  return i + 1;
}

inline size_t IndexCalculator::left(size_t i) const {
  assert((i % x_) >= 1);
  return i - 1;
}

// Returns false if x*y would overflow.
inline bool IndexCalculator::check_overflow(size_t x, size_t y) {
  return y == 0 || x*y/y == x;
}
