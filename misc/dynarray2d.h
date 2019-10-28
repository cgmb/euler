#pragma once

#include <assert.h>

#include <vector>

#include "math/vec2zu.h"

template <class T>
class dynarray2d : private std::vector<T> {
  using base = std::vector<T>;
public:
  using base::value_type;
  using base::size_type;
  using base::difference_type;
  using base::reference;
  using base::const_reference;
  using base::pointer;
  using base::const_pointer;
  using base::iterator;
  using base::const_iterator;
  using base::reverse_iterator;
  using base::const_reverse_iterator;

public:
  dynarray2d();
  dynarray2d(size_t width, size_t height);
  explicit dynarray2d(const vec2zu& dimensions);

  T& operator[](const vec2zu& index);
  const T& operator[](const vec2zu& index) const;
  using base::operator[];

  size_t width() const;
  size_t height() const;
  vec2zu dimensions() const;

  size_t flatten(const vec2zu& index) const;

  void resize(size_t x, size_t y);
  void resize(const vec2zu& dimensions);

  using base::data;
  using base::size;
  using base::clear;

  using base::begin;
  using base::end;
  using base::cbegin;
  using base::cend;
  using base::rbegin;
  using base::rend;
  using base::crbegin;
  using base::crend;

  void fill(const T& value);

private:
  bool is_valid(const vec2zu& index) const;

private:
  size_t width_;
  size_t height_;
};

template <class T>
dynarray2d<T>::dynarray2d()
  : base()
  , width_(0)
  , height_(0)
{}

template <class T>
dynarray2d<T>::dynarray2d(size_t width, size_t height)
  : base(width * height)
  , width_(width)
  , height_(height)
{
  assert(height == 0 || width*height/height == width); // overflow check
}

template <class T>
dynarray2d<T>::dynarray2d(const vec2zu& dimensions)
  : dynarray2d(dimensions.x, dimensions.y)
{}

template <class T>
T& dynarray2d<T>::operator[](const vec2zu& index) {
  assert(is_valid(index));
  return operator[](flatten(index));
}

template <class T>
const T& dynarray2d<T>::operator[](const vec2zu& index) const {
  assert(is_valid(index));
  return operator[](flatten(index));
}

template <class T>
size_t dynarray2d<T>::width() const {
  return width_;
}

template <class T>
size_t dynarray2d<T>::height() const {
  return height_;
}

template <class T>
vec2zu dynarray2d<T>::dimensions() const {
  return vec2zu(width_, height_);
}

template <class T>
size_t dynarray2d<T>::flatten(const vec2zu& index) const {
  return index.y*width_ + index.x;
}

template <class T>
bool dynarray2d<T>::is_valid(const vec2zu& index) const {
  return index.x < width_
      && index.y < height_;
}

template <class T>
void dynarray2d<T>::resize(size_t x, size_t y) {
  width_ = x;
  height_ = y;
  base::resize(x*y);
}

template <class T>
void dynarray2d<T>::resize(const vec2zu& dimensions) {
  resize(dimensions.x, dimensions.y);
}

template <class T>
void dynarray2d<T>::fill(const T& value) {
  const size_t n = size();
  for (size_t i = 0; i < n; ++i) {
    (*this)[i] = value;
  }
}
