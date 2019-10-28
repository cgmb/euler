#pragma once

#include <assert.h>

#include "math/vec2zu.h"
#include "misc/dynarray2d.h"
#include "index_calculator.h"

/* Cell Centre Grid

   0   1   2
   ╎   ╎   ╎
   ├───┼───┼╌ 2
   │ x │ x │
   ├───┼───┼╌ 1
   │ x │ x │
   └───┴───┴╌ 0
*/
class PBoolGridAdapter {
public:
  PBoolGridAdapter(const IndexCalculator& ic, const std::vector<bool>& arr);

  bool operator[](const vec2zu& index) const;
  bool at(const vec2zu& index) const;

  size_t width() const;
  size_t height() const;
  vec2zu dimensions() const;

private:
  const IndexCalculator& ic_;
  const std::vector<bool>& arr_;
};

inline PBoolGridAdapter::PBoolGridAdapter(
  const IndexCalculator& ic, const std::vector<bool>& arr)
  : ic_(ic)
  , arr_(arr)
{}

inline bool PBoolGridAdapter::operator[](const vec2zu& index) const {
  return arr_[ic_.to1d(index)];
}

inline bool PBoolGridAdapter::at(const vec2zu& index) const {
  return arr_[ic_.to1d(index)];
}

inline size_t PBoolGridAdapter::width() const {
  return ic_.width();
}

inline size_t PBoolGridAdapter::height() const {
  return ic_.height();
}

inline vec2zu PBoolGridAdapter::dimensions() const {
  return ic_.dimensions();
}

/* Horizontal Grid

   0   1
   ╎   ╎
   ├───┼─╌ 2
   x   x
   ├───┼─╌ 1
   x   x
   └───┴─╌ 0
*/
class UBoolGridAdapter {
public:
  explicit UBoolGridAdapter(const PBoolGridAdapter& pgrid);

  bool operator[](const vec2zu& index) const;
  bool at(const vec2zu& index) const;

  size_t width() const;
  size_t height() const;
  vec2zu dimensions() const;

private:
  const PBoolGridAdapter& pgrid_;
};

inline UBoolGridAdapter::UBoolGridAdapter(const PBoolGridAdapter& pgrid)
  : pgrid_(pgrid)
{}

inline bool UBoolGridAdapter::operator[](const vec2zu& index) const {
  assert(index.x > 0);
  vec2zu left = index - vec2zu(1,0);
  vec2zu right = index;
  return pgrid_[right] || pgrid_[left];
}

inline bool UBoolGridAdapter::at(const vec2zu& index) const {
  vec2zu left = index - vec2zu(1,0);
  vec2zu right = index;
  if (index.x == 0) {
    return pgrid_[right];
  } else if (index.x == pgrid_.width()) {
    return pgrid_[left];
  } else {
    return pgrid_[right] || pgrid_[left];
  }
}

inline size_t UBoolGridAdapter::width() const {
  return pgrid_.width() + 1;
}

inline size_t UBoolGridAdapter::height() const {
  return pgrid_.height();
}

inline vec2zu UBoolGridAdapter::dimensions() const {
  return vec2zu(width(), height());
}

/* Vertical Grid

   0   1   2
   ╎   ╎   ╎
   ├─x─┼─x─┼╌ 1
   │   │   │
   └─x─┴─x─┴╌ 0
*/
class VBoolGridAdapter {
public:
  explicit VBoolGridAdapter(const PBoolGridAdapter& pgrid);

  bool operator[](const vec2zu& index) const;
  bool at(const vec2zu& index) const;

  size_t width() const;
  size_t height() const;
  vec2zu dimensions() const;

private:
  const PBoolGridAdapter& pgrid_;
};

inline VBoolGridAdapter::VBoolGridAdapter(const PBoolGridAdapter& pgrid)
  : pgrid_(pgrid)
{}

inline bool VBoolGridAdapter::operator[](const vec2zu& index) const {
  assert(index.y > 0);
  vec2zu up = index;
  vec2zu down = index - vec2zu(0,1);
  return pgrid_[up] || pgrid_[down];
}

inline bool VBoolGridAdapter::at(const vec2zu& index) const {
  vec2zu up = index;
  vec2zu down = index - vec2zu(0,1);
  if (index.y == 0) {
    return pgrid_[up];
  } else if (index.y == pgrid_.height()) {
    return pgrid_[down];
  } else {
    return pgrid_[up] || pgrid_[down];
  }
}

inline size_t VBoolGridAdapter::width() const {
  return pgrid_.width();
}

inline size_t VBoolGridAdapter::height() const {
  return pgrid_.height() + 1;
}

inline vec2zu VBoolGridAdapter::dimensions() const {
  return vec2zu(width(), height());
}
