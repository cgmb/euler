#pragma once

#include <stddef.h>

#include "math/vec2f.h"
#include "math/vec2zu.h"
#include "misc/mixed_fraction2f.h"

struct SplitIndex2f {
  explicit SplitIndex2f(const vec2f& value);
  explicit SplitIndex2f(const MixedFraction2f& value);

  vec2zu whole;
  vec2f fractional;
};

inline SplitIndex2f::SplitIndex2f(const vec2f& value)
  : SplitIndex2f(MixedFraction2f(value))
{}

inline SplitIndex2f::SplitIndex2f(const MixedFraction2f& value)
  : whole(size_t(value.whole.x), size_t(value.whole.y))
  , fractional(value.fractional)
{}
