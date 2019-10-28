#pragma once

#include "math/vec2f.h"

struct MixedFraction2f {
  explicit MixedFraction2f(const vec2f& value);

  vec2f whole;
  vec2f fractional;
};

inline MixedFraction2f::MixedFraction2f(const vec2f& value) {
  fractional.x = modff(value.x, &whole.x);
  fractional.y = modff(value.y, &whole.y);
}
