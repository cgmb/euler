#pragma once

// [u]int8_t is often a typedef for [unsigned] char.
// Unfortunately, that typically results in treating them as exceptions to
// the strict aliasing rule, inhibiting optimizations anywhere they could
// potentially be modified.
//
// See GCC Bug 66110:
// https://gcc.gnu.org/bugzilla/show_bug.cgi?id=66110

enum u8 : uint8_t {};

inline u8& operator++(u8& x) {
  x = u8(x + 1);
  return x;
}

inline u8 operator++(u8& x, int) {
  u8 result = x;
  x = u8(x + 1);
  return result;
}


enum i8 : int8_t {};

inline i8& operator++(i8& x) {
  x = i8(x + 1);
  return x;
}

inline i8 operator++(i8& x, int) {
  i8 result = x;
  x = i8(x + 1);
  return result;
}
