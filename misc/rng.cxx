#include "misc/rng.h"

// State must point to a non-zero value.
// Based on "An Experimental Exploration of Marsaglia's xorshift Generators, Scrambled" (2016)
static uint64_t xorshift64star(uint64_t* state) {
  uint64_t x = *state;
  x ^= x >> 12; // a
  x ^= x << 25; // b
  x ^= x >> 27; // c
  *state = x;
  return x * 0x2545F4914F6CDD1D;
}

// State must point to a non-zero value.
//
// The high bits are of better quality than the low bits,
// so this is a better PRNG than xorshift64star.
uint32_t xorshift64_32star(uint64_t* state) {
  return xorshift64star(state) >> 32;
}
