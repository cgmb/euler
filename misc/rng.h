#ifndef RNG_H
#define RNG_H

#include <stdint.h>

uint32_t xorshift64_32star(uint64_t* state);

#endif
