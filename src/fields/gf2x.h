#ifndef GF2X_H
#define GF2X_H

/**
 * @file gf2x.h
 * @brief Header file for gf2x.c
 */

#include <stdint.h>

#include "../lib/shake_prng.h"

void vect_mul(uint64_t *o, const uint32_t *v1, const uint64_t *v2, const uint16_t weight);
void safe_mul(uint64_t *o, uint64_t *mask, uint32_t *a1, const uint64_t *a2, const uint16_t weight);


#endif
