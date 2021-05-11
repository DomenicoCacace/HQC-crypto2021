#ifndef GF2X_H
#define GF2X_H

/**
 * @file gf2x.h
 * @brief Header file for gf2x.c
 */

#include "rng.h"
#include <stdint.h>

void vect_mul(uint64_t *o, const uint32_t *v1, const uint64_t *v2, const uint16_t weight, AES_XOF_struct *ctx);

#endif
