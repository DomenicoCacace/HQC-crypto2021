#include <stdint.h>
#include "../common/parameters.h"
#include "../common/vector.h"
#include <string.h>

#pragma once

#define MASKS (MASK_LVL)


#if MASKS == 1
    typedef struct shares_t {
        uint64_t s0[VEC_N_SIZE_64];
    } shares_t;
#elif MASKS == 2
    typedef struct {
        uint64_t s0[VEC_N_SIZE_64];
        uint64_t s1[VEC_N_SIZE_64];
    } shares_t;
#elif MASKS == 3
typedef struct shares_t {
        uint64_t s0[VEC_N_SIZE_64];
        uint64_t s1[VEC_N_SIZE_64];
        uint64_t s2[VEC_N_SIZE_64];
    } shares_t;
#elif MASKS == 4
    typedef struct shares_t {
        uint64_t s0[VEC_N_SIZE_64];
        uint64_t s1[VEC_N_SIZE_64];
        uint64_t s2[VEC_N_SIZE_64];
        uint64_t s3[VEC_N_SIZE_64];
    } shares_t;
#else
#error TOO_MANY_SHARES
#endif

void shares_resize(shares_t *shares, const uint64_t *in);

void shares_add(shares_t *o, shares_t *a, shares_t *b);

static inline void shares_init(shares_t *x) {
    memset(x, 0x00, VEC_N_SIZE_BYTES * MASKS);
}
static inline void shares_reduce(uint64_t *o, shares_t *shares) {
#if MASKS == 1
    memcpy(o, shares->s0, VEC_N_SIZE_BYTES);
#elif MASKS == 2
    for(int i = 0; i < VEC_N_SIZE_64; i++)
        o[i] = shares->s0[i] ^ shares->s1[i];
#elif MASKS == 3
    for(int i = 0; i < VEC_N_SIZE_64; i++)
        o[i] = shares->s0[i] ^ shares->s1[i] ^ shares->s2[i];
#elif MASKS == 4
    for(int i = 0; i < VEC_N_SIZE_64; i++)
        o[i] = shares->s0[i] ^ shares->s1[i] ^ shares->s2[i] ^ shares->s3[i];
#endif
}
