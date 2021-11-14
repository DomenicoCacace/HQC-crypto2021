#include "shares.h"
#include <string.h>

void shares_resize(shares_t *shares, const uint64_t *in) {
    uint64_t mask[VEC_N_SIZE_64];
    shake_prng((uint8_t *)mask, VEC_N_SIZE_BYTES);
#if MASKS == 1
    memcpy(shares->s0, in, VEC_N1N2_SIZE_BYTES);
#elif MASKS == 2
    memcpy(shares->s0, in, VEC_N1N2_SIZE_BYTES/2);
    memcpy(shares->s1+VEC_N1N2_SIZE_64/2, in+VEC_N1N2_SIZE_64/2, VEC_N_SIZE_BYTES-VEC_N1N2_SIZE_BYTES/2);
#elif MASKS == 3
    memcpy(shares->s0, in, VEC_N1N2_SIZE_BYTES/3);
    memcpy(shares->s1+(VEC_N1N2_SIZE_64/3)*1, in+(VEC_N1N2_SIZE_64/3)*1, VEC_N1N2_SIZE_BYTES/3);
    memcpy(shares->s2+(VEC_N1N2_SIZE_64/3)*2, in+(VEC_N1N2_SIZE_64/3)*2, VEC_N1N2_SIZE_BYTES - (VEC_N1N2_SIZE_BYTES/3)*2);
#elif MASKS == 4
    memcpy(shares->s0, in, VEC_N1N2_SIZE_BYTES/4);
    memcpy(shares->s1+VEC_N1N2_SIZE_64/4*1, in+VEC_N1N2_SIZE_64/4*1, VEC_N1N2_SIZE_BYTES/4);
    memcpy(shares->s2+VEC_N1N2_SIZE_64/4*2, in+VEC_N1N2_SIZE_64/4*2, VEC_N1N2_SIZE_BYTES/4);
    memcpy(shares->s3+VEC_N1N2_SIZE_64/4*3, in+VEC_N1N2_SIZE_64/4*3, VEC_N1N2_SIZE_BYTES - VEC_N1N2_SIZE_BYTES/4*3);
#endif
}

void shares_add(shares_t *o, shares_t *a, shares_t *b) {
    uint64_t mask[VEC_N_SIZE_64];
    shake_prng((uint8_t *)mask, VEC_N_SIZE_BYTES);
#if MASKS == 1
    for(int i = 0; i < VEC_N_SIZE_64; i++)
        o->s0[i] = a->s0[i] ^ b->s0[i];
#elif MASKS == 2
    for(int i = 0; i < VEC_N_SIZE_64; i++) {
        o->s0[i] = a->s0[i] ^ (b->s0[i] ^ mask[i]);
        o->s1[i] = (a->s1[i] ^ mask[i]) ^ b->s1[i];
    }
#elif MASKS == 3
    for(int i = 0; i < VEC_N_SIZE_64/2; i++) {
        o->s0[i] = a->s0[i] ^ (b->s0[i] ^ mask[i]);
        o->s1[i] = a->s1[i] ^ b->s1[i];
        o->s2[i] = a->s2[i] ^ (b->s2[i] ^ mask[i]);
    }
    for(int i = VEC_N_SIZE_64/2; i < VEC_N_SIZE_64; i++) {
        o->s0[i] = a->s0[i] ^ (b->s0[i] ^ mask[i]);
        o->s1[i] = a->s1[i] ^ (b->s1[i] ^ mask[i]);
        o->s2[i] = a->s2[i] ^ b->s2[i];
    }
#elif MASKS == 4
    for(int i = 0; i < VEC_N_SIZE_64/2; i++) {
        o->s0[i] = a->s0[i] ^ (b->s0[i] ^ mask[i]);
        o->s1[i] = a->s1[i] ^ b->s1[i];
        o->s2[i] = a->s2[i] ^ (b->s2[i] ^ mask[i]);
        o->s3[i] = a->s3[i] ^ b->s3[i];
    }
    for(int i = VEC_N_SIZE_64/2; i < VEC_N_SIZE_64; i++) {
        o->s0[i] = a->s0[i] ^ b->s0[i];
        o->s1[i] = (a->s1[i] ^ mask[i]) ^ b->s1[i];
        o->s2[i] = a->s2[i] ^ b->s2[i];
        o->s3[i] = (a->s3[i] ^ mask[i]) ^ b->s3[i];
    }
#endif
}
