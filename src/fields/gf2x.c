/**
 * @file gf2x.c
 * @brief Implementation of multiplication of two polynomials
 */

#include <stdint.h>
#include <string.h>

#include "../common/parameters.h"
#include "../common/vector.h"
#include "gf2x.h"

#define TABLE 16
#define WORD 64

static void reduce(uint64_t *o, const uint64_t *a);
static void fast_convolution_mult(uint64_t *o, const uint32_t *a1, const uint64_t *a2, uint16_t weight, uint16_t size);


/**
 * @brief Compute o(x) = a(x) mod \f$ X^n - 1\f$
 *
 * This function computes the modular reduction of the polynomial a(x)
 *
 * @param[in] a Pointer to the polynomial a(x)
 * @param[out] o Pointer to the result
 */
static void reduce(uint64_t *o, const uint64_t *a) {
    size_t i;
    uint64_t r;
    uint64_t carry;

    for (i = 0; i < VEC_N_SIZE_64; i++) {
        r = a[i + VEC_N_SIZE_64 - 1] >> (PARAM_N & 0x3F);
        carry = (uint64_t ) (a[i + VEC_N_SIZE_64 ] << (64 - (PARAM_N & 0x3F)));
        o[i] = a[i] ^ r ^ carry;
    }
    o[VEC_N_SIZE_64 - 1] &= RED_MASK;
}


/**
 * @brief computes product of the polynomial a1(x) with the sparse polynomial a2; the dense polynomial
 * has always length VEC_N_SIZE_64/2
 *
 *  o(x) = a1(x)a2(x)
 *
 * @param[out] o Pointer to the result
 * @param[in] a1 Pointer to the sparse polynomial a2 (list of degrees of the monomials which appear in a2)
 * @param[in] a2 Pointer to the polynomial a1(x)
 * @param[in] weight Hamming wifht of the sparse polynomial a2
 */
static void fast_convolution_mult(uint64_t *o, const uint32_t *a1, const uint64_t *a2, const uint16_t weight, const uint16_t size){
    uint64_t carry;
    uint64_t tmp;
    uint64_t table[TABLE * (size + 1)];


    memcpy(table, a2, size*sizeof(uint64_t));
    table[size] = 0x0UL;

    for (size_t i = 1; i < TABLE; i++) {
        carry = 0x0UL;

        for (size_t j = 0; j < size; j++) {
            table[i*(size+1)+j] = (a2[j] << i) ^ carry;
            carry = (a2[j] >> ((WORD - i)));
        }
        table[i*(size+1)+size] = carry;
    }

    for (size_t i = 0; i < weight; i++) {
        uint16_t *res_16 = (uint16_t *) o+(a1[i] >> 4);

        for (size_t j = 0; j < size + 1; j++) {
            tmp = (uint64_t) res_16[0] | ((uint64_t) (res_16[1])) << 16 |
                           (uint64_t) (res_16[2]) << 32 | ((uint64_t) (res_16[3])) << 48;
            tmp ^= table[((a1[i] & 0xf) * (size + 1))+j];
            memcpy(res_16, &tmp, 8);
            res_16 += 4;
        }
    }
}


/**
 * @brief Multiply two polynomials modulo \f$ X^n - 1\f$.
 *
 * This functions multiplies a sparse polynomial <b>a1</b> (of Hamming weight equal to <b>weight</b>)
 * and a dense polynomial <b>a2</b>. The multiplication is done modulo \f$ X^n - 1\f$.
 *
 * @param[out] o Pointer to the result
 * @param[in] a1 Pointer to the sparse polynomial
 * @param[in] a2 Pointer to the dense polynomial
 * @param[in] weight Integer that is the weigt of the sparse polynomial
 * @param[in] ctx Pointer to the randomness context
 */
void vect_mul(uint64_t *o, const uint32_t *a1, const uint64_t *a2, uint16_t weight) {
    uint64_t tmp[(VEC_N_SIZE_64 << 1) + 1] = {0};

    fast_convolution_mult(tmp, a1, a2, weight, VEC_N_SIZE_64);
    reduce(o, tmp);
}


/**
 * @brief Sets the given vector null
 *
 * @param vec the raw_temp vector to nullify
 */
static inline
void reset(uint64_t *vec) {
    memset(vec, 0x00, ((VEC_N_SIZE_64<<1)+1)*8);
}


/**
 * @brief Multiply two polynomials modulo \f$ X^n - 1\f$, with masking
 *
 * This functions multiplies a sparse polynomial <b>a1</b> (of Hamming weight equal to <b>weight</b>)
 * and a dense polynomial <b>a2</b>; he multiplication is done modulo \f$ X^n - 1\f$. This function also
 * implements masking to avoid information leakage
 *
 * @param[out] o Pointer to the result
 * @param[in] a1 Pointer to the sparse polynomial
 * @param[in] a2 Pointer to the dense polynomial
 * @param[in] weight Integer that is the weight of the sparse polynomial
 */
void safe_mul(shares_t *o, const uint32_t *a1, const uint64_t *a2, uint16_t weight) {
#ifdef VERBOSE
    printf("\nsparse_in: ");
    for(int i=0;i<PARAM_OMEGA;i++) printf("%x ", a1[i]);
#endif
    uint64_t temp1[VEC_N_SIZE_64] = {0};
    uint64_t temp2[VEC_N_SIZE_64] = {0};
    uint64_t raw_temp[(VEC_N_SIZE_64 << 1) + 1];

    uint64_t s[VEC_N_SIZE_64] = {0};
    uint64_t s1[VEC_N_SIZE_64] = {0};

    seedexpander_state mask_seedexpander;
    uint8_t seed[SEED_BYTES];

#if MASKS == 1
    memset(raw_temp, 0x00, (VEC_N_SIZE_64*2+1)*8);
    fast_convolution_mult(raw_temp+(0*(VEC_N_SIZE_64/1)),
                          a1+(0*(weight/1)), a2+(0*(VEC_N_SIZE_64/1)),
                          weight - (weight/1)*0, VEC_N_SIZE_64 - (VEC_N_SIZE_64/1)*0);
    reduce(o->s0, raw_temp);
#elif MASKS == 2
    memset(raw_temp, 0x00, (VEC_N_SIZE_64*2+1)*8);
    fast_convolution_mult(raw_temp+(0*(VEC_N_SIZE_64/2)),
                          a1+(0*(weight/2)), a2+(0*(VEC_N_SIZE_64/2)),
                          weight/2, VEC_N_SIZE_64/2);
    reduce(o->s0, raw_temp);
    memset(raw_temp, 0x00, (VEC_N_SIZE_64*2+1)*8);
    fast_convolution_mult(raw_temp+(1*(VEC_N_SIZE_64/2)),
                          a1+(1*(weight/2)), a2+(1*(VEC_N_SIZE_64/2)),
                          weight - (weight/2)*1, VEC_N_SIZE_64 - (VEC_N_SIZE_64/2)*1);
    reduce(o->s1, raw_temp);
#elif MASKS == 3
    memset(raw_temp, 0x00, (VEC_N_SIZE_64*2+1)*8);
    fast_convolution_mult(raw_temp+(0*(VEC_N_SIZE_64/3)),
                          a1+(0*(weight/3)), a2+(0*(VEC_N_SIZE_64/3)),
                          weight/3, VEC_N_SIZE_64/3);
    reduce(o->s0, raw_temp);
    memset(raw_temp, 0x00, (VEC_N_SIZE_64*2+1)*8);
    fast_convolution_mult(raw_temp+(1*(VEC_N_SIZE_64/3)),
                          a1+(1*(weight/3)), a2+(1*(VEC_N_SIZE_64/3)),
                          weight/3, VEC_N_SIZE_64/3);
    reduce(o->s1, raw_temp);
    memset(raw_temp, 0x00, (VEC_N_SIZE_64*2+1)*8);
    fast_convolution_mult(raw_temp+(2*(VEC_N_SIZE_64/3)),
                          a1+(2*(weight/3)), a2+(2*(VEC_N_SIZE_64/3)),
                          weight - (weight/3)*2, VEC_N_SIZE_64 - (VEC_N_SIZE_64/3)*2);
    reduce(o->s2, raw_temp);
#elif MASKS == 4
    memset(raw_temp, 0x00, (VEC_N_SIZE_64*2+1)*8);
    fast_convolution_mult(raw_temp+(0*(VEC_N_SIZE_64/4)),
                          a1+(0*(weight/4)), a2+(0*(VEC_N_SIZE_64/4)),
                          weight/4, VEC_N_SIZE_64/4);
    reduce(o->s0, raw_temp);
    memset(raw_temp, 0x00, (VEC_N_SIZE_64*2+1)*8);
    fast_convolution_mult(raw_temp+(1*(VEC_N_SIZE_64/4)),
                          a1+(1*(weight/4)), a2+(1*(VEC_N_SIZE_64/4)),
                          weight/4, VEC_N_SIZE_64/4);
    reduce(o->s1, raw_temp);
    memset(raw_temp, 0x00, (VEC_N_SIZE_64*2+1)*8);
    fast_convolution_mult(raw_temp+(2*(VEC_N_SIZE_64/4)),
                          a1+(2*(weight/4)), a2+(2*(VEC_N_SIZE_64/4)),
                          weight/4, VEC_N_SIZE_64/4);
    reduce(o->s2, raw_temp);
    memset(raw_temp, 0x00, (VEC_N_SIZE_64*2+1)*8);
    fast_convolution_mult(raw_temp+(3*(VEC_N_SIZE_64/4)),
                          a1+(3*(weight/4)), a2+(3*(VEC_N_SIZE_64/4)),
                          weight - (weight/4)*3, VEC_N_SIZE_64 - (VEC_N_SIZE_64/4)*3);
    reduce(o->s3, raw_temp);
#elif MASKS == 5
    memset(raw_temp, 0x00, (VEC_N_SIZE_64*2+1)*8);
    fast_convolution_mult(raw_temp+(0*(VEC_N_SIZE_64/5)),
                          a1+(0*(weight/5)), a2+(0*(VEC_N_SIZE_64/5)),
                          weight/5, VEC_N_SIZE_64/5);
    reduce(o->s0, raw_temp);
    memset(raw_temp, 0x00, (VEC_N_SIZE_64*2+1)*8);
    fast_convolution_mult(raw_temp+(1*(VEC_N_SIZE_64/5)),
                          a1+(1*(weight/5)), a2+(1*(VEC_N_SIZE_64/5)),
                          weight/5, VEC_N_SIZE_64/5);
    reduce(o->s1, raw_temp);
    memset(raw_temp, 0x00, (VEC_N_SIZE_64*2+1)*8);
    fast_convolution_mult(raw_temp+(2*(VEC_N_SIZE_64/5)),
                          a1+(2*(weight/5)), a2+(2*(VEC_N_SIZE_64/5)),
                          weight/5, VEC_N_SIZE_64/5);
    reduce(o->s2, raw_temp);
    memset(raw_temp, 0x00, (VEC_N_SIZE_64*2+1)*8);
    fast_convolution_mult(raw_temp+(3*(VEC_N_SIZE_64/5)),
                          a1+(3*(weight/5)), a2+(3*(VEC_N_SIZE_64/5)),
                          weight/5, VEC_N_SIZE_64/5);
    reduce(o->s3, raw_temp);
    memset(raw_temp, 0x00, (VEC_N_SIZE_64*2+1)*8);
    fast_convolution_mult(raw_temp+(4*(VEC_N_SIZE_64/5)),
                          a1+(4*(weight/5)), a2+(4*(VEC_N_SIZE_64/5)),
                          weight - (weight/5)*4, VEC_N_SIZE_64 - (VEC_N_SIZE_64/5)*4);
    reduce(o->s4, raw_temp);
#elif MASKS == 6
    memset(raw_temp, 0x00, (VEC_N_SIZE_64*2+1)*8);
    fast_convolution_mult(raw_temp+(0*(VEC_N_SIZE_64/6)),
                          a1+(0*(weight/6)), a2+(0*(VEC_N_SIZE_64/6)),
                          weight/6, VEC_N_SIZE_64/6);
    reduce(o->s0, raw_temp);
    memset(raw_temp, 0x00, (VEC_N_SIZE_64*2+1)*8);
    fast_convolution_mult(raw_temp+(1*(VEC_N_SIZE_64/6)),
                          a1+(1*(weight/6)), a2+(1*(VEC_N_SIZE_64/6)),
                          weight/6, VEC_N_SIZE_64/6);
    reduce(o->s1, raw_temp);
    memset(raw_temp, 0x00, (VEC_N_SIZE_64*2+1)*8);
    fast_convolution_mult(raw_temp+(2*(VEC_N_SIZE_64/6)),
                          a1+(2*(weight/6)), a2+(2*(VEC_N_SIZE_64/6)),
                          weight/6, VEC_N_SIZE_64/6);
    reduce(o->s2, raw_temp);
    memset(raw_temp, 0x00, (VEC_N_SIZE_64*2+1)*8);
    fast_convolution_mult(raw_temp+(3*(VEC_N_SIZE_64/6)),
                          a1+(3*(weight/6)), a2+(3*(VEC_N_SIZE_64/6)),
                          weight/6, VEC_N_SIZE_64/6);
    reduce(o->s3, raw_temp);
    memset(raw_temp, 0x00, (VEC_N_SIZE_64*2+1)*8);
    fast_convolution_mult(raw_temp+(4*(VEC_N_SIZE_64/6)),
                          a1+(4*(weight/6)), a2+(4*(VEC_N_SIZE_64/6)),
                          weight/6, VEC_N_SIZE_64/6);
    reduce(o->s4, raw_temp);
    memset(raw_temp, 0x00, (VEC_N_SIZE_64*2+1)*8);
    fast_convolution_mult(raw_temp+(5*(VEC_N_SIZE_64/6)),
                          a1+(5*(weight/6)), a2+(5*(VEC_N_SIZE_64/6)),
                          weight - (weight/6)*5, VEC_N_SIZE_64 - (VEC_N_SIZE_64/6)*5);
    reduce(o->s5, raw_temp);
#endif

// PART 2

#if MASKS == 1
    // nothing, no masking applied
#elif MASKS == 2
    shake_prng(seed, SEED_BYTES);
    seedexpander_init(&mask_seedexpander, seed, SEED_BYTES);
    vect_set_random_fixed_weight(&mask_seedexpander, s, weight);
    memset(raw_temp, 0x00, (VEC_N_SIZE_64*2+1)*8);
    fast_convolution_mult(raw_temp+(1*(VEC_N_SIZE_64/2)),
                          a1+(0*(weight/2)), a2+(1*(VEC_N_SIZE_64/2)),
                          weight/2, VEC_N_SIZE_64 - (VEC_N_SIZE_64/2)*1);
    reduce(temp1, raw_temp);
    vect_add(temp1, temp1, s, VEC_N_SIZE_64);
    memset(raw_temp, 0x00, (VEC_N_SIZE_64*2+1)*8);
    fast_convolution_mult(raw_temp+(0*(VEC_N_SIZE_64/2)),
                          a1+(1*(weight/2)), a2+(0*(VEC_N_SIZE_64/2)),
                          weight - (weight/2)*1, VEC_N_SIZE_64/2);
    reduce(temp2, raw_temp);
    vect_add(s1, temp1, temp2, VEC_N_SIZE_64);
    vect_add(o->s0, o->s0, s, VEC_N_SIZE_64);
    vect_add(o->s1, o->s1, s1, VEC_N_SIZE_64);

#elif MASKS == 3
    shake_prng(seed, SEED_BYTES);
    seedexpander_init(&mask_seedexpander, seed, SEED_BYTES);
    vect_set_random_fixed_weight(&mask_seedexpander, s, weight);
    memset(raw_temp, 0x00, (VEC_N_SIZE_64*2+1)*8);
    fast_convolution_mult(raw_temp+(1*(VEC_N_SIZE_64/3)),
                          a1+(0*(weight/3)), a2+(1*(VEC_N_SIZE_64/3)),
                          weight/3, VEC_N_SIZE_64/3);
    reduce(temp1, raw_temp);
    vect_add(temp1, temp1, s, VEC_N_SIZE_64);
    memset(raw_temp, 0x00, (VEC_N_SIZE_64*2+1)*8);
    fast_convolution_mult(raw_temp+(0*(VEC_N_SIZE_64/3)),
                          a1+(1*(weight/3)), a2+(0*(VEC_N_SIZE_64/3)),
                          weight/3, VEC_N_SIZE_64/3);
    reduce(temp2, raw_temp);
    vect_add(s1, temp1, temp2, VEC_N_SIZE_64);
    vect_add(o->s0, o->s0, s, VEC_N_SIZE_64);
    vect_add(o->s1, o->s1, s1, VEC_N_SIZE_64);

    shake_prng(seed, SEED_BYTES);
    seedexpander_init(&mask_seedexpander, seed, SEED_BYTES);
    vect_set_random_fixed_weight(&mask_seedexpander, s, weight);
    memset(raw_temp, 0x00, (VEC_N_SIZE_64*2+1)*8);
    fast_convolution_mult(raw_temp+(2*(VEC_N_SIZE_64/3)),
                          a1+(0*(weight/3)), a2+(2*(VEC_N_SIZE_64/3)),
                          weight/3, VEC_N_SIZE_64 - (VEC_N_SIZE_64/3)*2);
    reduce(temp1, raw_temp);
    vect_add(temp1, temp1, s, VEC_N_SIZE_64);
    memset(raw_temp, 0x00, (VEC_N_SIZE_64*2+1)*8);
    fast_convolution_mult(raw_temp+(0*(VEC_N_SIZE_64/3)),
                          a1+(2*(weight/3)), a2+(0*(VEC_N_SIZE_64/3)),
                          weight - (weight/3)*2, VEC_N_SIZE_64/3);
    reduce(temp2, raw_temp);
    vect_add(s1, temp1, temp2, VEC_N_SIZE_64);
    vect_add(o->s0, o->s0, s, VEC_N_SIZE_64);
    vect_add(o->s2, o->s2, s1, VEC_N_SIZE_64);

    shake_prng(seed, SEED_BYTES);
    seedexpander_init(&mask_seedexpander, seed, SEED_BYTES);
    vect_set_random_fixed_weight(&mask_seedexpander, s, weight);
    memset(raw_temp, 0x00, (VEC_N_SIZE_64*2+1)*8);
    fast_convolution_mult(raw_temp+(2*(VEC_N_SIZE_64/3)),
                          a1+(1*(weight/3)), a2+(2*(VEC_N_SIZE_64/3)),
                          weight/3, VEC_N_SIZE_64 - (VEC_N_SIZE_64/3)*2);
    reduce(temp1, raw_temp);
    vect_add(temp1, temp1, s, VEC_N_SIZE_64);
    memset(raw_temp, 0x00, (VEC_N_SIZE_64*2+1)*8);
    fast_convolution_mult(raw_temp+(1*(VEC_N_SIZE_64/3)),
                          a1+(2*(weight/3)), a2+(1*(VEC_N_SIZE_64/3)),
                          weight - (weight/3)*2, VEC_N_SIZE_64/3);
    reduce(temp2, raw_temp);
    vect_add(s1, temp1, temp2, VEC_N_SIZE_64);
    vect_add(o->s1, o->s1, s, VEC_N_SIZE_64);
    vect_add(o->s2, o->s2, s1, VEC_N_SIZE_64);
#elif MASKS == 4
    shake_prng(seed, SEED_BYTES);
    seedexpander_init(&mask_seedexpander, seed, SEED_BYTES);
    vect_set_random_fixed_weight(&mask_seedexpander, s, weight);
    memset(raw_temp, 0x00, (VEC_N_SIZE_64*2+1)*8);
    fast_convolution_mult(raw_temp+(1*(VEC_N_SIZE_64/4)),
                          a1+(0*(weight/4)), a2+(1*(VEC_N_SIZE_64/4)),
                          weight/4, VEC_N_SIZE_64/4);
    reduce(temp1, raw_temp);
    vect_add(temp1, temp1, s, VEC_N_SIZE_64);
    memset(raw_temp, 0x00, (VEC_N_SIZE_64*2+1)*8);
    fast_convolution_mult(raw_temp+(0*(VEC_N_SIZE_64/4)),
                          a1+(1*(weight/4)), a2+(0*(VEC_N_SIZE_64/4)),
                          weight/4, VEC_N_SIZE_64/4);
    reduce(temp2, raw_temp);
    vect_add(s1, temp1, temp2, VEC_N_SIZE_64);
    vect_add(o->s0, o->s0, s, VEC_N_SIZE_64);
    vect_add(o->s1, o->s1, s1, VEC_N_SIZE_64);

    shake_prng(seed, SEED_BYTES);
    seedexpander_init(&mask_seedexpander, seed, SEED_BYTES);
    vect_set_random_fixed_weight(&mask_seedexpander, s, weight);
    memset(raw_temp, 0x00, (VEC_N_SIZE_64*2+1)*8);
    fast_convolution_mult(raw_temp+(2*(VEC_N_SIZE_64/4)),
                          a1+(0*(weight/4)), a2+(2*(VEC_N_SIZE_64/4)),
                          weight/4, VEC_N_SIZE_64/4);
    reduce(temp1, raw_temp);
    vect_add(temp1, temp1, s, VEC_N_SIZE_64);
    memset(raw_temp, 0x00, (VEC_N_SIZE_64*2+1)*8);
    fast_convolution_mult(raw_temp+(0*(VEC_N_SIZE_64/4)),
                          a1+(2*(weight/4)), a2+(0*(VEC_N_SIZE_64/4)),
                          weight/4, VEC_N_SIZE_64/4);
    reduce(temp2, raw_temp);
    vect_add(s1, temp1, temp2, VEC_N_SIZE_64);
    vect_add(o->s0, o->s0, s, VEC_N_SIZE_64);
    vect_add(o->s2, o->s2, s1, VEC_N_SIZE_64);

    shake_prng(seed, SEED_BYTES);
    seedexpander_init(&mask_seedexpander, seed, SEED_BYTES);
    vect_set_random_fixed_weight(&mask_seedexpander, s, weight);
    memset(raw_temp, 0x00, (VEC_N_SIZE_64*2+1)*8);
    fast_convolution_mult(raw_temp+(3*(VEC_N_SIZE_64/4)),
                          a1+(0*(weight/4)), a2+(3*(VEC_N_SIZE_64/4)),
                          weight/4, VEC_N_SIZE_64 - (VEC_N_SIZE_64/4)*3);
    reduce(temp1, raw_temp);
    vect_add(temp1, temp1, s, VEC_N_SIZE_64);
    memset(raw_temp, 0x00, (VEC_N_SIZE_64*2+1)*8);
    fast_convolution_mult(raw_temp+(0*(VEC_N_SIZE_64/4)),
                          a1+(3*(weight/4)), a2+(0*(VEC_N_SIZE_64/4)),
                          weight - (weight/4)*3, VEC_N_SIZE_64/4);
    reduce(temp2, raw_temp);
    vect_add(s1, temp1, temp2, VEC_N_SIZE_64);
    vect_add(o->s0, o->s0, s, VEC_N_SIZE_64);
    vect_add(o->s3, o->s3, s1, VEC_N_SIZE_64);

    shake_prng(seed, SEED_BYTES);
    seedexpander_init(&mask_seedexpander, seed, SEED_BYTES);
    vect_set_random_fixed_weight(&mask_seedexpander, s, weight);
    memset(raw_temp, 0x00, (VEC_N_SIZE_64*2+1)*8);
    fast_convolution_mult(raw_temp+(2*(VEC_N_SIZE_64/4)),
                          a1+(1*(weight/4)), a2+(2*(VEC_N_SIZE_64/4)),
                          weight/4, VEC_N_SIZE_64/4);
    reduce(temp1, raw_temp);
    vect_add(temp1, temp1, s, VEC_N_SIZE_64);
    memset(raw_temp, 0x00, (VEC_N_SIZE_64*2+1)*8);
    fast_convolution_mult(raw_temp+(1*(VEC_N_SIZE_64/4)),
                          a1+(2*(weight/4)), a2+(1*(VEC_N_SIZE_64/4)),
                          weight/4, VEC_N_SIZE_64/4);
    reduce(temp2, raw_temp);
    vect_add(s1, temp1, temp2, VEC_N_SIZE_64);
    vect_add(o->s1, o->s1, s, VEC_N_SIZE_64);
    vect_add(o->s2, o->s2, s1, VEC_N_SIZE_64);

    shake_prng(seed, SEED_BYTES);
    seedexpander_init(&mask_seedexpander, seed, SEED_BYTES);
    vect_set_random_fixed_weight(&mask_seedexpander, s, weight);
    memset(raw_temp, 0x00, (VEC_N_SIZE_64*2+1)*8);
    fast_convolution_mult(raw_temp+(3*(VEC_N_SIZE_64/4)),
                          a1+(1*(weight/4)), a2+(3*(VEC_N_SIZE_64/4)),
                          weight/4, VEC_N_SIZE_64 - (VEC_N_SIZE_64/4)*3);
    reduce(temp1, raw_temp);
    vect_add(temp1, temp1, s, VEC_N_SIZE_64);
    memset(raw_temp, 0x00, (VEC_N_SIZE_64*2+1)*8);
    fast_convolution_mult(raw_temp+(1*(VEC_N_SIZE_64/4)),
                          a1+(3*(weight/4)), a2+(1*(VEC_N_SIZE_64/4)),
                          weight - (weight/4)*3, VEC_N_SIZE_64/4);
    reduce(temp2, raw_temp);
    vect_add(s1, temp1, temp2, VEC_N_SIZE_64);
    vect_add(o->s1, o->s1, s, VEC_N_SIZE_64);
    vect_add(o->s3, o->s3, s1, VEC_N_SIZE_64);

    shake_prng(seed, SEED_BYTES);
    seedexpander_init(&mask_seedexpander, seed, SEED_BYTES);
    vect_set_random_fixed_weight(&mask_seedexpander, s, weight);
    memset(raw_temp, 0x00, (VEC_N_SIZE_64*2+1)*8);
    fast_convolution_mult(raw_temp+(3*(VEC_N_SIZE_64/4)),
                          a1+(2*(weight/4)), a2+(3*(VEC_N_SIZE_64/4)),
                          weight/4, VEC_N_SIZE_64 - (VEC_N_SIZE_64/4)*3);
    reduce(temp1, raw_temp);
    vect_add(temp1, temp1, s, VEC_N_SIZE_64);
    memset(raw_temp, 0x00, (VEC_N_SIZE_64*2+1)*8);
    fast_convolution_mult(raw_temp+(2*(VEC_N_SIZE_64/4)),
                          a1+(3*(weight/4)), a2+(2*(VEC_N_SIZE_64/4)),
                          weight - (weight/4)*3, VEC_N_SIZE_64/4);
    reduce(temp2, raw_temp);
    vect_add(s1, temp1, temp2, VEC_N_SIZE_64);
    vect_add(o->s2, o->s2, s, VEC_N_SIZE_64);
    vect_add(o->s3, o->s3, s1, VEC_N_SIZE_64);
#elif MASKS == 5
    shake_prng(seed, SEED_BYTES);
    seedexpander_init(&mask_seedexpander, seed, SEED_BYTES);
    vect_set_random_fixed_weight(&mask_seedexpander, s, weight);
    memset(raw_temp, 0x00, (VEC_N_SIZE_64*2+1)*8);
    fast_convolution_mult(raw_temp+(1*(VEC_N_SIZE_64/5)),
                          a1+(0*(weight/5)), a2+(1*(VEC_N_SIZE_64/5)),
                          weight/5, VEC_N_SIZE_64/5);
    reduce(temp1, raw_temp);
    vect_add(temp1, temp1, s, VEC_N_SIZE_64);
    memset(raw_temp, 0x00, (VEC_N_SIZE_64*2+1)*8);
    fast_convolution_mult(raw_temp+(0*(VEC_N_SIZE_64/5)),
                          a1+(1*(weight/5)), a2+(0*(VEC_N_SIZE_64/5)),
                          weight/5, VEC_N_SIZE_64/5);
    reduce(temp2, raw_temp);
    vect_add(s1, temp1, temp2, VEC_N_SIZE_64);
    vect_add(o->s0, o->s0, s, VEC_N_SIZE_64);
    vect_add(o->s1, o->s1, s1, VEC_N_SIZE_64);

    shake_prng(seed, SEED_BYTES);
    seedexpander_init(&mask_seedexpander, seed, SEED_BYTES);
    vect_set_random_fixed_weight(&mask_seedexpander, s, weight);
    memset(raw_temp, 0x00, (VEC_N_SIZE_64*2+1)*8);
    fast_convolution_mult(raw_temp+(2*(VEC_N_SIZE_64/5)),
                          a1+(0*(weight/5)), a2+(2*(VEC_N_SIZE_64/5)),
                          weight/5, VEC_N_SIZE_64/5);
    reduce(temp1, raw_temp);
    vect_add(temp1, temp1, s, VEC_N_SIZE_64);
    memset(raw_temp, 0x00, (VEC_N_SIZE_64*2+1)*8);
    fast_convolution_mult(raw_temp+(0*(VEC_N_SIZE_64/5)),
                          a1+(2*(weight/5)), a2+(0*(VEC_N_SIZE_64/5)),
                          weight/5, VEC_N_SIZE_64/5);
    reduce(temp2, raw_temp);
    vect_add(s1, temp1, temp2, VEC_N_SIZE_64);
    vect_add(o->s0, o->s0, s, VEC_N_SIZE_64);
    vect_add(o->s2, o->s2, s1, VEC_N_SIZE_64);

    shake_prng(seed, SEED_BYTES);
    seedexpander_init(&mask_seedexpander, seed, SEED_BYTES);
    vect_set_random_fixed_weight(&mask_seedexpander, s, weight);
    memset(raw_temp, 0x00, (VEC_N_SIZE_64*2+1)*8);
    fast_convolution_mult(raw_temp+(3*(VEC_N_SIZE_64/5)),
                          a1+(0*(weight/5)), a2+(3*(VEC_N_SIZE_64/5)),
                          weight/5, VEC_N_SIZE_64/5);
    reduce(temp1, raw_temp);
    vect_add(temp1, temp1, s, VEC_N_SIZE_64);
    memset(raw_temp, 0x00, (VEC_N_SIZE_64*2+1)*8);
    fast_convolution_mult(raw_temp+(0*(VEC_N_SIZE_64/5)),
                          a1+(3*(weight/5)), a2+(0*(VEC_N_SIZE_64/5)),
                          weight/5, VEC_N_SIZE_64/5);
    reduce(temp2, raw_temp);
    vect_add(s1, temp1, temp2, VEC_N_SIZE_64);
    vect_add(o->s0, o->s0, s, VEC_N_SIZE_64);
    vect_add(o->s3, o->s3, s1, VEC_N_SIZE_64);

    shake_prng(seed, SEED_BYTES);
    seedexpander_init(&mask_seedexpander, seed, SEED_BYTES);
    vect_set_random_fixed_weight(&mask_seedexpander, s, weight);
    memset(raw_temp, 0x00, (VEC_N_SIZE_64*2+1)*8);
    fast_convolution_mult(raw_temp+(4*(VEC_N_SIZE_64/5)),
                          a1+(0*(weight/5)), a2+(4*(VEC_N_SIZE_64/5)),
                          weight/5, VEC_N_SIZE_64 - (VEC_N_SIZE_64/5)*4);
    reduce(temp1, raw_temp);
    vect_add(temp1, temp1, s, VEC_N_SIZE_64);
    memset(raw_temp, 0x00, (VEC_N_SIZE_64*2+1)*8);
    fast_convolution_mult(raw_temp+(0*(VEC_N_SIZE_64/5)),
                          a1+(4*(weight/5)), a2+(0*(VEC_N_SIZE_64/5)),
                          weight - (weight/5)*4, VEC_N_SIZE_64/5);
    reduce(temp2, raw_temp);
    vect_add(s1, temp1, temp2, VEC_N_SIZE_64);
    vect_add(o->s0, o->s0, s, VEC_N_SIZE_64);
    vect_add(o->s4, o->s4, s1, VEC_N_SIZE_64);

    shake_prng(seed, SEED_BYTES);
    seedexpander_init(&mask_seedexpander, seed, SEED_BYTES);
    vect_set_random_fixed_weight(&mask_seedexpander, s, weight);
    memset(raw_temp, 0x00, (VEC_N_SIZE_64*2+1)*8);
    fast_convolution_mult(raw_temp+(2*(VEC_N_SIZE_64/5)),
                          a1+(1*(weight/5)), a2+(2*(VEC_N_SIZE_64/5)),
                          weight/5, VEC_N_SIZE_64/5);
    reduce(temp1, raw_temp);
    vect_add(temp1, temp1, s, VEC_N_SIZE_64);
    memset(raw_temp, 0x00, (VEC_N_SIZE_64*2+1)*8);
    fast_convolution_mult(raw_temp+(1*(VEC_N_SIZE_64/5)),
                          a1+(2*(weight/5)), a2+(1*(VEC_N_SIZE_64/5)),
                          weight/5, VEC_N_SIZE_64/5);
    reduce(temp2, raw_temp);
    vect_add(s1, temp1, temp2, VEC_N_SIZE_64);
    vect_add(o->s1, o->s1, s, VEC_N_SIZE_64);
    vect_add(o->s2, o->s2, s1, VEC_N_SIZE_64);

    shake_prng(seed, SEED_BYTES);
    seedexpander_init(&mask_seedexpander, seed, SEED_BYTES);
    vect_set_random_fixed_weight(&mask_seedexpander, s, weight);
    memset(raw_temp, 0x00, (VEC_N_SIZE_64*2+1)*8);
    fast_convolution_mult(raw_temp+(3*(VEC_N_SIZE_64/5)),
                          a1+(1*(weight/5)), a2+(3*(VEC_N_SIZE_64/5)),
                          weight/5, VEC_N_SIZE_64/5);
    reduce(temp1, raw_temp);
    vect_add(temp1, temp1, s, VEC_N_SIZE_64);
    memset(raw_temp, 0x00, (VEC_N_SIZE_64*2+1)*8);
    fast_convolution_mult(raw_temp+(1*(VEC_N_SIZE_64/5)),
                          a1+(3*(weight/5)), a2+(1*(VEC_N_SIZE_64/5)),
                          weight/5, VEC_N_SIZE_64/5);
    reduce(temp2, raw_temp);
    vect_add(s1, temp1, temp2, VEC_N_SIZE_64);
    vect_add(o->s1, o->s1, s, VEC_N_SIZE_64);
    vect_add(o->s3, o->s3, s1, VEC_N_SIZE_64);

    shake_prng(seed, SEED_BYTES);
    seedexpander_init(&mask_seedexpander, seed, SEED_BYTES);
    vect_set_random_fixed_weight(&mask_seedexpander, s, weight);
    memset(raw_temp, 0x00, (VEC_N_SIZE_64*2+1)*8);
    fast_convolution_mult(raw_temp+(4*(VEC_N_SIZE_64/5)),
                          a1+(1*(weight/5)), a2+(4*(VEC_N_SIZE_64/5)),
                          weight/5, VEC_N_SIZE_64 - (VEC_N_SIZE_64/5)*4);
    reduce(temp1, raw_temp);
    vect_add(temp1, temp1, s, VEC_N_SIZE_64);
    memset(raw_temp, 0x00, (VEC_N_SIZE_64*2+1)*8);
    fast_convolution_mult(raw_temp+(1*(VEC_N_SIZE_64/5)),
                          a1+(4*(weight/5)), a2+(1*(VEC_N_SIZE_64/5)),
                          weight - (weight/5)*4, VEC_N_SIZE_64/5);
    reduce(temp2, raw_temp);
    vect_add(s1, temp1, temp2, VEC_N_SIZE_64);
    vect_add(o->s1, o->s1, s, VEC_N_SIZE_64);
    vect_add(o->s4, o->s4, s1, VEC_N_SIZE_64);

    shake_prng(seed, SEED_BYTES);
    seedexpander_init(&mask_seedexpander, seed, SEED_BYTES);
    vect_set_random_fixed_weight(&mask_seedexpander, s, weight);
    memset(raw_temp, 0x00, (VEC_N_SIZE_64*2+1)*8);
    fast_convolution_mult(raw_temp+(3*(VEC_N_SIZE_64/5)),
                          a1+(2*(weight/5)), a2+(3*(VEC_N_SIZE_64/5)),
                          weight/5, VEC_N_SIZE_64/5);
    reduce(temp1, raw_temp);
    vect_add(temp1, temp1, s, VEC_N_SIZE_64);
    memset(raw_temp, 0x00, (VEC_N_SIZE_64*2+1)*8);
    fast_convolution_mult(raw_temp+(2*(VEC_N_SIZE_64/5)),
                          a1+(3*(weight/5)), a2+(2*(VEC_N_SIZE_64/5)),
                          weight/5, VEC_N_SIZE_64/5);
    reduce(temp2, raw_temp);
    vect_add(s1, temp1, temp2, VEC_N_SIZE_64);
    vect_add(o->s2, o->s2, s, VEC_N_SIZE_64);
    vect_add(o->s3, o->s3, s1, VEC_N_SIZE_64);

    shake_prng(seed, SEED_BYTES);
    seedexpander_init(&mask_seedexpander, seed, SEED_BYTES);
    vect_set_random_fixed_weight(&mask_seedexpander, s, weight);
    memset(raw_temp, 0x00, (VEC_N_SIZE_64*2+1)*8);
    fast_convolution_mult(raw_temp+(4*(VEC_N_SIZE_64/5)),
                          a1+(2*(weight/5)), a2+(4*(VEC_N_SIZE_64/5)),
                          weight/5, VEC_N_SIZE_64 - (VEC_N_SIZE_64/5)*4);
    reduce(temp1, raw_temp);
    vect_add(temp1, temp1, s, VEC_N_SIZE_64);
    memset(raw_temp, 0x00, (VEC_N_SIZE_64*2+1)*8);
    fast_convolution_mult(raw_temp+(2*(VEC_N_SIZE_64/5)),
                          a1+(4*(weight/5)), a2+(2*(VEC_N_SIZE_64/5)),
                          weight - (weight/5)*4, VEC_N_SIZE_64/5);
    reduce(temp2, raw_temp);
    vect_add(s1, temp1, temp2, VEC_N_SIZE_64);
    vect_add(o->s2, o->s2, s, VEC_N_SIZE_64);
    vect_add(o->s4, o->s4, s1, VEC_N_SIZE_64);

    shake_prng(seed, SEED_BYTES);
    seedexpander_init(&mask_seedexpander, seed, SEED_BYTES);
    vect_set_random_fixed_weight(&mask_seedexpander, s, weight);
    memset(raw_temp, 0x00, (VEC_N_SIZE_64*2+1)*8);
    fast_convolution_mult(raw_temp+(4*(VEC_N_SIZE_64/5)),
                          a1+(3*(weight/5)), a2+(4*(VEC_N_SIZE_64/5)),
                          weight/5, VEC_N_SIZE_64 - (VEC_N_SIZE_64/5)*4);
    reduce(temp1, raw_temp);
    vect_add(temp1, temp1, s, VEC_N_SIZE_64);
    memset(raw_temp, 0x00, (VEC_N_SIZE_64*2+1)*8);
    fast_convolution_mult(raw_temp+(3*(VEC_N_SIZE_64/5)),
                          a1+(4*(weight/5)), a2+(3*(VEC_N_SIZE_64/5)),
                          weight - (weight/5)*4, VEC_N_SIZE_64/5);
    reduce(temp2, raw_temp);
    vect_add(s1, temp1, temp2, VEC_N_SIZE_64);
    vect_add(o->s3, o->s3, s, VEC_N_SIZE_64);
    vect_add(o->s4, o->s4, s1, VEC_N_SIZE_64);
#elif MASKS == 6
    shake_prng(seed, SEED_BYTES);
    seedexpander_init(&mask_seedexpander, seed, SEED_BYTES);
    vect_set_random_fixed_weight(&mask_seedexpander, s, weight);
    memset(raw_temp, 0x00, (VEC_N_SIZE_64*2+1)*8);
    fast_convolution_mult(raw_temp+(1*(VEC_N_SIZE_64/6)),
                          a1+(0*(weight/6)), a2+(1*(VEC_N_SIZE_64/6)),
                          weight/6, VEC_N_SIZE_64/6);
    reduce(temp1, raw_temp);
    vect_add(temp1, temp1, s, VEC_N_SIZE_64);
    memset(raw_temp, 0x00, (VEC_N_SIZE_64*2+1)*8);
    fast_convolution_mult(raw_temp+(0*(VEC_N_SIZE_64/6)),
                          a1+(1*(weight/6)), a2+(0*(VEC_N_SIZE_64/6)),
                          weight/6, VEC_N_SIZE_64/6);
    reduce(temp2, raw_temp);
    vect_add(s1, temp1, temp2, VEC_N_SIZE_64);
    vect_add(o->s0, o->s0, s, VEC_N_SIZE_64);
    vect_add(o->s1, o->s1, s1, VEC_N_SIZE_64);

    shake_prng(seed, SEED_BYTES);
    seedexpander_init(&mask_seedexpander, seed, SEED_BYTES);
    vect_set_random_fixed_weight(&mask_seedexpander, s, weight);
    memset(raw_temp, 0x00, (VEC_N_SIZE_64*2+1)*8);
    fast_convolution_mult(raw_temp+(2*(VEC_N_SIZE_64/6)),
                          a1+(0*(weight/6)), a2+(2*(VEC_N_SIZE_64/6)),
                          weight/6, VEC_N_SIZE_64/6);
    reduce(temp1, raw_temp);
    vect_add(temp1, temp1, s, VEC_N_SIZE_64);
    memset(raw_temp, 0x00, (VEC_N_SIZE_64*2+1)*8);
    fast_convolution_mult(raw_temp+(0*(VEC_N_SIZE_64/6)),
                          a1+(2*(weight/6)), a2+(0*(VEC_N_SIZE_64/6)),
                          weight/6, VEC_N_SIZE_64/6);
    reduce(temp2, raw_temp);
    vect_add(s1, temp1, temp2, VEC_N_SIZE_64);
    vect_add(o->s0, o->s0, s, VEC_N_SIZE_64);
    vect_add(o->s2, o->s2, s1, VEC_N_SIZE_64);

    shake_prng(seed, SEED_BYTES);
    seedexpander_init(&mask_seedexpander, seed, SEED_BYTES);
    vect_set_random_fixed_weight(&mask_seedexpander, s, weight);
    memset(raw_temp, 0x00, (VEC_N_SIZE_64*2+1)*8);
    fast_convolution_mult(raw_temp+(3*(VEC_N_SIZE_64/6)),
                          a1+(0*(weight/6)), a2+(3*(VEC_N_SIZE_64/6)),
                          weight/6, VEC_N_SIZE_64/6);
    reduce(temp1, raw_temp);
    vect_add(temp1, temp1, s, VEC_N_SIZE_64);
    memset(raw_temp, 0x00, (VEC_N_SIZE_64*2+1)*8);
    fast_convolution_mult(raw_temp+(0*(VEC_N_SIZE_64/6)),
                          a1+(3*(weight/6)), a2+(0*(VEC_N_SIZE_64/6)),
                          weight/6, VEC_N_SIZE_64/6);
    reduce(temp2, raw_temp);
    vect_add(s1, temp1, temp2, VEC_N_SIZE_64);
    vect_add(o->s0, o->s0, s, VEC_N_SIZE_64);
    vect_add(o->s3, o->s3, s1, VEC_N_SIZE_64);

    shake_prng(seed, SEED_BYTES);
    seedexpander_init(&mask_seedexpander, seed, SEED_BYTES);
    vect_set_random_fixed_weight(&mask_seedexpander, s, weight);
    memset(raw_temp, 0x00, (VEC_N_SIZE_64*2+1)*8);
    fast_convolution_mult(raw_temp+(4*(VEC_N_SIZE_64/6)),
                          a1+(0*(weight/6)), a2+(4*(VEC_N_SIZE_64/6)),
                          weight/6, VEC_N_SIZE_64/6);
    reduce(temp1, raw_temp);
    vect_add(temp1, temp1, s, VEC_N_SIZE_64);
    memset(raw_temp, 0x00, (VEC_N_SIZE_64*2+1)*8);
    fast_convolution_mult(raw_temp+(0*(VEC_N_SIZE_64/6)),
                          a1+(4*(weight/6)), a2+(0*(VEC_N_SIZE_64/6)),
                          weight/6, VEC_N_SIZE_64/6);
    reduce(temp2, raw_temp);
    vect_add(s1, temp1, temp2, VEC_N_SIZE_64);
    vect_add(o->s0, o->s0, s, VEC_N_SIZE_64);
    vect_add(o->s4, o->s4, s1, VEC_N_SIZE_64);

    shake_prng(seed, SEED_BYTES);
    seedexpander_init(&mask_seedexpander, seed, SEED_BYTES);
    vect_set_random_fixed_weight(&mask_seedexpander, s, weight);
    memset(raw_temp, 0x00, (VEC_N_SIZE_64*2+1)*8);
    fast_convolution_mult(raw_temp+(5*(VEC_N_SIZE_64/6)),
                          a1+(0*(weight/6)), a2+(5*(VEC_N_SIZE_64/6)),
                          weight/6, VEC_N_SIZE_64 - (VEC_N_SIZE_64/6)*5);
    reduce(temp1, raw_temp);
    vect_add(temp1, temp1, s, VEC_N_SIZE_64);
    memset(raw_temp, 0x00, (VEC_N_SIZE_64*2+1)*8);
    fast_convolution_mult(raw_temp+(0*(VEC_N_SIZE_64/6)),
                          a1+(5*(weight/6)), a2+(0*(VEC_N_SIZE_64/6)),
                          weight - (weight/6)*5, VEC_N_SIZE_64/6);
    reduce(temp2, raw_temp);
    vect_add(s1, temp1, temp2, VEC_N_SIZE_64);
    vect_add(o->s0, o->s0, s, VEC_N_SIZE_64);
    vect_add(o->s5, o->s5, s1, VEC_N_SIZE_64);

    shake_prng(seed, SEED_BYTES);
    seedexpander_init(&mask_seedexpander, seed, SEED_BYTES);
    vect_set_random_fixed_weight(&mask_seedexpander, s, weight);
    memset(raw_temp, 0x00, (VEC_N_SIZE_64*2+1)*8);
    fast_convolution_mult(raw_temp+(2*(VEC_N_SIZE_64/6)),
                          a1+(1*(weight/6)), a2+(2*(VEC_N_SIZE_64/6)),
                          weight/6, VEC_N_SIZE_64/6);
    reduce(temp1, raw_temp);
    vect_add(temp1, temp1, s, VEC_N_SIZE_64);
    memset(raw_temp, 0x00, (VEC_N_SIZE_64*2+1)*8);
    fast_convolution_mult(raw_temp+(1*(VEC_N_SIZE_64/6)),
                          a1+(2*(weight/6)), a2+(1*(VEC_N_SIZE_64/6)),
                          weight/6, VEC_N_SIZE_64/6);
    reduce(temp2, raw_temp);
    vect_add(s1, temp1, temp2, VEC_N_SIZE_64);
    vect_add(o->s1, o->s1, s, VEC_N_SIZE_64);
    vect_add(o->s2, o->s2, s1, VEC_N_SIZE_64);

    shake_prng(seed, SEED_BYTES);
    seedexpander_init(&mask_seedexpander, seed, SEED_BYTES);
    vect_set_random_fixed_weight(&mask_seedexpander, s, weight);
    memset(raw_temp, 0x00, (VEC_N_SIZE_64*2+1)*8);
    fast_convolution_mult(raw_temp+(3*(VEC_N_SIZE_64/6)),
                          a1+(1*(weight/6)), a2+(3*(VEC_N_SIZE_64/6)),
                          weight/6, VEC_N_SIZE_64/6);
    reduce(temp1, raw_temp);
    vect_add(temp1, temp1, s, VEC_N_SIZE_64);
    memset(raw_temp, 0x00, (VEC_N_SIZE_64*2+1)*8);
    fast_convolution_mult(raw_temp+(1*(VEC_N_SIZE_64/6)),
                          a1+(3*(weight/6)), a2+(1*(VEC_N_SIZE_64/6)),
                          weight/6, VEC_N_SIZE_64/6);
    reduce(temp2, raw_temp);
    vect_add(s1, temp1, temp2, VEC_N_SIZE_64);
    vect_add(o->s1, o->s1, s, VEC_N_SIZE_64);
    vect_add(o->s3, o->s3, s1, VEC_N_SIZE_64);

    shake_prng(seed, SEED_BYTES);
    seedexpander_init(&mask_seedexpander, seed, SEED_BYTES);
    vect_set_random_fixed_weight(&mask_seedexpander, s, weight);
    memset(raw_temp, 0x00, (VEC_N_SIZE_64*2+1)*8);
    fast_convolution_mult(raw_temp+(4*(VEC_N_SIZE_64/6)),
                          a1+(1*(weight/6)), a2+(4*(VEC_N_SIZE_64/6)),
                          weight/6, VEC_N_SIZE_64/6);
    reduce(temp1, raw_temp);
    vect_add(temp1, temp1, s, VEC_N_SIZE_64);
    memset(raw_temp, 0x00, (VEC_N_SIZE_64*2+1)*8);
    fast_convolution_mult(raw_temp+(1*(VEC_N_SIZE_64/6)),
                          a1+(4*(weight/6)), a2+(1*(VEC_N_SIZE_64/6)),
                          weight/6, VEC_N_SIZE_64/6);
    reduce(temp2, raw_temp);
    vect_add(s1, temp1, temp2, VEC_N_SIZE_64);
    vect_add(o->s1, o->s1, s, VEC_N_SIZE_64);
    vect_add(o->s4, o->s4, s1, VEC_N_SIZE_64);

    shake_prng(seed, SEED_BYTES);
    seedexpander_init(&mask_seedexpander, seed, SEED_BYTES);
    vect_set_random_fixed_weight(&mask_seedexpander, s, weight);
    memset(raw_temp, 0x00, (VEC_N_SIZE_64*2+1)*8);
    fast_convolution_mult(raw_temp+(5*(VEC_N_SIZE_64/6)),
                          a1+(1*(weight/6)), a2+(5*(VEC_N_SIZE_64/6)),
                          weight/6, VEC_N_SIZE_64 - (VEC_N_SIZE_64/6)*5);
    reduce(temp1, raw_temp);
    vect_add(temp1, temp1, s, VEC_N_SIZE_64);
    memset(raw_temp, 0x00, (VEC_N_SIZE_64*2+1)*8);
    fast_convolution_mult(raw_temp+(1*(VEC_N_SIZE_64/6)),
                          a1+(5*(weight/6)), a2+(1*(VEC_N_SIZE_64/6)),
                          weight - (weight/6)*5, VEC_N_SIZE_64/6);
    reduce(temp2, raw_temp);
    vect_add(s1, temp1, temp2, VEC_N_SIZE_64);
    vect_add(o->s1, o->s1, s, VEC_N_SIZE_64);
    vect_add(o->s5, o->s5, s1, VEC_N_SIZE_64);

    shake_prng(seed, SEED_BYTES);
    seedexpander_init(&mask_seedexpander, seed, SEED_BYTES);
    vect_set_random_fixed_weight(&mask_seedexpander, s, weight);
    memset(raw_temp, 0x00, (VEC_N_SIZE_64*2+1)*8);
    fast_convolution_mult(raw_temp+(3*(VEC_N_SIZE_64/6)),
                          a1+(2*(weight/6)), a2+(3*(VEC_N_SIZE_64/6)),
                          weight/6, VEC_N_SIZE_64/6);
    reduce(temp1, raw_temp);
    vect_add(temp1, temp1, s, VEC_N_SIZE_64);
    memset(raw_temp, 0x00, (VEC_N_SIZE_64*2+1)*8);
    fast_convolution_mult(raw_temp+(2*(VEC_N_SIZE_64/6)),
                          a1+(3*(weight/6)), a2+(2*(VEC_N_SIZE_64/6)),
                          weight/6, VEC_N_SIZE_64/6);
    reduce(temp2, raw_temp);
    vect_add(s1, temp1, temp2, VEC_N_SIZE_64);
    vect_add(o->s2, o->s2, s, VEC_N_SIZE_64);
    vect_add(o->s3, o->s3, s1, VEC_N_SIZE_64);

    shake_prng(seed, SEED_BYTES);
    seedexpander_init(&mask_seedexpander, seed, SEED_BYTES);
    vect_set_random_fixed_weight(&mask_seedexpander, s, weight);
    memset(raw_temp, 0x00, (VEC_N_SIZE_64*2+1)*8);
    fast_convolution_mult(raw_temp+(4*(VEC_N_SIZE_64/6)),
                          a1+(2*(weight/6)), a2+(4*(VEC_N_SIZE_64/6)),
                          weight/6, VEC_N_SIZE_64/6);
    reduce(temp1, raw_temp);
    vect_add(temp1, temp1, s, VEC_N_SIZE_64);
    memset(raw_temp, 0x00, (VEC_N_SIZE_64*2+1)*8);
    fast_convolution_mult(raw_temp+(2*(VEC_N_SIZE_64/6)),
                          a1+(4*(weight/6)), a2+(2*(VEC_N_SIZE_64/6)),
                          weight/6, VEC_N_SIZE_64/6);
    reduce(temp2, raw_temp);
    vect_add(s1, temp1, temp2, VEC_N_SIZE_64);
    vect_add(o->s2, o->s2, s, VEC_N_SIZE_64);
    vect_add(o->s4, o->s4, s1, VEC_N_SIZE_64);

    shake_prng(seed, SEED_BYTES);
    seedexpander_init(&mask_seedexpander, seed, SEED_BYTES);
    vect_set_random_fixed_weight(&mask_seedexpander, s, weight);
    memset(raw_temp, 0x00, (VEC_N_SIZE_64*2+1)*8);
    fast_convolution_mult(raw_temp+(5*(VEC_N_SIZE_64/6)),
                          a1+(2*(weight/6)), a2+(5*(VEC_N_SIZE_64/6)),
                          weight/6, VEC_N_SIZE_64 - (VEC_N_SIZE_64/6)*5);
    reduce(temp1, raw_temp);
    vect_add(temp1, temp1, s, VEC_N_SIZE_64);
    memset(raw_temp, 0x00, (VEC_N_SIZE_64*2+1)*8);
    fast_convolution_mult(raw_temp+(2*(VEC_N_SIZE_64/6)),
                          a1+(5*(weight/6)), a2+(2*(VEC_N_SIZE_64/6)),
                          weight - (weight/6)*5, VEC_N_SIZE_64/6);
    reduce(temp2, raw_temp);
    vect_add(s1, temp1, temp2, VEC_N_SIZE_64);
    vect_add(o->s2, o->s2, s, VEC_N_SIZE_64);
    vect_add(o->s5, o->s5, s1, VEC_N_SIZE_64);

    shake_prng(seed, SEED_BYTES);
    seedexpander_init(&mask_seedexpander, seed, SEED_BYTES);
    vect_set_random_fixed_weight(&mask_seedexpander, s, weight);
    memset(raw_temp, 0x00, (VEC_N_SIZE_64*2+1)*8);
    fast_convolution_mult(raw_temp+(4*(VEC_N_SIZE_64/6)),
                          a1+(3*(weight/6)), a2+(4*(VEC_N_SIZE_64/6)),
                          weight/6, VEC_N_SIZE_64/6);
    reduce(temp1, raw_temp);
    vect_add(temp1, temp1, s, VEC_N_SIZE_64);
    memset(raw_temp, 0x00, (VEC_N_SIZE_64*2+1)*8);
    fast_convolution_mult(raw_temp+(3*(VEC_N_SIZE_64/6)),
                          a1+(4*(weight/6)), a2+(3*(VEC_N_SIZE_64/6)),
                          weight/6, VEC_N_SIZE_64/6);
    reduce(temp2, raw_temp);
    vect_add(s1, temp1, temp2, VEC_N_SIZE_64);
    vect_add(o->s3, o->s3, s, VEC_N_SIZE_64);
    vect_add(o->s4, o->s4, s1, VEC_N_SIZE_64);

    shake_prng(seed, SEED_BYTES);
    seedexpander_init(&mask_seedexpander, seed, SEED_BYTES);
    vect_set_random_fixed_weight(&mask_seedexpander, s, weight);
    memset(raw_temp, 0x00, (VEC_N_SIZE_64*2+1)*8);
    fast_convolution_mult(raw_temp+(5*(VEC_N_SIZE_64/6)),
                          a1+(3*(weight/6)), a2+(5*(VEC_N_SIZE_64/6)),
                          weight/6, VEC_N_SIZE_64 - (VEC_N_SIZE_64/6)*5);
    reduce(temp1, raw_temp);
    vect_add(temp1, temp1, s, VEC_N_SIZE_64);
    memset(raw_temp, 0x00, (VEC_N_SIZE_64*2+1)*8);
    fast_convolution_mult(raw_temp+(3*(VEC_N_SIZE_64/6)),
                          a1+(5*(weight/6)), a2+(3*(VEC_N_SIZE_64/6)),
                          weight - (weight/6)*5, VEC_N_SIZE_64/6);
    reduce(temp2, raw_temp);
    vect_add(s1, temp1, temp2, VEC_N_SIZE_64);
    vect_add(o->s3, o->s3, s, VEC_N_SIZE_64);
    vect_add(o->s5, o->s5, s1, VEC_N_SIZE_64);

    shake_prng(seed, SEED_BYTES);
    seedexpander_init(&mask_seedexpander, seed, SEED_BYTES);
    vect_set_random_fixed_weight(&mask_seedexpander, s, weight);
    memset(raw_temp, 0x00, (VEC_N_SIZE_64*2+1)*8);
    fast_convolution_mult(raw_temp+(5*(VEC_N_SIZE_64/6)),
                          a1+(4*(weight/6)), a2+(5*(VEC_N_SIZE_64/6)),
                          weight/6, VEC_N_SIZE_64 - (VEC_N_SIZE_64/6)*5);
    reduce(temp1, raw_temp);
    vect_add(temp1, temp1, s, VEC_N_SIZE_64);
    memset(raw_temp, 0x00, (VEC_N_SIZE_64*2+1)*8);
    fast_convolution_mult(raw_temp+(4*(VEC_N_SIZE_64/6)),
                          a1+(5*(weight/6)), a2+(4*(VEC_N_SIZE_64/6)),
                          weight - (weight/6)*5, VEC_N_SIZE_64/6);
    reduce(temp2, raw_temp);
    vect_add(s1, temp1, temp2, VEC_N_SIZE_64);
    vect_add(o->s4, o->s4, s, VEC_N_SIZE_64);
    vect_add(o->s5, o->s5, s1, VEC_N_SIZE_64);
#endif
}
