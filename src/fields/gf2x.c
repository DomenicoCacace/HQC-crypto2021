/**
 * @file gf2x.c
 * @brief Implementation of multiplication of two polynomials
 */

#include <stdint.h>
#include <string.h>

#include "../common/parameters.h"
#include "gf2x.h"

#define TABLE 16
#define WORD 64

static inline void swap(uint16_t *tab, uint16_t elt1, uint16_t elt2);
static void reduce(uint64_t *o, const uint64_t *a);
static void fast_convolution_mult(uint64_t *o, const uint32_t *a1, const uint64_t *a2, const uint16_t weight, seedexpander_state *ctx);

/**
 * @brief swap two elements in a table
 *
 * This function exchanges tab[elt1] with tab[elt2]
 *
 * @param[in] tab Pointer to the table
 * @param[in] elt1 Index of the first element
 * @param[in] elt2 Index of the second element
 */
static inline void swap(uint16_t *tab, uint16_t elt1, uint16_t elt2) {
    uint16_t tmp = tab[elt1];

    tab[elt1] = tab[elt2];
    tab[elt2] = tmp;
}



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
 * @brief computes product of the polynomial a1(x) with the sparse polynomial a2
 *
 *  o(x) = a1(x)a2(x)
 *
 * @param[out] o Pointer to the result
 * @param[in] a1 Pointer to the sparse polynomial a2 (list of degrees of the monomials which appear in a2)
 * @param[in] a2 Pointer to the polynomial a1(x)
 * @param[in] weight Hamming wifht of the sparse polynomial a2
 * @param[in] ctx Pointer to a seed expander used to randomize the multiplication process
 */
static void fast_convolution_mult(uint64_t *o, const uint32_t *a1, const uint64_t *a2, const uint16_t weight, seedexpander_state *ctx){
    uint64_t carry;
    int32_t dec, s;
    uint64_t table[TABLE * (VEC_N_SIZE_64 + 1)];
    uint16_t permuted_table[TABLE];
    uint16_t permutation_table[TABLE];
    uint16_t permuted_sparse_vect[PARAM_OMEGA_E];
    uint16_t permutation_sparse_vect[PARAM_OMEGA_E];

    for (size_t i = 0; i < TABLE; i++) {
        permuted_table[i] = (uint16_t) i;
    }

    seedexpander(ctx, (uint8_t *) permutation_table, TABLE << 1);

    for (size_t i = 0; i < TABLE - 1; i++) {
        swap(permuted_table + i, 0, permutation_table[i] % (TABLE - i));
    }

    uint64_t *pt = table + (permuted_table[0] * (VEC_N_SIZE_64 + 1));

    for (size_t i = 0 ; i < VEC_N_SIZE_64 ; i++) {
        pt[i] = a2[i];
    }

    pt[VEC_N_SIZE_64] = 0x0UL;

    for (size_t i = 1; i < TABLE; i++) {
        carry = 0x0UL;
        int32_t idx = permuted_table[i] * (VEC_N_SIZE_64 + 1);
        uint64_t *pt = table + idx;

        for (size_t j = 0; j < VEC_N_SIZE_64; j++) {
            pt[j] = (a2[j] << i) ^ carry;
            carry = (a2[j] >> ((WORD - i)));
        }

        pt[VEC_N_SIZE_64] = carry;
    }

    for (size_t i = 0; i < weight; i++) {
        permuted_sparse_vect[i] = (uint16_t) i;
    }

    seedexpander(ctx, (uint8_t *) permutation_sparse_vect, weight << 1);

    for (int32_t i = 0; i < (weight - 1); i++) {
        swap(permuted_sparse_vect + i, 0, permutation_sparse_vect[i] % (weight - i));
    }

    for (size_t i = 0; i < weight; i++) {
        carry = 0x0UL;
        dec = a1[permuted_sparse_vect[i]] & 0xf;
        s = a1[permuted_sparse_vect[i]] >> 4;

        uint16_t *res_16 = (uint16_t *) o;
        res_16 += s;
        uint64_t *pt = table + (permuted_table[dec] * (VEC_N_SIZE_64 + 1));

        for (size_t j = 0; j < VEC_N_SIZE_64 + 1; j++) {
            uint64_t tmp = (uint64_t) res_16[0] | ((uint64_t) (res_16[1])) << 16 |
                (uint64_t) (res_16[2]) << 32 | ((uint64_t) (res_16[3])) << 48;
            tmp ^= pt[j];
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
void vect_mul(uint64_t *o, const uint32_t *a1, const uint64_t *a2, const uint16_t weight, seedexpander_state *ctx) {
    uint64_t tmp[(VEC_N_SIZE_64 << 1) + 1] = {0};

    fast_convolution_mult(tmp, a1, a2, weight, ctx);
    reduce(o, tmp);
}
