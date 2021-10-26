#include "../common/api.h"
#include "../common/parameters.h"
#include "../common/vector.h"
#include "../fields/gf2x.h"
#include "../lib/shake_ds.h"
#include "board_config.h"
#include <stdint.h>
#include <string.h>
#include "timing_stats.h"
#include "../hqc/hqc.h"



int main() {
#ifdef CROSSCOMPILE
    setup();
    timer_init();
#endif
    const int ITERATIONS = 1000;

    unsigned char pk[PUBLIC_KEY_BYTES];
    unsigned char sk[SECRET_KEY_BYTES];
    uint64_t m[VEC_K_SIZE_64] = {0};
    uint8_t theta[SHAKE256_512_BYTES] = {0};
    uint64_t u[VEC_N_SIZE_64] = {0};
    uint64_t v[VEC_N1N2_SIZE_64] = {0};
    shake256incctx shake256state;

    uint32_t y[PARAM_OMEGA] = {0};
    seedexpander_state sk_seedexpander;

    uint8_t sk_seed[SEED_BYTES] = {0};
    uint64_t mulres[VEC_N_SIZE_64] = {0};
    uint64_t mulmask[VEC_N_SIZE_64] = {0};


    // timers declaration
    uint32_t start, end;
    welford_t enc_timer, dec_timer, mul_timer;


    // initialize timers
    welford_init(&enc_timer);
    welford_init(&dec_timer);
    welford_init(&mul_timer);


#ifdef CROSSCOMPILE
    ledOn();
#endif
    for(int i = 0; i < ITERATIONS; i++) {
        hqc_pke_keygen(pk, sk);

        vect_set_random_from_prng(m);
        shake256_512_ds(&shake256state, theta, (uint8_t*) m, VEC_K_SIZE_BYTES, G_FCT_DOMAIN);
        start = rdtsc();
        hqc_pke_encrypt(u, v, m, theta, pk);
        end = rdtsc();
        welford_update(&enc_timer, ((long double) (end - start)));

        memcpy(sk_seed, sk, SEED_BYTES);
        seedexpander_init(&sk_seedexpander, sk_seed, SEED_BYTES);
        vect_set_random_fixed_weight_by_coordinates(&sk_seedexpander, y, PARAM_OMEGA);
        start = rdtsc();
        safe_mul(mulres, mulmask, y, u, PARAM_OMEGA);
        end = rdtsc();
        welford_update(&mul_timer, ((long double) (end - start)));

        start = rdtsc();
        hqc_pke_decrypt(m, u, v, sk);
        end = rdtsc();
        welford_update(&dec_timer, ((long double) (end - start)));
    }


#ifdef DEBUG
    printf("\r\nEncryption \r\n");
    welford_print(enc_timer);
    printf("\r\nDecryption \r\n");
    welford_print(dec_timer);
    printf("\r\nMultiplication \r\n");
    welford_print(mul_timer);
#endif

#ifdef CROSSCOMPILE
    ledOff();
    printf("\r\nDONE\r\n");
#endif

}