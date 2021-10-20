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
    unsigned char pk_0[PUBLIC_KEY_BYTES];
    unsigned char sk_0[SECRET_KEY_BYTES];
    uint64_t m[VEC_K_SIZE_64] = {0};
    uint64_t m_0[VEC_K_SIZE_64] = {0};
    uint8_t theta[SHAKE256_512_BYTES] = {0};
    uint8_t theta_0[SHAKE256_512_BYTES] = {0};
    uint64_t u[VEC_N_SIZE_64] = {0};
    uint64_t v[VEC_N1N2_SIZE_64] = {0};
    shake256incctx shake256state;

    uint32_t y[PARAM_OMEGA] = {0};
    uint32_t y_0[PARAM_OMEGA] = {0};
    seedexpander_state sk_seedexpander;
    seedexpander_state perm_seedexpander;
    uint8_t sk_seed[SEED_BYTES] = {0};
    uint8_t perm_seed[SEED_BYTES] = {0};
    uint64_t mulres[VEC_N_SIZE_64] = {0};

    // timers declaration
    uint32_t start, end;
    welford_t enc_timer, dec_timer, mul_timer;
    welford_t enc_timer_0, dec_timer_0, mul_timer_0;


    // initialize timers
    welford_init(&enc_timer_0);
    welford_init(&dec_timer_0);
    welford_init(&mul_timer_0);

    welford_init(&enc_timer);
    welford_init(&dec_timer);
    welford_init(&mul_timer);

    hqc_pke_keygen(pk_0, sk_0);

    vect_set_random_from_prng(m_0);
    shake256_512_ds(&shake256state, theta_0, (uint8_t*) m_0, VEC_K_SIZE_BYTES, G_FCT_DOMAIN);

    memcpy(sk_seed, sk_0, SEED_BYTES);
    seedexpander_init(&sk_seedexpander, sk_seed, SEED_BYTES);
    vect_set_random_fixed_weight_by_coordinates(&sk_seedexpander, y_0, PARAM_OMEGA);


#ifdef CROSSCOMPILE
    ledOn();
#endif
    for(int i = 0; i < ITERATIONS; i++) {
        hqc_pke_keygen(pk, sk);

        start = rdtsc();
        hqc_pke_encrypt(u, v, m_0, theta_0, pk_0);
        end = rdtsc();
        welford_update(&enc_timer_0, ((long double) (end - start)));

        shake_prng(perm_seed, SEED_BYTES);
        seedexpander_init(&perm_seedexpander, perm_seed, SEED_BYTES);
        start = rdtsc();
        vect_mul(mulres, y_0, u, PARAM_OMEGA, &perm_seedexpander);
        end = rdtsc();
        welford_update(&mul_timer_0, ((long double) (end - start)));

        start = rdtsc();
        hqc_pke_decrypt(m_0, u, v, sk_0);
        end = rdtsc();
        welford_update(&dec_timer_0, ((long double) (end - start)));


        vect_set_random_from_prng(m);
        shake256_512_ds(&shake256state, theta, (uint8_t*) m, VEC_K_SIZE_BYTES, G_FCT_DOMAIN);
        start = rdtsc();
        hqc_pke_encrypt(u, v, m, theta, pk);
        end = rdtsc();
        welford_update(&enc_timer, ((long double) (end - start)));

        memcpy(sk_seed, sk, SEED_BYTES);
        seedexpander_init(&sk_seedexpander, sk_seed, SEED_BYTES);
        vect_set_random_fixed_weight_by_coordinates(&sk_seedexpander, y, PARAM_OMEGA);
        shake_prng(perm_seed, SEED_BYTES);
        seedexpander_init(&perm_seedexpander, perm_seed, SEED_BYTES);
        start = rdtsc();
        vect_mul(mulres, y, u, PARAM_OMEGA, &perm_seedexpander);
        end = rdtsc();
        welford_update(&mul_timer, ((long double) (end - start)));

        start = rdtsc();
        hqc_pke_decrypt(m, u, v, sk);
        end = rdtsc();
        welford_update(&dec_timer, ((long double) (end - start)));

    }


#ifdef DEBUG
    printf("\r\nEncryption \r\n");
    printf("%lf", welch_t_statistic(enc_timer_0, enc_timer));
    printf("\r\nDecryption \r\n");
    printf("%lf", welch_t_statistic(dec_timer_0, dec_timer));
    printf("\r\nMultiplication \r\n");
    printf("%lf", welch_t_statistic(mul_timer_0, mul_timer));
#endif

#ifdef CROSSCOMPILE
    ledOff();
    printf("\r\nDONE\r\n");
#endif

}