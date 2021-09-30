#include "../common/api.h"
#include "../common/parameters.h"
#include "board_config.h"
#include <stdint.h>
#include "timing_stats.h"

uint32_t mul_start, mul_end;

int main() {
#ifdef CROSSCOMPILE
    setup();
    timer_init();
#endif
    const int ITERATIONS = 1000;
    unsigned char pk[PUBLIC_KEY_BYTES];
    unsigned char sk[SECRET_KEY_BYTES];
    unsigned char ct[CIPHERTEXT_BYTES];
    unsigned char key1[SHARED_SECRET_BYTES];
    unsigned char key2[SHARED_SECRET_BYTES];

    // timers declaration
    uint32_t enc_start, enc_end;
    uint32_t dec_start, dec_end;
    welford_t enc_timer, dec_timer, mul_timer;
    welford_t enc_timer_0, dec_timer_0, mul_timer_0;


    // initialize timers
    welford_init(&enc_timer_0);
    welford_init(&dec_timer_0);
    welford_init(&mul_timer_0);

    welford_init(&enc_timer);
    welford_init(&dec_timer);
    welford_init(&mul_timer);


#ifdef CROSSCOMPILE
    ledOn();
#endif
    for(int i = 0; i < ITERATIONS; i++) {
        crypto_kem_keypair(pk, sk);

        enc_start = rdtsc();
        crypto_kem_enc(ct, key1, pk);
        enc_end = rdtsc();

        dec_start = rdtsc();
        crypto_kem_dec(key2, ct, sk);
        dec_end = rdtsc();

        welford_update(&enc_timer, ((long double) (enc_end - enc_start)));
        welford_update(&dec_timer, ((long double) (dec_end - dec_start)));
        welford_update(&mul_timer, ((long double) (mul_end - mul_start)));

    }

    // generate only once
    crypto_kem_keypair(pk, sk);

    for(int i = 0; i < ITERATIONS; i++) {

        enc_start = rdtsc();
        crypto_kem_enc_const(ct, key1, pk);
        enc_end = rdtsc();

        dec_start = rdtsc();
        crypto_kem_dec(key2, ct, sk);
        dec_end = rdtsc();

        welford_update(&enc_timer_0, ((long double)(enc_end - enc_start)));
        welford_update(&dec_timer_0, ((long double)(dec_end - dec_start)));
        welford_update(&mul_timer_0, ((long double)(mul_end - mul_start)));
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
