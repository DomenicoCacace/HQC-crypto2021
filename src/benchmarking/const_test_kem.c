#include "../common/api.h"
#include "../common/parameters.h"
#include "board_config.h"
#include <stdint.h>
#include "timing_stats.h"



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
    unsigned char ct[CIPHERTEXT_BYTES];
    unsigned char key1[SHARED_SECRET_BYTES];
    unsigned char key2[SHARED_SECRET_BYTES];


    // timers declaration
    uint32_t start, end;
    welford_t encaps_timer, decaps_timer;
    welford_t encaps_timer_0, decaps_timer_0;


    // initialize timers
    welford_init(&encaps_timer_0);
    welford_init(&decaps_timer_0);

    welford_init(&encaps_timer);
    welford_init(&decaps_timer);

    // same as hqc.pke_keygen
    crypto_kem_keypair(pk_0, sk_0);


#ifdef CROSSCOMPILE
    ledOn();
#endif
    for(int i = 0; i < ITERATIONS; i++) {

        start = rdtsc();
        crypto_kem_enc(ct, key1, pk_0);
        end = rdtsc();
        welford_update(&encaps_timer_0, ((long double) (end - start)));

        start = rdtsc();
        crypto_kem_dec(key2, ct, sk_0);
        end = rdtsc();
        welford_update(&decaps_timer_0, ((long double) (end - start)));

        crypto_kem_keypair(pk, sk);

        start = rdtsc();
        crypto_kem_enc(ct, key1, pk);
        end = rdtsc();
        welford_update(&encaps_timer, ((long double) (end - start)));

        start = rdtsc();
        crypto_kem_dec(key2, ct, sk);
        end = rdtsc();
        welford_update(&decaps_timer, ((long double) (end - start)));
    }


#ifdef DEBUG
    printf("\r\nEncapsulation \r\n");
    printf("%lf", welch_t_statistic(encaps_timer_0, encaps_timer));
    printf("\r\nDecapsulation \r\n");
    printf("%lf", welch_t_statistic(decaps_timer_0, decaps_timer));
#endif

#ifdef CROSSCOMPILE
    ledOff();
    printf("\r\nDONE\r\n");
#endif

}
