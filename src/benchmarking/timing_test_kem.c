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
    unsigned char ct[CIPHERTEXT_BYTES];
    unsigned char key1[SHARED_SECRET_BYTES];
    unsigned char key2[SHARED_SECRET_BYTES];


    uint32_t start, end;
    welford_t enc_timer, dec_timer;

    welford_init(&enc_timer);
    welford_init(&dec_timer);

#ifdef CROSSCOMPILE
    ledOn();
#endif
    for(int i = 0; i < ITERATIONS; i++) {
        crypto_kem_keypair(pk, sk);

        start = rdtsc();
        crypto_kem_enc(ct, key1, pk);
        end = rdtsc();
        welford_update(&enc_timer, ((long double)(end - start)));


        start = rdtsc();
        crypto_kem_dec(key2, ct, sk);
        end = rdtsc();
        welford_update(&dec_timer, ((long double)(end - start)));
    }

#ifdef DEBUG
    printf("\r\nEncapsulation \r\n");
    welford_print(enc_timer);
    printf("\r\nDecapsulation \r\n");
    welford_print(dec_timer);
#endif

#ifdef CROSSCOMPILE
    ledOff();
    printf("\r\nDONE\r\n");
#endif
}
