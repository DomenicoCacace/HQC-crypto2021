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
    const int ITERATIONS = 10;
    unsigned char pk[PUBLIC_KEY_BYTES];
    unsigned char sk[SECRET_KEY_BYTES];
    unsigned char ct[CIPHERTEXT_BYTES];
    unsigned char key1[SHARED_SECRET_BYTES];
    unsigned char key2[SHARED_SECRET_BYTES];


    uint32_t gen_start, gen_end;
    uint32_t enc_start, enc_end;
    uint32_t dec_start, dec_end;
    welford_t gen_timer, enc_timer, dec_timer;

    welford_init(&gen_timer);
    welford_init(&enc_timer);
    welford_init(&dec_timer);

    // heat cache to minimize variance
    crypto_kem_keypair(pk, sk);
    crypto_kem_enc(ct, key1, pk);
    crypto_kem_dec(key2, ct, sk);

#ifdef CROSSCOMPILE
    ledOn();
#endif
    for(int i = 0; i < ITERATIONS; i++) {
        gen_start = rdtsc();
        crypto_kem_keypair(pk, sk);
        gen_end = rdtsc();

        enc_start = rdtsc();
        crypto_kem_enc(ct, key1, pk);
        enc_end = rdtsc();

        dec_start = rdtsc();
        crypto_kem_dec(key2, ct, sk);
        dec_end = rdtsc();

        welford_update(&gen_timer, ((long double)(gen_end - gen_start)));
        welford_update(&enc_timer, ((long double)(enc_end - enc_start)));
        welford_update(&dec_timer, ((long double)(dec_end - dec_start)));
    }

#ifdef DEBUG
    printf("\r\nKey Generation \r\n");
    welford_print(gen_timer);
    printf("\r\nEncryption \r\n");
    welford_print(enc_timer);
    printf("\r\nDecryption \r\n");
    welford_print(dec_timer);
#endif

#ifdef CROSSCOMPILE
    ledOff();
    printf("\r\nDONE\r\n");
#endif

}