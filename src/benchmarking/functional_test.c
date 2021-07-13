#include <string.h>
#include "../common/api.h"
#include "../common/parameters.h"
#include "board_config.h"

int main() {
#ifdef CROSSCOMPILE
    setup();
#endif

    const int ITERATIONS = 1000;
    unsigned char pk[PUBLIC_KEY_BYTES];
    unsigned char sk[SECRET_KEY_BYTES];
    unsigned char ct[CIPHERTEXT_BYTES];
    unsigned char key1[SHARED_SECRET_BYTES];
    unsigned char key2[SHARED_SECRET_BYTES];

    int passed = 0;

#ifdef DEBUG
    printf("Running functional test - %d iterations \r\n\r\n", ITERATIONS);
#endif

      for(int i = 0; i < ITERATIONS; i++) {

        crypto_kem_keypair(pk, sk);
        crypto_kem_enc(ct, key1, pk);
        crypto_kem_dec(key2, ct, sk);

#ifdef DEBUG
        printf("Round %d\r\n", i+1);
        printf("secret1: ");
        for(int j = 0 ; j < SHARED_SECRET_BYTES ; ++j) printf("%x", key1[j]);

        printf("\r\nsecret2: ");
        for(int j = 0 ; j < SHARED_SECRET_BYTES ; ++j) printf("%x", key2[j]);
        printf("\r\n\r\n");
#endif

        if (memcmp(key1, key2, SHARED_SECRET_BYTES) == 0) {
            passed++;
        }
    }
#ifdef CROSSCOMPILE
    blink(50, 10);
#endif
      
    printf("%d\r\n", passed-ITERATIONS);
    return passed-ITERATIONS;
}