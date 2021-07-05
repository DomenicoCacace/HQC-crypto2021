/**
 * @file api.h
 * @brief NIST KEM API used by the HQC_KEM IND-CCA2 scheme
 */

#ifndef API_H
#define API_H


#if SECURITY_LEVEL == 128
    #define CRYPTO_ALGNAME                      "HQC-128"
    #define CRYPTO_SECRETKEYBYTES               2289
    #define CRYPTO_PUBLICKEYBYTES               2249
    #define CRYPTO_BYTES                        64
    #define CRYPTO_CIPHERTEXTBYTES              4481
#elif SECURITY_LEVEL == 192
    #define CRYPTO_ALGNAME                      "HQC-192"
    #define CRYPTO_SECRETKEYBYTES               4562
    #define CRYPTO_PUBLICKEYBYTES               4522
    #define CRYPTO_BYTES                        64
    #define CRYPTO_CIPHERTEXTBYTES              9026
#elif SECURITY_LEVEL == 256
    #define CRYPTO_ALGNAME                      "HQC-256"
    #define CRYPTO_SECRETKEYBYTES               7285
    #define CRYPTO_PUBLICKEYBYTES               7245
    #define CRYPTO_BYTES                        64
    #define CRYPTO_CIPHERTEXTBYTES              14469
#else
    #error INVALID SECURITY LEVEL
#endif 



// As a technicality, the public key is appended to the secret key in order to respect the NIST API.
// Without this constraint, CRYPTO_SECRETKEYBYTES would be defined as 32

int crypto_kem_keypair(unsigned char* pk, unsigned char* sk);
int crypto_kem_enc(unsigned char* ct, unsigned char* ss, const unsigned char* pk);
int crypto_kem_dec(unsigned char* ss, const unsigned char* ct, const unsigned char* sk);

#endif
