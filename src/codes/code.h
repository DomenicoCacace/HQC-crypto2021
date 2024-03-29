#ifndef CODE_H
#define CODE_H

/**
 * @file code.h
 * Header file of code.c
 */

#include <stddef.h>
#include <stdint.h>

#include "../common/parameters.h"

void code_encode(uint64_t *codeword, const uint64_t *message);
void code_decode(uint64_t *message, const uint64_t *vector);

#endif
