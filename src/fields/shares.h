#define MASKS (1 << MASK_LVL)

typedef struct {
    uint64_t share[MASKS][VEC_N_SIZE_64];
} shares_t;

void shares_reduce(shares_t shares, uint64_t *o, uint64_t *mask) {
    for(int i = 0; i < MASKS; i+=2) {
        vect_add(o, o, shares.share[i], VEC_N_SIZE_64);
        vect_add(mask, mask, shares.share[i+1], VEC_N_SIZE_64);

    }
}

uint16_t shares_size(int i, uint16_t weight) {
    if(i < MASKS-1)
        return weight/MASKS;
    return weight - weight/MASKS*i;
}

