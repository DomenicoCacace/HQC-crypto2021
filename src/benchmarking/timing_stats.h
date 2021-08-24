#pragma once
#include <math.h>

#ifdef CROSSCOMPILE

#define  ARM_CM_DEMCR      (*(uint32_t *)0xE000EDFC)
#define  ARM_CM_DWT_CTRL   (*(uint32_t *)0xE0001000)
#define  ARM_CM_DWT_CYCCNT (*(uint32_t *)0xE0001004)

void timer_init(void) {
    if (ARM_CM_DWT_CTRL != 0) {                  // See if DWT is available
        ARM_CM_DEMCR      |= 1 << 24;            // Set bit 24
        ARM_CM_DWT_CYCCNT  = 0;                
        ARM_CM_DWT_CTRL   |= 1 << 0;             // Set bit 0
    }
}

static inline 
uint32_t rdtsc(void) {
    return ARM_CM_DWT_CYCCNT;
}
#else
static inline
uint32_t rdtsc(void) {
    unsigned long long result;
    __asm__ __volatile__(
        "rdtscp;"
        "shl $32, %%rdx;"
        "or %%rdx, %%rax"
        : "=a"(result)
        :
        : "%rcx", "%rdx");
    return result;
}
#endif

typedef struct {
    long double mean;
    long double M2;
    long count;
} welford_t;

static inline
void welford_init(welford_t* state) {
    state->mean = 0.0;
    state->M2 = 0.0;
    state->count = 0;
    return;
}

static inline
void welford_update(welford_t* state, long double sample) {
    long double delta, delta2;
    state->count = state->count + 1;
    delta = sample - state->mean;
    state->mean += delta / (long double)(state->count);
    delta2 = sample - state->mean;
    state->M2 += delta * delta2;
}

static inline
void welford_print(const welford_t state) {
    printf("%.2Lf, %.2Lf\r\n",
           state.mean,
           sqrtl(state.M2/(long double)(state.count-1)));
}

static inline
double welch_t_statistic(const welford_t state1,
                         const welford_t state2) {
    long double num, den, var1, var2;
    var1 = state1.M2/(long double)(state1.count-1);
    var2 = state2.M2/(long double)(state2.count-1);

    num = state1.mean - state2.mean;
    den = sqrtl(var1/(long double) state1.count + var2/(long double) state2.count );

    return num/den;
}
