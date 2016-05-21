#include <math.h>
#include <string.h>
#include <strings.h>
#include <stdlib.h>
#include <complex.h>

#include <stdint.h>
#include <stdio.h>
#include <stdbool.h>
#include <assert.h>
#include <unistd.h>
#include <inttypes.h>

//Calculate the value of n_prb
/*
 * Calculate the value of number of physical RBs n_prb
 * input- N_sc , n_ul_rb for normal and extended CP
 * output- N_PRB for extended and normal CP;
 */
uint32_t N_prb(uint32_t *n_prb, uint32_t SRS_UL.n_ul_rb, uint32_t N_sc)
{
    uint32_t i;
    uint32_t k;
    /* resource elements in frequency domain 0.....n_ul_rb*N_sc*/
    k = SRS_UL.n_ul_rb * N_sc;
    n_prb = floor(k/N_sc);
}

