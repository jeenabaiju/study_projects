
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






void calc_prs_c(const uint32_t c_init, const uint32_t len, uint8_t* prs_c){
	uint32_t N_c = 1600;
	uint32_t x1, x2;
	x1 = 0x54D21B24; // x1 sequence moved forward
	x2 = c_init;

	uint32_t next_x1 = 0;
	uint32_t next_x2 = 0;

	// Move x2 forward
	uint32_t i;
	for(i = 0; i < N_c - 31; i++){
		next_x2 = ((x2 >> 3) ^ (x2 >> 2) ^ (x2 >> 1) ^ x2) & 0x1;

        x2 = (x2 >> 1) | (next_x2 << 30);
	}

	// Calculate prs_c
	for(i = 0; i < len; i++){

		next_x1 = ((x1 >> 3) ^ x1) & 0x1;
		next_x2 = ((x2 >> 3) ^ (x2 >> 2) ^ (x2 >> 1) ^ x2) & 0x1;

		*(prs_c++) = next_x1 ^ next_x2;

		x1 = (x1 >> 1) | (next_x1 << 30);
		x2 = (x2 >> 1) | (next_x2 << 30);
	}

}
uint32_t Sequence_Hopping(const uint16_t cell_ID, uint32_t delta_ss, uint32_t M_sc, uint32_t N_sc, uint32_t ns, uint32_t sequence_hopping, uint32_t group_hopping)
{
    /*Sequence hopping only applies for reference-signals of length M_sc_RS ≥ 6*N_sc_RB .For M sc RS < 6 N sc RB ,  v within u is v = 0 .*/
    uint32_t v;
    if ( M_sc >= 6*N_sc)
    {
        uint32_t len = ns;
        uint8_t n_prs[len];
        uint32_t c_init = ((cell_ID / 30) << 5) + (((cell_ID % 30) + delta_ss) % 30);
        calc_prs_c( c_init, len, &n_prs); /*generate_pseudo random sequence*/
        if(sequence_hopping ==1)
        {
            if(group_hopping == 1)
            {
                v = 0;
            }
            else
            {
                v=n_prs[ns];
                printf("n_prs = %d\n ",n_prs[ns]);
            }
        }
    }
    else
    {
        v = 0;
    }
return v;
}
uint32_t Group_hopping_f_gh(uint32_t *f_gh, const uint16_t cell_ID, uint32_t ns, uint32_t group_hopping)
{
    uint32_t len = ( 8 * ns + 7 );
    uint8_t n_prs[len];
    uint32_t c_init = floor( cell_ID / 30 );
    /*generate_pseudo random sequence /** Computes n_prs values as defined in 5.5.2.1.1 of 36.211 */
    calc_prs_c(c_init, len, &n_prs);
    /** Computes sequence-group pattern f_gh according to 5.5.1.3 of 36.211 */
    f_gh[ns] = 0;
    int i ;
    if(group_hopping)
    {
        for (i = 0; i < 8; i++)
        {
            f_gh[ns] += (((uint32_t) n_prs[8 * ns + i]) << i) % 30;/*

                         i=7
            summation of [c( 8*ns+ i )⋅ 2^i ]mod 30 if group hopping is enabled
                         i=0                                                      */
        }
    }
    else
    {
        f_gh[ns] = 0;
    }
}
/*
 * Calculate f_ss
 * input-  cellID
 * output - f_ss value
 */
uint32_t get_f_ss(const uint16_t cell_ID, uint32_t N_ID_PUCCH, uint32_t N_ID_PUSCH, const uint16_t PUCCH_ID, const uint16_t PUSCH_ID)
{
    uint32_t f_ss;
    const uint16_t n_ID = Get_cellID( N_ID_PUCCH, N_ID_PUSCH, cell_ID, PUCCH_ID, PUSCH_ID);
    f_ss = n_ID % 30;
    return f_ss;
}
/*
 * Calculate Group Hopping(u)u=(f_gh(ns)+f_ss)%30;
 * input- number of slots, group hopping -enabled/disabled, cellID
 * output - u value
 */
uint32_t Get_u_value(const uint16_t cell_ID, uint32_t ns, uint32_t group_hopping, uint32_t N_ID_PUCCH, uint32_t N_ID_PUSCH, const uint16_t PUCCH_ID, const uint16_t PUSCH_ID)
{
    uint32_t f_ss;
    uint32_t u;
    uint32_t f_gh;
    f_ss = get_f_ss( cell_ID, N_ID_PUCCH, N_ID_PUSCH, PUCCH_ID, PUSCH_ID);
    Group_hopping_f_gh(&f_gh, cell_ID, ns, group_hopping);
    u = ( f_gh + f_ss ) % 30;
    return u;
}
/* Calculate N_CellID
* input - cell_ID, Configurations for N_ID_PUCCH,N_ID_PUSCH
* output- N_ID
*/
const uint16_t Get_cellID(uint32_t N_ID_PUCCH, uint32_t N_ID_PUSCH, const uint16_t cell_ID, const uint16_t PUCCH_ID, const uint16_t PUSCH_ID)
{
    const uint16_t n_ID;
    if (N_ID_PUCCH)
    {
        n_ID= PUCCH_ID;
    }
    else if (N_ID_PUSCH)
    {
        n_ID= PUCCH_ID;
    }
    else
    {
        n_ID= cell_ID;// for SRS also//
    }
    return n_ID;
}


void main()
{

const uint16_t cell_ID = 100;
uint32_t delta_ss = 0; // delta_ss Δ ss ∈ { 0 , 1 ,..., 29 }
uint32_t group_hopping = 1;
uint32_t sequence_hopping = 0;
uint32_t N_sc = 12;
uint32_t M_sc = 6*N_sc;
uint32_t ns = 20, v, u;
uint32_t f_gh;
uint32_t N_ID_PUCCH = 0;
uint32_t N_ID_PUSCH = 0;
const uint16_t PUCCH_ID = 0;
const uint16_t PUSCH_ID = 0;
v = Sequence_Hopping( cell_ID, delta_ss, M_sc, N_sc, ns, sequence_hopping, group_hopping);

u = Get_u_value( cell_ID, ns, group_hopping, N_ID_PUCCH, N_ID_PUSCH,  PUCCH_ID, PUSCH_ID);
printf("u = %d v = %d\n ",u,v);
}
