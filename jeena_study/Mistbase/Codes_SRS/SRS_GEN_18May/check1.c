#include <unistd.h>
#include <inttypes.h>
#include <math.h>
#include <string.h>
#include <strings.h>
#include <stdlib.h>
#include <complex.h>
#include <stdint.h>
#include <stdio.h>
#include <stdbool.h>
#include <assert.h>


/* //Table 5.2.3-1: Resource block parameters according to 3GPP 36.211 5.2.3*/
uint32_t N_sc = 12;
uint32_t N_symbol_ul[2] = {7,6};
uint32_t m_srs_b[4][4][8] = {{
                        /* m_srs for 6<=n_rb<=40. Table 5.5.3.2-1 */
                        {36, 32, 24, 20, 16, 12, 8, 4},
                        {12, 16,  4,  4,  4,  4, 4, 4},
                        { 4,  8,  4,  4,  4,  4, 4, 4},
                        { 4,  4,  4,  4,  4,  4, 4, 4}},
                        {
                        /* m_srs for 40<n_rb<60. Table 5.5.3.2-2 */
                        {48, 48, 40, 36, 32, 24, 20, 16},
                        {24, 16, 20, 12, 16,  4,  4,  4},
                        {12,  8,  4,  4,  8,  4,  4,  4},
                        { 4,  4,  4,  4,  4,  4,  4,  4}},
                        {
                        /* m_srs for 60<n_rb<80. Table 5.5.3.2-3 */
                        {72, 64, 60, 48, 48, 40, 36, 32},
                        {24, 32, 20, 24, 16, 20, 12, 16},
                        {12, 16,  4, 12,  8,  4,  4,  8},
                        { 4,  4,  4,  4,  4,  4,  4,  4}},

                        {
                        /* m_srs for 80<n_rb<110. Table 5.5.3.2-4 */
                        {96, 96, 80, 72, 64, 60, 48, 48},
                        {48, 32, 40, 24, 32, 20, 24, 16},
                        {24, 16, 20, 12, 16,  4, 12,  8},
                        { 4,  4,  4,  4,  4,  4,  4,  4}}};

/*
/*
 * Calculate the value of M_sc
 * input-  NULRB,Bandwidth B and bw_configuration
 * output -M_sc m_srs
 */
uint32_t srsbwtable_idx(uint32_t n_ul_rb)
{
    if (n_ul_rb <= 40)
    {
        return 0;
    }
    else if (n_ul_rb <= 60)
    {
        return 1;
    }
    else if (n_ul_rb <= 80)
    {
        return 2;
    }
    else
    {
        return 3;
    }
}
uint32_t Get_Msc_values(uint32_t bw_cfg, uint32_t B, uint32_t n_ul_rb, uint32_t N_sc)
{
   
   return  m_srs_b[srsbwtable_idx(n_ul_rb)][B][bw_cfg] * N_sc / 2;/* According to 3GPP 36.211 5.5.3.2*/
}
/*
 * Calculate the value of N_zc
 * input- M_sc 
 * output -N_zc largest prime number less than M_sc
 */
uint32_t Get_Nzc(uint32_t M_sc)
{   
    
 /* get largest prime no N_zc<M_sc */
    int i,j;
    int prime;
    uint32_t N_zc;
    for(i = 1; i < M_sc; i++)
    {
        for(j = 2; j < i; j++) if(i % j == 0) break;
        if(j == i)
            prime = j;/* find prime numbers*/
    }
    N_zc = prime;
    return N_zc ;
}
/*
 * Calculate q value
 * input- u,v
 * output -q value to find Zadoff chu seq
 */
 static float get_qvalue(uint32_t u, uint32_t v, uint32_t N_zc)
{
    float q;
    float q_bar;
    float n_zc = (float) N_zc;
    q_bar = n_zc * (u + 1) / 31;
    q = ( floor ( q_bar + 0.5) ) + ( v * ( pow( (-1) , ( floor( 2 * q_bar) ) ) ) );
    return q;
}
/* Calculate N_CellID
* input - cell_ID, Configurations for N_ID_PUCCH,N_ID_PUSCH
* output- N_ID
*/
const uint16_t Get_cellID(uint32_t N_ID_PUCCH, uint32_t N_ID_PUSCH, const uint16_t cell_ID, const uint16_t PUCCH_ID, const uint16_t PUSCH_ID)
{
    const uint16_t N_ID;
    if (N_ID_PUCCH == 1)
    {
        const uint16_t N_ID= PUCCH_ID;
        return  N_ID;
    }
    else if (N_ID_PUSCH == 1)
    {
        const uint16_t N_ID= PUSCH_ID;
        return  N_ID;
    }
    else
    {
        const uint16_t N_ID= cell_ID;// for SRS also//
        return  N_ID;
    }
}


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
uint32_t get_f_gh(const uint16_t cell_ID, uint32_t ns, uint32_t group_hopping, uint32_t N_ID_PUCCH, uint32_t N_ID_PUSCH, const uint16_t PUCCH_ID, const uint16_t PUSCH_ID)
{
    uint32_t len = ( 8 * ns + 7 );
    uint8_t n_prs[len];
    uint32_t count = 0;
    uint32_t f_gh;
    const uint16_t n_ID = Get_cellID( N_ID_PUCCH, N_ID_PUSCH, cell_ID, PUCCH_ID, PUSCH_ID);
    uint32_t c_init = floor( n_ID / 30 );
    /*generate_pseudo random sequence /** Computes n_prs values as defined in 5.5.2.1.1 of 36.211 */
    calc_prs_c(c_init, len, n_prs);
    /** Computes sequence-group pattern f_gh according to 5.5.1.3 of 36.211 */
    int i ;
    if(group_hopping == 1)  
    {   

            for (i = 0; i < 8; i++)
            {
                count += (((uint32_t) n_prs[8 * ns + i]) << i) % 30;
            }
            uint32_t f_gh = count;
            return f_gh;
    }
    else
    {
       uint32_t f_gh = 0;
       return f_gh;
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
    uint32_t f_ss, f_gh;
    uint32_t u;
    f_ss = get_f_ss( cell_ID, N_ID_PUCCH, N_ID_PUSCH,PUCCH_ID, PUSCH_ID);
    f_gh = get_f_gh(cell_ID, ns, group_hopping,N_ID_PUCCH, N_ID_PUSCH, PUCCH_ID, PUSCH_ID);
    u = ( f_gh + f_ss ) % 30;
    return u;
}
uint32_t Get_v_value(const uint16_t cell_ID, uint32_t delta_ss, uint32_t M_sc, uint32_t N_sc, uint32_t ns, uint32_t sequence_hopping, uint32_t group_hopping,  uint32_t N_ID_PUCCH, uint32_t N_ID_PUSCH, const uint16_t PUCCH_ID, const uint16_t PUSCH_ID)
{
    /*Sequence hopping only applies for reference-signals of length M_sc_RS ≥ 6*N_sc_RB .For M sc RS < 6 N sc RB ,  v within u is v = 0 .*/
    uint32_t v;
    if ( M_sc >= 6*N_sc)
    {
        uint32_t len = ns;
        int i;
        uint8_t prs[len];
        uint8_t n_prs[len];
        const uint16_t n_ID = Get_cellID( N_ID_PUCCH, N_ID_PUSCH, cell_ID, PUCCH_ID, PUSCH_ID);
        uint32_t c_init = ((n_ID / 30) << 5) + (((n_ID % 30) + delta_ss) % 30);
        calc_prs_c( c_init, len, n_prs); /*generate_pseudo random sequence*/
        if ((sequence_hopping == 1) && (group_hopping == 0))
        { 
 
		     v = n_prs[ns];
                     return v;
		
        }
        else 
        {
            v = 0;
            return v;
        }
    }
    else 
    {
        v = 0;
        return v;
    }

}
// Calculate argument for Qth root Zadoff-Chu sequence according to 3GPP 36.211 5.5.1.1 to find R_ZC
/*
 * Calculate the value of qth root for ZC
 * input- q value , N_zc
 * output root_q for exponential calculation
 */
static void Root_q_arg(float *root_q_arg, uint32_t M_sc, uint32_t u, uint32_t v, uint32_t n_ul_rb, uint32_t N_zc)
{
    int m;
    float n_zc = (float) N_zc;
    float q = get_qvalue(u,v,N_sc);
    for (m = 0; m < N_zc; m++)
    {
        /* argument of x_q(m) according to 3GPP 36.211 5.5.1.1*/
        root_q_arg[m] = (-1 * M_PI * q * m * (m + 1)) / n_zc;
    }
}
void main()
{
uint32_t CP = 1;
uint32_t n_ul_rb = 6;
const uint16_t cell_ID = 100;
uint32_t delta_ss = 0; // delta_ss Δ ss ∈ { 0 , 1 ,..., 29 }
uint32_t group_hopping = 0;//u
uint32_t sequence_hopping = 1;//v
uint32_t N_sc = 12;
uint32_t B = 2;
//uint32_t M_sc = 3*N_sc;
uint32_t ns = 20;
//uint32_t v, u;
uint32_t bw_cfg = 0;//u
uint32_t f_ss, f_gh;
uint32_t N_ID_PUCCH = 0;
uint32_t N_ID_PUSCH = 0;
const uint16_t PUCCH_ID = 499;
const uint16_t PUSCH_ID = 466;
uint32_t M_sc = Get_Msc_values(bw_cfg, B, n_ul_rb, N_sc);
uint32_t N_zc = Get_Nzc( M_sc);
uint32_t fgh = get_f_gh(cell_ID, ns, group_hopping, N_ID_PUCCH, N_ID_PUSCH, PUCCH_ID, PUSCH_ID);
f_ss = get_f_ss(cell_ID, N_ID_PUCCH, N_ID_PUSCH, PUCCH_ID, PUSCH_ID);
printf("fgh = %d f_ss = %d \n",fgh,f_ss);
uint32_t u = Get_u_value(cell_ID, ns, group_hopping, N_ID_PUCCH, N_ID_PUSCH, PUCCH_ID, PUSCH_ID);
uint32_t v =  Get_v_value(cell_ID, delta_ss, M_sc, N_sc, ns, sequence_hopping, group_hopping, N_ID_PUCCH, N_ID_PUSCH, PUCCH_ID, PUSCH_ID);
float root_q_arg[N_zc];
float r_uv;
float q = get_qvalue(u, v, N_zc);
printf("u = %d v = %d q = %f\n",u,v,q);
 Root_q_arg(&root_q_arg, M_sc, u, v, n_ul_rb, N_zc);
/*int i;
for (i = 0; i < N_zc ; i++)
{
r_uv[i] = *root_q_arg[i];
}
//compute_r_uv_arg(r_uv, n_ul_rb, ns, N_sc, N_zc, sequence_hopping, group_hopping, u, v);*/
}
