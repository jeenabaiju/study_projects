
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
/* //Table 5.2.3-1: Resource block parameters according to 3GPP 36.211 5.2.3*/
uint32_t N_sc = 12;
uint32_t N_symbol_ul[2] = {7,6};
/* Table 5.5.3.3-1: Frame structure type 1 sounding reference signal subframe configuration. */
uint32_t T_sfc[15] = {1, 2, 2, 5, 5, 5, 5, 5, 5, 10, 10, 10, 10, 10, 10};   /*Configuration Period*/
uint32_t Delta_sfc1[7] = {0, 0, 1, 0, 1, 2, 3}; /*Transmission offset*/
uint32_t Delta_sfc2[4] = {0, 1, 2, 3};

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

/*Generating  φ (n ) for M_sc_RS = N_sc_RB . Table 5.5.1.2-1 */
int  Phi_M_sc_12[30][12]={{-1,1,3,-3,3,3,1,1,3,1,-3,3},
{1,1,3,3,3,-1,1,-3,-3,1,-3,3},
{1,1,-3,-3,-3,-1,-3,-3,1,-3,1,-1},
{-1,1,1,1,1,-1,-3,-3,1,-3,3,-1},
{-1,3,1,-1,1,-1,-3,-1,1,-1,1,3},
{1,-3,3,-1,-1,1,1,-1,-1,3,-3,1},
{-1,3,-3,-3,-3,3,1,-1,3,3,-3,1},
{-3,-1,-1,-1,1,-3,3,-1,1,-3,3,1},
{1,-3,3,1,-1,-1,-1,1,1,3,-1,1},
{1,-3,-1,3,3,-1,-3,1,1,1,1,1},
{-1,3,-1,1,1,-3,-3,-1,-3,-3,3,-1},
{3,1,-1,-1,3,3,-3,1,3,1,3,3},
{1,-3,1,1,-3,1,1,1,-3,-3,-3,1},
{3,3,-3,3,-3,1,1,3,-1,-3,3,3},
{-3,1,-1,-3,-1,3,1,3,3,3,-1,1},
{3,-1,1,-3,-1,-1,1,1,3,1,-1,-3},
{1,3,1,-1,1,3,3,3,-1,-1,3,-1},
{-3,1,1,3,-3,3,-3,-3,3,1,3,-1},
{-3,3,1,1,-3,1,-3,-3,-1,-1,1,-3},
{-1,3,1,3,1,-1,-1,3,-3,-1,-3,-1},
{-1,-3,1,1,1,1,3,1,-1,1,-3,-1},
{-1,3,-1,1,-3,-3,-3,-3,-3,1,-1,-3},
{1,1,-3,-3,-3,-3,-1,3,-3,1,-3,3},
{1,1,-1,-3,-1,-3,1,-1,1,3,-1,1},
{1,1,3,1,3,3,-1,1,-1,-3,-3,1},
{1,-3,3,3,1,3,3,1,-3,-1,-1,3},
{1,3,-3,-3,3,-3,1,-1,-1,3,-1,-3},
 {-3,-1,-3,-1,-3,3,1,-1,1,3,-3,-3},
 {-1,3,-3,3,-1,3,3,-3,3,3,-1,-1},
 {3,-3,-3,-1,-1,-3,-1,3,-3,3,1,-1}};
/*Generating  φ (n ) for M sc RS =2* N sc RB . Table 5.5.1.2-2 */
int Phi_M_sc_24[30][24]={{-1,3,1,-3,3,-1,1,3,-3,3,1,3,-3,3,1,1,-1,1,3,-3,3,-3,-1,-3},
{-3,3,-3,-3,-3,1,-3,-3,3,-1,1,1,1,3,1,-1,3,-3,-3,1,3,1,1,-3},
{3,-1,3,3,1,1,-3,3,3,3,3,1,-1,3,-1,1,1,-1,-3,-1,-1,1,3,3},
{-1,-3,1,1,3,-3,1,1,-3,-1,-1,1,3,1,3,1,-1,3,1,1,-3,-1,-3,-1},
{-1,-1,-1,-3,-3,-1,1,1,3,3,-1,3,-1,1,-1,-3,1,-1,-3,-3,1,-3,-1,-1},
{-3,1,1,3,-1,1,3,1,-3,1,-3,1,1,-1,-1,3,-1,-3,3,-3,-3,-3,1,1},
{1,1,-1,-1,3,-3,-3,3,-3,1,-1,-1,1,-1,1,1,-1,-3,-1,1,-1,3,-1,-3},
{-3,3,3,-1,-1,-3,-1,3,1,3,1,3,1,1,-1,3,1,-1,1,3,-3,-1,-1,1},
{-3,1,3,-3,1,-1,-3,3,-3,3,-1,-1,-1,-1,1,-3,-3,-3,1,-3,-3,-3,1,-3},
{1,1,-3,3,3,-1,-3,-1,3,-3,3,3,3,-1,1,1,-3,1,-1,1,1,-3,1,1},
{-1,1,-3,-3,3,-1,3,-1,-1,-3,-3,-3,-1,-3,-3,1,-1,1,3,3,-1,1,-1,3},
{1,3,3,-3,-3,1,3,1,-1,-3,-3,-3,3,3,-3,3,3,-1,-3,3,-1,1,-3,1},
{1,3,3,1,1,1,-1,-1,1,-3,3,-1,1,1,-3,3,3,-1,-3,3,-3,-1,-3,-1},
{3,-1,-1,-1,-1,-3,-1,3,3,1,-1,1,3,3,3,-1,1,1,-3,1,3,-1,-3,3},
{-3,-3, 3, 1, 3, 1, -3, 3, 1, 3, 1, 1, 3, 3, -1, -1, -3, 1, -3, -1, 3, 1, 1, 3},
{-1, -1, 1, -3, 1, 3, -3,1 ,-1 ,-3 ,-1 ,3 ,1 ,3 ,1 ,-1 ,-3 ,-3 ,-1 ,-1 ,-3 ,-3 ,-3, -1},
{-1 ,-3 ,3 ,-1 ,-1 ,-1 ,-1, 1, 1, -3, 3, 1 ,3, 3, 1 ,-1 ,1 ,-3 ,1 ,-3, 1, 1 ,-3, -1},
{1, 3 ,-1 ,3 ,3 ,-1, -3, 1 ,-1 ,-3 ,3 ,3, 3, -1, 1, 1, 3, -1, -3, -1, 3 ,-1, -1, -1},
{1 ,1 ,1 ,1 ,1 ,-1 ,3 ,-1 ,-3 ,1 ,1 ,3 ,-3, 1 ,-3, -1 ,1 ,1 ,-3 ,-3 ,3 ,1 ,1, -3},
{1 ,3, 3, 1 ,-1 ,-3 ,3 ,-1 ,3 ,3 ,3 ,-3, 1, -1, 1, -1, -3, -1, 1 ,3 ,-1 ,3 ,-3 ,-3},
{-1, -3 ,3 ,-3 ,-3 ,-3 ,-1 ,-1 ,-3 ,-1 ,-3 ,3 ,1 ,3 ,-3 ,-1 ,3 ,-1 ,1 ,-1 ,3 ,-3 ,1 ,-1},
{-3 ,-3, 1 ,1 ,-1, 1 ,-1, 1, -1 ,3, 1, -3 ,-1 ,1 ,-1, 1 ,-1 ,-1 ,3 ,3 ,-3 ,-1,1, -3},
{-3 ,-1 ,-3 ,3 ,1, -1, -3, -1, -3, -3 ,3 ,-3, 3, -3, -1, 1, 3, 1, -3, 1 ,3 ,3, -1, -3},
{-1 ,-1 ,-1 ,-1 ,3 ,3 ,3 ,1 ,3 ,3 ,-3, 1, 3 ,-1 ,3 ,-1 ,3 ,3 ,-3 ,3 ,1 ,-1, 3 ,3},
{1 ,-1 ,3, 3, -1, -3, 3,-3, -1, -1, 3 ,-1, 3, -1, -1, 1, 1, 1, 1, -1, -1, -3, -1, 3},
{1, -1, 1 ,-1 ,3 ,-1 ,3 ,1 ,1 ,-1 ,-1 ,-3 ,1 ,1 ,-3 ,1 ,3 ,-3 ,1 ,1 ,-3 ,-3 ,-1 ,-1},
{-3 ,-1, 1, 3, 1, 1, -3, -1, -1, -3, 3 ,-3, 3, 1 ,-3, 3, -3, 1, -1, 1, -3, 1 ,1, 1},
{-1, -3 ,3, 3 ,1 ,1 ,3 ,-1 ,-3 ,-1 ,-1, -1, 3, 1, -3, -3, -1, 3, -3, -1, -3, -1, -3, -1},
{-1, -3, -1 ,-1 ,1, -3, -1, -1, 1 ,-1, -3, 1, 1, -3, 1, -3, -3, 3, 1, 1, -1, 3, -1, -1},
{1, 1 ,-1 ,-1 ,-3 ,-1, 3, -1, 3, -1, 1, 3, 1, -1, 3, 1, 3, -3, -3, 1, -1, -1, 1, 3}};
/********************************************************************************************************************/

/***********************************************************************************************************************/
/*
/*
 * To Calculate the value of M_sc
 * input-  NULRB(uplink BW),Bandwidth B and bw_configuration
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

uint32_t  Get_Msc_values(uint32_t bw_cfg, uint32_t B, uint32_t n_ul_rb, uint32_t N_sc)
{
   uint32_t M_sc;
   M_sc = m_srs_b[srsbwtable_idx(n_ul_rb)][B][bw_cfg] * N_sc / 2;/* According to 3GPP 36.211 5.5.3.2*/
   return M_sc;
}
/***********************************************************************************************************************/
/*/*
 * Calculate the value of N_zc
 * input- bw_cfg,Uplink BW(n_ul_rb),N_sc, B_srs.
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


void calc_prs_c(const uint32_t c_init, const uint32_t len, uint8_t* prs_c)
{
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
uint32_t Get_v_value(const uint16_t cell_ID, uint32_t delta_ss, uint32_t M_sc, uint32_t N_sc, uint32_t ns, uint32_t sequence_hopping, uint32_t group_hopping,uint32_t N_ID_PUCCH, uint32_t N_ID_PUSCH, const uint16_t PUCCH_ID, const uint16_t  PUSCH_ID)
{
        /*Sequence hopping only applies for reference-signals of length M_sc_RS ≥ 6*N_sc_RB .For M sc RS < 6 N sc RB ,  v within u is v = 0 .*/
        uint32_t v;
        if ( M_sc >= 6*N_sc)
        {
            uint32_t len = ns;
            int i;
            uint8_t n_prs[len];
	        const uint16_t n_ID = Get_cellID( N_ID_PUCCH, N_ID_PUSCH, cell_ID,PUCCH_ID, PUSCH_ID);

            uint32_t  c_init = ((n_ID / 30) << 5) + (((n_ID  % 30) + delta_ss) % 30);
            calc_prs_c( c_init, len, n_prs); /*generate_pseudo random sequence*/
            if ((sequence_hopping == 1) && (group_hopping == 0))
            {
		       v = n_prs[ns];
			   printf( "v= %d \n", v);
		       return v;
            }
        }
        else
        {
            v = 0;
            return v;
        } 
}
uint32_t get_f_gh(const uint16_t cell_ID, uint32_t ns, uint32_t group_hopping, uint32_t N_ID_PUCCH, uint32_t N_ID_PUSCH, const uint16_t PUCCH_ID, const uint16_t PUSCH_ID)
{
    uint32_t len = ( 8 * ns + 7 );// Maximum length
    uint8_t n_prs[len];
    uint32_t count = 0;
    uint32_t f_gh;
	
    const uint16_t n_ID = Get_cellID( N_ID_PUCCH, N_ID_PUSCH,cell_ID, PUCCH_ID, PUSCH_ID);
    uint32_t c_init = floor( n_ID / 30 );
    /*generate_pseudo random sequence /** Computes n_prs values as defined in 5.5.2.1.1 of 36.211 */
    calc_prs_c(c_init, len, n_prs);
    /** Computes sequence-group pattern f_gh according to 5.5.1.3 of 36.211 */
    f_gh = 0;
        int i;
	
    if(group_hopping == 1)
    {
        for (i = 0; i < 8; i++)
        {
            count += (((uint32_t) n_prs[8 * ns + i]) << i) % 30;
        }
        uint32_t f_gh = count;
        return f_gh;

        /*             i=7
        summation of [c( 8*ns+ i )⋅ 2^i ]mod 30 if group hopping is enabled
                       i=0                                                      */
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
uint32_t Get_u_value(const uint16_t cell_ID, uint32_t ns, uint32_t  group_hopping, uint32_t  N_ID_PUCCH, uint32_t N_ID_PUSCH, const uint16_t PUCCH_ID, const uint16_t PUSCH_ID)
{
    uint32_t f_ss, f_gh;
    uint32_t u;
    f_ss = get_f_ss( cell_ID, N_ID_PUCCH, N_ID_PUSCH,PUCCH_ID, PUSCH_ID);
    if(group_hopping == 1)
    {
        f_gh = get_f_gh(cell_ID, ns,group_hopping,N_ID_PUCCH, N_ID_PUSCH, PUCCH_ID, PUSCH_ID);
        u = ( f_gh + f_ss ) % 30;
        return u;
    }
    else
    {
        const uint16_t n_ID = Get_cellID(N_ID_PUCCH, N_ID_PUSCH, cell_ID, PUCCH_ID, PUSCH_ID);
        u = n_ID % 30;
        return u;
    }
}

/*
 * Calculate Alpha value for SRS according to 5.5.3.1 of 36.211
 * input- N_Tx, K_Tc, n_srs_cs from higher layers
 * output - alpha
 */
static float alpha_p(uint32_t N_Tx, uint32_t Cyclic_shift, uint32_t K_Tc)
{
    uint32_t n_srs_max , p;
    uint32_t n_srs_p;
    float alpha;
    if ( K_Tc == 2 )
    {
        n_srs_max = 8;/* for all SRS signals*/
    }
    else
    {
        n_srs_max = 12;
    }
    p = N_Tx - 1;
    n_srs_p = ( Cyclic_shift + ( ( n_srs_max * p) / N_Tx ) ) % n_srs_max;
    alpha = 2 * M_PI * n_srs_p / n_srs_max;
    return alpha;
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
/*************************************************************************************************/
// Calculate argument for Qth root Zadoff-Chu sequence according to 3GPP 36.211 5.5.1.1 to find R_ZC
/*
 * Calculate the value of qth root for ZC
 * input- q value , N_zc
 * output root_q for exponential calculation
 */
static void Root_q_arg(float *root_q_arg,uint32_t u, uint32_t v, uint32_t n_ul_rb, uint32_t N_zc)
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
/***********************************************************************************************/
/*
 * Calculate the  argument value for exponential calculation of Base sequence for M_sc=N_sc
 * input-  u ,N_sc=no.ofsubcarriers
 * output -Base Sequence value for M_sc=N_sc
 */
static void Seq_Msc12_Exp(float *Seq_Nsc_exp, uint8_t u, uint32_t N_sc)
{
    int i;
    for (i = 0; i < N_sc; i++)
    {
        Seq_Nsc_exp[i] = Phi_M_sc_12[u][i] * M_PI / 4;
    }
}
/**************************************************************************************************/
 /*
 * Calculate the argument value for exponential calculation of Base sequence for M_sc=2*N_sc
 * input- u ,N_sc=no.ofsubcarriers
 * output -Base Sequence value for M_sc=2*N_sc
 */
static void Seq_Msc24_Exp(float *Seq_2Nsc_exp, uint8_t u, uint32_t N_sc)
{

    int i;
    for (i = 0 ;i < 2*N_sc; i++)
    {
        Seq_2Nsc_exp[i] = Phi_M_sc_24[u][i] * M_PI / 4;

    }
}
static void r_uv_mprb(float *r_uv, float *r_uv_xq,uint32_t M_sc, uint32_t N_zc)
{
    int n;
    for (n = 0; n < M_sc; n++)
    {
        r_uv[n] = r_uv_xq[n % N_zc];
    }
}
/***********************************************************************************************/
// Calculate arguments foe exponential in r_uv(n)
/*
 * Calculate the value of arguments for N_PRB=1,2,<=3
 * input- N_PRB,Grp no, Seq No, other values needed
 * output -argument values to find sequence */
/* Computes argument of r_u_v signal */
static void compute_r_uv_arg(float *r_uv, uint32_t n_ul_rb, uint32_t ns, uint32_t M_sc, uint32_t N_zc, uint32_t sequence_hopping, uint32_t group_hopping,  uint32_t u, uint32_t v)
{
    // get u

    if (n_ul_rb == 1)
    {
        Seq_Msc12_Exp(r_uv, u, N_sc);
    }
    else if (n_ul_rb == 2)
    {
        Seq_Msc24_Exp(r_uv, u,N_sc);
    }
    else
    {
        float r_uv_temp[N_zc];
        //float r_uv[M_sc];
        float r_uv_xq[M_sc];
		int i;
        //calculate sequence number v
        Root_q_arg(r_uv_temp, n_ul_rb, u, v, N_zc);  
		printf(" r_uv_temp[%d] = %f \n ",N_zc-2, r_uv_temp[N_zc-2]);
        r_uv_mprb(r_uv, r_uv_temp, M_sc, N_zc);
    }
}

/* Generate SRS signal as defined in Section 5.5.3.1 */
/*int srs_gen(struct srs_sequence *seq, uint32_t N_sc,uint32_t n_ul_rb uint32_t N_zc, uint32_t sf_idx, uint32_t group_hopping, uint32_t sequence_hopping, const uint16_t cell_ID, uint32_t delta_ss, uint32_t ns, uint32_t N_Tx, uint32_t Cyclic_shift, uint32_t K_Tc)*/
int srs_gen(double complex *r_srs, uint32_t N_sc, uint32_t n_ul_rb,uint32_t sf_idx, uint32_t group_hopping, uint32_t sequence_hopping, const uint16_t cell_ID, uint32_t delta_ss, uint32_t ns, uint32_t N_Tx, uint32_t Cyclic_shift, uint32_t K_Tc, uint32_t bw_cfg,uint32_t B,uint32_t N_ID_PUCCH, uint32_t N_ID_PUSCH, const uint16_t PUCCH_ID, const uint16_t PUSCH_ID)
{
    int n ;
    uint32_t nslot;
    float Seq_Nsc_exp[N_sc];
    float Seq_2Nsc_exp[N_sc];
    uint32_t M_sc = Get_Msc_values( bw_cfg, B, n_ul_rb, N_sc);
    uint32_t N_zc = Get_Nzc( M_sc);
    float r_uv[M_sc];
    float root_q_arg[N_zc];
    float r_uv_xq[M_sc];
    uint32_t f_gh = get_f_gh(cell_ID, ns, group_hopping, N_ID_PUCCH, N_ID_PUSCH, PUCCH_ID, PUSCH_ID);
    uint32_t f_ss = get_f_ss(cell_ID, N_ID_PUCCH, N_ID_PUSCH, PUCCH_ID, PUSCH_ID);
    uint32_t u = Get_u_value( cell_ID, ns, group_hopping, N_ID_PUCCH, N_ID_PUSCH, PUCCH_ID, PUSCH_ID);
    uint32_t v = Get_v_value( cell_ID, delta_ss, M_sc, N_sc, ns, sequence_hopping, group_hopping, N_ID_PUCCH, N_ID_PUSCH, PUCCH_ID, PUSCH_ID);
	printf(" uvalue = %d \n", u);
    float q = get_qvalue(u, v, N_zc);
	printf(" q value = %f \n", q);
    for (nslot = 2*sf_idx; nslot < 2*(sf_idx + 1); nslot++)
    {
        float alpha = alpha_p(N_Tx, Cyclic_shift, K_Tc);//n_srs
        compute_r_uv_arg(r_uv, n_ul_rb,ns,M_sc, N_zc,sequence_hopping,group_hopping,u,v);
        // Do complex exponential and adjust amplitude
        for ( n = 0; n < M_sc; n++)
        {
            //seq->r_srs[nslot][n] = cexpf( I * ( r_uv[n] + ( alpha * n ) ) );
            r_srs[(nslot % 2)*M_sc+n] = cexpf( I * ( r_uv[n] + ( alpha * n ) ) );
        }
    }
}

void main()
{

const uint16_t cell_ID = 100;
uint32_t delta_ss = 0; // delta_ss Δ ss ∈ { 0 , 1 ,..., 29 }
uint32_t group_hopping = 0;
uint32_t sequence_hopping = 1;
uint32_t N_sc = 12;
uint32_t ns = 20 ;
uint32_t f_gh;
uint32_t N_Tx = 1;
uint32_t sf_idx = 0;
uint32_t K_Tc = 2;
uint32_t bw_cfg = 0;
uint32_t B = 1;
uint32_t Cyclic_shift = 0;
uint32_t N_ID_PUCCH = 0;
uint32_t N_ID_PUSCH = 0;
const uint16_t PUCCH_ID = 499;
const uint16_t PUSCH_ID = 466;
int n;
uint32_t M_sc;
uint32_t n_ul_rb = 6;
M_sc = Get_Msc_values( bw_cfg, B, n_ul_rb, N_sc);
printf("M_sc = %d\n", M_sc);
uint32_t v = Get_v_value( cell_ID, delta_ss, M_sc, N_sc, ns, sequence_hopping, group_hopping, N_ID_PUCCH, N_ID_PUSCH, PUCCH_ID, PUSCH_ID);uint32_t u = Get_u_value( cell_ID, ns, group_hopping, N_ID_PUCCH, N_ID_PUSCH, PUCCH_ID, PUSCH_ID);
printf("Vrx= %d \n ",v);
double complex r_srs[N_sc * n_ul_rb*2];
srs_gen(r_srs, N_sc, n_ul_rb, sf_idx, group_hopping, sequence_hopping, cell_ID, delta_ss,ns, N_Tx, Cyclic_shift, K_Tc, bw_cfg, B, N_ID_PUCCH, N_ID_PUSCH,  PUCCH_ID, PUSCH_ID);
for ( n = 0; n < N_sc * n_ul_rb; n++)
        {
            //seq->r_srs[nslot][n] = cexpf( I * ( r_uv[n] + ( alpha * n ) ) );
            printf("r_srs = %f + i %f \n ",creal (r_srs[n]),cimag( r_srs[n]));
			//printf("Size %d",printf("r_srs = %f + i %f \n ",creal (r_srs[n]),cimag( r_srs[n]));)
        }

}
