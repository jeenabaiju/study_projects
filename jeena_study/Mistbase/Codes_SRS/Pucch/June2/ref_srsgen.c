/*
 * ref_srsgen.c
 *
 *  Created on: Jun 3, 2016
 *      Author: mistbasejeena
 */


/*
 *SRS_Gen.c
 *
 *  Created on: May12, 2016
 *      Author: Jeena
 */

//Calculate the value of n_prb
//Calculate the value of n_REs
// Calculate prime
//Calculate largest prime less than M_sc: N_zc values
//Calculate f_ss and f_gh
// Calculate Sequence Hopping
//Calculate Group Hopping

// Generate phi sequence
//Calculate q value
// calculate ZC_rv  sequence
//Calculate cyclic shift, alpha
//Calculate N_ID for SRS, PUCCH and PUSCH;
//Compute and generate pseudo random sequence
//Find // Base sequence =Nsc and 2Nsc
// Base sequence 3Nsc greater
//Calculate SRS seq
// find PUCCH seq
// place SRS in subframe
// TDD configuration





// srs-configIdx - bw-cfg, B, n_srs_cs(Cyclic shift), SRS_UL.sf_idx,
//SRS_UL.delta_ss, SRS_UL.group_hopping,sequence_hopping,uint8_t
//SRS_UL.n_ul_rb=6, HoppingBandwidth, freqdomainposition etc;
/**********************************************************************/
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
#include <ref_table.h>

struct SRS_UL{
           const uint16_t cell_ID;
           uint32_t  B;//B_srs={0,1,2,3} UE specific
           uint32_t bw_cfg; //C_srs takes {0,1,2,3,4,5,6,7} cell specific
           uint32_t srsSubframeConfig;
           uint32_t n_ul_rb;// must be 6 or greater in value
           uint32_t CyclicPrefixLength;
           uint32_t CP;// Normal -0, Extended-1
           uint32_t Duplex_Mode;// FDD-0,TDD-1
           uint32_t HoppingBandwidth; //b_hop={0,1,2,3}
           uint32_t freqDomainPosition;// n_RRC
           uint32_t freqDomainPosition_ap;// n_RRC
           uint32_t nf; //System frame number SFN
           uint32_t Config_idx;// I_srs {0,..... 644}
           uint32_t OffsetIdx;// {0,1}
           uint32_t K_Tc;//Transmission_comb =2 for SRS
           uint32_t transmissionComb;// kbar_tc = {0,1};
           uint32_t Cyclic_shift;// n_srs_cs
           uint32_t Cyclic_shift_ap;// n_srs_cs
           uint32_t N_Tx; //{1,2,4}
           uint32_t ns;//Slot number within a radio frame
           uint32_t N_sp;////No.of DL2UL switchpoints within a radioframe = 5ms
           uint32_t n_hf;// {0,1} first half frame 0 and second half 1 for UpPTS
           uint32_t n_RA;
           uint32_t sf_idx;
           uint32_t sequence_hopping;// enable=1 and disbale=0
           uint32_t group_hopping;// enable=1 and disbale=0
           uint32_t delta_ss; //delta_ss = { 0 , 1 ,..., 29 }
           uint32_t n_ID_PUCCH;// Configured=1 and Not Configured=0
           uint32_t n_ID_PUSCH;// Configured=1 and Not Configured=0
           uint32_t srsMaxUpPTS;// Cell specific
	   uint16_t PUCCH_ID;  // Configured  0 or 1
	   uint16_t PUSCH_ID;  // Configured  0 or 1
       };

/*****************************************************************************/
/*FUNCTIONS FOR SRS GENERATION*/
/******************************************************************************/
/*
 * Calculate the value of number of physical RBs n_prb
 * input- N_sc , .n_ul_rb for normal and extended CP
 * output- N_PRB for extended and normal CP;
 */
/******************************************************************************/
uint32_t N_prb(uint32_t n_ul_rb, uint32_t N_sc)
{
    uint32_t k;
    /* resource elements in frequency domain 0.....n_ul_rb*N_sc*/
    k = n_ul_rb * N_sc;
    uint32_t n_prb = floor(k / N_sc);
    return n_prb;// 6 PRBS in REl-13
}
/******************************************************************************/
/*
 * Calculate the value of number of resource elements
 * input- N_sc , N_symbol_ul for normal and extended CP
 * output -N_PRB for extended and normal CP;
*/
/******************************************************************************/
uint32_t SRS_NRE(struct SRS_UL *srs_ul)
{
    uint32_t n_RE;
    if (srs_ul->CP)
    {
        n_RE = N_sc * N_symbol_ul[0];/*if Normal CP=0  and n_re=(12*7=84)*/
    }
    else
    {
        n_RE = N_sc * N_symbol_ul[1];/*Extended CP=1 and n_re(12*6=72)*/

    }
    return n_RE;
}
/******************************************************************************/
/*
 * To Calculate the power
 * input-  a ,b
 * output -pow
 */
 /*****************************************************************************/

uint32_t power(uint32_t a , uint32_t b)
{
 int i,pow1;
 pow1=1;
 for(i=1;i<=b;i++)
 {
  pow1=pow1*a;
 }
 return pow1;
}
/******************************************************************************/
/*
 * To Calculate the value of M_sc
 * input-  NULRB(uplink BW),Bandwidth B and bw_configuration
 * output -M_sc m_srs
 */
 /*****************************************************************************/
uint32_t srsbwtable_idx(struct SRS_UL *srs_ul)
{
    if (srs_ul->n_ul_rb <= 40)
    {
        return 0;
    }
    else if (srs_ul->n_ul_rb <= 60)
    {
        return 1;
    }
    else if (srs_ul->n_ul_rb <= 80)
    {
        return 2;
    }
    else
    {
        return 3;
    }
}

uint32_t  Get_Msc_values(struct SRS_UL *srs_ul, uint32_t N_sc)
{
   uint32_t M_sc;
   /* According to 3GPP 36.211version 13 5.5.3.2 */
   M_sc = m_srs_b[srsbwtable_idx(srs_ul)][srs_ul->B][srs_ul->bw_cfg] * N_sc / 2;
   return M_sc;
}
/*****************************************************************************/
/*/*
 * Calculate the value of N_zc
 * input- bw_cfg,Uplink BW(n_ul_rb),N_sc, B_srs.
 * output -N_zc largest prime number less than M_sc
 */
 /*****************************************************************************/
uint32_t Get_Nzc(uint32_t M_sc)
{

 /* get largest prime no N_zc<M_sc */
    int i,j;
    uint32_t N_zc;
    N_zc = 0;
    for(i = 1; i < M_sc; i++)
    {
        for(j = 2; j < i; j++) if(i % j == 0) break;
        if(j == i)
            N_zc = j;/* find prime numbers*/
    }
    return N_zc ;
}
/******************************************************************************/
/*
* Calculate N_CellID
* input - SRS_UL.cell_ID, Configurations for n_ID_PUCCH,n_ID_PUSCH
* output- N_ID
*/
/******************************************************************************/
const uint16_t Get_cellID(struct SRS_UL *srs_ul)
{
    if (srs_ul->n_ID_PUCCH == 1)
    {
        const uint16_t N_ID= srs_ul->PUCCH_ID;
        return  N_ID;
    }
    else if (srs_ul->n_ID_PUSCH == 1)
    {
        const uint16_t N_ID= srs_ul->PUSCH_ID;
        return  N_ID;
    }
    else
    {
        const uint16_t N_ID= srs_ul->cell_ID;// for SRS also//
        return  N_ID;
    }
}
/****************************************************************************/
/* //Compute and generate pseudo random sequence
* input - c_init, len
* output- n_prs values

 * written by Henrik
 * */
/******************************************************************************/
void calc_prs_c(const uint32_t c_init, const uint32_t len, uint8_t* n_prs)
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

	// Calculate n_prs
	for(i = 0; i < len; i++){

		next_x1 = ((x1 >> 3) ^ x1) & 0x1;
		next_x2 = ((x2 >> 3) ^ (x2 >> 2) ^ (x2 >> 1) ^ x2) & 0x1;

		*(n_prs++) = next_x1 ^ next_x2;

		x1 = (x1 >> 1) | (next_x1 << 30);
		x2 = (x2 >> 1) | (next_x2 << 30);
	}

}
/******************************************************************************/
/*
 * Calculate Sequence Hopping(v)
 * input- CellID,SRTue 24 May 2016 10:18:08 AM CEST S_UL.delta_ss,
 * output-Calculate v
 */
/******************************************************************************/
uint32_t Get_v_value(struct SRS_UL *srs_ul, uint32_t M_sc, uint32_t N_sc)
{
    /*Sequence hopping only applies for reference-signals of
      length M_sc_RS < 6*N_sc_RB .
    * For M sc RS < 6 N sc RB ,  v within u is v = 0 .*/
   uint32_t v;
   v = 0;
   if ( M_sc >= 6*N_sc)
   {
      uint32_t len = srs_ul->ns+1;
      uint8_t n_prs[len];
      const uint16_t n_ID = Get_cellID( srs_ul);
      uint32_t temp;
      temp = ((n_ID / 30) << 5) + (n_ID % 30);
      uint32_t  c_init = (temp + srs_ul->delta_ss) % 30;
      calc_prs_c( c_init, len, n_prs); /*generate_pseudo random sequence*/
      if ((srs_ul->sequence_hopping == 1) && (srs_ul->group_hopping == 0))
      {
	  v = n_prs[len];
      }
   }
   else
   {
        v = 0;
   }
return v;
}
/*****************************************************************************/
/* Calculate f_gh
 * input- number of slots, group hopping -enabled/disabled, cellID
 * output - f_gh value
 */
/******************************************************************************/
uint32_t get_f_gh(struct SRS_UL *srs_ul)
{
    uint32_t len = ( 8 * srs_ul->ns + 8 );// Maximum length
    uint8_t n_prs[len];
    uint32_t count = 0;
    uint32_t f_gh;
    const uint16_t n_ID = Get_cellID( srs_ul);
    uint32_t c_init = floor( n_ID / 30 );
    /*generate_pseudo random sequence.
      Computes n_prs values as defined in 5.5.2.1.1 of 36.211 */
    calc_prs_c(c_init, len, n_prs);
    /** Computes sequence-group pattern f_gh according to 5.5.1.3 of 36.211 */
    f_gh = 0;
    int i;
	if(srs_ul->group_hopping == 1)
        {
           for (i = 0; i < 8; i++)
           {
               count += (((uint32_t) n_prs[8 * srs_ul->ns + i]) << i);
           }
           f_gh = count % 30;
           return f_gh;

        /*             i=7
       summation of [c( 8*srs_ul->ns+ i )⋅ 2^i ]mod 30,if  group hopping is
       enabled           i=0                                                 */
        }
        else
        {
          f_gh = 0;
          return f_gh;
        }
}
/******************************************************************************/
/*
 * Calculate f_ss
 * input-  cellID
 * output - f_ss value
 */
 /*****************************************************************************/
uint32_t get_f_ss(struct SRS_UL *srs_ul)
{
    uint32_t f_ss;
    const uint16_t n_ID = Get_cellID( srs_ul);
    f_ss = n_ID % 30;
    return f_ss;
}
/******************************************************************************/
/*
 * Calculate Group Hopping(u)
 * u=(f_gh(ns)+f_ss)%30;
 * input- number of slots, group hopping -enabled/disabled, cellID
 * output - u value
 */
/******************************************************************************/
static uint32_t Get_u_value(struct SRS_UL *srs_ul)
{
    uint32_t f_ss, f_gh;
    uint32_t u;
    f_ss = get_f_ss( srs_ul);
    f_gh = get_f_gh(srs_ul);
    u = ( f_gh + f_ss ) % 30;
    return u;
}
/******************************************************************************/
/*
 * Calculate Alpha value for SRS according to 5.5.3.1 of 36.211
 * input- N_Tx, K_Tc, n_srs_cs from higher layers
 * output - alpha
 */
 /*****************************************************************************/
static float alpha_p(struct SRS_UL *srs_ul)
{
    uint32_t n_srs_max , p;
    uint32_t n_srs_p;
    float alpha;
    if ( srs_ul->K_Tc == 2 )
    {
        n_srs_max = 8;/* for all SRS signals*/
    }
    else
    {
        n_srs_max = 12;
    }
    p = srs_ul->N_Tx - 1;
    n_srs_p = (srs_ul->Cyclic_shift + ((n_srs_max * p) / srs_ul->N_Tx)) % n_srs_max;
    alpha = 2 * M_PI * n_srs_p / n_srs_max;
    return alpha;
}
/******************************************************************************/
/*
 * Calculate q value
 * input- u,v
 * output -q value to find Zadoff chu seq
 */
 /*****************************************************************************/
 static float get_qvalue(uint32_t u, uint32_t v, uint32_t N_zc)
{
    float q;
    float q_bar;
    uint32_t temp, temp1;
    float n_zc = (float) N_zc;
    q_bar = n_zc * (u + 1) / 31;
    uint32_t pow1;
    temp = q_bar + 0.5;
    temp1 = 2 * q_bar;
    pow1 = power(-1, temp1);
    q = temp + (v * pow1);
    return q;
}
/******************************************************************************/
/*Calculate argument for Qth root Zadoff-Chu sequence
  according to 3GPP 36.211 5.5.1.1 to find R_ZC
 * Calculate the value of qth root for ZC
 * input- q value , N_zc
 * output root_q for exponential calculation
 */
 /*****************************************************************************/
static void Rootq_arg(float *rootq_arg, uint32_t u, uint32_t v, uint32_t N_zc)
{
    int m;
    float n_zc = (float) N_zc;
    float q = get_qvalue(u,v,N_zc);
    for (m = 0; m < N_zc; m++)
    {
        /* argument of x_q(m) according to 3GPP 36.211 5.5.1.1*/
        rootq_arg[m] = (-1.0 * M_PI * q * m * (m + 1)) / n_zc;
		printf(" M_PI= %f, n_zc = %f, q = %f, m= %d \n",M_PI,n_zc,q,m);
		printf(" rootq_arg[%d] = %f\n",m,rootq_arg[m]);
    }
}
/******************************************************************************/
/*
 * Calculate the  argument value for exponential calculation of
    Base sequence for M_sc=N_sc
 * input-  u ,N_sc=no.ofsubcarriers
 * output -Base Sequence value for M_sc=N_sc
 */
 /*****************************************************************************/
static void Seq_Msc12_Exp(float *Seq_Nsc_exp, uint8_t u, uint32_t N_sc)
{
    int i;
    for (i = 0; i < N_sc; i++)
    {
        Seq_Nsc_exp[i] = Phi_12[u][i] * M_PI / 4;
    }
}
/******************************************************************************/
 /*
 * Calculate the argument value for exponential calculation of
   Base sequence for M_sc=2*N_sc
 * input- u ,N_sc=no.ofsubcarriers
 * output -Base Sequence value for M_sc=2*N_sc
 */
 /*****************************************************************************/
static void Seq_Msc24_Exp(float *Seq_2Nsc_exp, uint8_t u, uint32_t N_sc)
{
    int i;
    for (i = 0 ;i < 2*N_sc; i++)
    {
        Seq_2Nsc_exp[i] = Phi_24[u][i] * M_PI / 4;
		printf(" Seq_2Nsc_exp[%d] = %f\n",i,Seq_2Nsc_exp[i]);

    }
}
// appending the sequences to fill M_sc
static void r_uv_mprb(float *r_uv, float *r_uv_xq,uint32_t M_sc, uint32_t N_zc)
{
    int n;
    for (n = 0; n < M_sc; n++)
    {
        r_uv[n] = r_uv_xq[n % N_zc];
    }
}
/******************************************************************************/
// Calculate
/*
 * Calculate the value of arguments for arguments foe exponential in r_uv(n)
 * input- N_PRB,Grp no, Seq No, seq hopping group hopping, u, v
 * output -argument values to find sequence */
/* Computes argument of r_u_v signal */
/*****************************************************************************/
static void r_uv_arg(float *r_uv, uint32_t M_sc, uint32_t N_zc, uint32_t u, uint32_t v)
{
    // get u

    if (M_sc == N_sc)
    {
        Seq_Msc12_Exp(r_uv, u, N_sc);
    }
    else if (M_sc == 2*N_sc)
    {
        Seq_Msc24_Exp(r_uv, u,N_sc);
    }
    else
    {
        float r_uv_temp[N_zc];
        //calculate sequence number v
        Rootq_arg(r_uv_temp, u, v, N_zc);
        r_uv_mprb(r_uv, r_uv_temp, M_sc, N_zc);
    }
}
/******************************************************************************/
/* Generate SRS signal as defined in Section 5.5.3.1
 * input - N_sc, n_ul_rb, sf_idx, group hopping, seq hopping,
    cell, PUCCH,PUSH,n_ID_PUCCH,n_ID_PUSCH, Ntx etc
 * output- SRS sequence
 * */
/******************************************************************************/
void srs_gen(double complex *r_srs,struct SRS_UL *srs_ul)
{
    int n ;
    uint32_t M_sc = Get_Msc_values( srs_ul, N_sc);
    uint32_t N_zc = Get_Nzc( M_sc);
    float r_uv[M_sc];
    uint32_t u = Get_u_value( srs_ul);
    uint32_t v = Get_v_value( srs_ul, M_sc, N_sc);
    float alpha = alpha_p(srs_ul);//n_srs
    r_uv_arg(r_uv, M_sc, N_zc, u, v);
    // Do complex exponential and adjust amplitude
    for ( n = 0; n < M_sc; n++)
    {
            r_srs[n] = cexpf( I * ( r_uv[n] + ( alpha * n ) ) );
    }
}
