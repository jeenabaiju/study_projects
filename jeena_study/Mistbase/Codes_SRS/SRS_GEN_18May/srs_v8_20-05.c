
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





// srs-configIdx - bw-cfg, B, n_srs_cs(Cyclic shift), sf_idx, delta_ss, group_hopping,sequence_hopping,uint8_t n_ul_rb=6, HoppingBandwidth, freqdomainposition etc;
/************************************************************************************************************************/
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

/************************************************************************************************************************/
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

/* Same tables for Nb */
uint32_t Nb[4][4][8] = {{
                        {1, 1, 1, 1, 1, 1, 1, 1},
                        {3, 2, 6, 5, 4, 3, 2, 1},
                        {3, 2, 1, 1, 1, 1, 1, 1},
                        {1, 2, 1, 1, 1, 1, 1, 1}},
                        {
                        {1, 1, 1, 1, 1, 1, 1, 1},
                        {2, 3, 2, 3, 2, 6, 5, 4},
                        {2, 2, 5, 3, 2, 1, 1, 1},
                        {3, 2, 1, 1, 2, 1, 1, 1}},
                        {
                        {1, 1, 1, 1, 1, 1, 1, 1},
                        {3, 2, 3, 2, 3, 2, 3, 2},
                        {2, 2, 5, 2, 2, 5, 3, 2},
                        {3, 4, 1, 3, 2, 1, 1, 2}},
                        {
                        {1, 1, 1, 1, 1, 1, 1, 1},
                        {2, 3, 2, 3, 2, 3, 2, 3},
                        {2, 2, 2, 2, 2, 5, 2, 2},
                        {6, 4, 5, 3, 4, 1, 3, 2}}};


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


/************************************************************************************************************************/
//CHECKED
/************************************************************************************************************************/

//Calculate the value of n_prb
/*
 * Calculate the value of number of physical RBs n_prb
 * input- N_sc , n_ul_rb for normal and extended CP
 * output- N_PRB for extended and normal CP;
 */
uint32_t N_prb(uint32_t SRS_UL.n_ul_rb, uint32_t N_sc)
{
    uint32_t i;
    uint32_t k;
    /* resource elements in frequency domain 0.....n_ul_rb*N_sc*/
    k = SRS_UL.n_ul_rb * N_sc;
    uint32_t n_prb = floor(k / N_sc);
    return n_prb;// 6 PRBS in REl-13
}
//checked
/*****************************************************************************************************************/

//Calculate the value of n_REs
/*
 * Calculate the value of number of resource elements
 * input- N_sc , N_symbol_ul for normal and extended CP
 * output -N_PRB for extended and normal CP;
*/
uint32_t SRS_NRE(SRS_UL.CP)
{
    uint32_t n_RE;
    if (SRS_UL.CP)
    {
        n_RE = N_sc * N_symbol_ul[0];/*if Normal CP=0  and n_re=(12*7=84)*/
        return n_RE;
    }
    else
    {
        n_RE = N_sc * N_symbol_ul[1];/*Extended CP=1 and n_re(12*6=72)*/
        return n_RE;
    }
    
}
//checked
/***********************************************************************************************************************/
/*
/*
 * To Calculate the value of M_sc
 * input-  NULRB(uplink BW),Bandwidth B and bw_configuration
 * output -M_sc m_srs
 */
uint32_t srsbwtable_idx(uint32_t SRS_UL.n_ul_rb)
{
    if (SRS_UL.n_ul_rb <= 40)
    {
        return 0;
    }
    else if (SRS_UL.n_ul_rb <= 60)
    {
        return 1;
    }
    else if (SRS_UL.n_ul_rb <= 80)
    {
        return 2;
    }
    else
    {
        return 3;
    }
}

void Get_Msc_values(uint32_t SRS_UL.bw_cfg, uint32_t SRS_UL.B, uint32_t SRS_UL.n_ul_rb, uint32_t N_sc)
{
   uint32_t M_sc;
   return  m_srs_b[srsbwtable_idx(SRS_UL.n_ul_rb)][SRS_UL.B][SRS_UL.bw_cfg] * N_sc / 2;/* According to 3GPP 36.211 5.5.3.2*/
}
/***********************************************************************************************************************/
/*
 * Calculate the value of prime numbers for M_sc
 * input- M_sc
 * output -largest prime_number less than M_sc
 */
uint32_t Prime_number(uint32_t M_sc)
{
    /* get largest prime no N_zc<M_sc */
    int i,j;
    int prime;
    int Largest_prime;
    /*printf("\nprime numbers of M_sc[%d] are:\n",M_sc);*/
    for(i = 1; i < M_sc; i++)
    {
        for(j = 2; j < i; j++) if(i % j == 0) break;
        if(j == i)
            prime = j;/* find prime numbers*/
    }
    Largest_prime = prime;
    /*printf("Largest primes %d\n",Largest_prime);*/
    return Largest_prime;
}
/********************************************************************************************************************/
/*
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
//CHECKED
/**********************************************************************************************************************/

/***********************************************************************************************/
/*  May 17th*/      //CHECKED
/***********************************************************************************************/
/* //Compute and generate pseudo random sequence
* input - c_init, len
* output- n_prs values
*/
void calc_n_prs_c(const uint32_t c_init, const uint32_t len, uint8_t* n_prs)
{
    uint32_t N_c = 1600;
    uint32_t x1, x2;
    x1 = 0x54D21B24; // x1 sequence moved forward
    x2 = c_init;
    uint32_t next_x1 = 0;
    uint32_t next_x2 = 0;

    // Move x2 forward
    uint32_t i;
    for(i = 0; i < N_c - 31; i++)
    {
        next_x2 = ((x2 >> 3) ^ (x2 >> 2) ^ (x2 >> 1) ^ x2) & 0x1;

        x2 = (x2 >> 1) | (next_x2 << 30);
    }

    // Calculate prs_c
    for(i = 0; i < len; i++)
    {
        next_x1 = ((x1 >> 3) ^ x1) & 0x1;
	next_x2 = ((x2 >> 3) ^ (x2 >> 2) ^ (x2 >> 1) ^ x2) & 0x1;

	*(n_prs++) = next_x1 ^ next_x2;

	x1 = (x1 >> 1) | (next_x1 << 30);
	x2 = (x2 >> 1) | (next_x2 << 30);
    }

}
/****************************************************************************************************************/
/*
 * Calculate Sequence Hopping(v)
 * input- CellID,Delta_ss,
 * output-Calculate v
 */
uint32_t Get_v_value(const uint16_t SRS_UL.cell_ID, uint32_t SRS_UL.delta_ss, uint32_t M_sc, uint32_t N_sc, uint32_t SRS_UL.ns, uint32_t SRS_UL.sequence_hopping, uint32_t SRS_UL.group_hopping, uint32_t SRS_UL.N_ID_PUCCH, uint32_t SRS_UL.N_ID_PUSCH, const uint16_t SRS_UL.PUCCH_ID, const uint16_t SRS_UL.PUSCH_ID)
{
    /*Sequence hopping only applies for reference-signals of length M_sc_RS ≥ 6*N_sc_RB .For M sc RS < 6 N sc RB ,  v within u is v = 0 .*/
    uint32_t v;
    if ( M_sc >= 6*N_sc)
    {
        uint32_t len = SRS_UL.ns;
        uint8_t n_prs[len];
        const uint16_t n_ID = Get_cellID( SRS_UL.N_ID_PUCCH, SRS_UL.N_ID_PUSCH, SRS_UL.cell_ID, SRS_UL.PUCCH_ID, SRS_UL.PUSCH_ID);
        c_init = ((n_ID / 30) << 5) + (((SRS_UL.cell_ID % 30) + SRS_UL.delta_ss) % 30);
        calc_prs_c( c_init, len, n_prs); /*generate_pseudo random sequence*/
        if(SRS_UL.sequence_hopping == 1)
        {
            if(SRS_UL.group_hopping == 1)
            {
                v = 0;
            }
            else
            {
                v=n_prs[ns];
            }
        }
    }
    else
    {
        v = 0;
    }
    return v;
}
/**************************************************************************************************/
/* Calculate N_CellID
* input - cell_ID, Configurations for SRS_UL.N_ID_PUSCH,SRS_UL.N_ID_PUSCH
* output- N_ID
*/
const uint16_t Get_cellID(uint32_t SRS_UL.N_ID_PUCCH, uint32_t SRS_UL.N_ID_PUSCH, const uint16_t SRS_UL.cell_ID, const uint16_t SRS_UL.PUCCH_ID, const uint16_t SRS_UL.PUSCH_ID)
{
    const uint16_t N_ID;
    if (SRS_UL.N_ID_PUCCH == 1)
    {
        const uint16_t N_ID= SRS_UL.PUCCH_ID;
        return  N_ID;
    }
    else if (SRS_UL.N_ID_PUSCH == 1)
    {
        const uint16_t N_ID= SRS_UL.PUSCH_ID;
        return  N_ID;
    }
    else
    {
        const uint16_t N_ID= SRS_UL.cell_ID;// for SRS also//
        return  N_ID;
    }
}

/****************************************************************************************************************/

/* Calculate f_gh
 * input- number of slots, group hopping -enabled/disabled, cellID
 * output - f_gh value
 */
uint32_t get_f_gh(const uint16_t SRS_UL.cell_ID, uint32_t SRS_UL.ns, uint32_t SRS_UL.group_hopping, uint32_t SRS_UL.N_ID_PUCCH, uint32_t SRS_UL.N_ID_PUSCH, const uint16_t SRS_UL.PUCCH_ID, const uint16_t SRS_UL.PUSCH_ID)
{
    uint32_t len = ( 8 * SRS_UL.ns + 7 );// Maximum length
    uint8_t n_prs[len];
    uint32_t count = 0;
    uint32_t f_gh;
    const uint16_t n_ID = Get_cellID( SRS_UL.N_ID_PUCCH, SRS_UL.N_ID_PUSCH, SRS_UL.cell_ID, SRS_UL.PUCCH_ID, SRS_UL.PUSCH_ID);
    uint32_t c_init = floor( n_ID / 30 );
        /*generate_pseudo random sequence /** Computes n_prs values as defined in 5.5.2.1.1 of 36.211 */
    calc_prs_c(c_init, len, n_prs);
    /** Computes sequence-group pattern f_gh according to 5.5.1.3 of 36.211 */
    int i;
    if(SRS_UL.group_hopping == 1)
    {
        for (i = 0; i < 8; i++)
        {
            count += (((uint32_t) n_prs[8 * SRS_UL.ns + i]) << i) % 30;
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
/******************************************/
/*
 * Calculate f_ss

 * input-  cellID
 * output - f_ss value
 */
uint32_t get_f_ss(const uint16_t SRS_UL.cell_ID, uint32_t SRS_UL.N_ID_PUCCH, uint32_t SRS_UL.N_ID_PUSCH, const uint16_t SRS_UL.PUCCH_ID, const uint16_t SRS_UL.PUSCH_ID)
{
    uint32_t f_ss;
    const uint16_t n_ID = Get_cellID( SRS_UL.N_ID_PUCCH, SRS_UL.N_ID_PUSCH, SRS_UL.cell_ID, SRS_UL.PUCCH_ID, SRS_UL.PUSCH_ID);
    f_ss = n_ID % 30;
    return f_ss;
}
/************************************************/
/*
 * Calculate Group Hopping(u)u=(f_gh(ns)+f_ss)%30;
 * input- number of slots, group hopping -enabled/disabled, cellID
 * output - u value
 */
uint32_t Get_u_value(const uint16_t SRS_UL.cell_ID, uint32_t SRS_UL.ns, uint32_t SRS_UL.group_hopping, uint32_t SRS_UL.N_ID_PUSCH, uint32_t SRS_UL.N_ID_PUSCH, const uint16_t SRS_UL.PUCCH_ID, const uint16_t SRS_UL.PUSCH_ID)
{
    uint32_t f_ss, f_gh;
    uint32_t u;
    f_ss = get_f_ss( SRS_UL.cell_ID, SRS_UL.N_ID_PUCCH, SRS_UL.N_ID_PUSCH,SRS_UL.PUCCH_ID, SRS_UL.PUSCH_ID);
    if(group_hopping == 1)
    {
        f_gh = get_f_gh(SRS_UL.cell_ID, SRS_UL.ns, SRS_UL.group_hopping,SRS_UL.N_ID_PUCCH, SRS_UL.N_ID_PUSCH, SRS_UL.PUCCH_ID, SRS_UL.PUSCH_ID);
        u = ( f_gh + f_ss ) % 30;
        return u;
    }
    else 
    {
        const uint16_t n_ID = Get_cellID( SRS_UL.N_ID_PUCCH, SRS_UL.N_ID_PUSCH, SRS_UL.cell_ID, SRS_UL.PUCCH_ID, SRS_UL.PUSCH_ID);
        u = n_ID % 30;
        return u;
    }
}
/**************************************************************************************************/
/*
 * Calculate Alpha value for SRS according to 5.5.3.1 of 36.211
 * input- N_Tx, SRS_UL.K_Tc, n_srs_cs from higher layers
 * output - alpha
 */
static float alpha_p(uint32_t SRS_UL.N_Tx, uint32_t SRS_UL.Cyclic_shift, uint32_t SRS_UL.K_Tc)
{
    uint32_t n_srs_max , p;
    uint32_t n_srs_p;
    float alpha;
    if ( SRS_UL.K_Tc == 2 )
    {
        n_srs_max = 8;/* for all SRS signals*/
    }
    else
    {
        n_srs_max = 12;
    }
    p = SRS_UL.N_Tx - 1;
    n_srs_p = ( SRS_UL.Cyclic_shift + ( ( n_srs_max * p) / SRS_UL.N_Tx ) ) % n_srs_max;
    alpha = 2 * M_PI * n_srs_p / n_srs_max;
    return alpha;
}

/***********************************************************************************************/
/*  May 16th*/   //..... CHECKED
/***********************************************************************************************/
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
static void Root_q_arg(float *root_q_arg, uint32_t M_sc, uint32_t u, uint32_t v, uint32_t SRS_UL.n_ul_rb, uint32_t N_zc)
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
        Seq_2Nsc_exp[i] = Phi_M_sc_24[u][i] * M_PI / 4));

    }
}
/***********************************************************************************************/
// Calculate arguments foe exponential in r_uv(n)
/*
 * Calculate the value of arguments for N_PRB=1,2,<=3
 * input- N_PRB,Grp no, Seq No, other values needed
 * output -argument values to find sequence */
/* Computes argument of r_u_v signal */
static void compute_r_uv_arg(float *r_uv, uint32_t SRS_UL.n_ul_rb, uint32_t SRS_UL.ns, uint32_t N_sc, uint32_t N_zc, uint32_t SRS_UL.sequence_hopping, uint32_t SRS_UL.group_hopping,  uint32_t u, uint32_t v)
{
    // get u

    if (n_ul_rb == 1)
    {
        Seq_Msc12_Exp(r_uv, u);
    }
    else if (n_ul_rb == 2)
    {
        Seq_Msc24_Exp(r_uv, u);
    }
    else
    {
        float r_uv_temp[Nzc];
        //calculate sequence number v
        v = Get_v_value( SRS_UL.cell_ID, SRS_UL.delta_ss, M_sc, N_sc, SRS_UL.ns, SRS_UL.sequence_hopping, SRS_UL.group_hopping);
        Root_q_arg(r_uv_temp, N_sc * SRS_UL.n_ul_rb, u, v, N_zc);
        r_uv_mprb(r_uv, r_uv_xq, M_sc, N_zc);
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
/*
struct {
            _Complex r_srs[20][100];
       } srs_sequence;
*/
/* Generate SRS signal as defined in Section 5.5.3.1 */
/*int srs_gen(struct srs_sequence *seq, uint32_t N_sc,uint32_t SRS_UL.n_ul_rb uint32_t N_zc, uint32_t SRS_UL.sf_idx, uint32_t SRS_UL.group_hopping, uint32_t SRS_UL.sequence_hopping, const uint16_t SRS_UL.cell_ID, uint32_t SRS_UL.delta_ss, uint32_t SRS_UL.ns, uint32_t SRS_UL.N_Tx, uint32_t SRS_UL.Cyclic_shift, uint32_t SRS_UL.K_Tc)*/
int srs_gen(_Complex *r_srs, uint32_t N_sc, uint32_t n_ul_rb,uint32_t sf_idx, uint32_t group_hopping, uint32_t sequence_hopping, const uint16_t cell_ID, uint32_t delta_ss, uint32_t ns, uint32_t N_Tx, uint32_t Cyclic_shift, uint32_t K_Tc, uint32_t bw_cfg,uint32_t B,uint32_t N_ID_PUCCH, uint32_t N_ID_PUSCH, const uint16_t PUCCH_ID, const uint16_t PUSCH_ID)
{
    int n ;
    uint32_t nslot;
    float root_q_arg[N_zc];
    float Seq_Nsc_exp[N_sc];
    float Seq_2Nsc_exp[N_sc];
    float r_uv[M_sc];
    float r_uv_xq[M_sc];

    uint32_t M_sc = Get_Msc_values( SRS_UL.bw_cfg, SRS_UL.B, SRS_UL.n_ul_rb, N_sc);
    uint32_t N_zc = Get_Nzc( M_sc);
    uint32_t fgh = get_f_gh(SRS_UL.cell_ID, SRS_UL.ns, SRS_UL.group_hopping, SRS_UL.N_ID_PUCCH, SRS_UL.N_ID_PUSCH, SRS_UL.PUCCH_ID, SRS_UL.PUSCH_ID);
    uint32_t f_ss = get_f_ss(SRS_UL.cell_ID, SRS_UL.N_ID_PUCCH, SRS_UL.N_ID_PUSCH, SRS_UL.PUCCH_ID, SRS_UL.PUSCH_ID);
    uint32_t u = Get_u_value( SRS_UL.cell_ID, SRS_UL.ns, SRS_UL.group_hopping, SRS_UL.N_ID_PUCCH, SRS_UL.N_ID_PUSCH, SRS_UL.PUCCH_ID, SRS_UL.PUSCH_ID);
    uint32_t v = Get_v_value( SRS_UL.cell_ID, SRS_UL.delta_ss, M_sc, N_sc, SRS_UL.ns, SRS_UL.sequence_hopping, SRS_UL.group_hopping, SRS_UL.N_ID_PUCCH, SRS_UL.N_ID_PUSCH, SRS_UL.PUCCH_ID, SRS_UL.PUSCH_ID);
    float q = get_qvalue(u, v, N_zc);
    for (nslot = 2*SRS_UL.sf_idx; nslot < 2*(SRS_UL.sf_idx + 1); nslot++)
    {
        float alpha = alpha_p(SRS_UL.N_Tx, SRS_UL.Cyclic_shift, SRS_UL.K_Tc);//n_srs
        compute_r_uv_arg(r_uv, SRS_UL.n_ul_rb, u, v, N_sc, N_zc, SRS_UL.ns);
        // Do complex exponential and adjust amplitude
        for ( n = 0; n < M_sc; n++)
        {
            //seq->r_srs[nslot][n] = cexpf( I * ( r_uv[n] + ( alpha * n ) ) );
            r_srs[(nslot % 2)*M_sc+n] = cexpf( I * ( r_uv[n] + ( alpha * n ) ) );
        }
    }
}


/****************************************************************************************************************/

/***********************************************************************************************/
/*  May 19th & 20th*/
/***********************************************************************************************/
uint32_t T_srs_FDD(uint32_t I_srs) {
  uint32_t T_srs; 
  /* This is Table 8.2-1 */
  if (I_srs < 2) {
    T_srs   = 2; 
  } else if (I_srs < 7) {
    T_srs   = 5; 
  } else if (I_srs < 17) {
    T_srs   = 10; 
  } else if (I_srs < 37) {
    T_srs   = 20; 
  } else if (I_srs < 77) {
    T_srs   = 40; 
  } else if (I_srs < 157) {
    T_srs   = 80; 
  } else if (I_srs < 317) {
    T_srs   = 160; 
  } else if (I_srs < 637) {
    T_srs   = 320; 
  } else {
    T_srs = 0; 
  }
  return T_srs; 
}
/*Table 8.2-2: UE Specific SRS Periodicity T_SRS for trigger type 0, TDD according to 3GPP TS 36.213 version 13.0.0 Release 13*/ 
uint32_t T_srs_TDD(uint32_t I_srs) {
  uint32_t T_srs; 
  /* This is Table 8.2-2 */
  if (I_srs < 10) {
    T_srs   = 2; 
  } else if (I_srs < 15) {
    T_srs   = 5; 
  } else if (I_srs < 25) {
    T_srs   = 10; 
  } else if (I_srs < 45) {
    T_srs   = 20; 
  } else if (I_srs < 85) {
    T_srs   = 40; 
  } else if (I_srs < 165) {
    T_srs   = 80; 
  } else if (I_srs < 325) {
    T_srs   = 160; 
  } else if (I_srs < 645) {
    T_srs   = 320; 
  } else {
    T_srs = 0; 
  }
  return T_srs; 
}


uint32_t T_offset_TDD(uint32_t I_srs) // I_srs is Config_Idx
{
  uint32_t T_offset; 
/*Table 8.2-2 of Subframe Offset Configuration T_offset for trigger type 0, TDD according to 3GPP TS 36.213 version 13.0.0 Release 13*/ 
  if (I_srs == 0) 
    {
        T_offset   = {0,1}; 
    }
    else if (I_srs == 1)
    {
        T_offset   = {0,2}; 
    }
    else if (I_srs == 2)
    {
        T_offset   = {1,2}; 
    }    
    else if (I_srs == 3)
    {
        T_offset   = {0,3}; 
    }
    else if (I_srs == 4)
    {
        T_offset   = {1,3}; 
    }
    else if (I_srs == 5)
    {
        T_offset   = {0,4}; 
    }
    else if (I_srs == 6)
    {
        T_offset   = {1,4};
    } 
    else if (I_srs == 7)
    {
        T_offset   = {2,3};
    }
    else if (I_srs == 8)
    {
        T_offset   = {2,4};
    }
    else if (I_srs == 9)
    {
        T_offset   = {3,4};
    }

    else if (I_srs < 15) 
    {
         T_offset   = I_srs - 10; 
    } 
    else if (I_srs < 25) 
    {
        T_offset   = I_srs - 15; 
    } 
    else if (I_srs < 45) 
    {
        T_offset   = I_srs - 25; 
    } 
    else if (I_srs < 85)     
    {
        T_offset   = I_srs - 45; 
    }   
    else if (I_srs < 165) 
    {
        T_offset   = I_srs - 85; 
    }
    else if (I_srs < 325) 
    {
        T_offset  = I_srs - 165; 
    } 
    else if (I_srs < 645) 
    {
        T_offset  = I_srs - 325; 
    } 
  else 
  {
      T_offset = 0; 
  }
  return T_offset; 
}

/*******************************************************/ 
/* /*Compute n_srs that counts the number of UE-specific SRS transmissions
* input - 
* output- 
*/
Get_n_srs (uint32_t SRS_UL.Config_idx, uint32_t SRS_UL.ns, uint32_t SRS_UL.nf, uint32_t SRS_UL.N_sp)
{
    uint32_t I_srs = SRS_UL.Config_idx;
    uint32_t T_srs_TDD = T_srs_TDD(I_srs);// srs periodicity T_srs
    uint32_t T_offset_TDD = T_offset_TDD(I_srs) 
    uint32_t n_srs;
    if( T_srs_TDD == 2)
    {
        //System frame numberSFN SRS_UL.nf
        //Number of downlink to uplink switch points within the radio frame N_sp
        n_srs = ((2 * SRS_UL.N_sp*SRS_UL.nf) + (2* (SRS_UL.N_sp- 1))* floor(SRS_UL.ns / 10)) + (floor(T_offset_TDD / I_srs - 325)); 
    }
    else
    {
        n_srs = floor(((SRS_UL.nf * 10) + floor(SRS_UL.ns / 2)) / T_srs_TDD);
    }
    return n_srs;
}  // Not checked
/*******************************************************/    
/* /*Compute Fb
* input - HoppingBandwidth,B_srs,freqDomainPosition,K_Tc,bw_cfg,n_ul_rb

* output- 
*/
uint32_t Get_Fb(uint32_t SRS_UL.n_ul_rb,uint32_t SRS_UL.B,uint32_t SRS_UL.bw_cfg, uint32_t SRS_UL.HoppingBandwidth,uint32_t n_srs)
{
    uint32_t b_hop = SRS_UL.HoppingBandwidth;
    uint32_t N_b = Nb[srsbwtable_idx(SRS_UL.n_ul_rb)][SRS_UL.B][SRS_UL.bw_cfg];
    uint32_t prod_NbDash ;
    prod_NbDash = 1;
    uint32_t prod_2NbDash ;
    prod_2NbDash = 1;
    uint32_t bdash;
    uint32_t Fb;
    for (bdash = b_hop+1; bdash <= (SRS_UL.B-1); bdash++) 
    {
        prod_NbDash *= Nb[srsbwtable_idx(SRS_UL.n_ul_rb)][bdash][SRS_UL.bw_cfg];
    }
    for (bdash = b_hop; bdash <= SRS_UL.B; bdash++) 
    {
        prod_2NbDash *= Nb[srsbwtable_idx(SRS_UL.n_ul_rb)][bdash][SRS_UL.bw_cfg];
    }
    if (N_b % 2 ==0)
    {
        Fb = (N_b / 2) * (((n_srs % prod_2NbDash) / (prod_NbDash)) + ((n_srs % prod_2NbDash) / (2 * prod_NbDash)));
    }
    else 
    { 
        Fb =(floor(N_b / 2)) * floor(n_srs / prod_NbDash);
    }  

    return Fb;
}

   //CHECKED
/*******************************************************/    
/*Find n_b values  */
/* input - HoppingBandwidth,B_srs,freqDomainPosition,K_Tc,bw_cfg,n_ul_rb
* output- n_b values
*/

uint32_t Get_nb(uint32_t SRS_UL.HoppingBandwidth, uint32_t SRS_UL.B, uint32_t SRS_UL.n_ul_rb, uint32_t SRS_UL.bw_cfg, uint32_t SRS_UL.freqDomainPosition,uint32_t n_srs)
{
    uint32_t n_b;
    if(SRS_UL.HoppingBandwidth <= SRS_UL.B)
    {
        n_b = (floor (4 *SRS_UL.freqDomainPosition)/ m_srs_b[srsbwtable_idx(SRS_UL.n_ul_rb)][SRS_UL.B][SRS_UL.bw_cfg])) % Nb[srsbwtable_idx(SRS_UL.n_ul_rb)][SRS_UL.B][SRS_UL.bw_cfg]);
    }
    else 
    {
        uint32_t Fb = Get_Fb(SRS_UL.n_ul_rb,SRS_UL.B,SRS_UL.bw_cfg,SRS_UL.HoppingBandwidth,n_srs);
        
        n_b = (Fb + (floor (4 *SRS_UL.freqDomainPosition)/ m_srs_b[srsbwtable_idx(SRS_UL.n_ul_rb)][SRS_UL.B][SRS_UL.bw_cfg]))) % Nb[srsbwtable_idx(SRS_UL.n_ul_rb)][SRS_UL.B][SRS_UL.bw_cfg]);
    }
     return n_b;
}


/***********************************************************************************************/
/*  May 18th*/
/***********************************************************************************************/
uint32_t Get_K_Tc_p(uint32_t SRS_UL.Cyclic_shift, uint32_t SRS_UL.N_Tx, uint32_t SRS_UL.K_Tc)
{
    uint32_t K_Tc_bar = SRS_UL.K_Tc;
    int p_index;
    uint32_t K_Tc_p;//K_Tc_p={0,1,...SRS_UL.K_Tc-1}
    static const unsigned n_SRS_cs = (1 << 4) | (1 << 5) | (1 << 6) | (1 << 7);
    static const unsigned p_hat = (1 << 1) | (1 << 3) ;
    /*find K_Tc_p*/
    for (p_index = 0; p_index <=3; p_index++)
    {
        if (((1 << SRS_UL.Cyclic_shift) & n_SRS_cs) && ((1 << p_index) & p_hat) && (SRS_UL.N_Tx == 4) ) // p_hat ands p relation table Table 5.2.1-1
        {
            K_Tc_p = 1 - K_Tc_bar;// for p_index = 1,3
        }
        else 
        {
            K_Tc_p = K_Tc_bar;// for p_index = 0,2
        }
    }
    return K_Tc_p;
}
/*******************************************************/

    /*For  normal UL subframes  find k_0_pbar  */
uint32_t get_k_0_pbar(uint32_t SRS_UL.bw_cfg, uint32_t N_sc, uint32_t SRS_UL.Cyclic_shift, uint32_t SRS_UL.N_Tx, uint32_t SRS_UL.K_Tc, uint32_t SRS_UL.n_ul_rb )
{
    uint32_t k_0_pbar;
    uint32_t K_Tc_p;//K_Tc_p={0,1,...SRS_UL.K_Tc-1}
    K_Tc_p = Get_K_Tc_p(SRS_UL.Cyclic_shift, SRS_UL.N_Tx, SRS_UL.K_Tc);
    if (SRS_UL.bw_cfg < 8)
    {
          k_0_pbar = (((floor(SRS_UL.n_ul_rb / 2)) - (m_srs_b[srsbwtable_idx(SRS_UL.n_ul_rb)][0][SRS_UL.bw_cfg] / 2 ))*N_sc ) + K_Tc_p;
    }

}


/*******************************************************/
/* /*Compute frequency-domain starting position k0_p -UE specific
* input - HoppingBandwidth,B_srs,freqDomainPosition,K_Tc,bw_cfg,n_ul_rb
* output- k_0_p values
*/

uint32_t Get_Freq_domain_start_k0p(uint32_t SRS_UL.HoppingBandwidth, uint32_t SRS_UL.B, uint32_t SRS_UL.n_ul_rb, uint32_t SRS_UL.bw_cfg, uint32_t SRS_UL.freqDomainPosition, uint32_t N_sc, uint32_t SRS_UL.K_Tc,, uint32_t SRS_UL.Config_idx, uint32_t SRS_UL.ns, uint32_t SRS_UL.nf, uint32_t SRS_UL.N_sp, uint32_t SRS_UL.Cyclic_shift, uint32_t SRS_UL.N_Tx)
{
    if (SRS_UL.bw_cfg < 8 && SRS_UL.B < 4 && SRS_UL.K_Tc < 2) 
    {
        int b;
        uint32_t k_0_p;
        uint32_t k_0_pbar = get_k_0_pbar( SRS_UL.bw_cfg, N_sc, SRS_UL.Cyclic_shift, SRS_UL.N_Tx, SRS_UL.K_Tc, SRS_UL.n_ul_rb);
        uint32_t n_srs = Get_n_srs(SRS_UL.Config_idx,SRS_UL.ns, SRS_UL.nf, SRS_UL.N_sp); 
        for (b = 0; b <= SRS_UL.B;  b++)
        {
	    uint32_t M_sc = Get_Msc_values( SRS_UL.bw_cfg, b, SRS_UL.n_ul_rb, N_sc);
	    uint32_t n_b = Get_nb(SRS_UL.HoppingBandwidth, b, SRS_UL.n_ul_rb, SRS_UL.bw_cfg, SRS_UL.freqDomainPosition,n_srs);
            k_0_p += k_0_pbar + (SRS_UL.K_Tc * M_sc * n_b);
        }
        return k_0_p;
    }
    return 0;
}
/****************************************************************************************************************************************************/




































/************************************************************/
// CHECK??????????????????????????????????????????????????????
// For UPpTS

uint32_t Get_k0p_UpPTS(uint32_t SRS_UL.HoppingBandwidth, uint32_t SRS_UL.B, uint32_t SRS_UL.n_ul_rb, uint32_t SRS_UL.bw_cfg, uint32_t SRS_UL.freqDomainPosition, uint32_t N_sc, uint32_t SRS_UL.K_Tc,, uint32_t SRS_UL.Config_idx, uint32_t SRS_UL.ns, uint32_t SRS_UL.nf, uint32_t SRS_UL.N_sp, uint32_t SRS_UL.Cyclic_shift, uint32_t SRS_UL.N_Tx)
{
    /*For  UpPTS  find k_0_pbar  */
    int b;
    uint32_t k_0_p;
    uint32_t k_0_pbar;
    uint32_t K_Tc_p = Get_K_Tc_p(SRS_UL.Cyclic_shift, SRS_UL.N_Tx, SRS_UL.K_Tc);
    if (((((SRS_UL.nf % 2) * (2 - SRS_UL.N_sp))+ n_hf) % 2 ) == 0)// find n_hf ??????????????????????????????
    {
       k_0_pbar = ((SRS_UL.n_ul_rb - (m_srs0_max[srsbwtable_idx(SRS_UL.n_ul_rb)][0][SRS_UL.bw_cfg]))*N_sc ) + K_Tc_p;
    }
    else
    {
       k_0_pbar = K_Tc_p ;
    }
    uint32_t n_srs = Get_n_srs(SRS_UL.Config_idx,SRS_UL.ns, SRS_UL.nf, SRS_UL.N_sp); 
    for (b = 0; b <= SRS_UL.B;  b++)
    {
	uint32_t M_sc = Get_Msc_values( SRS_UL.bw_cfg, b, SRS_UL.n_ul_rb, N_sc);
	uint32_t n_b = Get_nb(SRS_UL.HoppingBandwidth, b, SRS_UL.n_ul_rb, SRS_UL.bw_cfg, SRS_UL.freqDomainPosition,n_srs);
        k_0_p += k_0_pbar + (SRS_UL.K_Tc * M_sc * n_b);
    }
    return k_0_p;
}



/****************************************************************************************************************/
/* /*Compute m_srs_0[bw_cfg]

* input - bw_cfg
* output- m_srs_0 values   **************************** not required just a general calculation*/
uint32_t Get_m_srs_0( uint32_t SRS_UL.bw_cfg)
{
    int i;
    uint32_t c_srs = SRS_UL.bw_cfg; 
    uint32_t m_srs_0[8];
    for (i=0 ; i < 8; i++)
    {
        m_srs_0 [i] = *((uint32_t*)m_srs_b+ i);
    
    }
return m_srs_0[c_srs];
}
/************************************************************/
/* /*Compute m_srs_max[bw_cfg]

* input - bw_cfg
* output- m_srs_max values   **************************** not required just a general calculation*/


/****************************************************************************************************************/



struct {
           const uint16_t cell_ID;
           uint32_t  B;//B_srs={0,1,2,3} UE specific
           uint32_t bw_cfg; //C_srs is an element of {0,1,2,3,4,5,6,7}/ cell specific parameter
           uint32_t srsSubframeConfig;
           uint32_t n_ul_rb;// must be 6 or greater in value
           uint32_t CyclicPrefixLength;
           uint32_t CP;// Normal -0, Extended-1
           uint32_t Duplex_Mode;// FDD-0,TDD-1
           uint32_t HoppingBandwidth; //b_hop={0,1,2,3}
           uint32_t freqDomainPosition;// n_RRC
           uint32_t freqDomainPosition-ap;// n_RRC
           uint32_t nf; //System frame numberSFN 
           uint32_t Seq_no;
           uint32_t Config_idx;// I_srs {0,..... 644}
           uint32_t K_Tc;//Transmission_comb =2 for SRS
           uint32_t Cyclic_shift;// n_srs_cs
           uint32_t Cyclic_shift-ap;// n_srs_cs
           uint32_t N_Tx;
           uint32_t N_Tx;
           uint32_t ns;//Slot number within a radio frame
           uint32_t N_sp;////Number of downlink to uplink switch points within the radio frame = 5ms
           uint32_t sf_idx;
           uint32_t sequence_hopping;// enable=1 and disbale=0
           uint32_t group_hopping;// enable=1 and disbale=0
           uint32_t delta_ss;
           uint32_t N_ID_PUCCH;// Configured=1 and Not Configured=0
           uint32_t N_ID_PUSCH;// Configured=1 and Not Configured=0
           uint32_t srsMaxUpPTS;// Cell specific
       }SRS_UL;


