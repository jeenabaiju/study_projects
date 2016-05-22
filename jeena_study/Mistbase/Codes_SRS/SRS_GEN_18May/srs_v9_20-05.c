
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
//Find // Base sequence =SRS_UL.nsc and 2SRS_UL.nsc
// Base sequence 3SRS_UL.nsc greater
//Calculate SRS seq
// find PUCCH seq
// place SRS in subframe
// TDD configuration





// srs-configIdx - bw-cfg, B, n_srs_cs(Cyclic shift), SRS_UL.sf_idx, SRS_UL.delta_ss, SRS_UL.group_hopping,sequence_hopping,uint8_t SRS_UL.n_ul_rb=6, HoppingBandwidth, freqdomainposition etc;
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
/* //Table 5.2.3-1: Resource block parameters according to 3GPP 36.211 5.2.3*/
uint32_t N_sc = 12;
uint32_t N_symbol_ul[2] = {7,6};
/* Table 5.5.3.3-1: Frame structure type 1 sounding reference signal subframe configuration. */
uint32_t T_sfc[15] = {1, 2, 2, 5, 5, 5, 5, 5, 5, 10, 10, 10, 10, 10, 10};   /*Configuration Period*/
uint32_t Delta_sfc1[7] = {0, 0, 1, 0, 1, 2, 3}; /*TraSRS_UL.nsmission offset*/
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


/*Generating  phi (n ) for M_sc_RS = N_sc_RB . Table 5.5.1.2-1 */
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
/***********************************************************************************************************************/
/*FUNCTIOSRS_UL.ns FOR SRS GENERATION*/
/***********************************************************************************************************************/
/*
 * Calculate the value of number of physical RBs n_prb
 * input- N_sc , SRS_UL.n_ul_rb for normal and extended CP
 * output- N_PRB for extended and normal CP;
 */
/************************************************************************************************************************/
uint32_t N_prb(uint32_t SRS_UL.SRS_UL.n_ul_rb, uint32_t N_sc)
{
    uint32_t i;
    uint32_t k;
    /* resource elements in frequency domain 0.....SRS_UL.n_ul_rb*N_sc*/
    k = SRS_UL.SRS_UL.n_ul_rb * N_sc;
    uint32_t n_prb = floor(k / N_sc);
    return n_prb;// 6 PRBS in REl-13
}
/***********************************************************************************************/
/*
 * Calculate the value of number of resource elements
 * input- N_sc , N_symbol_ul for normal and extended CP
 * output -N_PRB for extended and normal CP;
*/
/************************************************************************************************/
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
/************************************************************************************************/
/*
 * To Calculate the value of M_sc
 * input-  NULRB(uplink BW),Bandwidth B and bw_configuration
 * output -M_sc m_srs
 */
 /************************************************************************************************/
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

uint32_t  Get_Msc_values(uint32_t SRS_UL.bw_cfg, uint32_t SRS_UL.B, uint32_t SRS_UL.n_ul_rb, uint32_t N_sc)
{
   uint32_t M_sc;
   M_sc = m_srs_b[srsbwtable_idx(SRS_UL.n_ul_rb)][SRS_UL.B][SRS_UL.bw_cfg] * N_sc / 2;/* According to 3GPP 36.211version 13 5.5.3.2 */
   return M_sc;
}
/************************************************************************************************/
/*/*
 * Calculate the value of N_zc
 * input- SRS_UL.bw_cfg,Uplink BW(SRS_UL.n_ul_rb),N_sc, B_srs.
 * output -N_zc largest prime number less than M_sc
 */
 /***********************************************************************************************/
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
/************************************************************************************************/
/*
/* Calculate N_CellID
* input - SRS_UL.cell_ID, ConfiguratioSRS_UL.ns for N_ID_PUCCH,N_ID_PUSCH
* output- N_ID
*/
/************************************************************************************************/
coSRS_UL.nst uint16_t Get_cellID(uint32_t SRS_UL.N_ID_PUCCH, uint32_t SRS_UL.N_ID_PUSCH, coSRS_UL.nst uint16_t SRS_UL.cell_ID, coSRS_UL.nst uint16_t SRS_UL.PUCCH_ID, coSRS_UL.nst uint16_t SRS_UL.PUSCH_ID)
{
    coSRS_UL.nst uint16_t N_ID;
    if (SRS_UL.N_ID_PUCCH == 1)
    {
        coSRS_UL.nst uint16_t N_ID= SRS_UL.PUCCH_ID;
        return  N_ID;
    }
    else if (SRS_UL.N_ID_PUSCH == 1)
    {
        coSRS_UL.nst uint16_t N_ID= SRS_UL.PUSCH_ID;
        return  N_ID;
    }
    else
    {
        coSRS_UL.nst uint16_t N_ID= SRS_UL.cell_ID;// for SRS also//
        return  N_ID;
    }
}
/***********************************************************************************************/
/* //Compute and generate pseudo random sequence
* input - c_init, len
* output- n_prs values

 * written by Henrik
 * */
/***********************************************************************************************/
void calc_prs_c(coSRS_UL.nst uint32_t c_init, coSRS_UL.nst uint32_t len, uint8_t* n_prs)
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
/***********************************************************************************************/
/*
 * Calculate Sequence Hopping(v)
 * input- CellID,SRS_UL.delta_ss,
 * output-Calculate v
 */
/***********************************************************************************************/
uint32_t Get_v_value(coSRS_UL.nst uint16_t SRS_UL.cell_ID, uint32_t SRS_UL.delta_ss, uint32_t M_sc, uint32_t N_sc, uint32_t SRS_UL.ns, uint32_t SRS_UL.sequence_hopping, uint32_t SRS_UL.group_hopping,uint32_t SRS_UL.N_ID_PUCCH, uint32_t SRS_UL.N_ID_PUSCH, coSRS_UL.nst uint16_t SRS_UL.PUCCH_ID, coSRS_UL.nst uint16_t  SRS_UL.PUSCH_ID)
{
        /*Sequence hopping only applies for reference-signals of length M_sc_RS < 6*N_sc_RB .
		 * For M sc RS < 6 N sc RB ,  v within u is v = 0 .*/
        uint32_t v;
        if ( M_sc >= 6*N_sc)
        {
            uint32_t len = SRS_UL.ns+1;
            int i;
            uint8_t n_prs[len];
	        coSRS_UL.nst uint16_t n_ID = Get_cellID( SRS_UL.N_ID_PUCCH, SRS_UL.N_ID_PUSCH, SRS_UL.cell_ID,SRS_UL.PUCCH_ID, SRS_UL.PUSCH_ID);

            uint32_t  c_init = ((n_ID / 30) << 5) + (((n_ID  % 30) + SRS_UL.delta_ss) % 30);
            calc_prs_c( c_init, len, n_prs); /*generate_pseudo random sequence*/
            if ((SRS_UL.sequence_hopping == 1) && (SRS_UL.group_hopping == 0))
            {
		       v = n_prs[len];
			   for (i=0; i<len ; i++)
			   { 
				  printf( "nprs[%d]= %d \n",i,n_prs[i]); 
               }
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
/***********************************************************************************************/
/* Calculate f_gh
 * input- number of slots, group hopping -enabled/disabled, cellID
 * output - f_gh value
 */
/***********************************************************************************************/
uint32_t get_f_gh(coSRS_UL.nst uint16_t SRS_UL.cell_ID, uint32_t SRS_UL.ns, uint32_t SRS_UL.group_hopping, uint32_t SRS_UL.N_ID_PUCCH, uint32_t SRS_UL.N_ID_PUSCH, coSRS_UL.nst uint16_t SRS_UL.PUCCH_ID, coSRS_UL.nst uint16_t SRS_UL.PUSCH_ID)
{
    uint32_t len = ( 8 * SRS_UL.ns + 8 );// Maximum length
    uint8_t n_prs[len];
    uint32_t count = 0;
    uint32_t f_gh;
	
    coSRS_UL.nst uint16_t n_ID = Get_cellID( SRS_UL.N_ID_PUCCH, SRS_UL.N_ID_PUSCH,SRS_UL.cell_ID, SRS_UL.PUCCH_ID, SRS_UL.PUSCH_ID);
    uint32_t c_init = floor( n_ID / 30 );
    /*generate_pseudo random sequence /** Computes n_prs values as defined in 5.5.2.1.1 of 36.211 */
    calc_prs_c(c_init, len, n_prs);
    /** Computes sequence-group pattern f_gh according to 5.5.1.3 of 36.211 */
    f_gh = 0;
    int i;
	if(SRS_UL.group_hopping == 1)
    {
        for (i = 0; i < 8; i++)
        {
            count += (((uint32_t) n_prs[8 * SRS_UL.ns + i]) << i);
        }
        uint32_t f_gh = count % 30;
        return f_gh;

        /*             i=7
        summation of [c( 8*SRS_UL.ns+ i )⋅ 2^i ]mod 30 if group hopping is enabled
                       i=0                                                      */
    }
    else
    {
       uint32_t f_gh = 0;
       return f_gh;
    }
}
/***********************************************************************************************/
/*
 * Calculate f_ss
 * input-  cellID
 * output - f_ss value
 */
 /***********************************************************************************************/
uint32_t get_f_ss(coSRS_UL.nst uint16_t SRS_UL.cell_ID, uint32_t SRS_UL.N_ID_PUCCH, uint32_t SRS_UL.N_ID_PUSCH, coSRS_UL.nst uint16_t SRS_UL.PUCCH_ID, coSRS_UL.nst uint16_t SRS_UL.PUSCH_ID)
{
    uint32_t f_ss;
    coSRS_UL.nst uint16_t n_ID = Get_cellID( SRS_UL.N_ID_PUCCH, SRS_UL.N_ID_PUSCH, SRS_UL.cell_ID, SRS_UL.PUCCH_ID, SRS_UL.PUSCH_ID);
    f_ss = n_ID % 30;
    return f_ss;
}
/***********************************************************************************************/
/*
 * Calculate Group Hopping(u)
 * u=(f_gh(SRS_UL.ns)+f_ss)%30;
 * input- number of slots, group hopping -enabled/disabled, cellID
 * output - u value
 */
/***********************************************************************************************/
static uint32_t Get_u_value(coSRS_UL.nst uint16_t SRS_UL.cell_ID, uint32_t SRS_UL.ns, uint32_t  SRS_UL.group_hopping, uint32_t  SRS_UL.N_ID_PUCCH, uint32_t SRS_UL.N_ID_PUSCH, coSRS_UL.nst uint16_t SRS_UL.PUCCH_ID, coSRS_UL.nst uint16_t SRS_UL.PUSCH_ID)
{
    uint32_t f_ss, f_gh;
    uint32_t u;
    f_ss = get_f_ss( SRS_UL.cell_ID, SRS_UL.N_ID_PUCCH, SRS_UL.N_ID_PUSCH,SRS_UL.PUCCH_ID, SRS_UL.PUSCH_ID);
    f_gh = get_f_gh(SRS_UL.cell_ID, SRS_UL.ns,SRS_UL.group_hopping,SRS_UL.N_ID_PUCCH, SRS_UL.N_ID_PUSCH, SRS_UL.PUCCH_ID, SRS_UL.PUSCH_ID);
    u = ( f_gh + f_ss ) % 30;
    return u;
}
/***********************************************************************************************/
/*
 * Calculate Alpha value for SRS according to 5.5.3.1 of 36.211
 * input- N_Tx, K_Tc, n_srs_cs from higher layers
 * output - alpha
 */
 /***********************************************************************************************/
static float alpha_p(uint32_t SRS_UL.N_Tx , uint32_t SRS_UL.Cyclic_shift, uint32_t SRS_UL.K_Tc)
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
/*
 * Calculate q value
 * input- u,v
 * output -q value to find Zadoff chu seq
 */
 /***********************************************************************************************/
 static float get_qvalue(uint32_t u, uint32_t v, uint32_t N_zc)
{
    float q;
    float q_bar;
    float n_zc = (float) N_zc;
    q_bar = n_zc * (u + 1) / 31;
    q = ( floor ( q_bar + 0.5) ) + ( v * ( pow( (-1) , ( floor( 2 * q_bar) ) ) ) );
	printf("u= %d, v=%d \n", u,v);
    return q;
}
/*************************************************************************************************/
// Calculate argument for Qth root Zadoff-Chu sequence according to 3GPP 36.211 5.5.1.1 to find R_ZC
/*
 * Calculate the value of qth root for ZC
 * input- q value , N_zc
 * output root_q for exponential calculation
 */
 /***********************************************************************************************/
static void Root_q_arg(float *root_q_arg, uint32_t SRS_UL.n_ul_rb, uint32_t u, uint32_t v,  uint32_t N_zc)
{
    int m;
    float n_zc = (float) N_zc;
    float q = get_qvalue(u,v,N_zc);
    for (m = 0; m < N_zc; m++)
    {
        /* argument of x_q(m) according to 3GPP 36.211 5.5.1.1*/
        root_q_arg[m] = (-1.0 * M_PI * q * m * (m + 1)) / n_zc;
		printf(" M_PI= %f, n_zc = %f, q = %f, m= %d \n",M_PI,n_zc,q,m);
		printf(" root_q_arg[%d] = %f\n",m,root_q_arg[m]);
    }
}
/***********************************************************************************************/
/*
 * Calculate the  argument value for exponential calculation of Base sequence for M_sc=N_sc
 * input-  u ,N_sc=no.ofsubcarriers
 * output -Base Sequence value for M_sc=N_sc
 */
 /***********************************************************************************************/
static void Seq_Msc12_Exp(float *Seq_SRS_UL.nsc_exp, uint8_t u, uint32_t N_sc)
{
    int i;
    for (i = 0; i < N_sc; i++)
    {
        Seq_SRS_UL.nsc_exp[i] = Phi_M_sc_12[u][i] * M_PI / 4;
    }
}
/**************************************************************************************************/
 /*
 * Calculate the argument value for exponential calculation of Base sequence for M_sc=2*N_sc
 * input- u ,N_sc=no.ofsubcarriers
 * output -Base Sequence value for M_sc=2*N_sc
 */
 /***********************************************************************************************/
static void Seq_Msc24_Exp(float *Seq_2SRS_UL.nsc_exp, uint8_t u, uint32_t N_sc)
{

    int i;
    for (i = 0 ;i < 2*N_sc; i++)
    {
        Seq_2SRS_UL.nsc_exp[i] = Phi_M_sc_24[u][i] * M_PI / 4;
		printf(" Seq_2SRS_UL.nsc_exp[%d] = %f\n",i,Seq_2SRS_UL.nsc_exp[i]);

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
/***********************************************************************************************/
// Calculate 
/*
 * Calculate the value of arguments for arguments foe exponential in r_uv(n)
 * input- N_PRB,Grp no, Seq No, seq hopping group hopping, u, v
 * output -argument values to find sequence */
/* Computes argument of r_u_v signal */
/***********************************************************************************************/
static void compute_r_uv_arg(float *r_uv, uint32_t SRS_UL.n_ul_rb, uint32_t M_sc, uint32_t N_zc, uint32_t SRS_UL.sequence_hopping, uint32_t SRS_UL.group_hopping,  uint32_t u, uint32_t v)
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
		int i;
        //calculate sequence number v
        Root_q_arg(r_uv_temp, SRS_UL.n_ul_rb, u, v, N_zc);  
		printf(" r_uv_temp[%d] = %f \n ",N_zc-2, r_uv_temp[N_zc-2]);
		printf("u1= %d, v1=%d \n", u,v);
        r_uv_mprb(r_uv, r_uv_temp, M_sc, N_zc);
    }
}
/***********************************************************************************************/
/* Generate SRS signal as defined in Section 5.5.3.1 
 * input - N_sc, SRS_UL.n_ul_rb, SRS_UL.sf_idx, group hopping, seq hopping, cell, PUCCH,PUSH, SRS_UL.N_ID_PUCCH,N_ID_PUSCH, Ntx etc
 * output- SRS sequence 
 * */
/***********************************************************************************************/
int srs_gen(double complex *r_srs, uint32_t N_sc, uint32_t SRS_UL.n_ul_rb,uint32_t SRS_UL.sf_idx, uint32_t SRS_UL.group_hopping, uint32_t SRS_UL.sequence_hopping, coSRS_UL.nst uint16_t SRS_UL.cell_ID, uint32_t SRS_UL.delta_ss, uint32_t SRS_UL.ns, uint32_t SRS_UL.N_Tx , uint32_t SRS_UL.Cyclic_shift, uint32_t SRS_UL.K_Tc, uint32_t SRS_UL.bw_cfg,uint32_t SRS_UL.B,uint32_t SRS_UL.N_ID_PUCCH, uint32_t SRS_UL.N_ID_PUSCH, coSRS_UL.nst uint16_t SRS_UL.PUCCH_ID, coSRS_UL.nst uint16_t SRS_UL.PUSCH_ID)
{
    int n ;
    uint32_t SRS_UL.nslot;
    float Seq_SRS_UL.nsc_exp[N_sc];
    float Seq_2SRS_UL.nsc_exp[N_sc];
    uint32_t M_sc = Get_Msc_values( SRS_UL.bw_cfg, SRS_UL.B, SRS_UL.n_ul_rb, N_sc);
    uint32_t N_zc = Get_Nzc( M_sc);
    float r_uv[M_sc];
    float root_q_arg[N_zc];
    float r_uv_xq[M_sc];
    uint32_t f_gh = get_f_gh(SRS_UL.cell_ID, SRS_UL.ns, SRS_UL.group_hopping, SRS_UL.N_ID_PUCCH, SRS_UL.N_ID_PUSCH, SRS_UL.PUCCH_ID, SRS_UL.PUSCH_ID);
    uint32_t f_ss = get_f_ss(SRS_UL.cell_ID, SRS_UL.N_ID_PUCCH, SRS_UL.N_ID_PUSCH, SRS_UL.PUCCH_ID, SRS_UL.PUSCH_ID);
    uint32_t u = Get_u_value( SRS_UL.cell_ID, SRS_UL.ns, SRS_UL.group_hopping, SRS_UL.N_ID_PUCCH, SRS_UL.N_ID_PUSCH, SRS_UL.PUCCH_ID, SRS_UL.PUSCH_ID);
    uint32_t v = Get_v_value( SRS_UL.cell_ID, SRS_UL.delta_ss, M_sc, N_sc, SRS_UL.ns, SRS_UL.sequence_hopping, SRS_UL.group_hopping, SRS_UL.N_ID_PUCCH, SRS_UL.N_ID_PUSCH, SRS_UL.PUCCH_ID, SRS_UL.PUSCH_ID);
	printf(" uvalue = %d \n", u);
    float q = get_qvalue(u, v, N_zc);
	printf(" q value = %f \n", q);
    float alpha = alpha_p(SRS_UL.N_Tx , SRS_UL.Cyclic_shift, SRS_UL.K_Tc);//n_srs
    compute_r_uv_arg(r_uv, SRS_UL.n_ul_rb,M_sc, N_zc,SRS_UL.sequence_hopping,SRS_UL.group_hopping,u,v);
	printf("u2= %d, v2=%d \n", u,v);
    // Do complex exponential and adjust amplitude
    for ( n = 0; n < M_sc; n++)
    {
            r_srs[n] = cexpf( I * ( r_uv[n] + ( alpha * n ) ) );
    }
}
/***************************************************************************************************************************/
struct {
           const uint16_t cell_ID;
           uint32_t  B;//B_srs={0,1,2,3} UE specific
           uint32_t bw_cfg; //C_srs is an element of {0,1,2,3,4,5,6,7}/ cell specific parameter
           uint32_t srsSubframeConfig;
           uint32_t n_ul_rb;// must be 6 or greater in value
           uint32_t CyclicPrefixLength;
           uint32_t CP;// Normal -0, Extended-1
           uint32_t Duplex_Mode;// FDD-0,TDD-1, Half Duplex -3
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
