/*
 * ref_srsgen.c
 *
 *  Created on: Jun 3, 2016
 *      Author: mistbasejeena
 */


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
//#include <ref_table.h>
//#include <ref_srsgen.h>

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
//#include <ref_table.h>

struct SRS_UL{
           uint32_t cell_ID;
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
	   uint32_t PUCCH_ID;  // Configured  0 or 1
	   uint32_t PUSCH_ID;  // Configured  0 or 1
       };
struct cell {

    uint32_t CP;
    uint32_t CellID;
    uint32_t n_ul_rb;// must be 6 or greater in value

};
#define SUCCESS                0
#define ERROR                  -1
#define ERROR_INVALID_INPUTS   -2

#define format_1                1
#define format_1a               2
#define format_1b               3
#define format_2                4
#define format_2a               5
#define format_2b               6
#define format_3                7
#define format_4                8
#define format_5                9

struct results {

    uint32_t n_cs_format1;
    uint32_t n_cs_format2;
    uint32_t n_cs_format3;

};

struct pucch_config {

    uint8_t CP;// Normal-1 and extended-0
    uint8_t CellID;
    uint8_t delta_ss;
    uint8_t N_cs_1; // No.of cyclic shift for PUCCH 1/1a/1b {0,1,...7}
    uint8_t N_RB_2; // N_RB_2 <= 0 denotes BW available for PUCCH 2/2a/2b
    uint8_t ns;// slot number
    uint32_t NSLOTS_X_FRAME;//number of slots per frame
    uint8_t N_ID_PUCCH;// cell ID if configured 1 else 0
    uint32_t N_UL_RB;
    uint8_t n_pucch_1;//{0,1,.. 2047}
    uint8_t n_pucch_2;
    uint8_t n_pucch_3;
    uint8_t n_pucch_4;
    uint8_t n_pucch_5;
    uint32_t sf_idx;
    uint32_t cyclicShift;
    uint32_t sequence_hopping;// enable=1 and disbale=0
    uint32_t group_hopping;// enable=1 and disbale=0
    uint16_t PUCCH_ID;  // Configured  0 or 1
    uint16_t PUSCH_ID;  // Configured  0 or 1
    uint32_t n_ID_PUCCH;// Configured=1 and Not Configured=0
    uint32_t n_ID_PUSCH;// Configured=1 and Not Configured=0
    uint32_t shortened; //Configured=1 and Not Configured=0
    uint32_t format;//{1,2,.....9}
};
/* //Table 5.2.3-1: Resource block parameters according to 3GPP 36.211 5.2.3*/
uint32_t N_sc = 12;
uint32_t N_symbol_ul[2] = {7,6};
/* Table 5.5.3.3-1: Frame structure type 1 sounding reference signal
   subframe configuration. */
           /*T -sfc Configuration Period*/
uint32_t T_sfc[15] = {1, 2, 2, 5, 5, 5, 5, 5, 5, 10, 10, 10, 10, 10, 10};
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


/*Generating  phi (n ) for M_sc_RS = N_sc_RB .
according to 3GPP TS 36.211 version 13.0.1 Release 13 Table 5.5.1.2-1 */
int  Phi_12[30][12]={{-1,1,3,-3,3,3,1,1,3,1,-3,3},
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

/*Generating  φ (n ) for M sc RS =2* N sc RB .
according to 3GPP TS 36.211 version 13.0.1 Release 13 Table 5.5.1.2-2 */
int Phi_24[30][24]={{-1,3,1,-3,3,-1,1,3,-3,3,1,3,-3,3,1,1,-1,1,3,-3,3,-3,-1,-3},
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
		   {-3,-3,3,1,3,1,-3,3,1,3,1,1,3,3,-1,-1,-3,1,-3,-1,3,1,1,3},
		   {-1,-1,1,-3,1,3,-3,1,-1,-3,-1,3,1,3,1,-1,-3,-3,-1,-1,-3,-3,-3,-1},
		   {-1,-3,3,-1,-1,-1,-1,1,1,-3,3,1,3,3,1,-1,1,-3,1,-3,1,1,-3,-1},
		   {1,3,-1,3,3,-1,-3,1,-1,-3,3,3,3,-1,1,1,3,-1,-3,-1,3,-1,-1,-1},
		   {1,1,1,1,1,-1,3,-1,-3,1,1,3,-3,1,-3,-1,1,1,-3,-3,3,1,1,-3},
		   {1 ,3, 3, 1,-1,-3,3,-1,3 ,3,3,-3,1,-1,1,-1,-3,-1,1,3,-1,3,-3,-3},
		   {-1,-3,3,-3,-3,-3,-1,-1,-3,-1,-3,3,1,3,-3,-1,3,-1,1,-1,3,-3,1,-1},
		   {-3,-3,1,1,-1,1,-1,1,-1,3,1,-3,-1,1,-1,1,-1,-1,3,3,-3,-1,1,-3},
		   {-3,-1,-3,3,1,-1,-3,-1,-3,-3,3,-3,3,-3,-1,1,3,1,-3,1,3,3,-1,-3},
		   {-1,-1,-1,-1,3,3,3,1,3,3,-3,1,3,-1,3,-1,3,3,-3,3,1,-1,3,3},
		   {1,-1,3,3,-1,-3,3,-3,-1,-1,3,-1,3,-1,-1,1,1,1,1,-1,-1,-3,-1,3},
		   {1,-1,1,-1,3,-1,3,1,1,-1,-1,-3,1,1,-3,1,3,-3,1,1,-3,-3,-1,-1},
		   {-3,-1,1,3,1,1,-3,-1,-1,-3,3,-3, 3,1,-3,3,-3,1,-1,1,-3,1,1,1},
		   {-1,-3,3,3,1,1,3,-1,-3,-1,-1,-1,3,1,-3,-3,-1,3,-3,-1,-3,-1,-3,-1},
		   {-1,-3,-1,-1,1,-3,-1,-1,1,-1,-3,1,1,-3,1,-3,-3,3,1,1,-1,3,-1,-1},
		   {1,1,-1,-1,-3,-1,3,-1,3,-1,1,3,1,-1,3,1,3,-3,-3,1,-1,-1,1,3}};
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
 /* To Calculate the value of index in table 5.5.3.2-1 to table 5.5.3.2-4*/
/******************************************************************************/
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
/******************************************************************************/
/*To Calculate the value of M_sc
 * input-  NULRB(uplink BW),Bandwidth B and bw_configuration
 * output -M_sc m_srs
 */
 /*****************************************************************************/
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

        /*             i=7
       summation of [c( 8*srs_ul->ns+ i )⋅ 2^i ]mod 30,if  group hopping is
       enabled           i=0                                                 */
        }
        else
        {
          f_gh = 0;
        }
    return f_gh;
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
		//printf(" M_PI= %f, n_zc = %f, q = %f, m= %d \n",M_PI,n_zc,q,m);
		//printf(" rootq_arg[%d] = %f\n",m,rootq_arg[m]);
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
		//printf(" Seq_2Nsc_exp[%d] = %f\n",i,Seq_2Nsc_exp[i]);

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
    float alpha = alpha_p(srs_ul);//n_srs_cs
    r_uv_arg(r_uv, M_sc, N_zc, u, v);
    // Do complex exponential and adjust amplitude
    for ( n = 0; n < M_sc; n++)
    {
            r_srs[n] = cexpf( I * ( r_uv[n] + ( alpha * n ) ) );
    }
}

/*************************************************************************************************************************************************/
//PUCCH
/*******************************************************/
// n_dmrs_2 table 5.5.2.1.1-1 from 36.211
uint32_t n_dmrs_2[4][8] = {{ 0, 6, 3, 4, 2, 8, 10, 9 },                 // for lambda = 0
		                   { 6, 0, 9, 10, 8, 2, 4, 3 },					// for lambda = 1
                           { 3, 9, 6, 7, 5, 11, 1, 0 },					// for lambda = 2
                           { 9, 3, 0, 1, 11, 5, 7, 6 }};				// for lambda = 3

// n_dmrs_1 table 5.5.2.1.1-2 from 36.211
uint32_t n_dmrs_1[8] = { 0, 2, 3, 4, 6, 8, 9, 10 }; // mapping of cyclic shift to n_dmrs values

/* Orthogonal sequences for PUCCH formats 1a, 1b and 1c. Table 5.5.2.2.1-2*/
float w_arg_pucch_format1_cpnorm[3][3] = {{1, 1, 1},
                                         {1, 2*M_PI/3, 4*M_PI/3},
                                         {1, 4*M_PI/3, 2*M_PI/3}};

float w_arg_pucch_format1_cpext[2][2]  = {{1, 1},
                                         {1, -1}};

/* Orthogonal sequences for PUCCH formats 2, 2a and 2b. Table 5.5.2.2.1-3  */

float w_arg_pucch_format2_cpnorm[2]  = {1, 1};
float w_arg_pucch_format2_cpext[1]   = {1};

/* Orthogonal sequences for PUCCH format 3. Table 5.5.2.2.1-3  */

float w_arg_pucch_format3_cpnorm[2]  = {1, 1};
float w_arg_pucch_format3_cpext[1]   = {1};

uint32_t pucch_dmrs_symbol_format1_cpnorm[3] = {2, 3, 4};
uint32_t pucch_dmrs_symbol_format1_cpext[2] = {2, 3};
uint32_t pucch_dmrs_symbol_format2_3_cpnorm[2] = {1, 5};
uint32_t pucch_dmrs_symbol_format2_3_cpext[1] = {3};
uint32_t pucch_dmrs_symbol_format2a_2b_cpnorm[2] = {1, 5};
uint32_t pucch_dmrs_symbol_format4_5_cpnorm[1] = {3};
uint32_t pucch_dmrs_symbol_format4_5_cpext[1] = {2};
/****************************************************************************/
//PUCCH FUNCTIONS
/****************************************************************************/
/* Number of PUCCH demodulation reference symbols per slot N_rs_pucch tABLE 5.5.2.2.1-1 36.211 */
static uint32_t get_N_rs_PUCCH(uint32_t format, uint32_t CP)
 {
	uint32_t n_rs = 0;
    static const unsigned n_RS_PUCCH1 = (1 << 1) | (1 << 2) | (1 << 3) ;
    static const unsigned n_RS_PUCCH2_3 = (1 << 4) | (1 << 7)  ;
    static const unsigned n_RS_PUCCH2a_2b = (1 << 5) | (1 << 6)  ;

     if (CP)
     {

         if ((1 << format) & n_RS_PUCCH1 )
         {
             n_rs =  3;
         }
         else if ((1 << format) & n_RS_PUCCH2_3 )
         {
        	 n_rs = 2;
         }
         else
         {
        	 if ((1 << format) & n_RS_PUCCH2a_2b)
             {
        		 n_rs = 2;
             }
         }
     }
     else
     {
         if ((1 << format) & n_RS_PUCCH1 )
         {
        	 n_rs = 2;
         }
         else if ((1 << format) & n_RS_PUCCH2_3 )
         {
        	 n_rs = 1;
         }
         else
         {
        	 if ((1 << format) & n_RS_PUCCH2a_2b)
             {
               return 0;
             }
         }
     }
     return n_rs;
 }
/****************************************************************************/
/* Table 5.5.2.2.2-1: Demodulation reference signal location for different PUCCH formats. 36.211 */
/*******************************************************/
static uint32_t get_pucch_dmrs_symbol(uint32_t m, uint32_t format, uint32_t CP)
{
	    static const unsigned n_RS_PUCCH1 = (1 << 1) | (1 << 2) | (1 << 3) ;
	    static const unsigned n_RS_PUCCH2_3 = (1 << 4) | (1 << 7)  ;
	    static const unsigned n_RS_PUCCH2a_2b = (1 << 5) | (1 << 6)  ;
	    static const unsigned n_RS_PUCCH4_5 = (1 << 8) | (1 << 9)  ;
        uint32_t l = 0;

	         if ((1 << format) & n_RS_PUCCH1 )
	         {
	        	 if (CP)   // CP==1
	        	 {
	        	   if (m < 4)
	        	   {
	        		  l = pucch_dmrs_symbol_format1_cpnorm[m];
	        	   }
	        	 }
	        	 else
	        	 {
	        	    if (m < 3)
	        	    {
	        	    	l = pucch_dmrs_symbol_format1_cpext[m];
	        	    }
	        	 }
	         }
	         else if ((1 << format) & n_RS_PUCCH2_3 )
	         {
	        	 if (CP)   // CP==1
	        	 {
	        	    if (m < 3)
	        	    {
	        	    	l =  pucch_dmrs_symbol_format2_3_cpnorm[m];
	        	    }
	        	 }
	             else
	             {
	        	   if (m < 2)
	        	   {
	        		   l =  pucch_dmrs_symbol_format2_3_cpext[m];
	        	   }
	             }
	         }
	         else if ((1 << format) & n_RS_PUCCH2a_2b)
	         {
	        	 if (CP)   // CP==1
	        	 {
	        	    if (m < 3)
	        	    {
	        	    	l =  pucch_dmrs_symbol_format2a_2b_cpnorm[m];
	        	    }
	        	 }
	        	 else
	        	 {
	        	    return 0;
	        	 }
	         }
	   	     else
	         {
	   	    	 if ((1 << format) & n_RS_PUCCH4_5)
	   	    	 {
	   	    		if (CP)   // CP==1
	   	    		{
	   	    			if(m < 2)
	   	    			{
	   	    				l =  pucch_dmrs_symbol_format4_5_cpnorm[m];
	   	    			}
	   	    		}
	   	    		else
	   	    		{
	   	    			if(m < 2)
	   	    			{
	   	    				l =  pucch_dmrs_symbol_format4_5_cpext[m];
	   	    			}
	   	    		}
	   	    	 }
	   	    	 else
	   	    	 {
	   	    		return 0;
	   	    	 }
	         }
	 return  l;
}
/****************************************************************************/
/* Modulated bits 20 and 21 for Formats 2a and 2b as in Table 5.4.2-1 in 36.211 */
/*******************************************************/
int pucch_format2a_2b_mod_symbol(uint32_t format, uint32_t bits[2], float complex *d_10)
{
  if (d_10)
  {
    if (format == format_2a)
    {
      *d_10 = bits[0]?-1.0:1.0;
      return SUCCESS;
    }
    else if (format == format_2b)
    {
      if (bits[0] == 0)
      {
        if (bits[1] == 0)
        {
          *d_10 = 1.0;
        }
        else
        {
          *d_10 = -I;
        }
      }
      else
      {
        if (bits[1] == 0)
        {
          *d_10 = I;
        }
        else
        {
          *d_10 = -1.0;
        }
      }
      return SUCCESS;
    }
    else
    {
      return ERROR;
    }
  }
  else
  {
    return ERROR;
  }
}
/****************************************************************************/
// cyclic shift n_cs_cell
/****************************************************************************/
int get_n_cs_cell(struct pucch_config *cfg, uint32_t n_cs_cell[cfg->NSLOTS_X_FRAME][3],struct SRS_UL *srs_ul,uint32_t ns)
{
    uint32_t len = ns+1;
    uint8_t n_prs[len];
	const uint16_t N_ID = Get_cellID(srs_ul);
	const uint32_t c_init = N_ID;
	calc_prs_c( c_init, len, n_prs);
	uint32_t CP_NSYMB = cfg->CP?3:2;
	uint32_t nslot;
	uint32_t l;
	int i;
	  for (nslot = 0; nslot < cfg->NSLOTS_X_FRAME; nslot++)
	  {
	    for (l = 0; l < CP_NSYMB; l++)
	    {
	      n_cs_cell[nslot][l] = 0;
	      for (i = 0;i < 8; i++)
	      {
	        n_cs_cell[nslot][l] += (((uint32_t) n_prs[8 * nslot * CP_NSYMB + 8 * l + i]) << i);
	      }
	    }
	  }
	  return SUCCESS;
}
/******************************************************************************************************************************/
/*PUCCH format1/1a/1b*/
/******************************************************************************************************************************/
uint32_t get_pucch_format1(struct pucch_config *cfg,uint32_t ns, uint32_t l, uint32_t n_cs_cell,uint32_t* n_oc)
{
	uint32_t temp;
	uint32_t temp1;
	int temp2;
    uint32_t n_prime_slot1;// for slot ns
    uint32_t n_prime_slot2;// for slot ns-1
    uint32_t temp3;
    uint32_t temp_n_cs , temp_1,temp_2;
    uint32_t h_p = 0;
    int d;
    int c ;
    c = cfg->CP?3:2;
    uint32_t Nprime;
    uint32_t n_oc_slot1;uint32_t n_oc_slot2;
    // Calculate Nprime
    Nprime = (cfg->n_pucch_1 < ((c * cfg->N_cs_1)/cfg->delta_ss ))?cfg->N_cs_1:N_sc;
     /*   if (((cfg->ns % 2) == 0)||((cfg->ns % 2) == 1))*/

    temp = (c * N_sc)/cfg->delta_ss;
    temp1 = (c * cfg->N_cs_1)/cfg->delta_ss;
    temp2 = (h_p % c)* Nprime;
    temp3 = temp2 / cfg->delta_ss;

    // Calculate nprime for ns and ns-1
    n_prime_slot1 = (cfg->n_pucch_1 < ((c * cfg->N_cs_1)/cfg->delta_ss ))?cfg->n_pucch_1:(cfg->n_pucch_1 - temp1) % temp;// for ns-1
    // Calculate h_p
    d = cfg->CP?2:0;
    h_p = (n_prime_slot1 + d) % ((c * Nprime)/cfg->delta_ss);
    n_prime_slot2 = (cfg->n_pucch_1 >= ((c * cfg->N_cs_1)/cfg->delta_ss ))?((c *(n_prime_slot1 + 1)) % (temp + 1) - 1):((h_p / c) + temp3);// for ns

    // Calculate n_oc
      n_oc_slot1 =  (n_prime_slot1 * cfg->delta_ss)/Nprime;// for slot ns-1
      n_oc_slot2 =  (n_prime_slot2 * cfg->delta_ss)/Nprime;// for slot ns

    // calculate n_cs
	if (cfg-> CP) // Normal CP
	{
		if ((ns % 2) == 0)// for slot ns-1
		{
			temp_1 = n_prime_slot1 * cfg->delta_ss;
			temp_2 = (int)n_oc_slot1 % (int)cfg->delta_ss;
		}
		else // for slot ns
		{
			temp_1 = n_prime_slot2 * cfg->delta_ss;
			temp_2 = (int) n_oc_slot2 % (int) cfg->delta_ss;
		}
	}
	else // extended CP
	{
		if ((ns % 2) == 0)
		{
			temp_1 = n_prime_slot2 * cfg->delta_ss;
			temp_2 = (int)n_oc_slot1 / 2;
		}
		else
		{
			temp_1 = n_prime_slot2 * cfg->delta_ss;
			temp_2 = (int) n_oc_slot2 / 2;
		}
	}
	temp_n_cs = (temp_1 + temp_2) % Nprime;
	uint32_t n_cs_format1;
	n_cs_format1 = (n_cs_cell + temp_n_cs) % N_sc;
	printf("n_cs1 = %d \n", n_cs_format1);
	/****************************************************************************/
	// ALPHA PUCCH format1/1a/1b
	/****************************************************************************/
    float alpha;
    *n_oc = ((ns % 2)== 0)?n_oc_slot1:n_oc_slot2;
    alpha = 2 * M_PI * n_cs_format1 / N_sc;
    return alpha;
}
/******************************************************************************************************************************/
/*format2/2a/2b */
/******************************************************************************************************************************/
uint32_t get_pucch_format2(struct pucch_config *cfg,uint32_t ns, uint32_t l,uint32_t n_cs_cell)
{

    uint32_t n_prime_slot1format2;
    uint32_t n_prime_slot2format2;
    uint32_t n_cs_format2;
    n_cs_format2 = 0;
    if (cfg->n_pucch_2 < ((cfg->N_RB_2 * N_sc) + ((cfg->N_cs_1/8) * (N_sc - cfg->N_cs_1 - 2))))
    {
    n_prime_slot1format2 = (cfg->n_pucch_2 < (N_sc * cfg->N_RB_2))?(cfg->n_pucch_2 % N_sc):((cfg->n_pucch_2 + cfg->N_cs_1 + 1) % N_sc);
    n_prime_slot2format2 = (cfg->n_pucch_2 < (N_sc * cfg->N_RB_2))?(((N_sc * (n_prime_slot1format2 + 1)) % (N_sc +1))-1):((N_sc - 2 - cfg->n_pucch_2 ) % N_sc);
	n_cs_format2 = ((ns % 2) == 0)?(n_cs_cell + n_prime_slot1format2) % N_sc:(n_cs_cell + n_prime_slot2format2) % N_sc;
    }
	/****************************************************************************/
	// ALPHA PUCCH format2/2a/2b
	/****************************************************************************/
    float alpha;
    alpha = 2 * M_PI * n_cs_format2 / N_sc;
    return alpha;
}
/******************************************************************************************************************************/
/*format 3 */
/******************************************************************************************************************************/
uint32_t get_pucch_format3(struct pucch_config *cfg,uint32_t ns,uint32_t l, uint32_t n_cs_cell)
{

    uint32_t n_prime_slot1format3;
    //uint32_t N_SF0_PUCCH = 0;
    uint32_t N_SF1_PUCCH = 0;
    //N_SF0_PUCCH = (cfg->shortened)?5:5;
    N_SF1_PUCCH = (cfg->shortened)?4:5;
    uint32_t n_oc0;
    uint32_t n_oc1;
    n_oc0 = cfg->n_pucch_3 % N_SF1_PUCCH;// for first slot
    n_oc1 = (N_SF1_PUCCH == 5)?((3 * n_oc0) % N_SF1_PUCCH):(n_oc0 % N_SF1_PUCCH);// for second slot
    switch(n_oc0)
    {
       case 0:
           n_prime_slot1format3 = 0;
           break;
       case 1:
           n_prime_slot1format3 = 3;
           break;
       case 2:
           n_prime_slot1format3 = 6;
           break;
       case 3:
           n_prime_slot1format3 =(N_SF1_PUCCH == 5)?8:9;
           break;
       case 4:
           if(N_SF1_PUCCH == 5)
           {
               n_prime_slot1format3 = 10;
           }
           else
           {
               return ERROR_INVALID_INPUTS;
           }
           break;
    }
    int n_prime_slot2format3;
    n_prime_slot2format3 = 0;
    uint32_t n_cs_format3;
	switch(n_oc1)
    {
	   case 0:
           n_prime_slot2format3 = 0;
           break;
       case 1:
           n_prime_slot2format3 = 3;
           break;
       case 2:
           n_prime_slot2format3 = 6;
           break;
       case 3:
           n_prime_slot2format3 =(N_SF1_PUCCH == 5)?8:9;
           break;
       case 4:
           if(N_SF1_PUCCH == 5)
           {
               n_prime_slot2format3 = 10;
           }
           break;
    }
	n_cs_format3 = ((ns % 2) == 0)?(n_cs_cell + n_prime_slot1format3) % N_sc:(n_cs_cell + n_prime_slot2format3) % N_sc;
	/****************************************************************************/
	// ALPHA PUCCH format 3
	/****************************************************************************/
	float alpha;
    alpha = 2 * M_PI * n_cs_format3 / N_sc;
    return alpha;
}
/****************************************************************************/
/* Generates DMRS for PUCCH formats 1/1a/1b,2/2a/2b,3 according to 5.5.2.2 in 36.211 */
/*******************************************************/
int pucch_dmrs_gen(uint32_t format,struct pucch_config *cfg,struct SRS_UL *srs_ul, uint32_t pucch_bits[2], float complex *r_pucch, struct cell *cell)
{
  int ret = ERROR_INVALID_INPUTS;
  uint32_t m;
  uint32_t ns;
  uint32_t n_rs;
  float complex z_m_1 = 1.0;// set for 1 default
  uint32_t M_sc;
  uint32_t mprime;
  M_sc = N_sc; //M_sc = 12
  if (format == format_2a || format == format_2b)
  {
	pucch_format2a_2b_mod_symbol(format,pucch_bits, &z_m_1);
  }

  for (ns = 2*cfg->sf_idx; ns < 2*(cfg->sf_idx+1);ns++)
  {
        //Calculate u value
        uint32_t u = Get_u_value(srs_ul);

        uint32_t n_oc;
        float Seq_Nsc_exp[N_sc];
        Seq_Msc12_Exp(Seq_Nsc_exp, u, N_sc); // M_sc = 12
        // Calculate n_rs
        n_rs=get_N_rs_PUCCH(format, cell->CP);
        for (m = 0; m < n_rs; m++)
         {
        	   uint32_t n_cs_cells[cfg->NSLOTS_X_FRAME][3];
        	   uint32_t l = get_pucch_dmrs_symbol(m, format, cell->CP);
        	   get_n_cs_cell(cfg, n_cs_cells,srs_ul,ns);
        	   uint32_t n_cs_cell = n_cs_cells[ns][l];
               // Add cyclic prefix alpha
               float alpha = 0.0;
               if (format == format_3)
               {
                 alpha = get_pucch_format3(cfg,ns,l,n_cs_cell);
               }
               if ((format >= format_2) && (format < format_3))
               {
                 alpha = get_pucch_format2(cfg,ns,l,n_cs_cell);
               }
               else
               {
                  alpha = get_pucch_format1(cfg, ns, l, n_cs_cell,&n_oc);
               }

              // Choose number of symbols and orthogonal sequence from Tables 5.5.2.2.1-1 to 5.5.2.2.1-3
             float *w=NULL;
             switch (format)
             {
                 case format_1 :
                 case format_1a :
                 case format_1b :
                      w = cfg->CP?w_arg_pucch_format1_cpnorm[n_oc]:w_arg_pucch_format1_cpext[n_oc];
                      break;
                 case format_2:
                      w = cfg->CP?w_arg_pucch_format2_cpnorm:w_arg_pucch_format2_cpext;
                      break;
                 case format_2a:
                 case format_2b:
                      w = w_arg_pucch_format2_cpnorm;// extended CP not applicable as per Table  5.5.2.2.1-1
                      break;
                 case format_3:
                	  w = cfg->CP?w_arg_pucch_format3_cpnorm:w_arg_pucch_format3_cpext;
             }
             float complex z_m = 1.0;
             if (m == 1)
             {
                z_m = z_m_1;// For format 2a & 2b, z(m) = d(10) when m=1 else z(m) = 1.
             }
             if (w)
             {
                 mprime= ns % 2;// m= 0,1 so it changes value just like for even and odd slot numbers
                 uint32_t n;
                 for (n = 0;n < N_sc; n++)
                 {
                     //r_uv_n[n] = cexpf( I * ( Seq_Nsc_exp[n] + ( alpha * n ) ) ) this is applied in the following with w(p). exp(A.B)=exp(A)+exp(B)
                     r_pucch[mprime*M_sc*n_rs+m*M_sc+n] = z_m*cexpf(I*(w[m]+Seq_Nsc_exp[n]+alpha*n));//  P=1
                 }
            }
            else
            {
              return ERROR;
            }
         }
    }
  return ret;
}
/****************************************************************************************************************************************/
// PUCCH format 4 and 5
/**********************************************/

/****************************************************************************/
// Alpha_lambda for PUCCH format 4 and 5
/****************************************************************************/
static float alpha_lambda(struct pucch_config *cfg,uint32_t format, uint32_t n_oc, uint32_t ns)
{
    float alpha;
    float n_cs_lambda;
    uint32_t n_dmrs_2 = 0;
    uint32_t nslot;
    uint32_t n_PN [cfg->NSLOTS_X_FRAME];
    int i;
    //find n_dmrs_2
    if (format == format_5)
    {
    	if (n_oc == 0)
    	{
    		n_dmrs_2 = 0;
    	}
    	else
    	{
    		if(n_oc ==1)
    		{
    			n_dmrs_2 = 6;
    		}
    	}
    }
    else
    {
    	 if (format == format_4)
    	 {
    		   n_dmrs_2 =0;
    	 }
    }
    uint8_t c_init = ((cfg->CellID / 30) << 5) + (((cfg->CellID % 30) + cfg->delta_ss) % 30);
    uint32_t len = cfg->NSLOTS_X_FRAME;
    uint8_t n_prs[len];
    calc_prs_c( c_init, len, n_prs);
    uint32_t CP_NSYMB = cfg->CP?3:2;
    for (nslot = 0;nslot < cfg->NSLOTS_X_FRAME;nslot++)
    {
      uint32_t n_pn = 0;
      for (i = 0; i < 8; i++)
      {
        n_pn += (n_prs[8 * CP_NSYMB * nslot + i] << i);
      }
      n_PN[nslot] = n_pn;
    }

    n_cs_lambda = (n_dmrs_1[cfg->cyclicShift] + n_dmrs_2 + n_PN[ns]) % 12;
    alpha = 2 * M_PI * n_cs_lambda / 12;
    return alpha;
}

/******************************************************************************/
/* * Calculate Sequence Hopping(v)Tue 24 May 2016 10:18:08 AM
 * input- CellID,SR CEST S_UL.delta_ss,
 * output-Calculate v */
/******************************************************************************/
uint32_t Get_v_value_pucch(struct pucch_config *cfg, uint32_t M_sc,struct SRS_UL *srs_ul)
{
    /*Sequence hopping only applies for reference-signals of
      length M_sc_RS < 6*N_sc_RB .
    * For M sc RS < 6 N sc RB ,  v within u is v = 0 .*/
   uint32_t v;
   v = 0;
   if ( M_sc >= 6*N_sc)
   {
      uint32_t len = cfg->ns+1;
      uint8_t n_prs[len];
      const uint16_t n_ID = Get_cellID(srs_ul);
      uint32_t temp;
      uint32_t temp1;
      temp = ((n_ID / 30) << 5);
      temp1 = (n_ID % 30) + cfg->delta_ss;
      uint32_t  c_init = temp + (temp1 % 30);
      calc_prs_c( c_init, len, n_prs); /*generate_pseudo random sequence*/
      if ((cfg->sequence_hopping == 1) && (cfg->group_hopping == 0))
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

/****************************************************************************/
// Generate DRMS for PUCCH format 4 and 5
/****************************************************************************/
void r_pucch_format4_5(struct pucch_config *cfg,uint32_t format ,uint32_t n_oc, uint32_t M_RB_PUCCH4, float complex *r_pucch,struct SRS_UL *srs_ul )
{
    int m;
    int n;
    uint32_t ns;
    uint32_t M_sc;
    uint32_t M_sc_PUCCH4 = 0;
    // Calculate M_sc
    if ((format_4) && ( M_RB_PUCCH4 <= cfg->N_UL_RB))// value of Msc Pucch is less or equal to nulrb
    {
        M_sc_PUCCH4 = M_RB_PUCCH4 * N_sc;
    }

    M_sc = (format_4)?M_sc_PUCCH4:N_sc;
    uint32_t N_zc = Get_Nzc(M_sc);
    // calculate drms signal
    for (ns = 2*cfg->sf_idx; ns < 2*(cfg->sf_idx+1);ns++)
    {
       // Get group hopping number u
       //Calculate u value
       uint32_t u = Get_u_value(srs_ul);
       //Calculate v value
       uint32_t v = Get_v_value_pucch(cfg,M_sc,srs_ul);

       float alpha = alpha_lambda(cfg,format,n_oc,ns);
       float r_uv[M_sc];
       float complex r_uv_4_5[M_sc];
       r_uv_arg(r_uv, M_sc, N_zc, u, v);
       for (m = 0 ; m < 2 ; m++)
       {
          for (n = 0; n < M_sc; n++)
           {
              r_uv_4_5[n] = cexpf( I * ( r_uv[n] + ( alpha * n ) ) );
              r_pucch[m * M_sc + n] = r_uv_4_5[n];
           }
        }
     }
}









void main ()
{

struct SRS_UL srs;

           srs.Duplex_Mode = 1;// FDD-0,TDD-1
           srs.CP = 1;// Normal -1, Extended-0
           srs.N_Tx =1; //{1,2,4}
           srs.nf = 0; //System frame number SFN
           srs.n_ul_rb = 6;// must be 6 or greater in value
           srs.sf_idx = 10; //Subframe number
           srs.Cyclic_shift = 0;// n_srs_cs
           srs.cell_ID = 1;
           srs.bw_cfg = 0; //C_srs takes {0,1,2,3,4,5,6,7} cell specific
           srs.B = 1;//B_srs={0,1,2,3} UE specific
           srs.sequence_hopping = 1;// enable=1 and disbale=0
           srs.group_hopping = 1;// enable=1 and disbale=0
           srs.Config_idx = 7;// I_srs {0,..... 644}
           srs.HoppingBandwidth = 1; //b_hop={0,1,2,3}


           srs.srsSubframeConfig = 7;
           srs.CyclicPrefixLength;
           srs.freqDomainPosition = 1 ;// n_RRC
           srs.freqDomainPosition_ap;// n_RRC

           srs.OffsetIdx = 0;// {0,1}
           srs.K_Tc = 2;//Transmission_comb =2 for SRS
           srs.transmissionComb = 1;// kbar_tc = {0,1};
           srs.Cyclic_shift_ap=0;// n_srs_cs
           srs.n_hf = 0;// {0,1} first half frame 0 and second half 1 for UpPTS
           srs.n_RA;
           srs.ns = 1;//Slot number within a radio frame
           srs.N_sp;////No.of DL2UL switchpoints within a radioframe = 5ms
           srs.delta_ss = 1; //delta_ss = { 0 , 1 ,..., 29 }
           srs.n_ID_PUCCH = 0;// Configured=1 and Not Configured=0
           srs.n_ID_PUSCH = 0;// Configured=1 and Not Configured=0
           srs.srsMaxUpPTS;// Cell specific
	       srs.PUCCH_ID =0;  //
	       srs.PUSCH_ID =0;  //


struct pucch_config pucch;
	       pucch.CP = 1;// Normal-0 and extended-1
	       pucch.CellID = 1 ;
	       pucch.delta_ss = 1;
	       pucch.N_cs_1 = 6; // No.of cyclic shift for PUCCH 1/1a/1b{0...7}
	       pucch.N_RB_2 = 0; // N_RB_2 <= 0 denotes BW available for PUCCH 2/2a/2b{0...98}
	       pucch.ns = 0;// slot number
	       pucch.NSLOTS_X_FRAME = 20;//number of slots per frame
	       pucch.N_UL_RB = 6;
	       pucch.n_pucch_1 = 1;//{0...2047}
	       pucch.n_pucch_2 =1;
	       pucch.n_pucch_3 =1;
	       pucch.n_pucch_4 =1;
	       pucch.sf_idx = 0;
	       pucch.cyclicShift = 1;
	       pucch.sequence_hopping =0;// enable=1 and disbale=0
	       pucch.group_hopping = 0;// enable=1 and disbale=0
	       pucch.PUCCH_ID = 0;  //
	       pucch.PUSCH_ID = 0;  //
	       pucch.n_ID_PUCCH = 0;// Configured=1 and Not Configured=0
	       pucch.n_ID_PUSCH = 0;// Configured=1 and Not Configured=0
	       pucch.shortened = 0; //Configured=1 and Not Configured=0
           pucch.format = 1; //{1,2.....9}

	       struct cell cells;

	           cells.CP = 1;
	           cells.CellID = 1;
int n;

uint32_t n_prb = N_prb(srs.n_ul_rb, N_sc);

uint32_t n_RE = SRS_NRE(&srs);

srsbwtable_idx(&srs);

uint32_t M_sc = Get_Msc_values(&srs, N_sc);

uint32_t N_zc = Get_Nzc( M_sc);

uint32_t u = Get_u_value( &srs);

_Complex r_srs[M_sc];
uint32_t v =  Get_v_value(&srs, M_sc, N_sc);

const uint16_t n_ID = Get_cellID( &srs);

float q = get_qvalue( u, v, N_zc);
float r_uv[M_sc];
r_uv_arg(r_uv,M_sc, N_zc,u,v);
float alpha = alpha_p(&srs);
printf ("q value = %f \n u value = %d \n v value = %d \n N_ID = %d \n Alpha = %f \n\n",q,u,v,n_ID,alpha);
srs_gen(r_srs,&srs);
for ( n = 0; n < M_sc; n++)
    {
	printf ("r_srs[%d] = %f + i%f \n\n", n, creal(r_srs[n]),cimag(r_srs[n]));
    }

//PUCCH
uint32_t n_rs = get_N_rs_PUCCH(pucch.format,pucch.CP);
uint32_t l =  get_pucch_dmrs_symbol(n_rs, pucch.format, pucch.CP);
uint32_t n_cs_cell[pucch.NSLOTS_X_FRAME][3];
get_n_cs_cell(&pucch,n_cs_cell,&srs,pucch.ns);
for (n = 0; n < pucch.NSLOTS_X_FRAME; n++)
{
  for (l = 0; l < 3; l++)
  {
    	printf ("n_cs_cell[%d][%d] = %d \n\n",n,l,n_cs_cell[n][l]);

  }
}
uint32_t n_oc;
float alpha1 = get_pucch_format1(&pucch, n_cs_cell[pucch.ns][l],pucch.ns, l,&n_oc);
float alpha2 = get_pucch_format2(&pucch,n_cs_cell[pucch.ns][l],pucch.ns,l);
float alpha3 = get_pucch_format3(&pucch,pucch.ns,l,n_cs_cell[pucch.ns][l]);
float alpha_lambda1 = alpha_lambda(&pucch,pucch.format,n_oc,pucch.ns);
printf ("alpha1 = %f \nalpha2 = %f \nalpha3= %f\nalphalambda = %f n_oc = %d \n\n",alpha1,alpha2,alpha3,alpha_lambda1,n_oc);
uint32_t pucch_bits[2] = {0,1};
float complex r_pucch[N_sc];
pucch_dmrs_gen(pucch.format,&pucch,&srs,pucch_bits,r_pucch, &cells);
for ( n = 0; n < N_sc; n++)
    {
	printf ("r_pucch[%d] = %f + i%f \n\n", n, creal(r_pucch[n]),cimag(r_pucch[n]));
    }

}

