/*
 * ref_ul_pucch_srs.c
 *
 *  Created on: May 31, 2016
 *      Author: mistbasejeena
 */

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
#define CP_NORM_NSYMB    7
#define CP_NORM_SF_NSYMB (2*CP_NORM_NSYMB)//subframe
#define CP_NORM_0_LEN    160
#define CP_NORM_LEN      144

#define CP_EXT_NSYMB     6
#define CP_EXT_SF_NSYMB  (2*CP_EXT_NSYMB)
#define CP_EXT_LEN       512
#define CP_EXT_7_5_LEN   1024

#define CP_NORM          1
#define CP_EXT           0
#define CP_ISNORM(cp) (cp==CP_NORM)
#define CP_ISEXT(cp) (cp==CP_EXT)
#define CP_NSYMB(cp) (CP_ISNORM(cp)?CP_NORM_NSYMB:CP_EXT_NSYMB)
#define RE_IDX(nof_prb, symbol_idx, sample_idx) ((symbol_idx)*(nof_prb)*(N_sc) + sample_idx)
/**********************************************************************/

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
           uint32_t n_RA;// number of format 4 PRACH
           uint32_t sf_idx;
           uint32_t sequence_hopping;// enable=1 and disbale=0
           uint32_t group_hopping;// enable=1 and disbale=0
           uint32_t delta_ss; //delta_ss = { 0 , 1 ,..., 29 }
           uint32_t n_ID_PUCCH;// Configured=1 and Not Configured=0
           uint32_t n_ID_PUSCH;// Configured=1 and Not Configured=0
           uint32_t srsMaxUpPTS;// Cell specific
	       uint32_t PUCCH_ID;  // PUCCH ID
           uint32_t PUSCH_ID;  // PUSCH ID
           uint32_t Beta_srs;
       };
struct cell {

    uint32_t CP;
    uint32_t CellID;
    uint32_t n_ul_rb;// must be 6 or greater in value

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

    uint32_t sf_idx;
    uint32_t cyclicShift;
    uint32_t sequence_hopping;// enable=1 and disbale=0
    uint32_t group_hopping;// enable=1 and disbale=0
    uint16_t PUCCH_ID;  // if Configured  PUCCH ID
    uint16_t PUSCH_ID;  // if Configured  PUCCH ID
    uint32_t n_ID_PUCCH;// Configured=1 and Not Configured=0
    uint32_t n_ID_PUSCH;// Configured=1 and Not Configured=0
    uint32_t shortened; //Configured=1 and Not Configured=0
    uint32_t format;//{1,2,.....9}
    uint8_t n_pucch_1;//{0,1,.. 2047}
    uint8_t n_pucch_2;
    uint8_t n_pucch_3;
    uint8_t n_pucch_4;
    uint8_t n_pucch_5;
    uint32_t M_RB_pucch4;//
    uint32_t n_oc5;//  {0 or 1 . this is only for format5}
    uint32_t  pucch_bits[2];
    uint32_t Beta_pucch;
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
    n_RE = (srs_ul->CP)?N_sc * N_symbol_ul[0]:N_sc * N_symbol_ul[1];
         /*if Normal CP=0  and n_re=(12*7=84)*/
         /*Extended CP=1 and n_re(12*6=72)*/
    return n_RE;
}
/******************************************************************************/
/* To Calculate the power
 * input-  a ,b
 * output -pow */
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
 * output -M_sc m_srs*/
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
/* Calculate N_CellID
* input - SRS_UL.cell_ID, Configurations for n_ID_PUCCH,n_ID_PUSCH
* output- N_ID*/
/******************************************************************************/
const uint16_t Get_cellID(struct SRS_UL *srs_ul)
{
    if (srs_ul->n_ID_PUCCH )
    {
        const uint16_t N_ID= srs_ul->PUCCH_ID;
        return  N_ID;
    }
    else if (srs_ul->n_ID_PUSCH)
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
 * written by Henrik*/
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
/* * Calculate Sequence Hopping(v)
 * input- CellID,SRTue 24 May 2016 10:18:08 AM CEST S_UL.delta_ss,
 * output-Calculate v */
/******************************************************************************/
uint32_t Get_v_value_srs(struct SRS_UL *srs_ul, uint32_t M_sc, uint32_t N_sc)
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
      uint32_t temp1;
      temp = ((n_ID / 30) << 5);
      temp1 = (n_ID % 30) + srs_ul->delta_ss;
      uint32_t  c_init = temp + (temp1 % 30);
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
/* * Calculate Sequence Hopping(v) for PUSCH
 * input- CellID,SRTue 24 May 2016 10:18:08 AM CEST S_UL.delta_ss,
 * output-Calculate v */
/******************************************************************************/
uint32_t Get_v_value_pusch(struct SRS_UL *srs_ul, uint32_t M_sc,uint32_t f_ss )
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
      uint32_t fss_pusch;
      fss_pusch = f_ss;
      temp = ((n_ID / 30) << 5);
      uint32_t  c_init = temp + fss_pusch;
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
 * output - f_gh value */
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
/* * Calculate f_ss
 * input-  cellID
 * output - f_ss value */
 /*****************************************************************************/
uint32_t get_f_ss(struct SRS_UL *srs_ul)
{
    uint32_t f_ss;
    const uint16_t n_ID = Get_cellID( srs_ul);
    f_ss = (srs_ul->n_ID_PUSCH == 0)?(srs_ul->cell_ID + srs_ul->delta_ss) % 30:n_ID % 30;
    return f_ss;
}
/******************************************************************************/
/*  Calculate Group Hopping(u)
 * u=(f_gh(ns)+f_ss)%30;
 * input- number of slots, group hopping -enabled/disabled, cellID
 * output - u value */
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
/* Calculate Alpha value for SRS according to 5.5.3.1 of 36.211
 * input- N_Tx, K_Tc, n_srs_cs from higher layers
 * output - alpha */
 /*****************************************************************************/
static float alpha_p(struct SRS_UL *srs_ul)
{
    uint32_t n_srs_max , p;
    uint32_t n_srs_p;
    float alpha;
    n_srs_max =  ( srs_ul->K_Tc == 2 )?8:12;
   /* for all SRS signals*/
    p = srs_ul->N_Tx - 1;
    n_srs_p = (srs_ul->Cyclic_shift + ((n_srs_max * p) / srs_ul->N_Tx)) % n_srs_max;
    alpha = 2 * M_PI * n_srs_p / n_srs_max;
    return alpha;
}
/******************************************************************************/
/*  Calculate q value
 * input- u,v
 * output -q value to find Zadoff chu seq */
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
 * output root_q for exponential calculation */
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
/* Calculate the  argument value for exponential calculation of
    Base sequence for M_sc=N_sc
 * input-  u ,N_sc=no.ofsubcarriers
 * output -Base Sequence value for M_sc=N_sc */
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
 /* Calculate the argument value for exponential calculation of
   Base sequence for M_sc=2*N_sc
 * input- u ,N_sc=no.ofsubcarriers
 * output -Base Sequence value for M_sc=2*N_sc */
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
/* Calculate the value of arguments for arguments foe exponential in r_uv(n)
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
 * output- SRS sequence*/
/******************************************************************************/
int srs_gen(float complex *r_srs,struct SRS_UL *srs_ul,float complex *a)
{
    int n ;
    uint32_t M_sc = Get_Msc_values( srs_ul, N_sc);
    uint32_t N_zc = Get_Nzc( M_sc);
    float r_uv[M_sc];
    uint32_t u = Get_u_value( srs_ul);
    uint32_t v = Get_v_value_srs( srs_ul, M_sc, N_sc);
    float alpha = alpha_p(srs_ul);//n_srs_cs
    r_uv_arg(r_uv, M_sc, N_zc, u, v);
    // Do complex exponential and adjust amplitude
    for ( n = 0; n < M_sc; n++)
    {
            r_srs[n] = cexpf( I * ( r_uv[n] + ( alpha * n ) ) );
            a[n] = srs_ul->Beta_srs * r_srs[n];
    }
return 0;
}
/******************************************************************************/
/* Mapping srs to physcal resources*/
/******************************************************************************/
int refsignal_srs_map(struct SRS_UL *srs_ul,float complex *r_srs, float complex *sf_symbols,uint32_t Fb, uint32_t N_b , uint32_t k_pbar,  uint32_t K_Tc_p)
{
  int ret = ERROR_INVALID_INPUTS;int i;
  if (r_srs && srs_ul)
    {
    ret = ERROR;
    uint32_t M_sc = Get_Msc_values(srs_ul,N_sc);
    uint32_t k0p;
    Get_Freq_domain_start_k0p(srs_ul,Fb, N_b ,k_pbar,K_Tc_p,&k0p);
    for (i=0;i<M_sc;i++)
    {
      sf_symbols[RE_IDX(srs_ul->n_ul_rb, 2*CP_NSYMB(srs_ul->CP)-1, k0p + 2*i)] = r_srs[i];
    }
    ret = SUCCESS;
  }
  return ret;
}

/******************************************************************************/
/* Generate T_srs for FDD according to 3GPP TS 36.213 version 13.0.0
   Release 13 Section 8.2
 * input - Config_idx
 * output- T_srs */
/******************************************************************************/
uint32_t T_srs_FDD(uint32_t Config_idx)
{
    uint32_t I_srs = Config_idx;
    uint32_t T_srs;
    /* This is Table 8.2-1 */
    if (I_srs < 2)
    {
        T_srs   = 2;
    }
    else if (I_srs < 7)
    {
        T_srs   = 5;
    }
    else if (I_srs < 17)
    {
        T_srs   = 10;
    }
    else if (I_srs < 37)
    {
        T_srs   = 20;
    }
    else if (I_srs < 77)
    {
        T_srs   = 40;
    }
    else if (I_srs < 157)
    {
        T_srs   = 80;
    }
    else if (I_srs < 317)
    {
        T_srs   = 160;
    }
    else if (I_srs < 637)
    {
        T_srs   = 320;
    }
    else
    {
        T_srs = 0;
    }
    return T_srs;
}
/******************************************************************************/
/* Generate T_srs for TDD according to 3GPP TS 36.213 version 13.0.0
    Release 13 Section 8.2
 * input - Config_idx
 * output- T_srs */
/******************************************************************************/

uint32_t T_srs_TDD(uint32_t Config_idx)
{
    uint32_t I_srs = Config_idx;
    uint32_t T_srs;
    /*Table 8.2-2: UE Specific SRS Periodicity T_SRS for trigger type 0,
      TDD according to 3GPP TS 36.213 version 13.0.0 Release 13*/
    if (I_srs < 10)
    {
      T_srs   = 2;
    }
    else if (I_srs < 15)
    {
      T_srs   = 5;
    }
    else if (I_srs < 25)
    {
      T_srs   = 10;
    }
    else if (I_srs < 45)
    {
      T_srs   = 20;
    }
    else if (I_srs < 85)
    {
      T_srs   = 40;
    }
    else if (I_srs < 165)
    {
      T_srs   = 80;
    }
    else if (I_srs < 325)
    {
      T_srs   = 160;
    }
    else if (I_srs < 645)
    {
      T_srs   = 320;
    }
    else
    {
      T_srs = 0;
    }
    return T_srs;
}
/******************************************************************************/
/* Compute T_srs
* input - Config_idx, Duplex_Mode
* output- T_srs*/
/******************************************************************************/
uint32_t T_srs_value(struct SRS_UL *srs_ul)
{
   uint32_t T_srs;
   T_srs =0;
   T_srs = (srs_ul->Duplex_Mode)? T_srs_FDD(srs_ul->Config_idx):T_srs_TDD(srs_ul->Config_idx) ;
return T_srs;
}
/******************************************************************************/
/* Compute n_srs that counts the number of UE-specific SRS transmissions
* input - Config_idx, ns,  nf, N_sp, T_srs, Duplex_Mode, T_offset,  T_offset_max
* output- n_srs*/
/******************************************************************************/
/*Table 8.2-2 of Subframe Offset Configuration T_offset for trigger type 0,
  TDD according to 3GPP TS 36.213 version 13.0.0 Release 13*/
uint32_t Offset(uint32_t T_srs, struct SRS_UL *srs_ul)
{
	uint32_t I_srs = srs_ul->Config_idx;// I_srs is Config_idx
        uint32_t T_offset;
	if (srs_ul->Duplex_Mode == 0) //  FDD mode
	{
	    T_offset = 0;
	}
	else if (srs_ul->Duplex_Mode == 1)
    {
        if (T_srs == 2)
	    {
            if ((I_srs == 0) &&  (srs_ul->OffsetIdx == 0))
            {
               T_offset   = 0;
            }
            else if ((I_srs == 0) &&  (srs_ul->OffsetIdx == 1))
            {
               T_offset  = 1;
            }
            else if ((I_srs == 1) &&  (srs_ul->OffsetIdx == 0))
            {
               T_offset   = 0;
            }
            else if ((I_srs == 1) &&  (srs_ul->OffsetIdx == 1))
            {
			   T_offset = 2;
            }
            else if ((I_srs == 2) &&  (srs_ul->OffsetIdx == 0))
            {
			   T_offset = 1;
            }
            else if ((I_srs == 2) &&  (srs_ul->OffsetIdx == 1))
            {
			   T_offset = 2;
            }
            else if ((I_srs == 3) &&  (srs_ul->OffsetIdx == 0))
            {
			   T_offset = 0;
            }
            else if ((I_srs == 3) &&  (srs_ul->OffsetIdx == 1))
            {
			   T_offset = 3;
            }
            else if ((I_srs == 4) &&  (srs_ul->OffsetIdx == 0))
            {
			   T_offset = 1;
            }
            else if ((I_srs == 4) &&  (srs_ul->OffsetIdx == 1))
            {
			   T_offset = 3;
            }
            else if ((I_srs == 5) &&  (srs_ul->OffsetIdx == 0))
            {
               T_offset   = 0;
            }
            else if ((I_srs == 5) &&  (srs_ul->OffsetIdx == 1))
            {
			   T_offset = 4;
            }
            else if ((I_srs == 6) &&  (srs_ul->OffsetIdx == 0))
            {
               T_offset   = 1;
            }
            else if ((I_srs == 6) &&  (srs_ul->OffsetIdx == 1))
            {
			   T_offset   = 4;
            }
            else if ((I_srs == 7) &&  (srs_ul->OffsetIdx == 0))
            {
               T_offset   = 2;
            }
            else if ((I_srs == 7) &&  (srs_ul->OffsetIdx == 1))
            {
			   T_offset   = 3;
            }
            else if ((I_srs == 8) &&  (srs_ul->OffsetIdx == 0))
            {
               T_offset   = 2;
            }
            else if ((I_srs == 8) &&  (srs_ul->OffsetIdx == 1))
            {
			   T_offset   = 4;
            }
            else if ((I_srs == 9) &&  (srs_ul->OffsetIdx == 0))
            {
               T_offset   = 3;
            }
            else if ((I_srs == 9) &&  (srs_ul->OffsetIdx == 1))
            {
			   T_offset   = 4;
            }
        }
        else
	    {
            if ( (I_srs >=10 ) && (I_srs < 15) )
            {
                T_offset   = I_srs - 10;
            }
            else if ((I_srs >=15) && (I_srs < 25) )
            {
                T_offset   = I_srs - 15;
            }
            else if ( ( I_srs >=25) && (I_srs < 45) )
            {
                T_offset   = I_srs - 25;
            }
            else if (  (I_srs >=45) && ( I_srs < 85)  )
            {
                T_offset   = I_srs - 45;
            }
            else if ((I_srs >=85) && (I_srs < 165) )
            {
                T_offset   = I_srs - 85;
            }
            else if ((I_srs >=165) && (I_srs < 325))
            {
                T_offset  = I_srs - 165;
            }
            else if ( (I_srs >=325) &&  (I_srs < 645) )
            {
                T_offset  = I_srs - 325;
            }
	    }
        }
        else
        {
            T_offset = 0;
        }

return T_offset;
}
/******************************************************************************/
/* Compute n_srs that counts the number of UE-specific SRS transmissions
* input - Config_idx, ns,  nf, N_sp, T_srs, Duplex_Mode, T_offset,  T_offset_max
* output- n_srs*/
/******************************************************************************/
float Get_n_srs(struct SRS_UL *srs_ul, uint32_t T_srs, uint32_t T_offset, uint32_t T_offset_max)
{
    float n_srs;
    n_srs = ((T_srs == 2) && (srs_ul->Duplex_Mode == 1))?((2 * srs_ul->N_sp*srs_ul->nf) + (2* (srs_ul->N_sp- 1))* floor(srs_ul->ns / 10)) + (floor(T_offset / T_offset_max)):floor(((srs_ul->nf * 10) + floor(srs_ul->ns / 2)) / T_srs);
         //System frame number SFN nf
        //Number of downlink to uplink switch points within the radio frame N_sp
    return n_srs;
}
/******************************************************************************/
/* Compute Fb
* input - HoppingBandwidth,B_srs,freqDomainPosition,K_Tc,bw_cfg,n_ul_rb,n_srs,N_b
* output- Fb*/
/******************************************************************************/

uint32_t Get_Fb(struct SRS_UL *srs_ul,uint32_t n_srs, uint32_t N_b)
{
    uint32_t b_hop = srs_ul->HoppingBandwidth;
    uint32_t prod_NbDash ;
    prod_NbDash = 1;
    uint32_t prod_2NbDash ;
    prod_2NbDash = 1;
    uint32_t bdash;
    uint32_t Fb;
    if(srs_ul->B > b_hop)
    {
       for (bdash = b_hop+1; bdash <= (srs_ul->B-1); bdash++)
       {
           prod_NbDash *= Nb[srsbwtable_idx(srs_ul)][bdash][srs_ul->bw_cfg];
       }
       for (bdash = b_hop; bdash <= srs_ul->B; bdash++)
       {
           prod_2NbDash *= Nb[srsbwtable_idx(srs_ul)][bdash][srs_ul->bw_cfg];
       }
       Fb = (N_b % 2 ==0)?(N_b / 2) * (((n_srs % prod_2NbDash) / (prod_NbDash)) + ((n_srs % prod_2NbDash) / (2 * prod_NbDash))):(floor(N_b / 2)) * floor(n_srs / prod_NbDash);
    }
       else
     {
          Fb = 0;
     }
     return Fb;
}
/*******************************************************/
/*Find n_b values  */
/* input - HoppingBandwidth,B_srs,freqDomainPosition,K_Tc,bw_cfg,n_ul_rb
* output- n_b values*/
/*******************************************************/
uint32_t Get_nb(struct SRS_UL *srs_ul, uint32_t Fb, uint32_t N_b)
{
    uint32_t n_b;
    uint32_t n_temp;
    uint32_t n_temp1;
    n_b = 0;
    int b;
    b = 0;
    uint32_t m_sc;
    uint32_t b_hop = srs_ul->HoppingBandwidth;
    if(b_hop < srs_ul->B)
    {
        for(b = 0; b < srs_ul->B ; b++)
        {
           if(b <= b_hop)
           {
	       m_sc = m_srs_b[srsbwtable_idx(srs_ul)][b][srs_ul->bw_cfg];
   	       n_temp = floor(4 *srs_ul->freqDomainPosition);
               n_temp1 = n_temp / m_sc ;
               n_b =n_temp1 % N_b;
           }
           else
           {
               // freq hopping not enabled
	       m_sc = m_srs_b[srsbwtable_idx(srs_ul)][b][srs_ul->bw_cfg];
   	       n_temp = floor(4 *srs_ul->freqDomainPosition);
               n_temp1 = n_temp / m_sc ;
               n_b = (Fb + n_temp1) % N_b;
           }
       }
    }
    else
    {
	       m_sc = m_srs_b[srsbwtable_idx(srs_ul)][b][srs_ul->bw_cfg];
   	       n_temp = floor(4 *srs_ul->freqDomainPosition);
               n_temp1 = n_temp / m_sc ;
               n_b =n_temp1 % N_b;
    }
    return n_b;
}

/*******************************************************/
/* Returns number of RB defined for the cell-specific SRS */
/*Compute m_srs_0[bw_cfg]
* input - bw_cfg
* output- m_srs_0 values   **************************** not required just a general calculation*/
/*******************************************************/
uint32_t Get_m_srs_0( uint32_t bw_cfg)
{
    int i;
    uint32_t c_srs = bw_cfg;
    uint32_t m_srs_0[8];
    for (i=0 ; i < 8; i++)
    {
        m_srs_0 [i] = *((uint32_t*)m_srs_b+ i);
	}
return m_srs_0[c_srs];
}
/*******************************************************/
/* Returns 1 if sf_idx is a valid subframe for SRS transmission according to subframe_config (cell-specific),
 * as defined in Section 5.5.3.3 of 36.211.
 * Returns 0 if no SRS shall be transmitted or -1 if error. */
/*******************************************************/
int srs_send_cs(struct SRS_UL *srs_ul)
{
  /*Table 5.5.3.3-1:Frame structure type 1*/
  uint32_t subframe_config;
  subframe_config = srs_ul->srsSubframeConfig;
  static const unsigned Delta_sfc = (1 << 0) | (1 << 1) | (1 << 2) | (1 << 3);
  static const unsigned Delta_sfc13 = (1 << 0) | (1 << 1) | (1 << 2) | (1 << 3) | (1 << 4) | (1 << 6) | (1 << 8);
  static const unsigned Delta_sfc14 = (1 << 0) | (1 << 1) | (1 << 2) | (1 << 3) | (1 << 4) | (1 << 5) | (1 << 6) | (1 << 8);
  if (subframe_config < 15 && srs_ul->sf_idx < 10)
  {
    //T_sfc[15] = {1, 2, 2, 5, 5, 5, 5, 5, 5, 10, 10, 10, 10, 10, 10};
    uint32_t t_sfc = T_sfc[subframe_config];
    if (subframe_config < 7)
    {
       //Delta_sfc1[7] = {0, 0, 1, 0, 1, 2, 3}
      if ((srs_ul->sf_idx % t_sfc)==Delta_sfc1[subframe_config])
      {
        return 1;
      }
      else
      {
        return 0;
      }
    }
    else if (subframe_config == 7)
    {
      if (((srs_ul->sf_idx % t_sfc)==0) || ((srs_ul->sf_idx % t_sfc)==1))
      {
        return 1;
      }
      else
      {
        return 0;
      }
    }
    else if (subframe_config == 8)
    {
      if (((srs_ul->sf_idx % t_sfc)==2) || ((srs_ul->sf_idx%t_sfc)==3))
      {
        return 1;
      }
      else
      {
        return 0;
      }
    } else if (subframe_config < 13) //Delta_sfc2[4] = {0, 1, 2, 3};
      {
      if ((1 << (srs_ul->sf_idx % t_sfc)) & Delta_sfc)
      {
        return 1;
      }
      else
      {
        return 0;
      }
    }
    else if (subframe_config == 13)
    {
      if ((1 << (srs_ul->sf_idx % t_sfc)) & Delta_sfc13)
      {
        return 1;
      }
      else
      {
        return 0;
      }
    }
    else if (subframe_config == 14)
    {
      if ((1 << (srs_ul->sf_idx % t_sfc)) & Delta_sfc14)
      {
        return 1;
      }
      else
      {
        return 0;
      }
    }
    else
    {
      return 0;
    }
  }
  else
  {
    return -1; // error
  }
}
/*******************************************************/
/*Compute m_srs_max
* input - bw_cfg,B,n_ul_rb,n_RA(PRACH format 4),msrs_0,srsMaxUpPTS
* output- m_srs_max values */
/*******************************************************/
uint32_t m_srs0_max(struct SRS_UL *srs_ul, uint32_t m_srs_0)
{
   int temp;
   int32_t msrs_max;
   int32_t m_srs_c[8];
   int count;
   count =0;
   uint32_t max_temp[8];
   if (srs_ul->srsMaxUpPTS == 1)
   {
	int i;
    temp = srs_ul->n_ul_rb - 6* srs_ul->n_RA;
	for (i = 0 ; i < 8; i++)
	{
	    m_srs_c[i] = m_srs_b[srsbwtable_idx(srs_ul)][0][i];
            if (m_srs_c[i] <= temp)
            {
               max_temp[0] = m_srs_c[i];
               if (max_temp[0] <= max_temp[count])
                {
                   max_temp[0] = max_temp[count];
                   msrs_max = max_temp[0];
                   /*c_srs = i;// c_srs value*/

                }
                count = count + 1;
	    }
            else
            {
               msrs_max = 0;
            }
        }
   }
   else
   {
	msrs_max = m_srs_0;
   }
return msrs_max;
}
/******************************************************************************/
//  Find K_Tc_bar from transmission combs
/******************************************************************************/
uint32_t ktc_value(uint32_t K_Tc, uint32_t *K_Tc_bar)
{
    int i;
    for (i = 0; i< K_Tc; i++)
    {
          K_Tc_bar[i] = i;
    }
return 0;
}
/******************************************************************************/
//  Find K_Tc_p from K_Tc_bar
/******************************************************************************/
uint32_t Get_K_Tc_p(struct SRS_UL *srs_ul)
{
    uint32_t K_Tc_bar[srs_ul->K_Tc];
    ktc_value(srs_ul->K_Tc, K_Tc_bar);
    uint32_t KTc_bar = K_Tc_bar[srs_ul->transmissionComb];// 0 or 1
    int p_index;
    uint32_t K_Tc_p;//K_Tc_p={0,1,...SRS_UL.K_Tc-1}
    static const unsigned n_SRS_cs = (1 << 4) | (1 << 5) | (1 << 6) | (1 << 7);
    static const unsigned p_hat = (1 << 1) | (1 << 3) ;
    /*find K_Tc_p*/
    for (p_index = 0; p_index <=3; p_index++)
    {
        // p_hat ands p relation table Table 5.2.1-1
    	 K_Tc_p = (((1 << srs_ul->Cyclic_shift) & n_SRS_cs) && ((1 << p_index) & p_hat) && (srs_ul->N_Tx == 4))?(1 - KTc_bar):KTc_bar;
        // 1 - KTc_bar for p_index = 1,3
        // KTc_bar for p_index = 0,2
    }
    return K_Tc_p;  // returns 0 or 1 for NTx=1
}
/*******************************************************/
//Find k_0_pbar with K_Tc_p
/*******************************************************/
    /*For  normal UL subframes  find k_0_pbar  */
uint32_t get_k_0_pbar(struct SRS_UL *srs_ul,uint32_t K_Tc_p,uint32_t m_srs_0)
{
    uint32_t k_pbar;
    k_pbar = 0;
    uint32_t K_0_temp;
    uint32_t m_sc,temp1;
    //K_Tc_p={0,1,...SRS_UL.K_Tc-1}
    if (srs_ul->bw_cfg < 8)
    {
	//m_sc = m_srs_b[srsbwtable_idx(n_ul_rb)][0][bw_cfg];
	m_sc= m_srs_0;
	temp1= floor(srs_ul->n_ul_rb / 2) ;
	K_0_temp = temp1-( m_sc / 2 );
        k_pbar =(K_0_temp *N_sc ) + K_Tc_p;
    }
	return k_pbar ;
}
/*******************************************************/
/*Compute frequency-domain starting position k0_p -UE specific
* input - HoppingBandwidth,B_srs,freqDomainPosition,K_Tc,bw_cfg,n_ul_rb, Fb,N_b, k_pbar
* output- k_0_p values*/
/*******************************************************/
uint32_t Get_Freq_domain_start_k0p(struct SRS_UL *srs_ul, uint32_t Fb, uint32_t N_b , uint32_t k_pbar,  uint32_t K_Tc_p,uint32_t *k_0_p)
{
    if (srs_ul->bw_cfg < 8 && srs_ul->B < 4 && K_Tc_p < 2)
    {
        int b;
        *k_0_p = k_pbar ;
        for (b = 0; b <= srs_ul->B;  b++)
        {
	    uint32_t M_sc = Get_Msc_values( srs_ul, N_sc);
	    uint32_t n_b = Get_nb(srs_ul,Fb,N_b);
	    k_0_p += (srs_ul->K_Tc * M_sc * n_b);
        }
    }
    return 0;
}
//*************************************************************************************************************************************************/
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
    switch (format)
    {
    case 1:// format1
    case 2:// format1a
    case 3:// format1b
     if (CP)
     {
         n_rs =  3;
     }
     else
     {
        n_rs = 2;
     }
     break;
    case 4:// format2
    case 5:// format3
     if (CP)
     {
      	 n_rs = 2;
     }
     else
     {
         n_rs = 1;
     }
     break;
    case 6:// format2a
    case 7:// format2b
         if (CP)
         {
          	 n_rs = 2;
         }
         break;
    }
     return n_rs;
 }
/****************************************************************************/
/* Table 5.5.2.2.2-1: Demodulation reference signal location for different PUCCH formats. 36.211 */
/*******************************************************/
 uint32_t get_pucch_dmrs_symbol(uint32_t format, uint32_t CP, uint32_t *loc)
{
  switch (format)
  {
	    case 1:
	    case 2:
	    case 3:

	        if (CP)
	        {
	        loc[0] = pucch_dmrs_symbol_format1_cpnorm[0];//location of ref symbols in pucch 1/1a/1b
            loc[1] = pucch_dmrs_symbol_format1_cpnorm[1];
            loc[2] = pucch_dmrs_symbol_format1_cpnorm[2];
            }
	        else
	        {
	        loc[0] = pucch_dmrs_symbol_format1_cpext[0];//location of ref symbols in pucch 1/1a/1b
            loc[1] = pucch_dmrs_symbol_format1_cpext[1];
            }
            break;
	    case 4:
	    case 7:
	    	if (CP)
	    	{
	    		loc[0] = pucch_dmrs_symbol_format2_3_cpnorm[0];
	    		loc[1] = pucch_dmrs_symbol_format2_3_cpnorm[1];
	  	    }
	  	    else
	  	    {
        		loc[0] =  pucch_dmrs_symbol_format2_3_cpext[0];
	    	}
	    	break;
       case 5:
       case 6:
	        if (CP)
	        {
	  	    	loc[0] = pucch_dmrs_symbol_format2a_2b_cpnorm[0];
	  	    	loc[1] = pucch_dmrs_symbol_format2a_2b_cpnorm[1];
	        }
	        else
            {
                return ERROR;
            }
	      break;
	    case 8:
	    case 9:
	    	if (CP)
	    	{
	    		loc[0] = pucch_dmrs_symbol_format4_5_cpnorm[0];
	    	}
	    	else
	    	{
	    		loc[0] = pucch_dmrs_symbol_format4_5_cpext[0];
	    	}
            break;
  }
  return 0;
}
 /****************************************************************************/
 // cyclic shift n_cs_cell
 /****************************************************************************/
 int get_n_cs_cell(struct pucch_config *cfg, uint32_t n_cs_cell[cfg->NSLOTS_X_FRAME][CP_NSYMB(cfg->CP)],struct SRS_UL *srs_ul,uint32_t n_rs)
 {
    uint32_t nslot,m;
 	uint32_t l[n_rs];
 	int i;
    uint32_t len = 8 * (cfg->NSLOTS_X_FRAME) * CP_NSYMB(cfg->CP) + 8 *CP_NSYMB(cfg->CP) + 7;
    uint8_t n_prs[len];
 	const uint16_t N_ID = Get_cellID(srs_ul);
 	const uint32_t c_init = N_ID;
 	calc_prs_c( c_init, len, n_prs);
    get_pucch_dmrs_symbol(cfg->format,cfg->CP,l);
   /* Generates n_cs_cell according to Sec 5.4 of 36.211 */
     for (nslot=0;nslot<cfg->NSLOTS_X_FRAME;nslot++)
     {
       for (m=0;m<n_rs ;m++)
       {
         n_cs_cell[nslot][l[m]] = 0;
         for ( i=0;i<8;i++)
         {
             n_cs_cell[nslot][l[m]] += n_prs[8*CP_NSYMB(cfg->CP) *nslot+8*l[m]+i]<<i;
         }
       }
     }
 	 return SUCCESS;
 }
/******************************************************************************************************************************/
/*PUCCH format1/1a/1b*/
/******************************************************************************************************************************/
uint32_t get_pucch_format1(struct pucch_config *cfg,struct SRS_UL *srs_ul,uint32_t* n_oc, uint32_t n_rs, float alpha[cfg->NSLOTS_X_FRAME][CP_NSYMB(cfg->CP)])
{
    uint32_t nslot,m, l[n_rs];
    uint32_t c ;
    c = cfg->CP?3:2;
    uint32_t n_cs_cell[cfg->NSLOTS_X_FRAME][CP_NSYMB(cfg->CP)];
   /* Generates n_cs_cell according to Sec 5.4 of 36.211 */
    get_n_cs_cell(cfg,n_cs_cell,srs_ul,n_rs);
    get_pucch_dmrs_symbol(cfg->format,cfg->CP,l);
      for ( nslot = 0; nslot < cfg->NSLOTS_X_FRAME; nslot++)
      {
         for (m= 0; m < n_rs  ; m++)
         {
           printf ("NCellCyclicShifts = [%d]\n ",n_cs_cell[nslot][l[m]]);//nprime[2]
    	 }
      }
	uint32_t temp;uint32_t temp1;    uint32_t temp3;    uint32_t temp_n_cs[2] , temp_1[2],temp_2[2]; int temp2;
    uint32_t h_p = 0;
    uint16_t d;
    uint32_t Nprime,nprime[2],n_oc1[2];

    // Calculate Nprime
    Nprime = (cfg->n_pucch_1 < ((c*cfg->N_cs_1)/cfg->delta_ss ))?cfg->N_cs_1:N_sc;
     /*   if (((cfg->ns % 2) == 0)||((cfg->ns % 2) == 1))*/

     for ( nslot = 0; nslot < cfg->NSLOTS_X_FRAME; nslot++)
     {
           if(nslot%2 == 0)
           {
              // Calculate nprime for ns and ns-1
              temp = (c*N_sc)/cfg->delta_ss;
              temp1 = (c*cfg->N_cs_1)/cfg->delta_ss;
              nprime[0] = (cfg->n_pucch_1 < temp1)?cfg->n_pucch_1:(cfg->n_pucch_1-temp1)%temp;// Nresource Idx for ns-1
              // Calculate h_p
              d = cfg->CP?2:0;
              h_p = (nprime[0]+d) %((c*Nprime)/cfg->delta_ss);
              temp2 = (h_p % c)*Nprime;
              temp3 = temp2/cfg->delta_ss;
              nprime[1] = (cfg->n_pucch_1 >= temp1)?((c *(nprime[0]+1)) % (temp+1)-1):((h_p/c)+temp3);// for ns

             // Calculate n_oc
           if (cfg->CP)
           {
               n_oc1[0] =(nprime[0]*cfg->delta_ss)/Nprime;// for slot ns-1
               n_oc1[1] = (nprime[1]*cfg->delta_ss)/Nprime;// for slot ns
           }
           else
           {
        	   n_oc1[0] =2*(nprime[0]*cfg->delta_ss)/Nprime;// for slot ns-1
        	   n_oc1[1] = 2*(nprime[1]*cfg->delta_ss)/Nprime;// for slot ns
           }
           //printf("\n OrthoSeq n_oc[%d  %d]  \n\n",n_oc1[nslot],n_oc1[nslot+1]);
           *n_oc = ((nslot % 2)== 0)?n_oc1[0]:n_oc1[1];
          }
     }
           uint32_t n_cs[cfg->NSLOTS_X_FRAME][CP_NSYMB(cfg->CP)];

           for ( nslot=0;nslot<cfg->NSLOTS_X_FRAME;nslot++)
           {
              for (m= 0; m < n_rs  ; m++)
              {
            	  n_cs[nslot][l[m]] = 0;
            	  if (cfg->CP)
            	  {
                      n_cs[nslot][l[m]] = (n_cs_cell[nslot][l[m]] + ((nprime[nslot]*cfg->delta_ss +(n_oc1[nslot]%cfg->delta_ss))%Nprime))%N_sc;
                      //printf ("n_cs_cell[%d][%d] = %d \n\n",nslot,l[m],n_cs_cell[nslot][l[m]]);
                     // printf ("n_cs[%d][%d] = %d \n\n",nslot,l[m],n_cs[nslot][l[m]]);
            	  }
            	  else
            	  {
            		  n_cs[nslot][l[m]] = (n_cs_cell[nslot][l[m]] + ((nprime[nslot]*cfg->delta_ss +(n_oc1[nslot]/2))%Nprime))%N_sc;
            	  }
              }
           }
	         /****************************************************************************/
	         // ALPHA PUCCH format1/1a/1b
	         /****************************************************************************/
          // float alpha[cfg->NSLOTS_X_FRAME][CP_NSYMB(cfg->CP)];
           for ( nslot=0;nslot<cfg->NSLOTS_X_FRAME;nslot++)
           {
              for (m= 0; m < n_rs  ; m++)
              {
            	  alpha[nslot][l[m]] = 0;
                  alpha[nslot][l[m]] = (2*M_PI*n_cs[nslot][l[m]])/N_sc;
                  //printf ("alpha[%d][%d] = %f \n\n",nslot,l[m],alpha[nslot][l[m]]);
              }
           }
}
/******************************************************************************************************************************/
/*format2/2a/2b */
/******************************************************************************************************************************/
uint32_t get_pucch_format2(struct pucch_config *cfg,struct SRS_UL *srs_ul,uint32_t n_rs,float alpha[cfg->NSLOTS_X_FRAME][CP_NSYMB(cfg->CP)])
{
	uint32_t nslot,l[n_rs],m;uint32_t nprime[2];
    uint32_t n_cs_cell[cfg->NSLOTS_X_FRAME][CP_NSYMB(cfg->CP)];
    get_pucch_dmrs_symbol(cfg->format,cfg->CP,l);
    uint32_t n_cs[cfg->NSLOTS_X_FRAME][CP_NSYMB(cfg->CP)];

   /* Generates n_cs_cell according to Sec 5.4 of 36.211 */
    get_n_cs_cell(cfg,n_cs_cell,srs_ul,n_rs);

    for ( nslot = 0; nslot < cfg->NSLOTS_X_FRAME; nslot++)
    {
       for (m = 0; m < n_rs  ; m++)
       {
          //printf ("NCellCyclicShifts = [%d]\n ",n_cs_cell[nslot][l[m]]);//nprime[2]
       }

    }

    for ( nslot = 0; nslot < cfg->NSLOTS_X_FRAME; nslot++)
    {
       for (m = 0; m < n_rs  ; m++)
       {
            nprime[0] = (cfg->n_pucch_2 < (N_sc * cfg->N_RB_2))?(cfg->n_pucch_2 % N_sc):((cfg->n_pucch_2 + cfg->N_cs_1 + 1) % N_sc);
            nprime[1] = (cfg->n_pucch_2 < (N_sc * cfg->N_RB_2))?(((N_sc * (nprime[0] + 1)) % (N_sc +1))-1):((N_sc - 2 - cfg->n_pucch_2 ) % N_sc);
            //printf("\n Resource Index nprime[%d  %d]  \n\n",nprime[0],nprime[1]);
            n_cs[nslot][l[m]] = 0;
	        n_cs[nslot][l[m]] = (n_cs_cell[nslot][l[m]] + nprime[nslot]) % N_sc;
	        //printf("\n n_cs[%d][%d]=%d  \n\n",nslot,l[m],n_cs[nslot][l[m]]);
       }
    }

	/****************************************************************************/
	// ALPHA PUCCH format2/2a/2b
	/****************************************************************************/
           for ( nslot=0;nslot<cfg->NSLOTS_X_FRAME;nslot++)
           {
              for (m= 0; m < n_rs  ; m++)
              {
            	  alpha[nslot][l[m]] = 0;
                  alpha[nslot][l[m]] = (2*M_PI*n_cs[nslot][l[m]])/N_sc;
                 //printf ("alpha[%d][%d] = %f \n\n",nslot,l[m],alpha[nslot][l[m]]);
              }
           }
    return 0;
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

/******************************************************************************************************************************/
/*format 3 */
/******************************************************************************************************************************/
uint32_t get_pucch_format3(struct pucch_config *cfg,struct SRS_UL *srs_ul,uint32_t n_rs,float alpha[cfg->NSLOTS_X_FRAME][CP_NSYMB(cfg->CP)])
{

    uint32_t nprime[2],i,m,l[n_rs],nslot;
    uint32_t n_cs[cfg->NSLOTS_X_FRAME][CP_NSYMB(cfg->CP)];
    uint32_t n_cs_cell[cfg->NSLOTS_X_FRAME][CP_NSYMB(cfg->CP)];
    get_pucch_dmrs_symbol(cfg->format,cfg->CP,l);
    /* Generates n_cs_cell according to Sec 5.4 of 36.211 */
    get_n_cs_cell(cfg,n_cs_cell,srs_ul,n_rs);
    //uint32_t N_SF0_PUCCH = 0;
    uint32_t N_SF1_PUCCH = 0;
    //N_SF0_PUCCH = (cfg->shortened)?5:5;
    N_SF1_PUCCH = (cfg->shortened)?4:5;
    uint32_t n_oc[2];
    n_oc[0]= cfg->n_pucch_3 % N_SF1_PUCCH;// for first slot
    n_oc[1] = (N_SF1_PUCCH == 5)?((3 * n_oc[0]) % N_SF1_PUCCH):(n_oc[0] % N_SF1_PUCCH);// for second slot
    //printf("\n OrthoSeq n_oc[%d  %d]  \n\n",n_oc[0],n_oc[1]);
    for (i=0;i<2;i++)
    {
       switch(n_oc[i])
    {
       case 0:
           nprime[0] = nprime[1] =0;
           break;
       case 1:
           nprime[0] = nprime[1] = 3;
           break;
       case 2:
           nprime[0] = nprime[1] = 6;
           break;
       case 3:
           nprime[0] = nprime[1] =(N_SF1_PUCCH == 5)?8:9;
           break;
       case 4:
           if(N_SF1_PUCCH == 5)
           {
               nprime[0] = 10;
           }
           else
           {
               return ERROR_INVALID_INPUTS;
           }
           break;
    }
    }
    //printf("\n Resource Index nprime[%d  %d]  \n\n",nprime[0],nprime[1]);
    for ( nslot = 0; nslot < cfg->NSLOTS_X_FRAME; nslot++)
    {
       for (m = 0; m < n_rs  ; m++)
       {
            n_cs[nslot][l[m]] = 0;
	        n_cs[nslot][l[m]] = (n_cs_cell[nslot][l[m]] + nprime[nslot]) % N_sc;
	        //printf("\n n_cs[%d][%d]=%d  \n\n",nslot,l[m],n_cs[nslot][l[m]]);
	        alpha[nslot][l[m]] = 0;
            alpha[nslot][l[m]] = (2*M_PI*n_cs[nslot][l[m]])/N_sc;
            //printf ("alpha[%d][%d] = %f \n\n",nslot,l[m],alpha[nslot][l[m]]);
       }
    }
return 0;
}
/******************************************************************************/
/* * Calculate Sequence Hopping(v)Tue 24 May 2016 10:18:08 AM
 * input- pucch config, srs config, ns
 * output-Calculate v */
/******************************************************************************/
uint32_t Get_v_value_pucch(struct pucch_config *cfg, uint32_t M_sc,struct SRS_UL *srs_ul, uint32_t ns)
{
    /*Sequence hopping only applies for reference-signals of
      length M_sc_RS < 6*N_sc_RB .
    * For M sc RS < 6 N sc RB ,  v within u is v = 0 .*/
   uint32_t v , i;
   v = 0;
   if ( M_sc >= 6*N_sc)
   {
      uint32_t len = M_sc;
      uint8_t n_prs[len];

      const uint16_t n_ID = Get_cellID(srs_ul);
      uint32_t temp;
      uint32_t temp1;
      temp = ((n_ID / 30) << 5);
      temp1 = (n_ID % 30) + cfg->delta_ss;
      uint32_t  c_init = temp + (temp1 % 30);
      calc_prs_c( c_init, len, n_prs); /*generate_pseudo random sequence*/
      for ( i=0;i<len;i++)
      {
             printf("n_prs[%d] = %d\n",i,n_prs[i]);
      }
      if ((cfg->sequence_hopping == 1) && (cfg->group_hopping == 0))
      {
	  v = n_prs[ns];
      }
   }
   else
   {
        v = 0;
   }
return v;
}
/****************************************************************************/
/* Generates DMRS for PUCCH formats 1/1a/1b,2/2a/2b,3 according to 5.5.2.2 in 36.211 */
/*******************************************************/
int pucch_dmrs_gen(uint32_t format,struct pucch_config *cfg,struct SRS_UL *srs_ul, struct cell *cell, uint32_t n_rs,float complex *r_pucch,uint32_t *l,float complex *r_pucch_a)
{
int ret = ERROR_INVALID_INPUTS;
  uint32_t m; uint32_t n_oc; uint32_t u; uint32_t n,tot;  uint32_t M_sc;  uint32_t mprime;
  float complex z_m_1 = 1.0;// set for 1 default
  float alpha[cfg->NSLOTS_X_FRAME][CP_NSYMB(cfg->CP)];
  float Seq_Nsc_exp[N_sc];
  float complex r_uv_12[N_sc];

  M_sc = N_sc; //M_sc = 12
  if (format == format_2a || format == format_2b)
  {
       pucch_format2a_2b_mod_symbol(format,cfg->pucch_bits, &z_m_1);
  }
  for (mprime = 0; mprime < 2; mprime++)
  {
       u = Get_u_value(srs_ul);// calculate u value for seq group hopping
	   Seq_Msc12_Exp(Seq_Nsc_exp, u, N_sc); // M_sc = 12 sequence
	   for (n = 0;n < N_sc; n++)
	   {
	 	  r_uv_12[n]= cexpf(I*(Seq_Nsc_exp[n]));
	   }
       for (m = 0; m < n_rs; m++)
       {
       // Add cyclic prefix alpha for format 1/1a/1b
          if (format < format_2)
          {
        	 get_pucch_format1(cfg,srs_ul,&n_oc,n_rs,alpha);
          }
          else if ((format >= format_2)&&(format < format_3))
          {
             get_pucch_format2(cfg,srs_ul,n_rs,alpha);
          }
          else if (format == format_3)
          {
              get_pucch_format3(cfg,srs_ul,n_rs,alpha);
          }

          // Choose number of symbols and orthogonal sequence from Tables 5.5.2.2.1-1 to 5.5.2.2.1-3
          float  *w=0;
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
                   break;
          }
          float complex z_m = 1.0;
          if (m == 1)
          {
        	 z_m = z_m_1;// For format 2a & 2b, z(m) = d(10) when m=1 else z(m) = 1.
          }
          if (w)
          {
            for (n = 0;n < N_sc; n++)
            {
              tot=(mprime*M_sc*n_rs+m*M_sc+n);
              r_pucch[tot] = z_m*w[m]*r_uv_12[n]*cexpf(I*alpha[mprime][l[m]]*n);//
        	 // printf("r_pucch[%d] = %.4f + i%.4f \n\n",tot,creal(r_pucch[tot] ),cimag(r_pucch[tot] ));
        	   r_pucch_a[tot] = cfg->Beta_pucch * r_pucch[tot] ;
            }
          }
          else
          {
              return ERROR;
          }
       }
  }
}
/****************************************************************************************************************************************/
// PUCCH format 4 and 5
/**********************************************/

/****************************************************************************/
// Alpha_lambda for PUCCH format 4 and 5
/****************************************************************************/
static float alpha_lambda(struct pucch_config *cfg,uint32_t format, uint32_t n_oc5,float *alpha)
{

    float n_cs_lambda[2];
    uint32_t n_dmrs_2 = 0;
    uint32_t nslot;
    uint32_t n_PN [cfg->NSLOTS_X_FRAME];
    int i;
    //find n_dmrs_2
    if (format == format_5)
    {
    	if (n_oc5 == 0)
    	{
    		n_dmrs_2 = 0;
    	}
    	else
    	{
    		if(n_oc5 ==1)
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
    uint32_t len = 8*CP_NSYMB(cfg->CP)*2+7;
    uint8_t n_prs[len];uint32_t n_pn ;
    calc_prs_c( c_init, len, n_prs);
    for (nslot = 0;nslot < cfg->NSLOTS_X_FRAME;nslot++)
    {
      n_pn = 0;
      for (i = 0; i < 8; i++)
      {
        n_pn += n_prs[8*CP_NSYMB(cfg->CP)*nslot+i] << i;
      }
      n_PN[nslot] = n_pn;
      n_cs_lambda[nslot] = (n_dmrs_1[cfg->cyclicShift] + n_dmrs_2 + n_PN[nslot]) % 12;
      alpha[nslot] = 2 * M_PI * n_cs_lambda[nslot] / 12;
    }
          //printf("n_PN =[%d %d] \n",n_PN[0],n_PN[1]);
          //printf("alphalambda =[%f %f] \n",alpha[0],alpha[1]);
    return 0;
}
/****************************************************************************/
// Generate DRMS for PUCCH format 4 and 5
/****************************************************************************/
void r_pucch_format4_5(struct pucch_config *cfg,struct SRS_UL *srs_ul,uint32_t format ,uint32_t n_oc5, uint32_t M_RB_pucch4, float complex *r_pucch,float complex *r_pucch_a)
{
    uint32_t u;uint32_t v; int n; uint32_t m, tot; uint32_t M_sc;uint32_t M_sc_PUCCH4;
    M_sc_PUCCH4 = 0;
    // Calculate M_sc
    if ((format == format_4) && ( M_RB_pucch4 <= cfg->N_UL_RB))// value of Msc Pucch is less or equal to nulrb
    {
        M_sc_PUCCH4 = M_RB_pucch4 * N_sc;
        M_sc = M_sc_PUCCH4;
    }
    else
    {
        M_sc = N_sc;
    }
    uint32_t N_zc = Get_Nzc(M_sc);
    float alpha[2];
    float r_uv[M_sc];
    float complex r_uv_4_5[M_sc];
    for (m = 0; m < 2;m++)
    {
       // Get group hopping number u
       //Calculate u value
       u = Get_u_value(srs_ul);
       //Calculate v value
       v = Get_v_value_pucch(cfg, M_sc,srs_ul,m);// calculate v value for seq index
       alpha_lambda(cfg,format,cfg->n_oc5,alpha);
       r_uv_arg(r_uv, M_sc, N_zc, u, v);
       for (n = 0; n < M_sc; n++)
           {
              r_uv_4_5[n] = cexpf( I * ( r_uv[n] + ( alpha[m]* n ) ) );
              tot = m * M_sc + n;
              r_pucch[tot] = r_uv_4_5[n];
             // printf("r_pucch[%d] = %.4f + i%.4f \n\n",tot,creal(r_pucch[tot] ),cimag(r_pucch[tot] ));
             // map to physical resources
             r_pucch_a[tot] = cfg->Beta_pucch * r_pucch[tot] ;

           }

    }
}
/****************************************************************************/
// Compute m according to Section 5.4.3 of 36.211
/****************************************************************************/
uint32_t get_pucch_mvalue(struct pucch_config *cfg, uint32_t format, uint32_t n_pucch,uint32_t CP)
{
  uint32_t m=0; int i;
  uint32_t N_SF0_PUCCH ;
  N_SF0_PUCCH = 5;
  uint32_t m1[cfg->M_RB_pucch4];
  switch (format)
  {
      case format_1:
      case format_1a:
      case format_1b:
           m = cfg->N_RB_2;
           // Find CP
           uint32_t c=(CP)?3:2;
           // Calculate m
           if (n_pucch >= c*cfg->N_cs_1 / cfg->delta_ss)
           {
             m = (n_pucch - c*cfg->N_cs_1 / cfg->delta_ss) / (c*N_sc/cfg->delta_ss) + cfg->N_RB_2+(uint32_t)ceilf((float) cfg->N_cs_1/8);
           }
           break;
     case format_2:
     case format_2a:
     case format_2b:
          m = n_pucch/N_sc;
          break;
     case format_3:
    	  m = n_pucch / N_SF0_PUCCH;
    	  break;
     case format_4:
         for( i =0;i < cfg->M_RB_pucch4; i++)
         {
            m1[i] = cfg->n_pucch_4 + i;
         }
            m = m1[i];
         break;
     case format_5:
    	  m = n_pucch;
    	  break;
  }
  return m;
}

/****************************************************************************/
// Mapping dmrs to physical resources according to Section 5.4.3 of 36.211
/****************************************************************************/
int refsignal_dmrs_pucch_map(struct pucch_config *cfg, struct cell *cell,uint32_t format, uint32_t n_pucch, float complex *r_pucch, float complex *output)
{
  int ret = ERROR_INVALID_INPUTS;
  if (cfg && output && r_pucch)
  {
    ret = ERROR;
    uint32_t nsymbols;
    nsymbols = cfg->CP?3:2;
    uint32_t nslot;
    uint32_t i;
    uint32_t n_rs;

    // Determine m
    uint32_t m = get_pucch_mvalue(cfg, format, n_pucch, cfg->CP);
    uint32_t n_PRB ;
    n_rs=get_N_rs_PUCCH(format, cfg->CP);
    uint32_t l[n_rs];
    for (nslot = 0 ; nslot < 2 ; nslot++)
    {
      // Determine n_prb as in 36.211 version 13.0 section 5.4.3
      if (format != format_4)
      {
          n_PRB = m/2;
          if ((m + nslot) % 2)
          {
             n_PRB = cfg->N_UL_RB - 1 - (m/2);
          }
      }
      else
      {
          if(nslot%2)
          {
              n_PRB = cfg->N_UL_RB - 1 - m;
          }
          else
          {
              n_PRB = m;
          }
      }

      for (i = 0;i < n_rs; i++)
      {
        get_pucch_dmrs_symbol(format, cfg->CP,l);
        memcpy(&output[RE_IDX(cfg->N_UL_RB, l[i] + nslot *nsymbols, n_PRB * N_sc)],&r_pucch[nslot * n_rs * N_sc + i * N_sc],N_sc*sizeof(float complex));
      }
    }

    ret = SUCCESS;
  }
  return ret;
}
