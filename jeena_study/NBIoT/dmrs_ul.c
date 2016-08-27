/*
 * dmrs_ul.c according to CR 0224 36.211v13.1.0
 *
 *  Created on: Jun 21, 2016
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

#define SUCC              		0
#define ERR                	   -1
#define ERROR_INVALID_INPUTS   -2
#define CP_NSYMB                7
#define N_sc_RB                12


struct NPUSCH {
           uint32_t sequence_hopping;/* enable=1 and disbale=0*/
           uint32_t group_hopping;/*/* enable=1 and disbale=0*/
           uint32_t groupAssignmentNPUSCH; /*//delta_ss = { 0 , 1 ,..., 29 }*/
           uint32_t ns;/*//Slot number within a radio frame*/
           uint32_t npusch_AllSymbols;/*/true -1 false -0*/
           uint32_t srsSubframeConfig;
           uint32_t I_RU;
           uint32_t N_CellID;
           uint32_t tone;
           uint32_t Subcarrierspacing;
           uint32_t threeTone_BaseSequence;
           uint32_t sixTone_BaseSequence;
           uint32_t twelveTone_BaseSequence;
           uint32_t NPUSCHformat; /*/ 0 if format 1, 1 if format 2*/
           uint32_t BaseSeq_set; /* 0 if not threeTone_BaseSequence,sixTone_BaseSequence, or twelveTone_BaseSequence not set, 1 if set*/
           uint32_t threeToneCyclicShift;/*/ {0,1,2}*/
           uint32_t sixToneCyclicShift;/*/ {0,1,,2,3}*/
           uint32_t Beta;
           uint32_t nf;
           uint32_t M_Rep_NPUSCH;
           uint32_t I_Rep;
           uint32_t ModScheme ;/*Modscheme = 0 if QPSK and 1 if BPSK*/
           };
void scramble_seq_gen(const uint32_t c_init, const uint32_t len,
                      uint8_t *prs_c)
{
    uint32_t N_c = 1600;
    uint32_t x1, x2;

    x1 = 0x54D21B24;            /* x1 sequence moved forward */
    x2 = c_init;

    uint32_t next_x1 = 0;
    uint32_t next_x2 = 0;

    /* Move x2 forward */
    uint32_t i;

    for (i = 0; i < N_c - 31; i++) {
        next_x2 = ((x2 >> 3) ^ (x2 >> 2) ^ (x2 >> 1) ^ x2) & 0x1;

        x2 = (x2 >> 1) | (next_x2 << 30);
    }

    /* Calculate prs_c */
    for (i = 0; i < len; i++) {

        next_x1 = ((x1 >> 3) ^ x1) & 0x1;
        next_x2 = ((x2 >> 3) ^ (x2 >> 2) ^ (x2 >> 1) ^ x2) & 0x1;

        *(prs_c++) = next_x1 ^ next_x2;

        x1 = (x1 >> 1) | (next_x1 << 30);
        x2 = (x2 >> 1) | (next_x2 << 30);
    }

}
//Table 16.5.1.1-3: Number of repetitions (N_rep_NPUSCH) for NPUSCH.

uint32_t N_rep[8] = {1,2,4,8,16,32,64,128};
           /* Table 10.1.4.1.1-1: Definition of   w(n)*/
int w_n[16][16] = {{1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1},
              {1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1},
              {1,1,-1,-1,1,1,-1,-1,1,1,-1,-1,1,1,-1,-1},
              {1,-1,-1,1,1,-1,-1,1,1,-1,-1,1,1,-1,-1,1},
              {1,1,1,1,-1,-1,-1,-1,1,1,1,1,-1,-1,-1,-1},
              {1,-1,1,-1,-1,1,-1,1,1,-1,1,-1,-1,1,-1,1},
              {1,1,-1,-1,-1,-1,1,1,1,1,-1,-1,-1,-1,1,1},
              {1,-1,-1,1,-1,1,1,-1,1,-1,-1,1,-1,1,1,-1},
              {1,1,1,1,1,1,1,1,-1,-1,-1,-1,-1,-1,-1,-1},
              {1,-1,1,-1,1,-1,1,-1,-1,1,-1,1,-1,1,-1,1},
              {1,1,-1,-1,1,1,-1,-1,-1,-1,1,1,-1,-1,1,1},
              {1,-1,-1,1,1,-1,-1,1,-1,1,1,-1,-1,1,1,-1},
              {1,1,1,1,-1,-1,-1,-1,-1,-1,-1,-1,1,1,1,1},
              {1,-1,1,-1,-1,1,-1,1,-1,1,-1,1,1,-1,1,-1},
              {1,1,-1,-1,-1,-1,1,1,-1,-1,1,1,1,1,-1,-1},
              {1,-1,-1,1,-1,1,1,-1,-1,1,1,-1,1,-1,-1,1}};
/* Orthogonal sequences for NPUSCH formats 2. Table 5.5.2.2.1-2*/
float w_arg_pusch_format2[3][3] = {{1, 1, 1},
                                   {1, 2*M_PI/3, 4*M_PI/3},
                                   {1, 4*M_PI/3, 2*M_PI/3}};
/**********************************NPUSCH ALPHA********************************************/
float alpha_3[3]= {0,2*M_PI/3, 4*M_PI/3};
float alpha_6[4]= {0,2*M_PI/6, 4*M_PI/6, 8*M_PI/6};
/*Generating  phi(n ) for N_sc_RU = 3 .
according to 3GPP TS 36.211 version 13.1.0 Release 13 Table 10.1.4.1.2-1*/
int  Phi_3[12][3]={{1,-3,-3},
                   {1,-3,-1},
                   {1,-3,3},
                   {1,-1,-1},
                   {1,-1,1},
                   {1,-1,3},
                   {1,1,-3},
                   {1,1-1},
                   {1,1,3},
                   {1,3,-1},
                   {1,3,1},
                   {1,3,3}};
/*Generating  phi(n ) for N_sc_RU = 6 .
according to 3GPP TS 36.211 version 13.1.0 Release 13 Table 10.1.4.1.2-2*/
int  Phi_6[14][6]= {{1,1,1,1,3,-3},
                    {1,1,3,1,-3,3},
                    {1,-1,-1,-1,1,-3},
                    {1,-1,3,-3,-1,-1},
                    {1,3,1,-1,-1,3},
                    {1,-3,-3,1,3,1},
                    {-1,-1,1,-3,-3,-1},
                    {-1,-1,-1,3,-3,-1},
                    {3,-1,1,-3,-3,3},
                    {3,-1,3,-3,-1,1},
                    {3,-3,3,-1,3,3},
                    {-3,1,3,1,-3,-1},
                    {-3,1,-3,3,-3,-1},
                    {-3,3,-3,1,1,-3}};
/*Generating  phi(n ) for N_sc_RU = N_sc_RB .
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
/******************************************************************************/
/* Table 10.1.4.2-1: Demodulation reference signal location for different NPUSCH formats. CR 36.211 */
/****************************************************************************/
uint32_t npusch_dmrs_symbol_format1[2] = {3,4};
uint32_t  npusch_dmrs_symbol_format2[6] = {0,1,2,2,3,4};
/******************************************************************************/
/* //Get Number of resource units N_RU
* input - I_RU
* output- N_RU values
*/
/*===========================================================================*/
uint32_t Get_N_RU(uint32_t I_RU)
{
	uint32_t N_RU;
	/*Table 16.5.1.1-2: Number of resource units (N_RU) for NPUSCH.*/
	if (I_RU < 6)
    		N_RU = I_RU +1;
	else if (I_RU == 6)
    		N_RU = 8;
	else
    		N_RU = 10;
	return N_RU;
}
/*===========================================================================*/
/* //Get N_sc_RU
* input - Subcarrierspacing = 0 if subcarrierspacing is 15kHz 1 if subcarrierspacing is 3.75kHz, tone - 0-singletone,1-multitone
* output- N_sc_RU and N_UL_slots values
*/
/*===========================================================================*/
uint32_t get_N_sc_RU(uint32_t Subcarrierspacing,uint32_t NPUSCHformat,uint32_t tone,uint32_t *N_UL_slots)
{
	uint32_t N_sc_RU=0 ;
	switch (NPUSCHformat){ /*NPUSCH case*/
	case 0:
        if(tone == 1){
        	N_sc_RU = 1;
        	*N_UL_slots =  16;
        }
        else if((Subcarrierspacing == 0) && (tone == 2)){
            N_sc_RU = 3;
            *N_UL_slots =  8;
        }
        else if((Subcarrierspacing == 0) && (tone == 3)){
           N_sc_RU = 6;
           *N_UL_slots =  4;
        }
        else if((Subcarrierspacing == 0) && (tone == 4)){
            N_sc_RU = 12;
            *N_UL_slots =  2;
        }
        else{
          if((Subcarrierspacing == 1) && (tone == 1)){
             N_sc_RU = 1;
             *N_UL_slots =  16;
          }
        }
        break;
    case 1:
             N_sc_RU = 1;
             *N_UL_slots =  4;
             break;
	}
	return N_sc_RU ;
}
/******************************************************************************/
/* //Get Nslots
* input - delta_f = 0 if subcarrierspacing is 15kHz 1 if subcarrierspacing is 3.75kHz
* output- M_identical values
*/
/****************************************************************************/
uint32_t Get_Nslots(uint32_t delta_f)
{
    uint32_t Nslots;
    Nslots = delta_f ? 1:2;
    return Nslots;
}
/******************************************************************************/
/* //Get N_seq_RU
* input - N_sc_RU
* output- N_seq_RU values
*/
/****************************************************************************/
uint32_t Get_N_seq_RU(uint32_t N_sc_RU)
{
	uint32_t N_seq_RU;
    if (N_sc_RU)
    	N_seq_RU = 16;
    else if (N_sc_RU == 3)
    	N_seq_RU = 12;
    else if (N_sc_RU == 6)
        N_seq_RU = 14;
    else if (N_sc_RU == 12)
        N_seq_RU = 30;
    return N_seq_RU;
}
/******************************************************************************/
/* //Compute M_identical
* input - N_sc_RU, M_rep_NPUSCH
* output- M_identical values
*/
/****************************************************************************/
uint32_t get_M_identical_NPUSCH(uint32_t N_sc_RU,uint32_t M_Rep_NPUSCH)
{
    uint32_t temp,M_identical;
    if(N_sc_RU > 1){
        temp = M_Rep_NPUSCH/2;
        if (temp < 4){
           M_identical = temp;
        }
        else{
            M_identical = 4;
    	}
    }
    else{
        M_identical = 1;
    }
	return M_identical;
}

/*****************************************************************************/
/* Calculate f_gh
 * input- slot number, group hopping -enabled/disabled, cellID
 * output - f_gh value */
/******************************************************************************/
uint32_t get_f_gh_NPUSCHformat1(struct NPUSCH *npusch,uint32_t N_seq_RU,uint32_t ns,uint32_t N_sc_RU)
{
	uint32_t len;
    uint32_t c_init;
    uint32_t count = 0;
    uint32_t f_gh;
    uint32_t ns_prime;
    c_init = (npusch->N_CellID/ N_seq_RU );
    if (N_sc_RU){
    	if (ns == 0){
            ns_prime = 0;
            len = ( 8 *ns_prime + 8 );// Maximum length
            uint8_t n_prs[len];
            /*generate_pseudo random sequence as defined in 5.5.2.1.1 of 36.211 */
            scramble_seq_gen(c_init, len, n_prs);
    	}
    }
    else if (N_sc_RU > 1){
    	if ((ns%2) == 0){
    	    len = (8*ns+8);// Maximum length
    	    uint8_t n_prs[len];
    	    /*generate_pseudo random sequence as defined in 5.5.2.1.1 of 36.211 */
    	    scramble_seq_gen(c_init, len, n_prs);
    	}
    }
    /** Computes sequence-group pattern f_gh according to 5.5.1.3 of 36.211 */
    f_gh = 0;
    int i;
    uint32_t n_prs[len];
	if(npusch->group_hopping){
    	for (i = 0; i < 8; i++){
        	count += (((uint32_t) n_prs[8 * ns_prime + i]) << i);
    	}
    	f_gh = count % N_seq_RU;
    	return f_gh;
        /*           i=7
       summation of [c( 8*ns+ i )⋅ 2^i ]mod N_seq_RU,if  group hopping is enabled
                    i=0                                                 */
	}
	else{
    	f_gh = 0;
    	return f_gh;
	}
}
/******************************************************************************/
/* * Calculate f_ss  * input-  cellID  * output - f_ss value */
 /*****************************************************************************/
uint32_t get_f_ss_NPUSCHformat1(struct NPUSCH *npusch,uint32_t N_seq_RU)
{
	uint32_t f_ss;
    f_ss = (npusch->N_CellID+ npusch->groupAssignmentNPUSCH) % N_seq_RU;
    return f_ss;
}
/******************************************************************************/
/* Group Hopping*/
/*  Calculate Group Hopping(u) u=(f_gh(ns)+f_ss)%30;
 * input- slot number, group hopping -enabled/disabled, cellID,NPUSCHformat,N_seq_RU
 * output - u value */
/******************************************************************************/
uint32_t Get_uValue_NPUSCHformat1(struct NPUSCH *npusch,uint32_t ns, uint32_t N_seq_RU)
{
        uint32_t u = 0;
        uint32_t N_UL_slots;
        uint32_t f_ss = get_f_ss_NPUSCHformat1(npusch,N_seq_RU);
        uint32_t N_sc_RU = get_N_sc_RU(npusch->Subcarrierspacing,npusch->NPUSCHformat,npusch->tone,&N_UL_slots);
        uint32_t f_gh = get_f_gh_NPUSCHformat1(npusch,N_seq_RU,ns,N_sc_RU );
        if (npusch->group_hopping)
        {
            if (npusch->NPUSCHformat==0)
            {
                u = (f_gh+ f_ss)% N_seq_RU;
            }
        }
        return u;
}
/******************************************************************************/
uint32_t Get_uvalue_ifnotsignalled(struct NPUSCH *npusch,uint32_t N_sc_RU)
{
     uint32_t u = 0;
     if(npusch->BaseSeq_set)
     {
         switch(N_sc_RU)
         {
           case 3:
                u = npusch->N_CellID % 12;
                break;
           case 6:
                u = npusch->N_CellID % 14;
                break;
           case 12:
                u = npusch->N_CellID % 30;
                break;
         }
     }
	 else{
     	switch(N_sc_RU) {
     	case 3:
        	u = npusch->threeTone_BaseSequence;
         	break;
     	case 6:
         	u = npusch->sixTone_BaseSequence;
            break;
     	case 12:
        	u = npusch->twelveTone_BaseSequence;
        	break;
         }
     }
	return u;
}
/******************************************************************************/
/* Demodulation reference signal */
/******************************************************************************/
void Get_phi_seq3(uint32_t N_sc_RU,uint32_t u, float complex *Seq_3_exp)
{
    uint32_t i;
    for (i = 0 ;i < N_sc_RU; i++)
    {
        Seq_3_exp[i] = cexpf(I*(Phi_3[u][i]*M_PI/4));
    }
}
/******************************************************************************/
void Get_phi_seq6(uint32_t N_sc_RU,uint32_t u, float complex *Seq_6_exp)
{
    uint32_t i;
    for (i = 0; i < N_sc_RU; i++)
    {
        Seq_6_exp[i] = cexpf(I*(Phi_6[u][i]*M_PI/4));
    }
}
/******************************************************************************/
void Get_phi_seq12(float complex *Seq_12_exp, uint32_t u, uint32_t N_sc_RU)
{
    uint32_t i;
    for (i = 0; i < N_sc_RU; i++)
    {
        Seq_12_exp[i] = cexpf(I*(Phi_12[u][i] * M_PI / 4));
    }
}
/******************************************************************************/
float get_alphavalue(struct NPUSCH *npusch,uint32_t N_sc_RU)
{
    uint32_t p;
    float alpha = -1;
    if ( N_sc_RU == 3 )
    {
    	 p = npusch->threeToneCyclicShift ;
    	 alpha = alpha_3[p];
    }else if (N_sc_RU == 6){
        p =npusch->sixToneCyclicShift;/* {0,1,,2,3}*/
        alpha = alpha_6[p];
     }else if (N_sc_RU == 12){
    	alpha = 0;
      }
    return alpha;
}
/******************************************************************************/
void Get_phi_sequence(uint32_t N_sc_RU,uint32_t u, float complex *Phi_Seq)
{
    if (N_sc_RU == 3)
    	Get_phi_seq3( N_sc_RU,u,Phi_Seq);
    else if (N_sc_RU == 6)
    	Get_phi_seq6( N_sc_RU,u,Phi_Seq);
    else
    	Get_phi_seq12(Phi_Seq,u,N_sc_RU);
}
/******************************************************************************/
/* Reference signal sequence for  N_sc_RU = 1*/
/******************************************************************************/
void Get_refsignal_r_u1(struct NPUSCH *npusch,uint32_t N_UL_slots,uint32_t I_Rep,uint32_t N_RU,uint32_t ns,uint32_t N_seq_RU,_Complex *r_u_bar)
{
    uint32_t n;
    uint32_t c_init = 35;
    uint32_t N_rep_NPUSCH;
    uint32_t u;
    N_rep_NPUSCH = N_rep[I_Rep];
    uint32_t len = N_rep_NPUSCH*N_UL_slots*N_RU ;
    uint8_t n_prs[len];
    float temp,temp1,temp2,temp3;
    u = npusch->N_CellID % 16;
    if (npusch->group_hopping)
    	u = Get_uValue_NPUSCHformat1(npusch,ns,N_seq_RU);
    temp1=0;
    for (n = 0; n < N_rep_NPUSCH*N_UL_slots*N_RU ; n++)
    {
        scramble_seq_gen( c_init, len, n_prs); /*generate_pseudo random sequence*/
    	temp = 1/sqrt(2);
    	temp3= 2*n_prs[n];
    	temp1= 1-temp3;
    	temp2 = w_n[u][n%16];
        r_u_bar[n] = (temp)*(1+I)*temp1*temp2;
    }
}
/******************************************************************************/
/* Reference signal sequence for  NPUSCHformat1 for N_sc_RU = 1*/
/******************************************************************************/
void Get_refsignal_NPUSCHformat1(struct NPUSCH *npusch,uint32_t ns,uint32_t N_UL_slots,uint32_t I_Rep,uint32_t N_RU,uint32_t N_seq_RU,_Complex *r_u_format1)
{
	uint32_t N_rep_NPUSCH,n;
    N_rep_NPUSCH = N_rep[I_Rep];
    _Complex r_u_bar[N_rep_NPUSCH*N_UL_slots*N_RU];
    Get_refsignal_r_u1(npusch,N_UL_slots,I_Rep,N_RU,ns,N_seq_RU,r_u_bar);
    if(npusch->NPUSCHformat==0){
    	for (n = 0; n < N_rep_NPUSCH*N_UL_slots*N_RU ; n++){
    		r_u_format1[n] = r_u_bar[n] ;
    	}
    }
}
/******************************************************************************/
/*Sequence index n_oc for  NPUSCHformat2 for N_sc_RU = 1*/
/******************************************************************************/
uint32_t Get_n_oc(uint32_t N_CellID,uint32_t ns,uint32_t *n_oc)
{
    int i;
    uint32_t c_init = N_CellID;
    uint32_t len = 8 * ns + 8;
    uint8_t n_prs[len];
    uint32_t temp;
    temp = 0;
    scramble_seq_gen( c_init, len, n_prs); /*generate_pseudo random sequence*/
    	  for (i = 0;i < 8; i++)
	      {
	        temp += (((uint32_t) n_prs[8 * ns + i]) << i);
	      }
    *n_oc = temp % 3;
    return 0;
}
/******************************************************************************/
/* Reference signal sequence for  NPUSCHformat2 for N_sc_RU = 1*/
/******************************************************************************/
void Get_refsignal_NPUSCHformat2(struct NPUSCH *npusch, uint32_t ns,uint32_t N_UL_slots,uint32_t N_RU,uint32_t N_seq_RU,uint32_t n_oc,uint32_t I_Rep,_Complex *r_u_format2)
{
    uint32_t n;
    int m;
    uint32_t N_rep_NPUSCH;
    float w_m;
    N_rep_NPUSCH = N_rep[I_Rep];
    _Complex r_u_bar[N_rep_NPUSCH*N_UL_slots*N_RU];
    Get_refsignal_r_u1(npusch,N_UL_slots,I_Rep,N_RU,ns,N_seq_RU,r_u_bar);
    if(npusch->NPUSCHformat){
    	for (m = 0; m < 3; m++){
    		w_m = w_arg_pusch_format2[n_oc][m];
    		for (n = 0; n < N_rep_NPUSCH*N_UL_slots*N_RU ; n++){
    			r_u_format2[3*n+m] = w_m*r_u_bar[n];
    		}
    	}
    }
}
/******************************************************************************/
/* Reference signal sequence for N_sc_RU > 1*/
/******************************************************************************/
void Get_refsignal_r_u_greaterthan1(struct NPUSCH *npusch,uint32_t N_sc_RU,uint32_t N_seq_RU,uint32_t ns,_Complex *r_u)
{
     float alpha;
     uint32_t n;
     uint32_t u;
     if(npusch->group_hopping)
     {
         u = Get_uValue_NPUSCHformat1(npusch,ns,N_seq_RU);
     }
     else
     {
         u = Get_uvalue_ifnotsignalled(npusch,N_sc_RU);
     }
     float complex Phi_Seq[N_sc_RU];
     for (n = 0; n < N_sc_RU ; n++)
     {
    	 Get_phi_sequence(N_sc_RU,u,Phi_Seq);
    	 alpha = get_alphavalue(npusch,N_sc_RU);
         r_u[n] = expf(I*alpha*n)* Phi_Seq[n];
     }
}
/******************************************************************************/
/* Table 10.1.4.2-1: Demodulation reference signal location for different NPUSCH formats. CR 36.211 */
/****************************************************************************/
uint32_t get_npusch_dmrs_symbollocation(uint32_t NPUSCHformat, uint32_t Subcarrierspacing, uint32_t *loc)
{
  switch (NPUSCHformat){    /*/NPUSCHformat 2*/
  case 0:
	  loc[0] = npusch_dmrs_symbol_format1[Subcarrierspacing];
	  break;
  case 1:
	  if (Subcarrierspacing) {                   /* 3.75kHz subcarrier spacing*/
	      loc[0] = npusch_dmrs_symbol_format2[0];/*location of ref symbols in 3.75kHz subcarrier spacing*/
	      loc[1] = npusch_dmrs_symbol_format2[1];
	      loc[2] = npusch_dmrs_symbol_format2[2];
       }else{
	      loc[0] = npusch_dmrs_symbol_format2[3];/*location of ref symbols in 15kHz subcarrier spacing*/
	      loc[1] = npusch_dmrs_symbol_format2[4];
	      loc[2] = npusch_dmrs_symbol_format2[5];
        }
        break;
  }
  return 0;
}
/****************************************************************************/
/* Number of PUSCH demodulation reference symbols  */
/****************************************************************************/
uint32_t get_N_rs_PUSCH(uint32_t format)
 {
	uint32_t n_rs = 0;
    switch (format){
    case 0:/*format1*/
         n_rs =  1;
         break;
    case 1:/* format2*/
      	 n_rs = 3;
           break;
    }
   return n_rs;
}
/****************************************************************************/
// Mapping dmrs to physical resources according to Section 10.1.4.2 of CR 0224 36.211v13.1.0
/****************************************************************************/
int dmrs_npusch_map(struct NPUSCH *npusch,_Complex Symbolmap[12][1120*npusch->M_Rep_NPUSCH])
{
  int ret = ERROR_INVALID_INPUTS;
  if (npusch)
  {
    ret = ERR;
    uint32_t nsymbols;
    uint32_t N_UL_slots;
    uint32_t N_sc_RU;
    uint32_t N_seq_RU;
    uint32_t N_RU;
    uint32_t N_rep_NPUSCH;
    N_rep_NPUSCH = N_rep[npusch->I_Rep];
	N_sc_RU=get_N_sc_RU(npusch->Subcarrierspacing,npusch->NPUSCHformat,npusch->tone,&N_UL_slots);
	N_seq_RU = Get_N_seq_RU(N_sc_RU);
	N_RU = Get_N_RU(npusch->I_RU);
    /* find the number of dmrs symbols*/
    nsymbols = npusch->NPUSCHformat?3:1;
    // find the number of slots to repeat according to CR 36.211 10.1.3.6
    uint32_t i;
    uint32_t k;
    uint32_t l;
    uint32_t n;
    uint32_t n_oc;
    uint32_t loc[nsymbols];
    // find the symbol location according to table 10.1.4.2-1 in CR 36.211 rel13 v13.1.0
    get_npusch_dmrs_symbollocation(npusch->NPUSCHformat, npusch->Subcarrierspacing,loc);
    // find the number of repetitions for Nslots
    uint32_t M_identical = get_M_identical_NPUSCH(N_sc_RU,npusch->M_Rep_NPUSCH);
    _Complex r_npusch[N_rep_NPUSCH*N_UL_slots*N_RU*3];
    double SymbolMap_real[N_UL_slots*N_RU*npusch->M_Rep_NPUSCH*CP_NSYMB];
    double SymbolMap_imag[N_UL_slots*N_RU*npusch->M_Rep_NPUSCH*CP_NSYMB];
    Get_n_oc(npusch->N_CellID,npusch->ns,&n_oc);
	if(N_sc_RU == 1){
		if (npusch->NPUSCHformat==0)
		  Get_refsignal_NPUSCHformat1(npusch,npusch->ns,N_UL_slots,npusch->I_Rep,N_RU,N_seq_RU,r_npusch);
		else
			Get_refsignal_NPUSCHformat2(npusch,npusch->ns,N_UL_slots,N_RU,N_seq_RU,n_oc,npusch->I_Rep,r_npusch);
		for (n=0;n<(N_rep_NPUSCH*N_UL_slots*N_RU);n++){
		  for (l=0;l<nsymbols;l++){
		    for (k=0;k<N_sc_RU;k++){
		    	Symbolmap[k][loc[l]+n*CP_NSYMB] = r_npusch[k+n*N_sc_RU] ;
		    }
		  }
		}
	}
	else if (N_sc_RU >1){
		Get_refsignal_r_u_greaterthan1(npusch,N_sc_RU,N_seq_RU,npusch->ns,r_npusch);
		if (M_identical>1){
			  for (i=0;i<N_UL_slots*N_RU*npusch->M_Rep_NPUSCH;i++){
		     		for (l=0;l<nsymbols;l++){
		     			for (k=0;k<N_sc_RU;k++){
		                	Symbolmap[k][loc[l]+i*CP_NSYMB] =  r_npusch[k] ;
		                	//printf(" Symbolmapping[%d][%d] = %.4f+i%.4f \n",k,loc[l]+i*CP_NSYMB,creal(Symbolmap[k][loc[l]+i*CP_NSYMB]),cimag(Symbolmap[k][loc[l]+i*CP_NSYMB]));
		     			}

		     		}
		     }
		}
		else if (M_identical){
		  for (i=0;i<N_UL_slots*N_RU*CP_NSYMB;i++) {
		    for (l=0;l<nsymbols;l++){
		    	for (k=0;k<N_sc_RU;k++)		    	{
		    		Symbolmap[k][loc[l]+i*CP_NSYMB] = r_npusch[k];
			    	printf(" Symbolmapping[%d][%d] = %.4f+i%.4f \n",k,loc[l]+i*CP_NSYMB,creal(Symbolmap[k][loc[l]+i*CP_NSYMB]),cimag(Symbolmap[k][loc[l]+i*CP_NSYMB]));
		     	}		    	
		    }
		  }
		}
	}
  ret = SUCC;
  }
  return ret;
}
void run_dmrs_uplink(struct NPUSCH *npusch,double *r_npusch_real, double *r_npusch_imag)
{
	uint32_t N_sc_RU;
	uint32_t N_UL_slots;
	uint32_t N_RU;
	uint32_t loc[3];
	uint32_t n_oc;
	uint32_t N_seq_RU;
 	uint32_t N_rep_NPUSCH;
	int i;
	int len;
	int n;
	N_rep_NPUSCH = N_rep[npusch->I_Rep];
	N_sc_RU=get_N_sc_RU(npusch->Subcarrierspacing,npusch->NPUSCHformat,npusch->tone,&N_UL_slots);
	N_seq_RU = Get_N_seq_RU(N_sc_RU);
	N_RU = Get_N_RU(npusch->I_RU);
	len = (N_sc_RU == 1)?(N_rep_NPUSCH*N_UL_slots*N_RU):N_sc_RU;
	_Complex r_npusch[len*3];
	get_npusch_dmrs_symbollocation(npusch->NPUSCHformat,npusch->Subcarrierspacing,loc);
	Get_n_oc(npusch->N_CellID,npusch->ns,&n_oc);
	if(N_sc_RU == 1){
    	if (npusch->NPUSCHformat==0){
    		Get_refsignal_NPUSCHformat1(npusch, npusch->ns, N_UL_slots, npusch->I_Rep, N_RU, N_seq_RU, r_npusch);
    		for(i=0;i<len;i++){
            	r_npusch_real[i]= creal(r_npusch[i]);
                r_npusch_imag[i]= cimag(r_npusch[i]);
              }
            }
    	else{
    		Get_refsignal_NPUSCHformat2(npusch,npusch->ns,N_UL_slots,N_RU,N_seq_RU,n_oc,npusch->I_Rep,r_npusch);
    		for (n=0;n<len*3;n++){
             	r_npusch_real[n]= creal(r_npusch[n]);
            	r_npusch_imag[n]= cimag(r_npusch[n]);
        	}
    	}
	}
	else{
    	Get_refsignal_r_u_greaterthan1(npusch,N_sc_RU,N_seq_RU,npusch->ns,r_npusch);
    	for(i=0;i<len;i++){
        	r_npusch_real[i]= creal(r_npusch[i]);
        	r_npusch_imag[i]= cimag(r_npusch[i]);
    	}
   	}
	/*struct NPUSCH nbiot;
	nbiot.sequence_hopping=0;// enable=1 and disable=0
	nbiot.group_hopping=0;// enable=1 and disable=0
	nbiot.groupAssignmentNPUSCH=0; //delta_ss = { 0 , 1 ,..., 29 }
	nbiot.ns=0;//Slot number within a radio frame
	nbiot.npusch_AllSymbols=1;//true -1 false -0
	nbiot.srsSubframeConfig=0;
	nbiot.I_RU=5;
	nbiot.N_CellID=1;
	nbiot.tone=1;
	nbiot.Subcarrierspacing=0;
	nbiot.threeTone_BaseSequence=1;
	nbiot.sixTone_BaseSequence=1;
	nbiot.twelveTone_BaseSequence=1;
	nbiot.NPUSCHformat=1; // 0 if format 1, 1 if format 2
	nbiot.BaseSeq_set=1; // 0 if not threeTone_BaseSequence,sixTone_BaseSequence, or twelveTone_BaseSequence not set, 1 if set
	nbiot.threeToneCyclicShift=0;// {0,1,2}
	nbiot.sixToneCyclicShift=0;// {0,1,,2,3}
    nbiot.Beta=1;
    nbiot.nf=0;
    nbiot.M_Rep_NPUSCH=4;
    nbiot.I_Rep=2;
    nbiot.ModScheme =0 ;//ModScheme = 0 if QPSK and 1 if BPSK*/
}

void main()
{
	struct NPUSCH nbiot;
	nbiot.sequence_hopping=0;// enable=1 and disable=0
	nbiot.group_hopping=0;// enable=1 and disable=0
	nbiot.groupAssignmentNPUSCH=0; //delta_ss = { 0 , 1 ,..., 29 }
	nbiot.ns=0;//Slot number within a radio frame
	nbiot.npusch_AllSymbols=1;//true -1 false -0
	nbiot.srsSubframeConfig=0;
	nbiot.I_RU=5;
	nbiot.N_CellID=1;
	nbiot.tone=3;
	nbiot.Subcarrierspacing=0;
	nbiot.threeTone_BaseSequence=1;
	nbiot.sixTone_BaseSequence=1;
	nbiot.twelveTone_BaseSequence=1;
	nbiot.NPUSCHformat=0; // 0 if format 1, 1 if format 2
	nbiot.BaseSeq_set=1; // 0 if not threeTone_BaseSequence,sixTone_BaseSequence, or twelveTone_BaseSequence not set, 1 if set
	nbiot.threeToneCyclicShift=0;// {0,1,2}
	nbiot.sixToneCyclicShift=0;// {0,1,,2,3}
    nbiot.Beta=1;
    nbiot.nf=0;
    nbiot.M_Rep_NPUSCH=4;
    nbiot.I_Rep=2;
    nbiot.ModScheme =0 ;//ModScheme = 0 if QPSK and 1 if BPSK*/
	uint32_t N_sc_RU,N_UL_slots,N_RU,loc[3],n_oc,N_seq_RU,n;
	uint32_t N_rep_NPUSCH,i,l,k;
	N_rep_NPUSCH = N_rep[nbiot.I_Rep];
	N_sc_RU=get_N_sc_RU(nbiot.Subcarrierspacing,nbiot.NPUSCHformat,nbiot.tone,&N_UL_slots);
	N_seq_RU = Get_N_seq_RU(N_sc_RU);
	N_RU = Get_N_RU(nbiot.I_RU);
	_Complex r_npusch[N_rep_NPUSCH*N_UL_slots*N_RU*3];
	printf(" N_UL_slots = %d \nN_sc_RU = %d \n N_rep_NPUSCH = %d \n N_seq_RU= %d \n N_RU = %d \n",N_UL_slots,N_sc_RU,N_rep_NPUSCH,N_seq_RU,N_RU);
	get_npusch_dmrs_symbollocation(nbiot.NPUSCHformat,nbiot.Subcarrierspacing,loc);
	uint32_t nsymbols = get_N_rs_PUSCH(nbiot.NPUSCHformat);
	Get_n_oc(nbiot.N_CellID,nbiot.ns,&n_oc);
	//uint32_t N_Uplinkslots = N_UL_slots;
	printf(" nsymbols = %d \n n_oc= %d \nN_Uplinkslots = %d\n",nsymbols,n_oc,N_UL_slots);
	_Complex Symbolmap[N_sc_RU][N_UL_slots*N_RU*nbiot.M_Rep_NPUSCH*CP_NSYMB*CP_NSYMB];
	if(N_sc_RU == 1){
		if (nbiot.NPUSCHformat==0){
			Get_refsignal_NPUSCHformat1(&nbiot,nbiot.ns,N_UL_slots,nbiot.I_Rep,N_RU,N_seq_RU,r_npusch);
			uint32_t n;
				for (n = 0; n < N_rep_NPUSCH*N_UL_slots*N_RU; n++){
					//printf(" r_npusch[%d] = %.4f+i%.4f \n",n,creal(r_npusch[n]),cimag(r_npusch[n]));
				}
		}
		else{
			Get_refsignal_NPUSCHformat2(&nbiot,nbiot.ns,N_UL_slots,N_RU,N_seq_RU,n_oc,nbiot.I_Rep,r_npusch);
			for (n = 0; n < N_rep_NPUSCH*N_UL_slots*N_RU; n++){
				//printf(" r_npusch[%d] = %.4f+i%.4f \n",n,creal(r_npusch[n]),cimag(r_npusch[n]));
		    }
		}
	}
	else{
		Get_refsignal_r_u_greaterthan1(&nbiot,N_sc_RU,N_seq_RU,nbiot.ns,r_npusch);
		for (n = 0; n < N_sc_RU; n++){
			 //printf(" r_npusch[%d] = %.4f+i%.4f \n",n,creal(r_npusch[n]),cimag(r_npusch[n]));
	     }
	}
	dmrs_npusch_map(&nbiot,Symbolmap);
	for (i=0;i<N_UL_slots*N_RU*nbiot.M_Rep_NPUSCH*CP_NSYMB;i++){
		for (l=0;l<nsymbols;l++){
			for (k=0;k<N_sc_RU;k++){
				//printf(" Symbolmapping[%d][%d] = %.4f+i%.4f \n",k,loc[l]+i,creal(Symbolmap[k][loc[l]+i]),cimag(Symbolmap[k][loc[l]+i]));

			}
		}
	}
}
