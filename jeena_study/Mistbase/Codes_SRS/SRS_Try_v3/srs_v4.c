
/*
 *SRS_Gen.c
 *
 *  Created on: May12, 2016
 *      Author: Jeena
 */
// Calculate prime
//Calculate largest prime less than M_sc: N_zc values
// Calculate Sequence Hopping
//Calculate Group Hopping

// Generate phi sequence
//Calculate q value
// calculate ZC sequence
//Calculate cyclic shift, alpha
//Calculate N_SRS

//Find // Base sequence =Nsc and 2Nsc
// Base sequence 3Nsc greater
//Calculate SRS seq
// find PUCCH seq
// place SRS in subframe
// TDD configuration





// srs-config - bw-cfg, B, n_srs, sf_idx, delta_ss, group_hopping,sequence_hopping,uint8_t n_ul_rb=6;
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
#include "scramble_util.h"
/************************************************************************************************************************/
/* //Table 5.2.3-1: Resource block parameters according to 3GPP 36.211 5.2.3*/
uint32_t N_sc=12;
uint32_t N_symbol_ul[2]={7,6};
/* Table 5.5.3.3-1: Frame structure type 1 sounding reference signal subframe configuration. */
uint32_t T_sfc[15] = {1, 2, 2, 5, 5, 5, 5, 5, 5, 10, 10, 10, 10, 10, 10};   /*Configuration Period*/
uint32_t Delta_sfc1[7] = {0, 0, 1, 0, 1, 2, 3}; /*Transmission offset*/
uint32_t Delta_sfc2[4] = {0, 1, 2, 3};


/************************************************************************************************************************/

//Calculate the value of n_prb
/*
 * Calculate the value of number of physical RBs n_prb
 * input- N_sc , n_ul_rb for normal and extended CP
 * output- N_PRB for extended and normal CP;
 */
uint32_t N_prb(uint32_t *n_prb, uint32_t n_ul_rb, uint32_t N_sc)
{
    uint32_t i;
    uint32_t k;
    /* resource elements in frequency domain 0.....n_ul_rb*N_sc*/
    k = n_ul_rb * N_sc;

    n_prb=floor(k/N_sc);
}



//Calculate the value of n_REs
/*
 * Calculate the value of number of resource elements
 * input- N_sc , N_symbol_ul for normal and extended CP
 * output -N_PRB for extended and normal CP;
*/
uint32_t SRS_NRE(){
    uint32_t n_RE[2];
if (srsLTE.CP='Normal'){
    n_RE[0]=N_sc*N_symbol_ul[0];/*Normal CP-(12*7)*/

}
else {
    n_RE[1]=N_sc*N_symbol_ul[1];/*Extended CP-(12*6)*/
}
return n_RE;
}
/**************************************************************************************************/
/*
/*
 * Calculate the v value of for M_sc
 * input-  NULRB, N_Sc subcarriers, Hopping
 * output -v[],M_sc[] for 1<=m<=5 and 6<=m<=NULRB,
 */

void Get_Msc_values(const char *Hopping, uint32_t *v, uint32_t *M_sc, uint32_t n_ul_rb, uint32_t N_sc)
{
    int m;
    for (m = 1; m <= n_ul_rb; m++)
    {
        M_sc[m-1] = m*N_sc;

        if (m>=1 && m<=5)
        {
            v[m-1] = 0;
        }
        else if (m>=6&&m<=n_ul_rb)
        {
            if (strcmp(Hopping,"Off") ==0)
            {
                v[m-1] = 0;
            }
            else
            {
                v[m-1]= 1;
            }
        }

    }
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
    for(i=1;i<M_sc;i++)
    {
        for(j=2;j<i;j++) if(i%j==0) break;
        if(j==i)
            prime=j;/* find prime numbers*/
    }
    Largest_prime=prime;
    /*printf("Largest primes %d\n",Largest_prime);*/
    return Largest_prime;
}
/********************************************************************************************************************/
/*
 * Calculate the value of largest prime number less than M_sc
 * input- M_sc for 1<=m<=5 and 6<=m<=NULRB, NULRB, N_Sc subcarriers
 * output -N_zc largest prime number less than M_sc
 */
uint32_t Get_Nzc( uint32_t n_ul_rb, uint32_t M_sc)
{

    return Prime_number(M_sc));

}

/**********************************************************************************************************************/

/*
 * Calculate q value
 * input- M_sc ,N_sc=no.ofsubcarriers,Grp_no,Seq_no
 * output -q value to find Zadoff chu seq
 */
 static float get_qvalue(uint32_t Grp_no, uint32_t Seq_no, uint32_t N_zc) {
  float q;
  float q_bar;
  float n_zc = (float) N_zc;
  q_bar= n_zc*(Grp_no+1)/31;
  q=(floor(q_bar+0.5))+(Seq_no*(pow((-1),(floor(2*q_bar)))));
  return q;
}
/*************************************************************************************************/
// Calculate argument for Qth root Zadoff-Chu sequence according to 3GPP 36.211 5.5.1.1 to find R_ZC
/*
 * Calculate the value of qth root for ZC
 * input- q value , N_zc
 * output root_q for exponential calculation
 */
static void Root_q_arg(float *root_q_arg, uint32_t M_sc, uint32_t Grp_no, uint32_t Seq_no, uint32_t n_ul_rb)
{
    float m;
    uint32_t N_zc = Get_Nzc(n_ul_rb, M_sc);
    float n_zc = (float) N_zc;
    float q = get_qvalue(Grp_no,Seq_no,N_zc);
    for (m = 0; m < N_zc; m++)
    {
        /* argument of x_q(m) according to 3GPP 36.211 5.5.1.1*/
        root_q_arg[m]= (-1 * M_PI * q * m * (m + 1)) / n_zc;
    }
}

/***********************************************************************************************/
/*
 * Calculate the  argument value for exponential calculation of Base sequence for M_sc=N_sc
 * input-  Grp. No ,N_sc=no.ofsubcarriers
 * output -Base Sequence value for M_sc=N_sc
 */
static void Seq_Msc12_Exp(float *Seq_Nsc_exp, uint8_t Grp_No, uint32_t N_sc)
{
    int i;
    for (i=0;i<N_sc;i++)
    {
        Seq_Nsc_exp[i]= Phi_M_sc_12[Grp_No][i] * M_PI / 4;
    }
}
 /*
 * Calculate the argument value for exponential calculation of Base sequence for M_sc=2*N_sc
 * input- Grp. No ,N_sc=no.ofsubcarriers
 * output -Base Sequence value for M_sc=2*N_sc
 */
static void Seq_Msc24_Exp(float *Seq_2Nsc_exp, uint8_t Grp_No, uint32_t N_sc)
{

    int i;
    for (i=0;i<2*N_sc;i++)
    {
        Seq_2Nsc_exp[i]= Phi_M_sc_24[Grp_No][i] * M_PI / 4));

    }
}
/***********************************************************************************************/

  // Calculate arguments foe exponential in r_uv(n)
/*
 * Calculate the value of arguments for N_PRB=1,2,<=3
 * input- N_PRB,Grp no, Seq No, other values needed
 * output -argument values to find sequence */
/* Computes argument of r_u_v signal */
static void compute_r_uv_arg(float *r_uv, uint32_t n_ul_rb, uint32_t ns, uint32_t N_sc, uint32_t N_zc, uint32_t seq_hopping, uint32_t group_hopping)
{
  // get u
  uint32_t Grp_no;

  Group_hopping_f_gh(&Grp_no, ns, const uint16_t cell_ID, group_hopping);/* fillllllll*/

  if (n_ul_rb == 1) {
    Seq_Msc12_Exp(r_uv, Grp_no);
  } else if (n_ul_rb == 2) {
    Seq_Msc24_Exp(r_uv, Grp_no);
  } else {
      float r_uv_temp[Nzc];
      //calculate sequence number v
      Root_q_arg(r_uv_temp, N_sc*n_ul_rb, Grp_no, Seq_no);

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

/***********************************************************************************************/


/*  May 16th*/
/***********************************************************************************************/
/*
 * Calculate Sequence Hopping(v)
 * input- M_sc,N_sc,
 * output-Calculate f_ss
 */
uint32_t Sequence_Hopping(const uint16_t cell_ID, uint32_t n_delta_ss){
for (uint32_t delta_ss=0;delta_ss<n_delta_ss;delta_ss++) {
     c_init = ((cell_ID / 30) << 5) + (((cell_ID % 30) + delta_ss) % 30);
     calc_prs_c( c_init, len, prs_c); /*generate_pseudo random sequence*/
     }
    v=prs_c[ns];
    }??????????????????????????????????????????????????????????
}
/**************************************************************************************************/
/*
 * Calculate Group Hopping(u)
 * input- M_sc ,N_sc=no.ofsubcarriers,Seq_no, N_UL_RB
 * output -
 */
 /* Calculate f_gh
u=(f_gh(ns)+f_ss)%30;*/
/** Computes sequence-group pattern f_gh according to 5.5.1.3 of 36.211 */
int Group_hopping_f_gh(uint32_t *Grp_no, const uint16_t cell_ID, uint32_t ns, uint32_t group_hopping)
{
    uint16_t n_ID= cell_ID;
    uint32_t len=(8*ns+7);
    uint8_t prs_c[len];
    uint32_t f_ss=n_ID%30;
    uint32_t c_init = floor(cell_ID/30);
    /*generate_pseudo random sequence /** Computes n_prs values as defined in 5.5.2.1.1 of 36.211 */
    calc_prs_c(c_init, len, prs_c);
    f_gh[ns] = 0;
    if(group_hopping)
    {
        for (int i = 0; i < 8; i++)
        {
            f_gh[ns] += (((uint32_t) prs_c[8 * ns + i]) << i)%30;/*

                         i=7
            summation of [c( 8*ns+ i )⋅ 2^i ]mod 30 if group hopping is enabled
                         i=0                                                      */
        }
    }
    else
    {
        f_gh[ns] = 0;
    }
    Grp_no=(f_gh[ns]+f_ss)%30;
}
/**************************************************************************************************/











































/**************************************************************************************************/
uint32_t srsbwtable_idx(uint32_t n_ul_rb) {
  if (n_ul_rb <= 40) {
    return 0;
  } else if (n_ul_rb <= 60) {
    return 1;
  } else if (n_ul_rb <= 80) {
    return 2;
  } else {
    return 3;
  }
}
uint32_t srs_M_sc()????????? {
  return m_srs_b[srsbwtable_idx(n_ul_rb)][.B][.bw_cfg]*2/2;/*SRSLTE_NRE=2 for Ntx=1*/
}


















/***********************************************************************************************/

struct {
    cf_t r_srs[20][100];
} srs_sequence;

/* Generate SRS signal as defined in Section 5.5.3.1 */
int srs_gen(struct srs_sequence *seq, uint32_t N_sc, uint32_t N_zc, sf_idx, uint32_t group_hopping, uint32_t sequence_Hopping)

    uint32_t M_sc = srs_M_sc(); ????????????

    for (uint32_t ns = 2*sf_idx; ns < 2*(sf_idx+1); ns++)
    {
        float alpha = 2*M_PI* n_srs/8;

        compute_r_uv_arg(&r_uv, n_ul_rb, Grp_no, Seq_no, N_sc, N_zc, ns);

        // Do complex exponential and adjust amplitude
        for (int n=0;n<M_sc;n++)
        {
            seq->r_srs[ns][n] = cexpf(I*(r_uv[n] + alpha*n));
        }
    }

}
