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
struct SRS_UL{
           const uint16_t cell_ID;
           uint32_t  B;//B_srs={0,1,2,3} UE specific
           uint32_t bw_cfg; //C_srs is an element of {0,1,2,3,4,5,6,7}/ cell specific parameter
           uint32_t srsSubframeConfig;
           uint32_t n_ul_rb;// must be 6 or greater in value
           uint32_t CyclicPrefixLength;
           uint32_t CP;// Normal -0, Extended-1
           uint32_t Duplex_Mode ;// FDD-0,TDD-1, Half Duplex -2
           uint32_t HoppingBandwidth; //b_hop={0,1,2,3}
           uint32_t freqDomainPosition;// n_RRC
           uint32_t freqDomainPosition_ap;// n_RRC
           uint32_t nf; //System frame number SFN 
           uint32_t Config_idx;// I_srs {0,..... 644}
           uint32_t  Offset_Idx ; // takes 0 or 1 
           uint32_t K_Tc;//Transmission_comb =2 for SRS
           uint32_t Cyclic_shift;// n_srs_cs
           uint32_t Cyclic_shift_ap;// n_srs_cs
           uint32_t N_Tx; //{1,2,4}
           uint32_t ns;//Slot number within a radio frame
           uint32_t N_sp;////Number of downlink to uplink switch points within the radio frame = 5ms
           uint32_t sf_idx;
           uint32_t sequence_hopping;// enable=1 and disbale=0
           uint32_t group_hopping;// enable=1 and disbale=0
           uint32_t delta_ss; //delta_ss = { 0 , 1 ,..., 29 }
           uint32_t N_ID_PUCCH;// Configured=1 and Not Configured=0
           uint32_t N_ID_PUSCH;// Configured=1 and Not Configured=0
           uint32_t srsMaxUpPTS;// Cell specific
           uint16_t PUCCH_ID;  // Configured  0 or 1
           uint16_t PUSCH_ID;  // Configured  0 or 1
       };
	   
	   
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
/******************************************************************************/


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
/***********************************************************************************************/
uint32_t T_srs_FDD(uint32_t Config_idx) 
{
    uint32_t I_srs;
    I_srs = Config_idx;
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

/*Table 8.2-2: UE Specific SRS Periodicity T_SRS for trigger type 0, TDD according to 3GPP TS 36.213 version 13.0.0 Release 13*/ 

uint32_t T_srs_TDD(uint32_t Config_idx) 
{
    uint32_t I_srs = Config_idx;
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
/******************************************************************************/ 
/* /*Compute T_srs that counts the number of UE-specific SRS transmissions for FDD and TDD mode
* input - Config_idx, Duplex_Mode (FDD -0 and TDD-1)
* output- T_srs
*/
/******************************************************************************/
uint32_t T_srs_value(uint32_t Config_idx,uint32_t Duplex_Mode) 
{
	uint32_t T_srs;
        if (Duplex_Mode == 0) //  FDD mode
	{
	    T_srs = T_srs_FDD(Config_idx) ;
	}
	else 
        {
        if (Duplex_Mode == 1 )// TDD mode
	{
            T_srs = T_srs_TDD(Config_idx) ;
        }
        }
return T_srs;
}


/*Table 8.2-2 of Subframe Offset Configuration T_offset for trigger type 0, TDD according to 3GPP TS 36.213 version 13.0.0 Release 13*/ 
uint32_t Offset_value(uint32_t T_srs,uint32_t Duplex_Mode, uint32_t Config_idx, uint32_t Offset_Idx, uint32_t *T_offset_max)
{
        uint32_t I_srs = Config_idx;// I_srs is Config_Idx
        uint32_t T_offset;
        if (Duplex_Mode == 0) //  FDD mode
	{
	    T_offset = 0;
	}
	else if (Duplex_Mode == 1)
        {
	    if (T_srs == 2)
	    {
		      if ((I_srs == 0) &&  (Offset_Idx == 0))
                      {
                          T_offset   = 0; 
			  *T_offset_max = 1;
		      }
		      else if ((I_srs == 0) &&  (Offset_Idx == 1))
		      {
		          T_offset  = 1;
                          *T_offset_max = 1;
                      }
	              else if ((I_srs == 1) &&  (Offset_Idx == 0))
	              {
                  	  T_offset   = 0; 
			  *T_offset_max = 2;
                      }
                      else if ((I_srs == 1) &&  (Offset_Idx == 1))
                      {
			  T_offset = 2;
			  *T_offset_max = 2;
		      }
		      else if ((I_srs == 2) &&  (Offset_Idx == 0))
		      {
			  T_offset = 1;
			  *T_offset_max = 2;
		      }
		      else if ((I_srs == 2) &&  (Offset_Idx == 1))
		      {
			  T_offset = 2;
			  *T_offset_max = 2;
		      }
                      else if ((I_srs == 3) &&  (Offset_Idx == 0))
		      {
			  T_offset = 0;
			  *T_offset_max = 3;
		      }
		      else if ((I_srs == 3) &&  (Offset_Idx == 1))
		      {
			  T_offset = 3;
			  *T_offset_max = 3;
		      }
		      else if ((I_srs == 4) &&  (Offset_Idx == 0))
		      {
			  T_offset = 1;
			  *T_offset_max = 3;
		      }
		      else if ((I_srs == 4) &&  (Offset_Idx == 1))
		      {
			  T_offset = 3;
			  *T_offset_max = 3;
		      }
              	      else if ((I_srs == 5) &&  (Offset_Idx == 0))
                      {
                          T_offset   = 0; 
			  *T_offset_max = 4;
                      }
		      else if ((I_srs == 5) &&  (Offset_Idx == 1))
		      {
			  T_offset = 4;
			  *T_offset_max = 4;
		      }
                      else if ((I_srs == 6) &&  (Offset_Idx == 0))
                      {
                          T_offset   = 1;
			  *T_offset_max = 4;
		      }
		      else if ((I_srs == 6) &&  (Offset_Idx == 1))
		      {
			  T_offset   = 4;
			  *T_offset_max = 4;
		      }
                      else if ((I_srs == 7) &&  (Offset_Idx == 0))
                      {
                          T_offset   = 2;
			  *T_offset_max = 3;
		      }
		      else if ((I_srs == 7) &&  (Offset_Idx == 1))
                      {
			  T_offset   = 3;
			  *T_offset_max = 3;
		      }
                      else if ((I_srs == 8) &&  (Offset_Idx == 0))
                      {
                          T_offset   = 2;
			  *T_offset_max = 4;
		      }
		      else if ((I_srs == 8) &&  (Offset_Idx == 1))
		      {
			  T_offset   = 4;
			  *T_offset_max = 4;
		      }
                      else if ((I_srs == 9) &&  (Offset_Idx == 0))
                      {
                          T_offset   = 3;
			  *T_offset_max = 4;
		      }
		      else if ((I_srs == 9) &&  (Offset_Idx == 1))
		      {
			  T_offset   = 4;
			  *T_offset_max = 4;
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
		      *T_offset_max = T_offset;
	    }
        }
        else 
        {
            T_offset = 0;
        }

        return T_offset;
}
/******************************************************************************/ 
/* /*Compute n_srs that counts the number of UE-specific SRS transmissions
* input - Config_idx, ns,  nf, N_sp, T_srs, Duplex_Mode, T_offset,  T_offset_max
* output- n_srs
*/
/******************************************************************************/
uint32_t Get_n_srs(uint32_t Config_idx, uint32_t ns, uint32_t nf, uint32_t N_sp, uint32_t T_srs, uint32_t Duplex_Mode, uint32_t T_offset, uint32_t T_offset_max)
{
    uint32_t n_srs;
    if( T_srs == 2)
    {
        //System frame numberSFN SRS_UL.nf
        //Number of downlink to uplink switch points within the radio frame N_sp
        n_srs = ((2 * N_sp*nf) + (2* (N_sp- 1))* floor(ns / 10)) + (floor(T_offset / T_offset_max)); 
    }
    else
    {
        n_srs = floor(((nf * 10) + floor(ns / 2)) / T_srs);
    }
    return n_srs;
} 
 // checked


/******************************************************************************/ 
/* /*Compute Fb
* input - HoppingBandwidth,B_srs,freqDomainPosition,K_Tc,bw_cfg,n_ul_rb,n_srs,N_b
* output- Fb
*/
/******************************************************************************/

uint32_t Get_Fb(uint32_t n_ul_rb,uint32_t B,uint32_t bw_cfg, uint32_t HoppingBandwidth,uint32_t n_srs, uint32_t N_b)
{
    uint32_t b_hop = HoppingBandwidth;
    uint32_t prod_NbDash ;
    prod_NbDash = 1;
    uint32_t prod_2NbDash ;
    prod_2NbDash = 1;
    uint32_t bdash;
    uint32_t Fb;
    if(B > b_hop)
    {
       for (bdash = b_hop+1; bdash <= (B-1); bdash++) 
       {
           prod_NbDash *= Nb[srsbwtable_idx(n_ul_rb)][bdash][bw_cfg];
       }
       for (bdash = b_hop; bdash <= B; bdash++) 
       {
           prod_2NbDash *= Nb[srsbwtable_idx(n_ul_rb)][bdash][bw_cfg];
       }
       if (N_b % 2 ==0)
       {
           Fb = (N_b / 2) * (((n_srs % prod_2NbDash) / (prod_NbDash)) + ((n_srs % prod_2NbDash) / (2 * prod_NbDash)));
       }
       else 
       { 
           Fb =(floor(N_b / 2)) * floor(n_srs / prod_NbDash);
       }  
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
* output- n_b values
*/

uint32_t Get_nb(uint32_t HoppingBandwidth, uint32_t B, uint32_t n_ul_rb, uint32_t bw_cfg, uint32_t freqDomainPosition,uint32_t n_srs, uint32_t Fb, uint32_t N_b)
{
    uint32_t n_b;
    uint32_t n_temp;
    int b;
    uint32_t b_hop = HoppingBandwidth;
    if(b_hop < B)
    {
        for(b = 0; b < B ; b++) 
        {
           if(b <= b_hop)
           {
               n_temp = (floor ((4 *freqDomainPosition)/ m_srs_b[srsbwtable_idx(n_ul_rb)][b][bw_cfg])) ;
               printf("N_temp = %d\n",n_temp);
               n_b =n_temp % N_b;
           }
           else 
           {
               // freq hopping not enabled
               n_temp = (floor ((4 *freqDomainPosition)/ m_srs_b[srsbwtable_idx(n_ul_rb)][b][bw_cfg])) ;
               n_b = (Fb + n_temp) % N_b;
           }
       }
    }
    else
    {
       n_temp = (floor ((4 *freqDomainPosition)/ m_srs_b[srsbwtable_idx(n_ul_rb)][b][bw_cfg])) ;
       printf("N_tempB = %d\n",n_temp);
       n_b =n_temp % N_b;
    }
    return n_b;
}

/***********************************************************************************************/
uint32_t Get_K_Tc_p(uint32_t Cyclic_shift, uint32_t N_Tx, uint32_t K_Tc)
{
    uint32_t K_Tc_bar = K_Tc;
    int p_index;
    uint32_t K_Tc_p;//K_Tc_p={0,1,...SRS_UL.K_Tc-1}
    static const unsigned n_SRS_cs = (1 << 4) | (1 << 5) | (1 << 6) | (1 << 7);
    static const unsigned p_hat = (1 << 1) | (1 << 3) ;
    /*find K_Tc_p*/
    for (p_index = 0; p_index <=3; p_index++)
    {
        if (((1 << Cyclic_shift) & n_SRS_cs) && ((1 << p_index) & p_hat) && (N_Tx == 4) ) // p_hat ands p relation table Table 5.2.1-1
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
uint32_t get_k_0_pbar(uint32_t bw_cfg, uint32_t N_sc, uint32_t n_ul_rb ,uint32_t K_Tc_p)
{
    uint32_t k_0_pbar;
    //K_Tc_p={0,1,...SRS_UL.K_Tc-1}
    if (bw_cfg < 8)
    {
          k_0_pbar = (((floor(n_ul_rb / 2)) - (m_srs_b[srsbwtable_idx(n_ul_rb)][0][bw_cfg] / 2 ))*N_sc ) + K_Tc_p;
    }

}



/* Returns number of RB defined for the cell-specific SRS */
uint32_t srslte_refsignal_srs_rb_L_cs(uint32_t bw_cfg, uint32_t n_ul_rb) {
  if (bw_cfg < 8) {
    return m_srs_b[srsbwtable_idx(n_ul_rb)][0][bw_cfg];
	printf (" M_SRS = %d\n",m_srs_b[srsbwtable_idx(n_ul_rb)][0][bw_cfg]);
  } 
  return 0; 
}
/* /*Compute m_srs_0[bw_cfg]

* input - bw_cfg
* output- m_srs_0 values   **************************** not required just a general calculation*/
uint32_t Get_m_srs_0( uint32_t bw_cfg)
{
    int i;
    uint32_t c_srs = bw_cfg; 
    uint32_t m_srs_0[8];
    for (i=0 ; i < 8; i++)
    {
        m_srs_0 [i] = *((uint32_t*)m_srs_b+ i);
        
    }
printf (" M_SRS1= %d\n",m_srs_0[c_srs]);
return m_srs_0[c_srs];
}
void main()
{
           const uint16_t cell_ID = 100;
           uint32_t  B = 0;//B_srs={0,1,2,3} UE specific
           uint32_t bw_cfg = 5; //C_srs is an element of {0,1,2,3,4,5,6,7}/ cell specific parameter
           uint32_t srsSubframeConfig;
           uint32_t n_ul_rb = 6;// must be 6 or greater in value
           uint32_t CyclicPrefixLength;
           uint32_t CP = 0;// Normal -0, Extended-1
           uint32_t Duplex_Mode = 1;// FDD-0,TDD-1
           uint32_t HoppingBandwidth = 1; //b_hop={0,1,2,3} if b_hop is > B, frequency hopping enabled
           uint32_t freqDomainPosition = 0;// n_RRC
           uint32_t freqDomainPosition_ap;// n_RRC if enabled b_hop is disabled.
           uint32_t nf = 0; //System frame number SFN 
           uint32_t Config_idx = 7;// I_srs {0,..... 644}
           uint32_t Offset_Idx = 0;// {0,1}
           uint32_t K_Tc = 2;//Transmission_comb =2 for SRS
           uint32_t Cyclic_shift = 0;// n_srs_cs
           uint32_t Cyclic_shift_ap;// n_srs_cs
           uint32_t N_Tx = 1; //{1,2,4}
           uint32_t ns = 20;//Slot number within a radio frame
           uint32_t N_sp = 5;////Number of downlink to uplink switch points within the radio frame = 5ms
           uint32_t sf_idx = 0;
           uint32_t sequence_hopping = 0;// enable=1 and disbale=0
           uint32_t group_hopping = 1;// enable=1 and disbale=0
           uint32_t delta_ss = 0; //delta_ss = { 0 , 1 ,..., 29 }
           uint32_t N_ID_PUCCH = 0;// Configured=1 and Not Configured=0
           uint32_t N_ID_PUSCH = 0;// Configured=1 and Not Configured=0
           uint32_t srsMaxUpPTS;// Cell specific
	   uint16_t PUCCH_ID = 499;  // Configured  0 or 1
	   uint16_t PUSCH_ID = 499;  // Configured  0 or 1

        uint32_t T_offset_max;
        uint32_t T_srs = T_srs_value(Config_idx,Duplex_Mode); 
	uint32_t T_offset = Offset_value(T_srs,Duplex_Mode, Config_idx,Offset_Idx,&T_offset_max);
	uint32_t n_srs = Get_n_srs( Config_idx, ns,  nf, N_sp, T_srs, Duplex_Mode, T_offset,  T_offset_max);
        uint32_t N_b = Nb[srsbwtable_idx(n_ul_rb)][B][bw_cfg];
        printf("N_b = %d\n",N_b);
        uint32_t Fb = Get_Fb(n_ul_rb,B,bw_cfg,HoppingBandwidth,n_srs,N_b);
        uint32_t n_b = Get_nb( HoppingBandwidth,B, n_ul_rb, bw_cfg, freqDomainPosition, n_srs, Fb,N_b);
		uint32_t K_Tc_p = Get_K_Tc_p(Cyclic_shift,  N_Tx, K_Tc);
    	printf("T_offset = %d ,T_offset_max = %d , T_srs = %d , n_srs = %d, Fb = %d, n_b = %d, K_Tc_p = %d\n", T_offset,T_offset_max,T_srs,n_srs,Fb,n_b,K_Tc_p);
		srslte_refsignal_srs_rb_L_cs(bw_cfg,  n_ul_rb);
		Get_m_srs_0(  bw_cfg);
}
