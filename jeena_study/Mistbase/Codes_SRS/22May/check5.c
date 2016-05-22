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
           uint32_t Duplex_Mode;// FDD-0,TDD-1, Half Duplex -2
           uint32_t HoppingBandwidth; //b_hop={0,1,2,3}
           uint32_t freqDomainPosition;// n_RRC
           uint32_t freqDomainPosition_ap;// n_RRC
           uint32_t nf; //System frame number SFN 
           uint32_t Config_idx;// I_srs {0,..... 644}
		   uint32_t  OffsetIdx; // takes 0 or 1 
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



/***********************************************************************************************/
uint32_t T_srs_FDD(uint32_t Config_Idx) 
{
  uint32_t I_srs = Config_Idx;
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
uint32_t T_srs_TDD(uint32_t Config_Idx) 
{
  uint32_t I_srs = Config_Idx;
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

 
/*Table 8.2-2 of Subframe Offset Configuration T_offset for trigger type 0, TDD according to 3GPP TS 36.213 version 13.0.0 Release 13*/ 
uint32_t Mode(uint32_t Duplex_Mode, uint32_t Config_Idx,  uint32_t OffsetIdx) 
{ 
	uint32_t I_srs = Config_Idx;// I_srs is Config_Idx
	uint32_t T_srs;
	uint32_t T_offset;
	if (Duplex_Mode == 0) //  FDD mode
	{
		 T_srs = T_srs_FDD(Config_Idx) ;
		 T_offset = 0;
	}
	else if (Duplex_Mode == 1 )
	{
		T_srs = T_srs_TDD(Config_Idx) ;
		if (T_srs == 2)
		{
		      if ((Config_Idx == 0)&&  (OffsetIdx == 0))
              {
                  T_offset   = 0; 
		      }
		      else 
		      {
		          T_offset  = 1;
              }
	          if ((Config_Idx == 1)&&  (OffsetIdx == 0))
	          {
                  T_offset   = 0; 
              }
              else 
              {
			      T_offset = 2;
			  }
		      if ((Config_Idx == 2)&&  (OffsetIdx == 0))
		      {
			      T_offset = 1;
		      }
		      else
		      {
			      T_offset = 2;
		      }
              if ((Config_Idx == 3)&&  (OffsetIdx == 0))
		      {
			      T_offset = 0;
		      }
		     else
		      {
			     T_offset = 3;
		      }
		      if ((Config_Idx == 4)&&  (OffsetIdx == 0))
		      {
			     T_offset = 1;
		      }
		      else
		      {
			      T_offset = 3;
		      }
              if ((Config_Idx == 5)&&  (OffsetIdx == 0))
              {
                 T_offset   = 0; 
              }
		      else 
		      {
			      T_offset = 4;
		      }
              if ((Config_Idx == 6)&&  (OffsetIdx == 0))
              {
                  T_offset   = 1;
		      }
		      else
		      {
			      T_offset   = 4;
			  }
              if ((Config_Idx == 7)&&  (OffsetIdx == 0))
              {
                 T_offset   = 2;
		      }
		      else 
              {
			      T_offset   = 3;
		      }
              if ((Config_Idx == 8)&&  (OffsetIdx == 0))
              {
                  T_offset   = 2;
		      }
		      else 
		      {
			      T_offset   = 4;
		      }
              if ((Config_Idx == 9)&&  (OffsetIdx == 0))
              {
                  T_offset   = 3;
		      }
		      else 
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
               else if ( (I_srs >=325) && () I_srs < 645) )
               {
                  T_offset  = I_srs - 325; 
               } 
	       }
	   }
      else
	  {
		      if (Duplex_Mode == 2)
               {
		            T_srs = T_srs_Half(Config_Idx) ;
					T_offset = 0; 
			   }  
               else	
				{
	                 T_offset = 0; 
				}
	  }
   return T_offset; 
}