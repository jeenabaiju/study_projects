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
#define format_1                 1
#define format_1a                 2
#define format_1b                3
#define format_2                 4
#define format_2a               5
#define format_2b               6
#define format_3                 7
#define format_4                8
#define format_5                9
struct pucch_format {

    uint8_t format1;
    uint8_t format2;
    uint8_t format3;
    uint8_t format1_a;
    uint8_t format1_b;
    uint8_t format2_a;
    uint8_t format2_b;
};
struct cell {

    uint8_t CP;
    uint8_t CellID;

};
struct pucch_config {

    uint8_t CP;// Normal-3 and extended-2
    uint8_t CellID;
    uint8_t delta_ss_pucch;
    uint8_t N_cs_1; // No.of cyclic shift for PUCCH 1/1a/1b
    uint8_t N_RB_2; // N_RB_2 <= 0 denotes BW available for PUCCH 2/2a/2b
    uint8_t ns;// number of slots
    uint8_t N_ID_PUCCH;// cell ID if configured 1 else 0
    uint8_t PUCCH_ID;
    uint8_t N_UL_symbol;
    uint8_t n_pucch_1;
    uint8_t n_pucch_2;
    uint8_t n_pucch_3;
    uint8_t n_pucch_4;
};
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

uint32_t pucch_dmrs_symbol_format1_cpnorm[3] = {2, 3, 4};
uint32_t pucch_dmrs_symbol_format1_cpext[2] = {2, 3};
uint32_t pucch_dmrs_symbol_format2_3_cpnorm[2] = {1, 5};
uint32_t pucch_dmrs_symbol_format2_3_cpext[1] = {3};
uint32_t pucch_dmrs_symbol_format2a_2b_cpnorm[2] = {1, 5};
uint32_t pucch_dmrs_symbol_format4_5_cpnorm[1] = {3};
uint32_t pucch_dmrs_symbol_format4_5_cpext[1] = {2};
/****************************************************************************/

/****************************************************************************/
//PUCCH
/****************************************************************************/
/* Number of PUCCH demodulation reference symbols per slot N_rs_pucch tABLE 5.5.2.2.1-1 36.211 */
static uint32_t get_N_rs_PUCCH(uint32_t format, uint32_t CP)
 {

    static const unsigned n_RS_PUCCH1 = (1 << 1) | (1 << 2) | (1 << 3) ;
    static const unsigned n_RS_PUCCH2_3 = (1 << 4) | (1 << 7)  ;
    static const unsigned n_RS_PUCCH2a_2b = (1 << 5) | (1 << 6)  ;
     if (CP)
     {
         if ((1 << format) & n_RS_PUCCH1 )
         {
             return 3;
         }
         else if ((1 << format) & n_RS_PUCCH2_3 )
         {
             return 2;
         }
         else
         {
        	 if ((1 << format) & n_RS_PUCCH2a_2b)
             {
                return 2;
             }
         }
     }
     else
     {
         if ((1 << format) & n_RS_PUCCH1 )
         {
             return 2;
         }
         else if ((1 << format) & n_RS_PUCCH2_3 )
         {
             return 1;
         }
         else
         {
        	 if ((1 << format) & n_RS_PUCCH2a_2b)
             {
               return 0;
             }
         }
     }
 }
/****************************************************************************/
     /* Table 5.5.2.2.2-1: Demodulation reference signal location for different PUCCH formats. 36.211 */
static uint32_t get_pucch_dmrs_symbol(uint32_t l, uint32_t format, uint32_t CP)
{
	    static const unsigned n_RS_PUCCH1 = (1 << 1) | (1 << 2) | (1 << 3) ;
	    static const unsigned n_RS_PUCCH2_3 = (1 << 4) | (1 << 7)  ;
	    static const unsigned n_RS_PUCCH2a_2b = (1 << 5) | (1 << 6)  ;
	    static const unsigned n_RS_PUCCH4_5 = (1 << 8) | (1 << 9)  ;


	         if ((1 << format) & n_RS_PUCCH1 )
	         {
	        	 if (CP)   // CP==1
	        	 {
	        	   if (l < 4)
	        	   {
	        	      return pucch_dmrs_symbol_format1_cpnorm[l];
	        	   }
	        	 }
	        	 else
	        	 {
	        	    if (l < 3)
	        	    {
	        	       return pucch_dmrs_symbol_format1_cpext[l];
	        	    }
	        	 }
	         }
	         else if ((1 << format) & n_RS_PUCCH2_3 )
	         {
	        	 if (CP)   // CP==1
	        	 {
	        	    if (l < 3)
	        	    {
	        	        return pucch_dmrs_symbol_format2_3cpnorm[l];
	        	    }
	        	 }
	             else
	             {
	        	   if (l < 2)
	        	   {
	        	      return pucch_dmrs_symbol_format2_3cpext[l];
	        	   }
	             }
	         }
	         else if ((1 << format) & n_RS_PUCCH2a_2b)
	         {
	        	 if (CP)   // CP==1
	        	 {
	        	    if (l < 3)
	        	    {
	        	 	   return pucch_dmrs_symbol_format2a_2bcpnorm[l];
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
	   	    			if(l < 2)
	   	    			{
	   	    			return pucch_dmrs_symbol_format4_5_cpnorm[l];
	   	    			}
	   	    		}
	   	    		else
	   	    		{
	   	    			if(l < 2)
	   	    			{
	   	    				return pucch_dmrs_symbol_format4_5_cpext[l];
	   	    			}
	   	    		}
	   	    	 }
	         }
}
/****************************************************************************/
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
/****************************************************************************/
/* Generates DMRS for PUCCH according to 5.5.2.2 in 36.211 */
int pucch_dmrs_gen(uint32_t format, uint32_t n_pucch, uint32_t sf_idx, uint8_t pucch_bits[2], float complex *r_pucch, struct cell *cell)
{
  int ret = ERROR_INVALID_INPUTS;

    uint32_t N_rs=get_N_rs(format, cell->CP);

    float complex z_m_1 = 1.0;
    for (uint32_t m=0;m<N_rs;m++)
    {
        uint32_t n_oc=0;
        uint32_t l = get_pucch_dmrs_symbol(m, format, cell->CP);
    }
        // Add cyclic prefix alpha
        float alpha = 0.0;
        if (format < format2)
        {
          alpha = pucch_alpha_format1(n_cs_cell, pucch_cfg, n_pucch, cell->CP, true, ns, l, &n_oc, NULL);
        }
        else
        {
          alpha = pucch_alpha_format2(n_cs_cell, pucch_cfg, n_pucch, ns, l);
        }
        /****************************************************************************/

  return ret;
}
int get_n_cs_ptilda(uint32_t n_cs_cell,uint32_t delta_ss_pucch,n_pucch,struct cell *cell,uint32_t l )
{
   uint32_t n_cs;
}
/****************************************************************************/
void get_ndash_ptilda(struct pucch_config *config, uint32_t c, uint32_t Ndash ,uint32_t *N_dash_slot1,uint32_t *N_dash_slot2)
{
    uint32_t temp, temp1;
    uint32_t N_dash_ptilda;
    temp1 = (c * config->N_cs_1)/config->delta_ss_pucch;
    temp = (c * N_sc)/config->delta_ss_pucch;
    uint32_t len = config->ns+1;
    uint8_t n_prs[len];
	const uint16_t N_ID = Get_cellID(srs_ul);
    const uint32_t c_init = N_ID;
 /*   if (((config->ns % 2) == 0)||((config->ns % 2) == 1))*/

       if (config->n_pucch_1 < ((c * config->N_cs_1)/config->delta_ss_pucch ))
       {
          N_dash_slot1= config->n_pucch_1;
       }
       else
       {
          N_dash_slot1 = (config->n_pucch_1 - temp1) % temp;
       }
       if (config->n_pucch_1 >= ((c * config->N_cs_1)/config->delta_ss_pucch ))
       {
	      calc_prs_c( c_init, len, n_prs);
          N_dash_slot2 = (n_prs[N_dash_slot1 + 1]) % (temp + 1) - 1 ;
       }
       else
       {
           uint32_t temp2, temp3;
           temp2 = (h_p % c)* Ndash;
           temp3 = temp2 / config->delta_ss_pucch;
           uint32_t h_p = get_h_p(config,c,N_dash_slot1, Ndash);
           N_dash_slot2 = (h_p / c) + temp3;
       }
 }
/****************************************************************************/
//get h_p
/****************************************************************************/
int get_h_p(struct pucch_config *config, uint32_t c, uint32_t N_dash_slot1,uint32_t Ndash)
{
   uint32_t h_p;
   if (config->CP)
   {
       d = 2;
   }
   else
   {
       d = 0;
   }
   h_p = (N_dash_slot1 + d) % ((c * Ndash)/config->delta_ss_pucch);
   return h_p;
}

/****************************************************************************/
//Ndash
/****************************************************************************/
int get_Ndash(struct pucch_config *config ,uint32_t c)
{
   uint32_t Ndash;
   if (config->n_pucch_1 < ((c * config->N_cs_1)/config->delta_ss_pucch ))
   {
       Ndash = config->N_cs_1;
   }
   else
    {
        Ndash = N_sc;
    }
    return Ndash;
}
/****************************************************************************/
// cyclic shift n_cs_cell
/****************************************************************************/
int get_n_cs_cell(n_pucch,struct cell *cell,uint32_t l,struct SRS_UL *srs_ul,struct pucch_config *config)
{
    uint32_t len = config->ns+1;
    uint8_t n_prs[len];
    uint32_t n_cs_cell;
    n_cs_cell = 0;
	const uint16_t N_ID = Get_cellID(srs_ul);
	const uint32_t c_init = N_ID;
	calc_prs_c( c_init, len, n_prs);
	for ( i = 0; i < 8 ; i++)
    {
	  n_cs_cell += (((uint32_t) n_prs[8 * config->ns * config->N_UL_symbol + 8* l + i]) << i);
    }
    return n_cs_cell;
}
/****************************************************************************/
int get_CP_pucch(struct pucch_config *config, )
{
    uint32_t c ;
    if (config->CP)
    {
         c = 3;// Normal CP
    }
    else
    {
         c = 2;// Extended CP
    }
return c;
}
/****************************************************************************/
void get_n_oc(struct pucch_config *config,uint32_t Ndash, uint32_t c, uint32_t *n_oc_slot1,  uint32_t *n_oc_slot2)
{
    uint32_t n_oc_slot1;
    uint32_t n_oc_slot2;
    uint32_t N_dash_slot1;// for slot ns
    uint32_t N_dash_slot2;// for slot ns-1
    get_ndash_ptilda(config,c,Ndash ,&N_dash_slot1,&N_dash_slot2);
    if (config->CP)// normal CP
    {
      n_oc_slot1 =  (&N_dash_slot1 * config->delta_ss_pucch)/Ndash;// for slot ns
      n_oc_slot2 =  (&N_dash_slot2 * config->delta_ss_pucch)/Ndash;// for slot ns-1
    }
    else  // extended CP
    {
      n_oc_slot1 = 2*((&N_dash_slot1  * config->delta_ss_pucch)/Ndash);
      n_oc_slot2 = 2*((&N_dash_slot2  * config->delta_ss_pucch)/Ndash);
    }
 }
/****************************************************************************/





/****************************************************************************/

