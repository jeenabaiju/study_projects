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
#define CP_NSYMB                7
uint32_t N_sc_RB = 12;
struct NPUSCH {
           uint32_t sequence_hopping;// enable=1 and disbale=0
           uint32_t group_hopping;// enable=1 and disbale=0
           uint32_t groupAssignmentNPUSCH; //delta_ss = { 0 , 1 ,..., 29 }
           uint32_t ns;//Slot number within a radio frame
           uint32_t npusch_AllSymbols;//true -1 false -0
           uint32_t srsSubframeConfig;
           uint32_t I_RU;
           uint32_t NPUSCHformat;
           uint32_t threeTone_CyclicShift;
           uint32_t sixTone_CyclicShift;
           };



           // Table 10.1.4.1.1-1: Definition of   w(n)

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
/******************************************************************************/
uint32_t Get_N_RU(uint32_t I_RU)
{
    uint32_t N_RU;
    //Table 16.5.1.1-2: Number of resource units (N_RU) for NPUSCH.
  if (I_RU < 6)
  {
     N_RU = I_RU +1;
  }
  else if (I_RU == 6)
  {
     N_RU = 8;
  }
  else
  {
     N_RU = 10;
  }
  return N_RU;
}
/******************************************************************************/
uint32_t Get_Nslots(delta_f)
{
    uint32_t Nslots;
    Nslots = delta_f ? 1:2;
    return Nslots;
}
/******************************************************************************/
uint32_t get_M_identical_NPUSCH(uint32_t N_sc_RU,uint32_t M_rep_NPUSCH)
{
    uint16_t temp,M_identical;
    if (N_sc_RU)
    {
        M_identical = 1;
    }
    else if(N_sc_RU > 2)
    {
        temp = M_rep_NPUSCH/2;
        if (temp < 4)
        {
           M_identical = temp;
        }
        else
        {
             M_identical = 4;
        }
     }
     return M_identical;
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
/* Group Hopping*/

uint32_t Get_uValue(uint32_t N_CellID,uint32_t group_hopping,uint32_t NPUSCHformat)
{
        uint32_t u;
        if (group_hopping)
        {
            if (NPUSCHformat)
        }


            u=N_CellID % 16;

}


/******************************************************************************/
/* Demodulation reference signal */
/******************************************************************************/

/* Reference signal sequence for  N_sc_RU = 1*/
uint32_t Get_r_u_bar(uint32_t N_CellID,uint32_t N_UL_slots,uint32_t N_rep_NPUSCH,uint32_t N_RU,uint32_t u,float complex *r_u_bar)
{
    int n;
    uint32_t c_init = 35;
    uint32_t len = 35;
    calc_prs_c( c_init, len, n_prs); /*generate_pseudo random sequence*/
    for (n = 0; n < N_rep_NPUSCH*N_UL_slots*N_RU ; n++)
    {
        r_u_bar(n) = 1/sqrt(2)*(1+M_PI)*(1-(2*n_prs[n]))*w[u][n%16];
    }
return 0;
}
/******************************************************************************/
uint32_t Get_r_u_format1(uint32_t N_CellID, uint32_t ns,uint32_t N_UL_slots,uint32_t N_rep_NPUSCH,uint32_t N_RU,uint32_t u)
{
    Get_r_u_bar(N_CellID,N_UL_slots,N_rep_NPUSCH, N_RU,u,r_u_bar);
    for (n = 0; n < N_rep_NPUSCH*N_UL_slots*N_RU ; n++)
    {
        r_u_format1(n) = r_u_bar(n) ;
    }
  return 0;
}
/******************************************************************************/

uint32_t Get_n_oc(uint32_t N_CellID,uint32_t ns,uint32_t *n_oc)
{
    int i;
    uint32_t c_init = N_CellID;
    uint32_t len = 8 * ns + 8;
    uint32_t temp;
    temp = 0;
    calc_prs_c( c_init, len, n_prs); /*generate_pseudo random sequence*/
    	  for (i = 0;i < 8; i++)
	      {
	        temp += (((uint32_t) n_prs[8 * ns + i]) << i);
	      }
    n_oc = temp % 3;
    return 0;
}
uint32_t Get_r_u_format2(uint32_t N_CellID, uint32_t ns,uint32_t N_UL_slots,uint32_t N_rep_NPUSCH,uint32_t N_RU,uint32_t u,uint32_t n_oc,float complex *r_u_format2)
{
    int n,m;
    uint32_t w_m;
    float complex r_u_bar[N_rep_NPUSCH*N_UL_slots*N_RU];
    Get_r_u_bar(N_CellID,N_UL_slots,N_rep_NPUSCH, N_RU,u,r_u_bar);
    for (m = 0; m < 3; m++)
    {
       for (n = 0; n < N_rep_NPUSCH*N_UL_slots*N_RU ; n++)
       {
         w_m = w_arg_pusch_format2[n_oc][m];
         r_u_format2[3*n+m] = w_m*r_u_bar[n];
       }
    }
    return 0;
}
/******************************************************************************/





