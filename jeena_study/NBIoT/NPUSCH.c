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
struct NPUSCH npusch{
           uint32_t sequence_hopping;// enable=1 and disbale=0
           uint32_t group_hopping;// enable=1 and disbale=0
           uint32_t delta_ss; //delta_ss = { 0 , 1 ,..., 29 }
           uint32_t ns;//Slot number within a radio frame
           uint32_t npusch-AllSymbols;//true -1 false -0
           uint32_t srsSubframeConfig;
           uint32_t I_RU;
           };

uint32_t Get_N_RU(uint32_t I_RU)
{
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
}
uint32_t Get_Nslots(delta_f)
{
    Nslots = delta_f ? 1:2;
}

uint32_t get_M_identical_NPUSCH(uint32_t N_sc_RU,uint32_t M_rep_NPUSCH)
{
    uint16_t temp;
    if (N_sc_RU)
    {
        M_identical = 1;
    }
    else if(N_sc_RU > 2)
    {
        temp = M_rep_NPUSCH/2;
        if temp < 4
        {
           M_identical = temp;
        }
           M_identical = 4;
    }
}


/* Demodulation reference signal*/


/* Reference signal sequence for  N_sc_RU = 1*/


// Table 10.1.4.1.1-1: Definition of   w(n)

w_n[16][16] = {{1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1},
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





