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
int get_n_cs_cell(uint32_t CP, uint32_t n_cs_cell[uint32_t NSLOTS_X_FRAME][3])
{
    uint32_t len = 7;
    uint8_t n_prs[len];
	const uint16_t N_ID = 1;
	const uint32_t c_init = N_ID;
	calc_prs_c( c_init, len, n_prs);
	uint32_t CP_NSYMB = CP?3:2;
	uint32_t nslot;
	uint32_t l;
	int i;
	  for (nslot = 0; nslot < NSLOTS_X_FRAME; nslot++)
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
void main ()
{
int nslot, l;
uint32_t NSLOTS_X_FRAME = 2;
uint32_t CP = 1;
uint32_t n_cs_cell[NSLOTS_X_FRAME ][3] = 0;
get_n_cs_cell(CP,n_cs_cell);
 for (nslot = 0; nslot < NSLOTS_X_FRAME; nslot++)
	  {
	    for (l = 0; l < CP_NSYMB; l++)
	    {
	    printf ("N_cs_cell[%d][%d] = %d", nslot, l, n_cs_cell[nslot][l]);
	    }
	    }
	    }
