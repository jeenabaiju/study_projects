
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
/* Calculate N_CellID
* input - cell_ID, Configurations for SRS_UL.N_ID_PUSCH,SRS_UL.N_ID_PUSCH
* output- N_ID
*/
const uint16_t Get_cellID(uint32_t N_ID_PUCCH, uint32_t N_ID_PUSCH, const uint16_t  cell_ID, const uint16_t  PUCCH_ID, const uint16_t  PUSCH_ID)
{
    const uint16_t  N_ID;
    if (N_ID_PUCCH == 1)
    {
        const uint16_t  N_ID = PUCCH_ID;
        return  N_ID;
    }
    
    else if (N_ID_PUSCH == 1)
        {
            const uint16_t  N_ID = PUSCH_ID;
            return  N_ID;
        }
    else
    {
        const uint16_t  N_ID = cell_ID;// for SRS also//
        printf("N_ID = %d \n",N_ID);
        return  N_ID;
    }
    
}
void calc_n_prs_c(const uint32_t c_init, const uint32_t len, uint8_t* n_prs)
{
    uint32_t N_c = 1600;
    uint32_t x1, x2;
    x1 = 0x54D21B24; // x1 sequence moved forward
    x2 = c_init;
    uint32_t next_x1 = 0;
    uint32_t next_x2 = 0;

    // Move x2 forward
    uint32_t i;
    for(i = 0; i < N_c - 31; i++)
    {
        next_x2 = ((x2 >> 3) ^ (x2 >> 2) ^ (x2 >> 1) ^ x2) & 0x1;

        x2 = (x2 >> 1) | (next_x2 << 30);
    }

    // Calculate prs_c
    for(i = 0; i < len; i++)
    {
        next_x1 = ((x1 >> 3) ^ x1) & 0x1;
	next_x2 = ((x2 >> 3) ^ (x2 >> 2) ^ (x2 >> 1) ^ x2) & 0x1;

	*(n_prs++) = next_x1 ^ next_x2;

	x1 = (x1 >> 1) | (next_x1 << 30);
	x2 = (x2 >> 1) | (next_x2 << 30);
    }

}

void main()
{

    uint32_t N_ID_PUCCH= 1;
    uint32_t N_ID_PUSCH= 0;
    const uint16_t PUSCH_ID= 400; 
    const uint16_t PUCCH_ID= 400;
    const uint16_t cell_ID = 300;
    const uint32_t len = 167;
    uint8_t n_prs [len];
    int i;
    const uint16_t n_ID = Get_cellID( N_ID_PUCCH, N_ID_PUSCH, cell_ID, PUCCH_ID, PUSCH_ID);
    const uint32_t c_init = floor( cell_ID / 30 );
    printf("c_init = %d \n",c_init);
    
    printf("n_ID = %d \n",n_ID);
    calc_n_prs_c(c_init, len, n_prs);
    for (i = 0; i< len; i++)
    {
        printf("N_PRS[%d] = %d \n",i,n_prs[i]);
    }
    
}

