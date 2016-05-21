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
uint32_t Get_K_Tc_p(uint32_t Cyclic_shift, uint32_t N_Tx, uint32_t K_Tc)
{
    uint32_t K_Tc_bar = K_Tc ;
    int p_index;
    uint32_t K_Tc_p;//K_Tc_p={0,1,...K_Tc-1}
    static const unsigned n_SRS_cs = (1 << 4) | (1 << 5) | (1 << 6) | (1 << 7);
    static const unsigned p_hat = (1 << 1) | (1 << 3) ;
    /*find K_Tc_p*/
    for (p_index = 0; p_index <=3; p_index++)
    {
        if (((1 << Cyclic_shift) & n_SRS_cs) && ((1 << p_index) & p_hat) && (N_Tx == 4) ) // p_hat ands p relation table Table 5.2.1-1
        {
            K_Tc_p = 1 - K_Tc_bar;
            printf("K_Tc_p = %d \n",K_Tc_p);

        }
        else 
        {
            K_Tc_bar = K_Tc;        
            K_Tc_p = K_Tc_bar;
        }
    }
    return K_Tc_p;
}
void main()
{
    uint32_t Cyclic_shift = 0;
    uint32_t K_Tc_p ;//K_Tc_p={0,1,...K_Tc-1}
    uint32_t K_Tc_bar;
    uint32_t N_Tx = 1;
    uint32_t K_Tc = 2;
    Get_K_Tc_p(Cyclic_shift, N_Tx, K_Tc);
    printf("K_Tc_p = %d \n",K_Tc_p);
}


