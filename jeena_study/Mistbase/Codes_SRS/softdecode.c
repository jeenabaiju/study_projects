/*
 * softdecode.c
 *
 *  Created on: Jun 7, 2016
 *      Author: mistbasejeena
 */
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>

#define  N_SYMB_TMP 14
#define  N_sc 72
#define MAX_CW 1024

static int32_t skip_indices[] = {
  0, 12, 24, 36, 48, 60, 72, 84, 96, 108, 120, 132, 582, 594, 606,
  618, 630, 642, 654, 666, 678, 690, 702, 714, 720, 722, 724, 726,
  728, 730, 732, 734, 736, 738, 740, 742, 744, 746, 748, 750, 752,
  754, 756, 758, 760, 762, 764, 766, 768, 770, 772, 774, 776, 778,
  780, 782, 784, 786, 788, 790, 792, 794, 796, 798, 800, 802, 804,
  806, 808, 810, 812, 814, 816, 818, 820, 822, 824, 826, 828, 830,
  832, 834, 836, 838, 840, 842, 844, 846, 848, 850, 852, 854, 856,
  858, 860, 862, 1008, 1020, 1032, 1044, 1056, 1068, 1080, 1092,
  1104, 1116, 1128, 1140, 1590, 1602, 1614, 1626, 1638, 1650, 1662,
  1674, 1686, 1698, 1710, 1722,
  -1 // end marker
};
inline int stuff_td_data16( uint16_t datavalue, uint8_t** td_data_access , uint32_t *stream, uint32_t * counter, uint16_t CW){

  uint8_t *relocate = *td_data_access;

  **td_data_access = datavalue; // write 4 nibbles at the time
  (*td_data_access)++;
  // Last of termination bits
  if(*counter>CW){
    *(*td_data_access) = 0;
    (*td_data_access)++;
    (*stream)++;
    if (*stream == 3) {
       return 0;
    }
    if (*counter > CW+2){ // relocate blocks to next stream
        relocate++;
        uint8_t first;
        first= *relocate;
        *relocate = 0;
        relocate +=3;
        *relocate =first;
        *td_data_access=relocate;
        *counter=2;
    }
    else {
      *counter = 0;
    }
  }
  else
    (*counter)+=4;
  return 1;
}
void decode_QPSK_soft(int16_t* data, uint16_t CW, uint8_t* td_data)
{
	/* 01  													  00
	 * x-----------------x----x----xx----x----x---------------x
	 * -sq1           -alpha -sq2   0    sq2 alpha            sq1
	    =-1/sqrt(2)                                         =1/sqrt(2)
	 */

    int32_t i;
    int32_t *next_skip_index = skip_indices;
    uint32_t cw = CW;
    uint32_t td_pos = 0;
    uint32_t stream_nbr = 0;
    int32_t alpha = 0x2D41; //  // 1/2* sqrt(2) bcz [ (1/sqrt(2)+0)/2 = 1/2* sqrt(2)]
    int32_t sq1 =  0x5A82;   //  1/sqrt(2)
    int32_t sq2 =  0x16A0;   //  1/(4*sqrt(2)) between alpha and 0
    int32_t subalphascaler = (21);  //scalefactor -alpha<x<alpha  [7*0xffff/sq1]
    int32_t superalphascaler = (21); //scalefactor x<-alpha;x>alpha [7*0xffff/sq1-alpha]// alpha =0
    if (CW > MAX_CW*2)
    {
      cw = MAX_CW*2;
    }
    for (i =0 ; i < N_sc * N_SYMB_TMP * 2; i += 2)
    {
      if (i == *next_skip_index)
      {
        next_skip_index++;
        data += 2;
        continue;
      }
      int32_t re, im;
      uint32_t result[2];  // the two bit patterns
      re = *data++;
      im = *data++;
      if (re > 0)
      {
        result[0] = 7 - (abs(re - sq1)*superalphascaler);   //handling z0
      }
      else
      {
        result[0] = -8 + (abs(re + sq1)*superalphascaler);
      }
      if (im > 0)
      {
        result[1]= 7- (abs(im - sq1)*superalphascaler);         //handling 0z
      }
      else
      {
        result[1]= -8 + (abs(im + sq1)*superalphascaler);
      }
      result[0] = result[0] >>2;
      result[1] = result[1] >>2;
       printf(" result[0] = %d \nresult[1] = %d \n",result[0],result[1]);
      uint16_t stuffer;
      stuffer = result[0] & 0x0F;
      stuffer |= result[1] <<2 & 0xF0;

      printf(" stuffer = %d \n ",stuffer);
     if (!stuff_td_data16(stuffer, &td_data, &stream_nbr, &td_pos, cw))
       break;
    }
}
void main()
{
	int16_t data[32]= {3,1,4,-1,4,-4,.1,-.1,1,0,2,7,4,-3,-7,-5,4,3,-2,-4,.1,-.1,1,0,2,7,4,-.5,0.7,-.4,3};
	uint8_t td_data;
	uint16_t CW = 1024;
	decode_QPSK_soft(data, CW, &td_data);

}

