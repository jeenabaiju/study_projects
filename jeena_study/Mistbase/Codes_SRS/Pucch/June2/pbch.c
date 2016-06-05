/**
 * pbch.c
 * Jeena
 */
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

// Common error codes
#define SUCCESS                 0
#define ERROR                  -1
#define INVALID_INPUTS         -2


// First you have to define properites of a eNodeB.
// NDLRB indicate System Bandwith in the unit of RBs.
//NDLRB 6 = 1.4 Mhz, NDLRB 15 = 3.0 Mhz, NDLRB 25 = 5.0 Mhz,
// CellRefP indicate number of downlink Antenna.CellRefP = 1 =>1tx antenna(SISO)
// NCellID indicate PCI (Physical Channel Identity) of the Cell
// NSubframe indicate the subframe number.
/****************************************************************************/
uint32_t tti_bch = 40; // 40ms
uint32_t NSCRB = 12; // for CP normal and extended cases
uint32_t NDLRB = 6;
struct PBCH{
	uint32_t PHICH_R_1_6;
	uint32_t PHICH_R_1_2;
	uint32_t PHICH_R_1;
	uint32_t PHICH_R_2;
	uint32_t phich_CP;
	int NDL_RB ;
	int Ng;
	int DL_bw_idx;
};

struct pbch_t{
       uint32_t nports; //N_Rx = nports =1
       uint32_t CP;
       };
struct cell{
       uint32_t CellID;
       uint32_t N_Rx; //N_Rx =1
       uint32_t CP;
       uint32_t NFrame;
       uint32_t MOD_SCHEME;// QPSK for PBCH
};
uint32_t BCH_PAYLOAD_LEN;

/*****************************************************************************/
bool pbch_exists(int NFrame, int nslot) 
{
  return (!(NFrame % 5) && nslot == 1);// confirms first 4 frame and slot 1 for eah frame
}

uint32_t PHICH_EXT = 1728;
uint32_t PHICH_NORM = 1920;
uint32_t MAX_PORTS = 1;
uint32_t MAX_LAYERS = 1;


 const uint8_t srslte_crc_mask[4][16] = {
    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, 
    { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 }, 
    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, 
    { 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1 } }; 

     
// CP_NORM 1920bits and CP_EXTD 1728 bits
//generate MIB bits.same MIB  will be sent in n, n+1, n+2 and n+3 System Frames.
/*generate physical layer symbols for BCH which  is PBCH symbol.
*PBCH symbol is complex and has NSymbols x N_Tx values*/
//split the generated symbols into 4 clusters, each of them has 240 symbols


// l - index for 4 OFDM symbols in slot 1 {0,1,2,3}
/****************************************************************************/
/* //Compute and generate pseudo random sequence
* input - c_init, len
* output- n_prs values

 * written by Henrik
 * */
/******************************************************************************/
void calc_pseudoseq(const uint32_t c_init, const uint32_t len, uint8_t* n_prs)
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
// RECEIVER
/****************************************************************************/
/*DESCRAMBLING*/
/******************************************************************************/
int pbch_descrambling(uint32_t CellID, const uint32_t  Mbits, uint32_t *demod_bits, uint32_t *bits_out, uint32_t NFrame)
{
  int i;
  uint8_t n_prs[Mbits];
  uint32_t c_init;
  if (NFrame % 4 == 0)
  {
	  c_init = CellID;

     for (i = 0; i < Mbits; i++)
     {
         calc_pseudoseq( c_init, Mbits, n_prs);// calculate pseudo random sequence
         bits_out[i]= (demod_bits[i] + n_prs[i]) % 2; // de-scrambled output
     }
  }
  else
  {
	  c_init = 0;
  }
  return c_init;
}
/******************************************************************************/
/*DEMODULATION*/
/******************************************************************************/
int pbch_softdemod(int Msymbol,float complex *rx_symbols, uint32_t *soft_bits, double N0)
{
	double factor = 2 * sqrt(2.0) / N0;
	double exp_pi4 ;
	int i;
	float complex temp;

	soft_bits[Msymbol] = 0;
	exp_pi4 = exp( cos(M_PI / 4) + sin(M_PI / 4));
	for (i = 0; i < Msymbol; i++)
	{
		temp = rx_symbols[i] * exp_pi4;
		soft_bits[(i << 1) + 1] = creal(temp) * factor;
		soft_bits[i << 1] = cimag(temp) * factor;
	}

return Msymbol;
}
/******************************************************************************/
/*DEMAPPING*/
/******************************************************************************/

int pbch_res_demapping(uint32_t NFrame, uint32_t N_Rx, float complex *data, float complex *y, uint32_t NDL_RB)
{
//int temp, kdash,i;
uint32_t k;
k = 0;

return k;
}
/******************************************************************************/

int predecoding_single_gen(float complex *y, float complex *h, float complex *x, int Msymbols, float noise_estimate)
{
  int i ;
  for (i = 0; i < Msymbols; i++)
  {
    x[i] = y[i]*conj(h[i])/(conj(h[i])*h[i]+noise_estimate);
  }
  return 0;
}
/******************************************************************************/
int predecoding_1port(float complex *symbols[0], float complex *ch_est[0], float complex *d, int Msymbols, float noise_estimate)
{
	/* ZF/MMSE SISO equalizer x=y(h'h+no)^(-1)h' (ZF if n0=0.0)*/
	predecoding_single_gen(symbols[0],ch_est[0], d, Msymbols, noise_estimate);
	return Msymbols;
}

/******************************************************************************/
/* Decodes the PBCH channel
 *The PBCH spans for 40 ms. This function is called every 10 ms. It tries to 
decode the MIB from the symbols of a subframe (1 ms).*/
/******************************************************************************/
int pbch_decode(struct cell *cell,float complex *slot1_symbols, float complex *ch_est_slot1[MAX_PORTS], float noise_estimate,struct pbch_t *t)
{
  int i;
  //_Complex *x[MAX_LAYERS];// only 1 layer here
  int ret = ERROR; // incase of error returns this
      /* Set pointers for layermapping & precoding */
      /* number of layers equals number of ports */
    for (i = 0; i < t->nports; i++)
    {
      if (ch_est_slot1[i] == NULL) 
      {
        return ERROR;
      } 
    }
     ///////* extract symbols */
     ////// switch to DL and hence processed in the DL section gets Msymbols from this

   /*   bch_symbols[0] = DL_input[0];// from up layer*/
   /*  int Msymbols = no_f_symbols;*/
    /* channel estimation*/
    // switch to hardware processing and gets ch_est from this
   /*    ch_est_slot1[0]= ch_est_output[0];*/
    /* no need for layer demapping  but does de-precoding*/

     /*  symbols[0] = bch_symbols[0];*/
  /*  int Msymbols = predecoding_1port(symbols[0], ch_est_slot1[0], &d, Msymbols, noise_estimate);*/

    /* demodulate symbols */
   /*int Mbits =  pbch_softdemod(Msymbol,d,demod_bits);*/

    /* de-scrambling*/
   /* pbch_descrambling(cell->CellID, Mbits,demod_bits, bits_out, cell->NFrame);*/
return ret;
}
/******************************************************************************/


int pbch_decode_frame(uint32_t NFrame)
{
if (NFrame % 4 == 0)
{
}
return 0;
}
/******************************************************************************/
/*CRC-CHECK*/
/******************************************************************************/
uint32_t pbch_crc_check(uint8_t bch_payload[BCH_PAYLOAD_LEN])
{



return 0;
}
/******************************************************************************/



//FINISHED
/******************************************************************************/
/**Packs bits.*/
uint32_t bit_pack(uint8_t **bits, int nbits)
{
    int i;
    uint32_t num=0;
    for(i=0; i<nbits; i++) 
    {
      num |= (uint32_t) (*bits)[i] << (nbits-i-1);// shifting the bits to find the value
    }
    *bits += nbits;
    return num;
}
/******************************************************************************/
/* MIB - 24 bit information with following information 
    3 bits - system bandwidth
    3 bits -PHICH information,
        1 bit - normal or extended PHICH
        2 bit - the PHICH Ng value
    8 bits-system frame number
    10 bits - reserved for future use (spare)
    Apart from the information in the payload, the MIB CRC also has the number 
    of tx antennas used by the eNodeB.
    The MIB CRC is scrambled or XORed with a antenna specific mask.*/
/******************************************************************************/

/** Unpacks MIB from PBCH message.
 * msg buffer must be 24 byte length at least
 */

void pbch_mib_unpack(uint8_t *msg, uint32_t *sfn, struct PBCH *t)
{
  int phich_resources;

   /* gets the DL bandwidth from first 3 bits*/
  t-> DL_bw_idx = bit_pack(&msg, 3);// gets the bw from the MIB
  switch (t->DL_bw_idx)
  {
      case 0:
    	  t-> NDL_RB = 6;
        break;
      case 1:
    	  t-> NDL_RB = 15;
        break;
      default:
    	  t-> NDL_RB = (t->DL_bw_idx - 1) * 25;
        break;
  }// return the value of NDLRB
  /*gets the CP from the next bit - PHICH configuration(3 bits)*/
  if (*msg) // 1 bit for NORMCP or EXTCP
  {
    t->phich_CP = PHICH_EXT;// if bit is 1 extended CP
  } 
  else 
  {
	  t->phich_CP = PHICH_NORM;// if bit is 0 , Normal CP
  }
  msg++;// return the value of CP

  phich_resources = bit_pack(&msg, 2);// PHICH Ng value
  /* according to 36.211 13.01 section 6.9 */
  switch (phich_resources) 
  {
      case 0:
    	  t->  Ng = t->PHICH_R_1_6;
        break;
      case 1:
    	  t-> Ng = t->PHICH_R_1_2;
        break;
      case 2:
    	  t->  Ng= t->PHICH_R_1;
        break;
      case 3:
    	  t-> Ng = t->PHICH_R_2;
        break;
  }// return the value of Ng
  if (sfn) 
  {
    *sfn = bit_pack(&msg, 8) << 2;    
  }// return the value of system frame number
}
/******************************************************************************/









