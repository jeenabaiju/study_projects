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

int pbch_softdemod(int Msymbol,float complex *rx_symbols, uint32_t *soft_bits, double N0)
{
	double factor = 2 * sqrt(2.0) / N0;
	float  exp_pi4 ;
	int i;
	double complex temp;

	soft_bits[Msymbol] = 0;
	exp_pi4 = exp( cos(M_PI / 4) + sin(M_PI / 4));
        printf("exp value is %f\n", exp_pi4 );
	for (i = 0; i < Msymbol; i++)
	{
		temp = rx_symbols[i] * exp_pi4;
		soft_bits[(i << 1) + 1] = creal(temp) * factor;
		soft_bits[i << 1] = cimag(temp) * factor;
	}

return Msymbol;
}
int pbch_softdemod1(int Msymbol,float complex *d, uint32_t *demod_bits)
{
   int i;
   int qpsklen;
   qpsklen = 2* Msymbol;
    float real[qpsklen];
    float imag[qpsklen];

    for (i = 0; i < Msymbol; i++)
    {
        real[i] = creal(d[i]);
        imag[i] = cimag(d[i]);
        if (real[i] >= 0.9)
        {
          demod_bits[2*i] = 0;
        }
        else if ((real[i] < 0.9) && (real[i] >= 0.6))
        {
          demod_bits[2*i] = 1;
        }
        else if ((real[i] < 0.6) && (real[i] >= 0.3))
        {
          demod_bits[2*i] = 2;
        }
        else if ((real[i] < 0.3) && (real[i] >= 0))
        {
          demod_bits[2*i] = 3;
        }
        else if ((real[i] < 0) && (real[i] >= -0.3))
        {
          demod_bits[2*i] = 4;
        }
        else if ((real[i] < -0.3) && (real[i] >= -0.6))
        {
          demod_bits[2*i] = 5;
        }
        else if ((real[i] < -0.6) && (real[i] >= -0.9))
        {
          demod_bits[2*i] = 6;
        }
        else if(real[i] < -0.9)
        {
        	demod_bits[2*i] = 7;
        }
        else
        {
            return -2;
        }
        if (imag[i] >= 0.9)
                {
                  demod_bits[2*i+1] = 0;
                }
                else if ((imag[i] > 0.9) && (imag[i] >= 0.6))
                {
                  demod_bits[2*i+1] = 1;
                }
                else if ((imag[i] > 0.6) && (imag[i] >= 0.3))
                {
                  demod_bits[2*i+1] = 2;
                }
                else if ((imag[i] > 0.3) && (imag[i] >= 0))
                {
                  demod_bits[2*i+1] = 3;
                }
                else if ((imag[i] < 0) && (imag[i] >= -0.3))
                {
                  demod_bits[2*i+1] = 4;
                }
                else if ((imag[i] < -0.3) && (imag[i] >= -0.6))
                {
                  demod_bits[2*i+1] = 5;
                }
                else if ((imag[i] < -0.6) && (imag[i] >= -0.9))
                {
                  demod_bits[2*i+1] = 6;
                }
                else if(imag[i] < -0.9)
                {
                	demod_bits[2*i+1] = 7;
                }
                else
                {
                    return -2;
                }
    }
return Msymbol;
}
void main()
{

int i,Msymbol = 6;

double N0 =1;
uint32_t soft_bits[12];
float complex rx_symbols[6] = { 1+1*I,0+-1.5*I, .1-1*I,3+3*I,-I,-0.1};
pbch_softdemod1(Msymbol,rx_symbols,soft_bits);
for(i = 0;i< 2*Msymbol; i++)
{
printf("SOFTBITS[%d] is %d\n\n", i,soft_bits[i] );
//printf("Symbols[%d] is %f + i %f\n\n", i,creal(rx_symbols[i]),cimag(rx_symbols[i]) );

}
}

