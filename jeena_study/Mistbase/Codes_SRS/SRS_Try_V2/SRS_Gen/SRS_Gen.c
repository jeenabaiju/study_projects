/*
 *SRS_Gen.c
 *
 *  Created on: May12, 2016
 *      Author: Jeena
 */

#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <strings.h>
#include <stdlib.h>
#include <stdbool.h>
#include <assert.h>
#include <math.h>
//#include "srslte/common/phy_common.h"
/*#include < srslte/srslte.h>
#include "srslte/common/sequence.h"
#include "srslte/utils/vector.h"
#include "srslte/utils/bit.h"*/
int32_t NULRB=6;
int32_t N_SC_RB=12;
void main()
{

if(N_Tx==1)
{
uint32_t m;
for (m=1;m<=NULRB;m++)
 {
   if (m>=1&&m<=5)
    {
	uint32_t M_SC_RS[NULRB];
       M_SC_RS(m)=m*N_SC_RB;
       v(m)=0;
//      N_zc_RS(m)=max(primes((M_SC_RS(m))));
    }
   else if(m>=6&&m<=NULRB)
         {
          M_SC_RS(m)=m*N_SC_RB;
//        N_zc_RS(m)=max(primes((M_SC_RS(m))));
          else if(strcmp(ue.Hopping,'Off')==1){
           v(m)=0;
	   }
          else
	  {
           v(m)=1;
          }
         }

 }
uint32_t k=1;
uint32_t u=1;
uint32_t l=1;
int Phi[30];
float root_q[30];
uint32_t M_zc[NULRB];
uint32_t kk=1;
float RootSeq_q[30,4];
 float SRSSeq=cell(1,4);
 /*SRSSeq_Nsc=cell(1,4);
 SRSSeq_2Nsc=cell(1,4);*/
float Alpha= 2*pi*(N_SRS/8);// cyclic shift
// double Qth_root[];
for (m=1;m<=NULRB;m++){
     if (M_SC_RS(m)>=3*N_SC_RB)
       {
       M_zc(k)=M_SC_RS(m);
       N_zc_RS(k)=max(primes(M_zc(k)));
        uint32_t Seq_Group;
       for (Seq_Group=0; Seq_Group<=29;Seq_Group++){
           Q_Bar((Seq_Group+1),k)=(N_zc_RS(k)*(Seq_Group+1))/31;
           RootSeq_q(Seq_Group+1,k)=fix(Q_Bar(Seq_Group+1,k)+0.5)+v(m)*(-1)^fix(2*Q_Bar(Seq_Group+1,k)); // root q value
           }
          k=k+1;
        int pi=3.14;
           for (kk=1; kk<=k-1;kk++)
	     {
               root_q=RootSeq_q(:,kk);
               int num_padd=(M_zc(kk)-N_zc_RS(kk));
               for (Seq_Group=0; Seq_Group<=29;Seq_Group++)
	   	   {
                   for m=0;m<=N_zc_RS(kk)-1;m++){
                   double Qth_root(m+1)= double exp((-1i*pi* root_q(Seq_Group+1)*m*(m+1))/N_zc_RS(kk));
                   }
//                   r= lteZadoffChuSeq(Qth_root(Seq_Group+1),N_zc_RS(kk));
//                    plot(abs(xcorr(r)./length(r)),'r-*')

                  for (j=1;j<=M_zc(kk);j++)
		  {
                     double X_q(j)=Qth_root(1+mod(j,N_zc_RS(kk)));
//                 plot(abs(xcorr(X_q)./length(X_q)),'r-*')
                        double Seq_SRS(j)=(exp(1i*Alpha*j))*X_q(j);
                  }
                 double SRSSeq{1,kk}=Seq_SRS;// SRSSeq{1,1}(:,1)
               }
            }

      }
int Seq_Indx_phi[NULRB];
int M_normal_NSC[NULRB];
int M_normal_2NSC[NULRB];
    else if (M_SC_RS(m)==N_SC_RB)
                {
               int M_normal_NSC(u)=M_SC_RS(m);
               int Seq_Indx_phi(u)=v(m);

                for (Seq_Group=0; Seq_Group<=29;Seq_Group++)
		  {
                   phii.  erate_Phi(Seq_Group,M_normal_NSC(u),N_SC_RB,NULRB);
                   Phi(Seq_Group+1,:)=phi;

                  for (n=1;n<=M_normal_NSC(u);n++)
           	  {
                          r(n)=exp((-1i*pi*phi(n))/4);% Base sequence r(n)
                  }

                   for (j=1;j<=M_normal_NSC(u);j++)
	           {
                        double Seq_SRS_NSC=(exp(1i*Alpha*j))*r;
                   }
                //SRSSeq_Nsc{1,u}=Seq_SRS_NSC;% SRSSeq{1,1}(:,1)

                 }
                u=u+1;
               }
          else if(M_SC_RS(m)==2*N_SC_RB)
                {
                   M_normal_2NSC(l)=M_SC_RS(m);
                     Seq_Indx_phi(l)=v(m);
double Seq_phi_2N_SC[30];
                    for (Seq_Group=0; Seq_Group<=29;Seq_Group++)
		      {
                       phi= Generate_Phi(Seq_Group,M_normal_2NSC(l),N_SC_RB,NULRB);
                     //double Seq_phi_2N_SC(Seq_Group+1,:)=phi;


                       for (n=1;n<=M_normal_2NSC(l);n++)
			{
                        float r(n)=exp((-1i*pi*phi(n))/4); // Base sequence r(n)
                        }
                    // SRS_2N_SC(Seq_Group+1,:)=r;

                        for (j=1;j<=M_normal_2NSC(l);j++)
			{
                          double Seq_SRS_2NSC=(exp(1i*Alpha*j))*r;
                        }
                       //SRSSeq_2Nsc{1,l}=Seq_SRS_2NSC;//SRSSeq{1,1}

                        }
                       l=l+1;


               }

}
}
};
