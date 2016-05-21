
#include <stdio.h>
#include <strings.h>
#include <stdlib.h>
#include <stdbool.h>
#include <assert.h>
#include <math.h>
void main()
{
    float q_hat;
    float bb,cc;
    int u=0;
    float n_sz= 71;
    q_hat = n_sz *(u + 1) / 31;
    bb=floor(q_hat+0.5);
    printf("bb=%d \n",bb);
    /*cc=(-1)^floor(2*q_hat);*/
    vv=floor(2*q_hat);
    printf("vv=%d \n",vv);
    cc=pow((-1),vv);
    printf("cc=%d \n",cc);
    q=bb+ Seq_No*cc;
    getch();
}
