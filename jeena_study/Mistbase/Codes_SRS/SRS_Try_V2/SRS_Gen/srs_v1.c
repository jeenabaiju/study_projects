//Table 5.2.3-1: Resource block parameters


uint32_t m_srs_b[4][4][8] = {{
                        /* m_srs for 6<n_rb<40. Table 5.5.3.2-1 */
                        {36, 32, 24, 20, 16, 12, 8, 4},
                        {12, 16,  4,  4,  4,  4, 4, 4},
                        { 4,  8,  4,  4,  4,  4, 4, 4},
                        { 4,  4,  4,  4,  4,  4, 4, 4}},
                        {
                        /* m_srs for 40<n_rb<60. Table 5.5.3.2-2 */
                        {48, 48, 40, 36, 32, 24, 20, 16},
                        {24, 16, 20, 12, 16,  4,  4,  4},
                        {12,  8,  4,  4,  8,  4,  4,  4},
                        { 4,  4,  4,  4,  4,  4,  4,  4}},
                        {
                        /* m_srs for 60<n_rb<80. Table 5.5.3.2-3 */
                        {72, 64, 60, 48, 48, 40, 36, 32},
                        {24, 32, 20, 24, 16, 20, 12, 16},
                        {12, 16,  4, 12,  8,  4,  4,  8},
                        { 4,  4,  4,  4,  4,  4,  4,  4}},

                        {
                        /* m_srs for 80<n_rb<110. Table 5.5.3.2-4 */
                        {96, 96, 80, 72, 64, 60, 48, 48},
                        {48, 32, 40, 24, 32, 20, 24, 16},
                        {24, 16, 20, 12, 16,  4, 12,  8},
                        { 4,  4,  4,  4,  4,  4,  4,  4}}};

/* Same tables for Nb */
uint32_t Nb[4][4][8] = {{
                        {1, 1, 1, 1, 1, 1, 1, 1},
                        {3, 2, 6, 5, 4, 3, 2, 1},
                        {3, 2, 1, 1, 1, 1, 1, 1},
                        {1, 2, 1, 1, 1, 1, 1, 1}},
                        {
                        {1, 1, 1, 1, 1, 1, 1, 1},
                        {2, 3, 2, 3, 2, 6, 5, 4},
                        {2, 2, 5, 3, 2, 1, 1, 1},
                        {3, 2, 1, 1, 2, 1, 1, 1}},
                        {
                        {1, 1, 1, 1, 1, 1, 1, 1},
                        {3, 2, 3, 2, 3, 2, 3, 2},
                        {2, 2, 5, 2, 2, 5, 3, 2},
                        {3, 4, 1, 3, 2, 1, 1, 2}},
                        {
                        {1, 1, 1, 1, 1, 1, 1, 1},
                        {2, 3, 2, 3, 2, 3, 2, 3},
                        {2, 2, 2, 2, 2, 5, 2, 2},
                        {6, 4, 5, 3, 4, 1, 3, 2}}};


/**************************************************************************************************/
// Calculate prime
/*
 * Calculate the value of largest prime number less than M_sc
 * input- M_sc 1<=m<=5,
 * output -N_zc largest prime number less than M_sc
 */
uint8_t Largestprime(uint8_t M_sc)
{

}
/*
 * Calculate Sequence Hopping
 * input- M_sc,N_sc,
 * output-
 */

/*
 * Calculate Group Hopping
 * input- M_sc ,N_sc=no.ofsubcarriers,Seq_no, N_UL_RB
 * output -phi for 0:11 and 0:23
 */
// Base sequence 3Nsc greater
// Generate phi sequence
/*
 * Calculate the value of phi for M_sc=N_sc and M_sc=2*N_sc
 * input- M_sc ,N_sc=no.ofsubcarriers,Seq_no, N_UL_RB
 * output -phi for 0:11 and 0:23
 */
 int8_t* Generate_Phi(uint8_t M_sc, uint8_t N_sc, uint8_t Seq_No, uint8_t N_ul_rb)
 {

 }
/*
 * Calculate q value
 * input- M_sc ,N_sc=no.ofsubcarriers,Grp_no,Seq_no
 * output -q value to find Zadoff chu seq
 */

 static uint32_t get_qvalue(uint32_t Grp_no, uint32_t Seq_no, uint32_t N_zc) {
  float q;
  float q_bar;
  float n_zc = (float) N_zc;
  q_bar= n_zc*(Grp_no+1)/31;
  q=(floor(q_bar+0.5))+(Seq_no*(pow((-1),(floor(2*q_bar)))));
  return (uint32_t) q;
}
/*************************************************************************************************/
// Calculate argument for Qth root Zadoff-Chu sequence as in according to 3GPP 36.211 5.5.1.1 to find R_ZC
/*
 * Calculate the value of qth root for ZC
 * input- q value , N_zc
 * output root_q for exponential calculation
 */
static void Root_q_arg(float *root_q_arg,uint32_t M_sc, uint32_t Grp_no, uint32_t Seq_no)
{
  float m;
  N_zc=Largestprime( M_sc);
  float n_zc = (float) N_zc;
  float q = get_qvalue(Grp_no,Seq_no,N_zc);
  for (m=0;m<N_zc;m++){
    root_q_arg[i]=-M_PI * q * m * (m + 1) / n_zc;/* argument of x_q(m) according to 3GPP 36.211 5.5.1.1
  }

/***********************************************************************************************/



  // Calculate Exponential to find R_ZC
/*
 * Calculate the value of arguments for NULRB=1,2,<=3
 * input- NULRB,Grp no, Seq No, other values needed
 * output -argument values to find sequence */
/* Computes argument of r_u_v signal */
static void compute_r_uv_arg(uint32_t n_ul_rb, uint32_t Grp_no, uint32_t Seq_no) {
  if (n_ul_rb == 1) {
    arg_phi_Nsc(phi_M_sc_12_arg, Grp_no);
  } else if (n_ul_rb == 2) {
    arg_phi_2Nsc(phi_M_sc_24, Grp_no);
  } else {
    Root_q_arg(root_q_arg, N_sc*n_ul_rb, Grp_no, Seq_no);

  }
}



   // Calculate Exponential to find R_ZC
/*
 * Calculate the value of Exponential for ZC
 * input- argument for exponential calculation
 * output -cexpo value*/
    // Do complex exponential
    root_q= exp(I* root_q_arg);
