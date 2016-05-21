void main ()
{


uint32_t n_prb = N_prb(SRS_UL.n_ul_rb, N_sc);

uint32_t n_RE = SRS_NRE(SRS_UL.CP);

srsbwtable_idx(SRS_UL.n_ul_rb);

uint32_t M_sc = Get_Msc_values(SRS_UL.bw_cfg, SRS_UL.B, SRS_UL.n_ul_rb, N_sc);

uint32_t N_zc = Get_Nzc( M_sc);

uint32_t u = Get_u_value( SRS_UL.cell_ID, SRS_UL.ns, SRS_UL.group_hopping, SRS_UL.N_ID_PUSCH, SRS_UL.N_ID_PUSCH, SRS_UL.PUCCH_ID, SRS_UL.PUSCH_ID);


uint32_t v =  Get_v_value(SRS_UL.cell_ID, SRS_UL.delta_ss, M_sc, N_sc, SRS_UL.ns, SRS_UL.sequence_hopping, SRS_UL.group_hopping, SRS_UL.N_ID_PUCCH, SRS_UL.N_ID_PUSCH, SRS_UL.PUCCH_ID, SRS_UL.PUSCH_ID);

const uint16_t n_ID = Get_cellID( SRS_UL.N_ID_PUCCH, SRS_UL.N_ID_PUSCH, SRS_UL.cell_ID, SRS_UL.PUCCH_ID, SRS_UL.PUSCH_ID);

float q = get_qvalue( u, v, N_zc);

compute_r_uv_arg(r_uv, SRS_UL.n_ul_rb, SRS_UL.ns, N_sc, N_zc, SRS_UL.sequence_hopping, SRS_UL.group_hopping,u,v);
float alpha = alpha_p(SRS_UL.N_Tx, SRS_UL.Cyclic_shift, SRS_UL.K_Tc);
}
