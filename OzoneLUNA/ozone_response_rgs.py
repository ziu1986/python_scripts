xu_rgs, xu_rgs_sigma = ratio(xu_gs_o3, xu_gs_cf, xu_gs_o3_sigma, xu_gs_cf_sigma)
pelle_rgs, pelle_rgs_sigma = ratio(pelle_gs_o3[2::3], pelle_gs_cf[2::3], pelle_gs_o3_sigma[2::3], pelle_gs_cf_sigma[2::3])
watanabe_rgs, watanabe_rgs_sigma = ratio(watanabe_gs_o3[watanabe_gs_o3.index.str.match("OO_")],
                                         watanabe_gs_o3[watanabe_gs_o3.index.str.match("CC_")],
                                         watanabe_gs_o3_sigma[watanabe_gs_o3_sigma.index.str.match("OO_")],
                                         watanabe_gs_o3_sigma[watanabe_gs_o3_sigma.index.str.match("CC_")])
watanabe_rgs_oc, watanabe_rgs_sigma_oc = ratio(watanabe_gs_o3[watanabe_gs_o3.index.str.match("OC_")],
                                               watanabe_gs_o3[watanabe_gs_o3.index.str.match("CC_")],
                                               watanabe_gs_o3_sigma[watanabe_gs_o3_sigma.index.str.match("OC_")],
                                               watanabe_gs_o3_sigma[watanabe_gs_o3_sigma.index.str.match("CC_")])
watanabe_rgs_co, watanabe_rgs_sigma_co = ratio(watanabe_gs_o3[watanabe_gs_o3.index.str.match("CO_")],
                                               watanabe_gs_o3[watanabe_gs_o3.index.str.match("CC_")],
                                               watanabe_gs_o3_sigma[watanabe_gs_o3_sigma.index.str.match("CO_")],
                                               watanabe_gs_o3_sigma[watanabe_gs_o3_sigma.index.str.match("CC_")])
pelle14_rgs, pelle14_rgs_sigma = ratio(pelle14_gs_o3[2::3], pelle14_gs_cf[2::3], pelle14_gs_o3_sigma[2::3], pelle14_gs_cf_sigma[2::3])

kinose_rgs, kinose_rgs_sigma = ratio(kinose_gs_o3[kinose_gs_o3.index.str.match("S1_")][1:],
                                     kinose_gs_o3[kinose_gs_o3.index.str.match("CF_")][1:],
                                     kinose_gs_o3_sigma[kinose_gs_o3_sigma.index.str.match("S1_")][1:],
                                     kinose_gs_o3_sigma[kinose_gs_o3_sigma.index.str.match("CF_")][1:])
kinose_rgs_s15, kinose_rgs_s15_sigma = ratio(kinose_gs_o3[kinose_gs_o3.index.str.match("S15_")][1:],
                                             kinose_gs_o3[kinose_gs_o3.index.str.match("CF_")][1:],
                                             kinose_gs_o3_sigma[kinose_gs_o3_sigma.index.str.match("S15_")][1:],
                                             kinose_gs_o3_sigma[kinose_gs_o3_sigma.index.str.match("CF_")][1:])

watanabe13_rgs_beech, watanabe13_rgs_beech_sigma = ratio(watanabe13_gs_o3[watanabe13_gs_o3.index.str.match("S_Beech_o3")],
                                     watanabe13_gs_o3[watanabe13_gs_o3.index.str.match("S_Beech_amb")],
                                     watanabe13_gs_o3_sigma[watanabe13_gs_o3_sigma.index.str.match("S_Beech_o3")],
                                     watanabe13_gs_o3_sigma[watanabe13_gs_o3_sigma.index.str.match("S_Beech_amb")])
watanabe13_rgs_oak, watanabe13_rgs_oak_sigma = ratio(watanabe13_gs_o3[watanabe13_gs_o3.index.str.match("S_Oak_o3")],
                                     watanabe13_gs_o3[watanabe13_gs_o3.index.str.match("S_Oak_amb")],
                                     watanabe13_gs_o3_sigma[watanabe13_gs_o3_sigma.index.str.match("S_Oak_o3")],
                                     watanabe13_gs_o3_sigma[watanabe13_gs_o3_sigma.index.str.match("S_Oak_amb")])

gao_rgs, gao_rgs_sigma = ratio(gao_gs_o3[1::2],
                               gao_gs_o3[0::2],
                               gao_gs_o3_sigma[1::2],
                               gao_gs_o3_sigma[0::2])


harmens_rgs_1, harmens_rgs_1_sigma = ratio(harmens_gs_o3[0::2][1::2],
                               harmens_gs_o3[0::2][0::2],
                               harmens_gs_o3_sigma[0::2][1::2],
                               harmens_gs_o3_sigma[0::2][0::2])

harmens_rgs_2, harmens_rgs_2_sigma = ratio(harmens_gs_o3[1::2][1::2],
                               harmens_gs_o3[1::2][0::2],
                               harmens_gs_o3_sigma[1::2][1::2],
                               harmens_gs_o3_sigma[1::2][0::2])

