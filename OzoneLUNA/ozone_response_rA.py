xu_rA, xu_rA_sigma = ratio(xu_A_o3, xu_A_cf, xu_A_o3_sigma, xu_A_cf_sigma)
pelle_rA, pelle_rA_sigma = ratio(pelle_A_o3[2::3], pelle_A_cf[2::3], pelle_A_o3_sigma[2::3], pelle_A_cf_sigma[2::3])
watanabe_rA, watanabe_rA_sigma = ratio(watanabe_A_o3[watanabe_A_o3.index.str.match("OO_")],
                                         watanabe_A_o3[watanabe_A_o3.index.str.match("CC_")],
                                         watanabe_A_o3_sigma[watanabe_A_o3_sigma.index.str.match("OO_")],
                                         watanabe_A_o3_sigma[watanabe_A_o3_sigma.index.str.match("CC_")])
watanabe_rA_oc, watanabe_rA_sigma_oc = ratio(watanabe_A_o3[watanabe_A_o3.index.str.match("OC_")],
                                               watanabe_A_o3[watanabe_A_o3.index.str.match("CC_")],
                                               watanabe_A_o3_sigma[watanabe_A_o3_sigma.index.str.match("OC_")],
                                               watanabe_A_o3_sigma[watanabe_A_o3_sigma.index.str.match("CC_")])
watanabe_rA_co, watanabe_rA_sigma_co = ratio(watanabe_A_o3[watanabe_A_o3.index.str.match("CO_")],
                                               watanabe_A_o3[watanabe_A_o3.index.str.match("CC_")],
                                               watanabe_A_o3_sigma[watanabe_A_o3_sigma.index.str.match("CO_")],
                                               watanabe_A_o3_sigma[watanabe_A_o3_sigma.index.str.match("CC_")])
pelle14_rA, pelle14_rA_sigma = ratio(pelle14_A_o3[2::3], pelle14_A_cf[2::3], pelle14_A_o3_sigma[2::3], pelle14_A_cf_sigma[2::3])

kinose_rA, kinose_rA_sigma = ratio(kinose_A_o3[kinose_A_o3.index.str.match("S1_")][1:],
                                     kinose_A_o3[kinose_A_o3.index.str.match("CF_")][1:],
                                     kinose_A_o3_sigma[kinose_A_o3_sigma.index.str.match("S1_")][1:],
                                     kinose_A_o3_sigma[kinose_A_o3_sigma.index.str.match("CF_")][1:])
kinose_rA_s15, kinose_rA_s15_sigma = ratio(kinose_A_o3[kinose_A_o3.index.str.match("S15_")][1:],
                                             kinose_A_o3[kinose_A_o3.index.str.match("CF_")][1:],
                                             kinose_A_o3_sigma[kinose_A_o3_sigma.index.str.match("S15_")][1:],
                                             kinose_A_o3_sigma[kinose_A_o3_sigma.index.str.match("CF_")][1:])

watanabe13_rA_beech, watanabe13_rA_beech_sigma = ratio(watanabe13_A_o3[watanabe13_A_o3.index.str.match("S_Beech_o3")],
                                     watanabe13_A_o3[watanabe13_A_o3.index.str.match("S_Beech_amb")],
                                     watanabe13_A_o3_sigma[watanabe13_A_o3_sigma.index.str.match("S_Beech_o3")],
                                     watanabe13_A_o3_sigma[watanabe13_A_o3_sigma.index.str.match("S_Beech_amb")])
watanabe13_rA_oak, watanabe13_rA_oak_sigma = ratio(watanabe13_A_o3[watanabe13_A_o3.index.str.match("S_Oak_o3")],
                                     watanabe13_A_o3[watanabe13_A_o3.index.str.match("S_Oak_amb")],
                                     watanabe13_A_o3_sigma[watanabe13_A_o3_sigma.index.str.match("S_Oak_o3")],
                                     watanabe13_A_o3_sigma[watanabe13_A_o3_sigma.index.str.match("S_Oak_amb")])

gao_rA, gao_rA_sigma = ratio(gao_A_o3[1::2],
                               gao_A_o3[0::2],
                               gao_A_o3_sigma[1::2],
                               gao_A_o3_sigma[0::2])


harmens_rA_1, harmens_rA_1_sigma = ratio(harmens_A_o3[0::2][1::2],
                               harmens_A_o3[0::2][0::2],
                               harmens_A_o3_sigma[0::2][1::2],
                               harmens_A_o3_sigma[0::2][0::2])

harmens_rA_2, harmens_rA_2_sigma = ratio(harmens_A_o3[1::2][1::2],
                               harmens_A_o3[1::2][0::2],
                               harmens_A_o3_sigma[1::2][1::2],
                               harmens_A_o3_sigma[1::2][0::2])

