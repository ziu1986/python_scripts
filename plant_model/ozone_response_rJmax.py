xu_rJmax, xu_rJmax_sigma = ratio(xu_Jmax_o3, xu_Jmax_cf, xu_Jmax_o3_sigma, xu_Jmax_cf_sigma)
pelle_rJmax, pelle_rJmax_sigma = ratio(pelle_Jmax_o3, pelle_Jmax_cf, pelle_Jmax_o3_sigma, pelle_Jmax_cf_sigma)
pelle14_rJmax, pelle14_rJmax_sigma = ratio(pelle14_Jmax_o3, pelle14_Jmax_cf, pelle14_Jmax_o3_sigma, pelle14_Jmax_cf_sigma)
kinose_rJmax, kinose_rJmax_sigma = ratio(kinose_Jmax_o3[kinose_Jmax_o3.index.str.match("S1_")][1:],
                                     kinose_Jmax_o3[kinose_Jmax_o3.index.str.match("CF_?")][1:],
                                     kinose_Jmax_o3_sigma[kinose_Jmax_o3_sigma.index.str.match("S1_")][1:],
                                     kinose_Jmax_o3_sigma[kinose_Jmax_o3_sigma.index.str.match("CF_?")][1:])
kinose_rJmax_s15, kinose_rJmax_s15_sigma = ratio(kinose_Jmax_o3[kinose_Jmax_o3.index.str.match("S15_")][1:],
                                             kinose_Jmax_o3[kinose_Jmax_o3.index.str.match("CF_?")][1:],
                                             kinose_Jmax_o3_sigma[kinose_Jmax_o3_sigma.index.str.match("S15_")][1:],
                                             kinose_Jmax_o3_sigma[kinose_Jmax_o3_sigma.index.str.match("CF_?")][1:])

watanabe13_rJmax_beech, watanabe13_rJmax_beech_sigma = ratio(watanabe13_Jmax_o3[watanabe13_Jmax_o3.index.str.match("S_Beech_o3")],
                                     watanabe13_Jmax_o3[watanabe13_Jmax_o3.index.str.match("S_Beech_amb")],
                                     watanabe13_Jmax_o3_sigma[watanabe13_Jmax_o3_sigma.index.str.match("S_Beech_o3")],
                                     watanabe13_Jmax_o3_sigma[watanabe13_Jmax_o3_sigma.index.str.match("S_Beech_amb")])
watanabe13_rJmax_oak, watanabe13_rJmax_oak_sigma = ratio(watanabe13_Jmax_o3[watanabe13_Jmax_o3.index.str.match("S_Oak_o3")],
                                     watanabe13_Jmax_o3[watanabe13_Jmax_o3.index.str.match("S_Oak_amb")],
                                     watanabe13_Jmax_o3_sigma[watanabe13_Jmax_o3_sigma.index.str.match("S_Oak_o3")],
                                     watanabe13_Jmax_o3_sigma[watanabe13_Jmax_o3_sigma.index.str.match("S_Oak_amb")])

gao_rJmax, gao_rJmax_sigma = ratio(gao_Jmax_o3[1::2],
                               gao_Jmax_o3[0::2],
                               gao_Jmax_o3_sigma[1::2],
                               gao_Jmax_o3_sigma[0::2])

harmens_rJmax_1, harmens_rJmax_1_sigma = ratio(harmens_Jmax_o3[0::2][1::2],
                               harmens_Jmax_o3[0::2][0::2],
                               harmens_Jmax_o3_sigma[0::2][1::2],
                               harmens_Jmax_o3_sigma[0::2][0::2])

harmens_rJmax_2, harmens_rJmax_2_sigma = ratio(harmens_Jmax_o3[1::2][1::2],
                               harmens_Jmax_o3[1::2][0::2],
                               harmens_Jmax_o3_sigma[1::2][1::2],
                               harmens_Jmax_o3_sigma[1::2][0::2])
