xu_rVcmax, xu_rVcmax_sigma = ratio(xu_Vcmax_o3, xu_Vcmax_cf, xu_Vcmax_o3_sigma, xu_Vcmax_cf_sigma)
pelle_rVcmax, pelle_rVcmax_sigma = ratio(pelle_Vcmax_o3, pelle_Vcmax_cf, pelle_Vcmax_o3_sigma, pelle_Vcmax_cf_sigma)
watanabe_rVcmax, watanabe_rVcmax_sigma = ratio(watanabe_Vcmax_o3[watanabe_Vcmax_o3.index.str.match("OO_?")], watanabe_Vcmax_o3[watanabe_Vcmax_o3.index.str.match("CC_?")], watanabe_Vcmax_o3_sigma[watanabe_Vcmax_o3_sigma.index.str.match("OO_?")], watanabe_Vcmax_o3_sigma[watanabe_Vcmax_o3_sigma.index.str.match("CC_?")])
watanabe_rVcmax_oc, watanabe_rVcmax_sigma_oc = ratio(watanabe_Vcmax_o3[watanabe_Vcmax_o3.index.str.match("OC_?")], watanabe_Vcmax_o3[watanabe_Vcmax_o3.index.str.match("CC_?")], watanabe_Vcmax_o3_sigma[watanabe_Vcmax_o3_sigma.index.str.match("OC_?")], watanabe_Vcmax_o3_sigma[watanabe_Vcmax_o3_sigma.index.str.match("CC_?")])
watanabe_rVcmax_co, watanabe_rVcmax_sigma_co = ratio(watanabe_Vcmax_o3[watanabe_Vcmax_o3.index.str.match("CO_?")], watanabe_Vcmax_o3[watanabe_Vcmax_o3.index.str.match("CC_?")], watanabe_Vcmax_o3_sigma[watanabe_Vcmax_o3_sigma.index.str.match("CO_?")], watanabe_Vcmax_o3_sigma[watanabe_Vcmax_o3_sigma.index.str.match("CC_?")])
pelle14_rVcmax, pelle14_rVcmax_sigma = ratio(pelle14_Vcmax_o3, pelle14_Vcmax_cf, pelle14_Vcmax_o3_sigma, pelle14_Vcmax_cf_sigma)
kinose_rVcmax, kinose_rVcmax_sigma = ratio(kinose_Vcmax_o3[kinose_Vcmax_o3.index.str.match("S1_")][1:],
                                     kinose_Vcmax_o3[kinose_Vcmax_o3.index.str.match("CF_?")][1:],
                                     kinose_Vcmax_o3_sigma[kinose_Vcmax_o3_sigma.index.str.match("S1_")][1:],
                                     kinose_Vcmax_o3_sigma[kinose_Vcmax_o3_sigma.index.str.match("CF_?")][1:])
kinose_rVcmax_s15, kinose_rVcmax_s15_sigma = ratio(kinose_Vcmax_o3[kinose_Vcmax_o3.index.str.match("S15_")][1:],
                                             kinose_Vcmax_o3[kinose_Vcmax_o3.index.str.match("CF_?")][1:],
                                             kinose_Vcmax_o3_sigma[kinose_Vcmax_o3_sigma.index.str.match("S15_")][1:],
                                             kinose_Vcmax_o3_sigma[kinose_Vcmax_o3_sigma.index.str.match("CF_?")][1:])

watanabe13_rVcmax_beech, watanabe13_rVcmax_beech_sigma = ratio(watanabe13_Vcmax_o3[watanabe13_Vcmax_o3.index.str.match("S_Beech_o3")],
                                     watanabe13_Vcmax_o3[watanabe13_Vcmax_o3.index.str.match("S_Beech_amb")],
                                     watanabe13_Vcmax_o3_sigma[watanabe13_Vcmax_o3_sigma.index.str.match("S_Beech_o3")],
                                     watanabe13_Vcmax_o3_sigma[watanabe13_Vcmax_o3_sigma.index.str.match("S_Beech_amb")])
watanabe13_rVcmax_oak, watanabe13_rVcmax_oak_sigma = ratio(watanabe13_Vcmax_o3[watanabe13_Vcmax_o3.index.str.match("S_Oak_o3")],
                                     watanabe13_Vcmax_o3[watanabe13_Vcmax_o3.index.str.match("S_Oak_amb")],
                                     watanabe13_Vcmax_o3_sigma[watanabe13_Vcmax_o3_sigma.index.str.match("S_Oak_o3")],
                                     watanabe13_Vcmax_o3_sigma[watanabe13_Vcmax_o3_sigma.index.str.match("S_Oak_amb")])


gao_rVcmax, gao_rVcmax_sigma = ratio(gao_Vcmax_o3[1::2],
                               gao_Vcmax_o3[0::2],
                               gao_Vcmax_o3_sigma[1::2],
                               gao_Vcmax_o3_sigma[0::2])

