xu_rRd, xu_rRd_sigma = ratio(xu_Rd_o3, xu_Rd_cf, xu_Rd_o3_sigma, xu_Rd_cf_sigma)
pelle_rRd, pelle_rRd_sigma = ratio(pelle_Rd_o3, pelle_Rd_cf, pelle_Rd_o3_sigma, pelle_Rd_cf_sigma)
pelle14_rRd, pelle14_rRd_sigma = ratio(pelle14_Rd_o3, pelle14_Rd_cf, pelle14_Rd_o3_sigma, pelle14_Rd_cf_sigma)
kinose_rRd, kinose_rRd_sigma = ratio(kinose_Rd_o3[kinose_Rd_o3.index.str.match("S1_")][1:],
                                     kinose_Rd_o3[kinose_Rd_o3.index.str.match("CF_?")][1:],
                                     kinose_Rd_o3_sigma[kinose_Rd_o3_sigma.index.str.match("S1_")][1:],
                                     kinose_Rd_o3_sigma[kinose_Rd_o3_sigma.index.str.match("CF_?")][1:])
kinose_rRd_s15, kinose_rRd_s15_sigma = ratio(kinose_Rd_o3[kinose_Rd_o3.index.str.match("S15_")][1:],
                                             kinose_Rd_o3[kinose_Rd_o3.index.str.match("CF_?")][1:],
                                             kinose_Rd_o3_sigma[kinose_Rd_o3_sigma.index.str.match("S15_")][1:],
                                             kinose_Rd_o3_sigma[kinose_Rd_o3_sigma.index.str.match("CF_?")][1:])
watanabe13_rRd_beech, watanabe13_rRd_beech_sigma = ratio(watanabe13_Rd_o3[watanabe13_Rd_o3.index.str.match("S_Beech_o3")],
                                     watanabe13_Rd_o3[watanabe13_Rd_o3.index.str.match("S_Beech_amb")],
                                     watanabe13_Rd_o3_sigma[watanabe13_Rd_o3_sigma.index.str.match("S_Beech_o3")],
                                     watanabe13_Rd_o3_sigma[watanabe13_Rd_o3_sigma.index.str.match("S_Beech_amb")])
watanabe13_rRd_oak, watanabe13_rRd_oak_sigma = ratio(watanabe13_Rd_o3[watanabe13_Rd_o3.index.str.match("S_Oak_o3")],
                                     watanabe13_Rd_o3[watanabe13_Rd_o3.index.str.match("S_Oak_amb")],
                                     watanabe13_Rd_o3_sigma[watanabe13_Rd_o3_sigma.index.str.match("S_Oak_o3")],
                                     watanabe13_Rd_o3_sigma[watanabe13_Rd_o3_sigma.index.str.match("S_Oak_amb")])
