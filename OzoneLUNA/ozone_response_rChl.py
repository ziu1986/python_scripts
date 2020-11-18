xu_rChl, xu_rChl_sigma = ratio(xu_Chl_o3, xu_Chl_cf, xu_Chl_o3_sigma, xu_Chl_cf_sigma)
pelle_rChl, pelle_rChl_sigma = ratio(pelle_Chl_o3, pelle_Chl_cf, pelle_Chl_o3_sigma, pelle_Chl_cf_sigma)
watanabe13_rChl_beech, watanabe13_rChl_beech_sigma = ratio(watanabe13_Chl_o3[watanabe13_Chl_o3.index.str.match("S_Beech_o3")],
                                     watanabe13_Chl_o3[watanabe13_Chl_o3.index.str.match("S_Beech_amb")],
                                     watanabe13_Chl_o3_sigma[watanabe13_Chl_o3_sigma.index.str.match("S_Beech_o3")],
                                     watanabe13_Chl_o3_sigma[watanabe13_Chl_o3_sigma.index.str.match("S_Beech_amb")])
watanabe13_rChl_oak, watanabe13_rChl_oak_sigma = ratio(watanabe13_Chl_o3[watanabe13_Chl_o3.index.str.match("S_Oak_o3")],
                                     watanabe13_Chl_o3[watanabe13_Chl_o3.index.str.match("S_Oak_amb")],
                                     watanabe13_Chl_o3_sigma[watanabe13_Chl_o3_sigma.index.str.match("S_Oak_o3")],
                                     watanabe13_Chl_o3_sigma[watanabe13_Chl_o3_sigma.index.str.match("S_Oak_amb")])

gao_rChl, gao_rChl_sigma = ratio(gao_Chl_o3[1::2],
                               gao_Chl_o3[0::2],
                               gao_Chl_o3_sigma[1::2],
                               gao_Chl_o3_sigma[0::2])
