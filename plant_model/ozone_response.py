import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from sample_from_norm import compute_cuo
from mytools.met_tools import print_all

def ratio(x1, x2, sigma_x1, sigma_x2):
    rg = x1.values/x2.values
    rg_sigma = np.sqrt((1/x2.values*sigma_x1.values)**2 + (x1.values/x2.values**2*sigma_x2.values)**2)
    return(rg, rg_sigma)

def cuo_watanabe(o3_mu, o3_sigma, gs_o3, gs_o3_sigma, o3_fumi, o3_days, **kwarg):
    exp = kwarg.pop("exp", 'CC')
    pcuo = []
    pcuo_std = []
    cuo_mean, cuo_std = compute_cuo(o3_mu['%s_1' % exp], o3_sigma['%s_1' % exp], gs_o3['%s_1' % exp], gs_o3_sigma['%s_1' % exp], int(o3_fumi['%s_1' % exp]), o3_days['%s_1' % exp])
    pcuo.append(cuo_mean)
    pcuo_std.append(cuo_std)  
    cuo_mean, cuo_std = compute_cuo(o3_mu['%s_1.5' % exp], o3_sigma['%s_1.5' % exp], (gs_o3['%s_1' % exp]+gs_o3['%s_2' % exp])*0.5, np.sqrt(gs_o3_sigma['%s_1' % exp]**2+gs_o3_sigma['%s_2' % exp]**2)*0.5, int(o3_fumi['%s_1.5' % exp]), o3_days['%s_1.5' % exp]-o3_days['%s_1' % exp])
    cuo_mean = cuo_mean+pcuo[-1]
    cuo_std = np.sqrt(cuo_std**2+pcuo_std[-1]**2)
    pcuo.append(cuo_mean)
    pcuo_std.append(cuo_std)
    cuo_mean, cuo_std = compute_cuo(o3_mu['%s_2' % exp], o3_sigma['%s_2' % exp], gs_o3['%s_2' % exp], gs_o3_sigma['%s_2' % exp], int(o3_fumi['%s_2' % exp]), (o3_days['%s_2' % exp]-o3_days['%s_1.5' % exp]))
    cuo_mean = cuo_mean+pcuo[-1]
    cuo_std = np.sqrt(cuo_std**2+pcuo_std[-1]**2)
    pcuo.append(cuo_mean)
    pcuo_std.append(cuo_std)
    cuo_mean, cuo_std = compute_cuo(o3_mu['%s_3' % exp], o3_sigma['%s_3' % exp], gs_o3['%s_3' % exp], gs_o3_sigma['%s_3' % exp], int(o3_fumi['%s_3' % exp]), (o3_days['%s_3' % exp]-o3_days['%s_2' % exp]))
    cuo_mean = cuo_mean+pcuo[-1]
    cuo_std = np.sqrt(cuo_std**2+pcuo_std[-1]**2)
    pcuo.append(cuo_mean)
    pcuo_std.append(cuo_std)
    cuo_mean, cuo_std = compute_cuo(o3_mu['%s_4' % exp], o3_sigma['%s_4' % exp], gs_o3['%s_4' % exp], gs_o3_sigma['%s_4' % exp], int(o3_fumi['%s_4' % exp]), (o3_days['%s_4' % exp]-o3_days['%s_3' % exp]))
    cuo_mean = cuo_mean+pcuo[-1]
    cuo_std = np.sqrt(cuo_std**2+pcuo_std[-1]**2)
    pcuo.append(cuo_mean)
    pcuo_std.append(cuo_std)

    pcuo = np.array(pcuo)
    pcuo_std = np.array(pcuo_std)

    return(pcuo,pcuo_std)


try:
    kinose_o3_mu
except NameError:
    execfile("ozone_response_read_data.py")
            
# Clean up
plt.close('all')
   

# Plot data
xu_pcuo = []
xu_pcuo_std = []
for i in range(4):
    cuo_mean, cuo_std = compute_cuo(xu_o3_mu[i], xu_o3_sigma[i], xu_gs_o3[i], xu_gs_o3_sigma[i], 10, xu_o3_days[i])
    xu_pcuo.append(cuo_mean)
    xu_pcuo_std.append(cuo_std)

pelle_pcuo = []
pelle_pcuo_std = []
for i in range(1,6):
    cuo_mean, cuo_std = compute_cuo(pelle_o3_mu[i], pelle_o3_sigma[i], (pelle_gs_o3[i]-pelle_gs_o3[i-1])*0.5+pelle_gs_o3[i-1], 0.5*np.sqrt(pelle_gs_o3_sigma[i]**2+pelle_gs_o3_sigma[i-1]**2), 5, pelle_o3_days[i]-pelle_o3_days[i-1])
    
    if i>1:
       cuo_mean = cuo_mean+pelle_pcuo[-1]
       cuo_std = np.sqrt(cuo_std**2+pelle_pcuo_std[-1]**2)
    pelle_pcuo.append(cuo_mean)
    pelle_pcuo_std.append(cuo_std)    

watanabe_pcuo_cf, watanabe_pcuo_std_cf = cuo_watanabe(watanabe_o3_mu, watanabe_o3_sigma, watanabe_gs_o3, watanabe_gs_o3_sigma, watanabe_o3_fumi, watanabe_o3_days, exp='CC')

watanabe_pcuo, watanabe_pcuo_std = cuo_watanabe(watanabe_o3_mu, watanabe_o3_sigma, watanabe_gs_o3, watanabe_gs_o3_sigma, watanabe_o3_fumi, watanabe_o3_days, exp='OO')
watanabe_pcuo_oc, watanabe_pcuo_std_oc = cuo_watanabe(watanabe_o3_mu, watanabe_o3_sigma, watanabe_gs_o3, watanabe_gs_o3_sigma, watanabe_o3_fumi, watanabe_o3_days, exp='OC')
watanabe_pcuo_co, watanabe_pcuo_std_co = cuo_watanabe(watanabe_o3_mu, watanabe_o3_sigma, watanabe_gs_o3, watanabe_gs_o3_sigma, watanabe_o3_fumi, watanabe_o3_days, exp='CO')

pelle14_pcuo = []
pelle14_pcuo_std = []
pelle14_pcuo1 = []
pelle14_pcuo1_std = []
pelle14_pcuo2 = []
pelle14_pcuo2_std = []
for i in range(1,6):
    cuo_mean, cuo_std = compute_cuo(pelle14_o3_mu[i], pelle14_o3_sigma[i], (pelle14_gs_o3[i]-pelle14_gs_o3[i-1])*0.5+pelle14_gs_o3[i-1], 0.5*np.sqrt(pelle14_gs_o3_sigma[i]**2+pelle14_gs_o3_sigma[i-1]**2), 5, pelle14_o3_days[i]-pelle14_o3_days[i-1])
    cuo1_mean, cuo1_std = compute_cuo(pelle14_o3_mu[i], pelle14_o3_sigma[i], pelle14_gs_o3[i], pelle14_gs_o3_sigma[i], 5, pelle14_o3_days[i]-pelle14_o3_days[i-1])
    cuo2_mean, cuo2_std = compute_cuo(pelle14_o3_mu[i], pelle14_o3_sigma[i], pelle14_gs_o3[i-1], pelle14_gs_o3_sigma[i-1], 5, pelle14_o3_days[i]-pelle14_o3_days[i-1])
    if i>1:
       cuo_mean = cuo_mean+pelle14_pcuo[-1]
       cuo_std = np.sqrt(cuo_std**2+pelle14_pcuo_std[-1]**2)
       cuo1_mean = cuo1_mean+pelle14_pcuo1[-1]
       cuo1_std = np.sqrt(cuo1_std**2+pelle14_pcuo1_std[-1]**2)
       cuo2_mean = cuo2_mean+pelle14_pcuo2[-1]
       cuo2_std = np.sqrt(cuo2_std**2+pelle14_pcuo2_std[-1]**2)
       
    pelle14_pcuo.append(cuo_mean)
    pelle14_pcuo_std.append(cuo_std)
    pelle14_pcuo1.append(cuo1_mean)
    pelle14_pcuo1_std.append(cuo1_std)
    pelle14_pcuo2.append(cuo2_mean)
    pelle14_pcuo2_std.append(cuo2_std)

kinose_pcuo = []
kinose_pcuo_std = []
for j in range(3):
    o3_mu = kinose_o3_mu[j::3]
    o3_sigma = kinose_o3_sigma[j::3]
    gs_o3 = kinose_gs_o3[j::3].interpolate()
    gs_o3_sigma = kinose_gs_o3_sigma[j::3].interpolate()
    o3_days = kinose_o3_days[j::3]
        
    for i in range(1,kinose_o3_mu[0::3].size):
        cuo_mean, cuo_std = compute_cuo(o3_mu[i], o3_sigma[i], (gs_o3[i]-gs_o3[i-1])*0.5+gs_o3[i-1], 0.5*np.sqrt(gs_o3_sigma[i]**2+gs_o3_sigma[i-1]**2), 12, o3_days[i]-o3_days[i-1])
        if i>1:
            cuo_mean = cuo_mean+kinose_pcuo[-1]
            cuo_std = np.sqrt(cuo_std**2+kinose_pcuo_std[-1]**2)
        kinose_pcuo.append(cuo_mean)
        kinose_pcuo_std.append(cuo_std)

kinose_pcuo = np.array(kinose_pcuo).reshape(3,len(kinose_pcuo)/3)
kinose_pcuo_std = np.array(kinose_pcuo_std).reshape(3,len(kinose_pcuo_std)/3)

xu_rgs, xu_rgs_sigma = ratio(xu_gs_o3, xu_gs_cf, xu_gs_o3_sigma, xu_gs_cf_sigma)
pelle_rgs, pelle_rgs_sigma = ratio(pelle_gs_o3[2::3], pelle_gs_cf[2::3], pelle_gs_o3_sigma[2::3], pelle_gs_cf_sigma[2::3])
watanabe_rgs, watanabe_rgs_sigma = ratio(watanabe_gs_o3[watanabe_gs_o3.index.str.match("OO_?")],
                                         watanabe_gs_o3[watanabe_gs_o3.index.str.match("CC_?")],
                                         watanabe_gs_o3_sigma[watanabe_gs_o3_sigma.index.str.match("OO_?")],
                                         watanabe_gs_o3_sigma[watanabe_gs_o3_sigma.index.str.match("CC_?")])
watanabe_rgs_oc, watanabe_rgs_sigma_oc = ratio(watanabe_gs_o3[watanabe_gs_o3.index.str.match("OC_?")],
                                               watanabe_gs_o3[watanabe_gs_o3.index.str.match("CC_?")],
                                               watanabe_gs_o3_sigma[watanabe_gs_o3_sigma.index.str.match("OC_?")],
                                               watanabe_gs_o3_sigma[watanabe_gs_o3_sigma.index.str.match("CC_?")])
watanabe_rgs_co, watanabe_rgs_sigma_co = ratio(watanabe_gs_o3[watanabe_gs_o3.index.str.match("CO_?")],
                                               watanabe_gs_o3[watanabe_gs_o3.index.str.match("CC_?")],
                                               watanabe_gs_o3_sigma[watanabe_gs_o3_sigma.index.str.match("CO_?")],
                                               watanabe_gs_o3_sigma[watanabe_gs_o3_sigma.index.str.match("CC_?")])
pelle14_rgs, pelle14_rgs_sigma = ratio(pelle14_gs_o3[2::3], pelle14_gs_cf[2::3], pelle14_gs_o3_sigma[2::3], pelle14_gs_cf_sigma[2::3])

kinose_rgs, kinose_rgs_sigma = ratio(kinose_gs_o3[kinose_gs_o3.index.str.match("S1_")][1:],
                                     kinose_gs_o3[kinose_gs_o3.index.str.match("CF_?")][1:],
                                     kinose_gs_o3_sigma[kinose_gs_o3_sigma.index.str.match("S1_")][1:],
                                     kinose_gs_o3_sigma[kinose_gs_o3_sigma.index.str.match("CF_?")][1:])
kinose_rgs_s15, kinose_rgs_s15_sigma = ratio(kinose_gs_o3[kinose_gs_o3.index.str.match("S15_")][1:],
                                             kinose_gs_o3[kinose_gs_o3.index.str.match("CF_?")][1:],
                                             kinose_gs_o3_sigma[kinose_gs_o3_sigma.index.str.match("S15_")][1:],
                                             kinose_gs_o3_sigma[kinose_gs_o3_sigma.index.str.match("CF_?")][1:])


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

xu_rChl, xu_rChl_sigma = ratio(xu_Chl_o3, xu_Chl_cf, xu_Chl_o3_sigma, xu_Chl_cf_sigma)
pelle_rChl, pelle_rChl_sigma = ratio(pelle_Chl_o3, pelle_Chl_cf, pelle_Chl_o3_sigma, pelle_Chl_cf_sigma)

# Plot Jcmax and Vmax
fig1 = plt.figure(1, figsize=(16,9))
fig1.canvas.set_window_title("ozone_response")
ax11 = plt.subplot(231)
ax11.errorbar(xu_pcuo, xu_rgs, xerr=xu_pcuo_std, yerr=xu_rgs_sigma, ls='None', marker='o', label='Xu2019')
ax11.errorbar(pelle_pcuo[1::3], pelle_rgs, xerr=pelle_pcuo_std[1::3], yerr=pelle_rgs_sigma, ls='None', marker='d', label='Pellegrini2011')
ax11.errorbar(watanabe_pcuo-watanabe_pcuo_cf, watanabe_rgs, xerr=np.sqrt(watanabe_pcuo_std**2+watanabe_pcuo_std_cf**2), yerr=watanabe_rgs_sigma, ls='None', marker='s', label='Watanabe2014')
ax11.errorbar(watanabe_pcuo-watanabe_pcuo_oc, watanabe_rgs_oc, xerr=np.sqrt(watanabe_pcuo_std**2+watanabe_pcuo_std_oc**2), yerr=watanabe_rgs_sigma_oc, ls='None', marker='s', color='black')
ax11.errorbar(watanabe_pcuo-watanabe_pcuo_co, watanabe_rgs_co, xerr=np.sqrt(watanabe_pcuo_std**2+watanabe_pcuo_std_co**2), yerr=watanabe_rgs_sigma_co, ls='None', marker='s', color='black')
ax11.errorbar(pelle14_pcuo[1::3], pelle14_rgs, xerr=pelle14_pcuo_std[1::3], yerr=pelle14_rgs_sigma, ls='None', marker='d', label='Pellegrini2014')
ax11.errorbar(np.take(kinose_pcuo, (2,4,7)), kinose_rgs, xerr=np.take(kinose_pcuo_std, (2,4,7)), yerr=kinose_rgs_sigma, ls='None', marker='^', label='Kinose2019')
ax11.errorbar(np.take(kinose_pcuo, (2,4,7)), kinose_rgs_s15, xerr=np.take(kinose_pcuo_std, (2,4,7)), yerr=kinose_rgs_s15_sigma, ls='None', marker='^', color=ax11.lines[-1].get_color())

ax12 = plt.subplot(232)
ax12.errorbar(xu_pcuo, xu_rJmax, xerr=xu_pcuo_std, yerr=xu_rJmax_sigma, ls='None', marker='o')
ax12.errorbar(pelle_pcuo[1::3], pelle_rJmax, xerr=pelle_pcuo_std[1::3], yerr=pelle_rJmax_sigma, ls='None', marker='d')
ax12.plot(-1,-1, ls='None', marker='s')
ax12.errorbar(pelle14_pcuo[1::3], pelle14_rJmax, xerr=pelle14_pcuo_std[1::3], yerr=pelle14_rJmax_sigma, ls='None', marker='d')
ax12.errorbar(np.take(kinose_pcuo, (2,4,7)), kinose_rJmax, xerr=np.take(kinose_pcuo_std, (2,4,7)), yerr=kinose_rJmax_sigma, ls='None', marker='^', label='Kinose2019')
ax12.errorbar(np.take(kinose_pcuo, (2,4,7)), kinose_rJmax_s15, xerr=np.take(kinose_pcuo_std, (2,4,7)), yerr=kinose_rJmax_s15_sigma, ls='None', marker='^', color=ax11.lines[-1].get_color())

ax13 = plt.subplot(233)
ax13.errorbar(xu_pcuo, xu_rVcmax, xerr=xu_pcuo_std, yerr=xu_rVcmax_sigma, ls='None', marker='o')
ax13.errorbar(pelle_pcuo[1::3], pelle_rVcmax, xerr=pelle_pcuo_std[1::3], yerr=pelle_rVcmax_sigma, ls='None', marker='d')
ax13.errorbar(watanabe_pcuo-watanabe_pcuo_cf, watanabe_rVcmax, xerr=np.sqrt(watanabe_pcuo_std**2+watanabe_pcuo_std_cf**2), yerr=watanabe_rVcmax_sigma, ls='None', marker='s')
ax13.errorbar(watanabe_pcuo-watanabe_pcuo_oc, watanabe_rVcmax, xerr=np.sqrt(watanabe_pcuo_std**2+watanabe_pcuo_std_oc**2), yerr=watanabe_rVcmax_sigma_oc, ls='None', marker='s', color='black')
ax13.errorbar(watanabe_pcuo-watanabe_pcuo_co, watanabe_rVcmax, xerr=np.sqrt(watanabe_pcuo_std**2+watanabe_pcuo_std_co**2), yerr=watanabe_rVcmax_sigma_co, ls='None', marker='s', color='black')
ax13.errorbar(pelle14_pcuo[1::3], pelle14_rVcmax, xerr=pelle14_pcuo_std[1::3], yerr=pelle14_rVcmax_sigma, ls='None', marker='d')
ax13.errorbar(np.take(kinose_pcuo, (2,4,7)), kinose_rVcmax, xerr=np.take(kinose_pcuo_std, (2,4,7)), yerr=kinose_rVcmax_sigma, ls='None', marker='^', label='Kinose2019')
ax13.errorbar(np.take(kinose_pcuo, (2,4,7)), kinose_rVcmax_s15, xerr=np.take(kinose_pcuo_std, (2,4,7)), yerr=kinose_rVcmax_s15_sigma, ls='None', marker='^', color=ax11.lines[-1].get_color())

ax14 = plt.subplot(234)
ax14.errorbar(xu_pcuo, xu_rRd, xerr=xu_pcuo_std, yerr=xu_rRd_sigma, ls='None', marker='o')
ax14.errorbar(pelle_pcuo[1::3], pelle_rRd, xerr=pelle_pcuo_std[1::3], yerr=pelle_rRd_sigma, ls='None', marker='d')
ax14.plot(-1,-1, ls='None', marker='s')
ax14.errorbar(pelle14_pcuo[1::3], pelle14_rRd, xerr=pelle14_pcuo_std[1::3], yerr=pelle14_rRd_sigma, ls='None', marker='d')
ax14.errorbar(np.take(kinose_pcuo, (2,4,7)), kinose_rRd, xerr=np.take(kinose_pcuo_std, (2,4,7)), yerr=kinose_rRd_sigma, ls='None', marker='^', label='Kinose2019')
ax14.errorbar(np.take(kinose_pcuo, (2,4,7)), kinose_rRd_s15, xerr=np.take(kinose_pcuo_std, (2,4,7)), yerr=kinose_rRd_s15_sigma, ls='None', marker='^', color=ax11.lines[-1].get_color())

ax15 = plt.subplot(235)
ax15.errorbar(xu_pcuo, xu_rChl, xerr=xu_pcuo_std, yerr=xu_rChl_sigma, ls='None', marker='o')
ax15.errorbar([0,]+pelle_pcuo, pelle_rChl, xerr=[0,]+pelle_pcuo_std, yerr=pelle_rChl_sigma, ls='None', marker='d')

ax14.set_xlabel("CUO (mmol $m^{-2}$)",x=1)
ax11.set_ylabel("$g_s^{O_3}/g_s^{CF}$")
ax12.set_ylabel("$J_{max}^{O_3}/J_{max}^{CF}$")
ax13.set_ylabel("$V_{cmax}^{O_3}/V_{cmax}^{CF}$")
ax14.set_ylabel("$R_{d}^{O_3}/R_{d}^{CF}$")
ax15.set_ylabel("$Chl_{a+b}^{O_3}/Chl_{a+b}^{CF}$")

ax11.legend(bbox_to_anchor=(3.05, -0.5), loc='lower right', borderaxespad=0.)

for ax in fig1.axes:
    ax.set_ylim(0.,3)
    ax.set_xlim(0.,80)

fig2 = plt.figure(2)
fig2.canvas.set_window_title("ozone_response_pellegrini2014_cuo")
ax21 = plt.subplot()
ax21.errorbar((4.36,9.31,14.41,16.57,33.8), pelle14_pcuo, yerr=pelle14_pcuo_std, color='black', label="interpol linear")
ax21.errorbar((4.36,9.31,14.41,16.57,33.8), pelle14_pcuo1, yerr=pelle14_pcuo1_std, color='red', label="endpoint gs")
ax21.errorbar((4.36,9.31,14.41,16.57,33.8), pelle14_pcuo2, yerr=pelle14_pcuo2_std, color='blue', label="startpoint gs")
ax21.set_xlabel("$CUO_{lom}$ (mmol $m^{-2}$)")
ax21.set_ylabel("$CUO_{tw}$ (mmol $m^{-2}$)")
ax21.set_xlim(0,50)
ax21.set_ylim(0,50)
ax21.plot(np.arange(0,50), np.arange(0,50), ls='--', color='black')

ax21.legend()

plt.show(block=False)

