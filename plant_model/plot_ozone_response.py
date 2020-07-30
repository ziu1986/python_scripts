# Clean up
plt.close('all')
# Plot data
# Stomatal conductance gs
fig1 = plt.figure(1, figsize=(18,9))
fig1.canvas.set_window_title("ozone_response")

ax11 = plt.subplot(231)
ax11.errorbar(xu_pcuo, xu_rgs,
              xerr=xu_pcuo_std, yerr=xu_rgs_sigma,
              ls='None', marker='o')
ax11.errorbar(pelle_pcuo[1::3], pelle_rgs,
              xerr=pelle_pcuo_std[1::3], yerr=pelle_rgs_sigma,
              ls='None', marker='d')
ax11.errorbar(watanabe_pcuo-watanabe_pcuo_cf, watanabe_rgs,
              xerr=np.sqrt(watanabe_pcuo_std**2+watanabe_pcuo_std_cf**2), yerr=watanabe_rgs_sigma,
              ls='None', marker='s')
ax11.errorbar(watanabe_pcuo-watanabe_pcuo_oc, watanabe_rgs_oc,
              xerr=np.sqrt(watanabe_pcuo_std**2+watanabe_pcuo_std_oc**2), yerr=watanabe_rgs_sigma_oc,
              ls='None', marker='s', color=ax11.lines[-1].get_color())
ax11.errorbar(watanabe_pcuo-watanabe_pcuo_co, watanabe_rgs_co,
              xerr=np.sqrt(watanabe_pcuo_std**2+watanabe_pcuo_std_co**2), yerr=watanabe_rgs_sigma_co,
              ls='None', marker='s', color=ax11.lines[-1].get_color())
ax11.errorbar(np.take(kinose_pcuo[1], (2,4,7))-np.take(kinose_pcuo[0], (2,4,7)), kinose_rgs,
              xerr=np.sqrt(np.take(kinose_pcuo_std[0], (2,4,7))**2+np.take(kinose_pcuo_std[1], (2,4,7))**2), yerr=kinose_rgs_sigma,
              ls='None', marker='^')
ax11.errorbar(pelle14_pcuo[1::3], pelle14_rgs,
              xerr=pelle14_pcuo_std[1::3], yerr=pelle14_rgs_sigma,
              ls='None', marker='d')
ax11.errorbar(np.take(kinose_pcuo[2], (2,4,7))-np.take(kinose_pcuo[0], (2,4,7)), kinose_rgs_s15,
              xerr=np.sqrt(np.take(kinose_pcuo_std[0], (2,4,7))**2+np.take(kinose_pcuo_std[2], (2,4,7))**2), yerr=kinose_rgs_s15_sigma,
              ls='None', marker='^', color=ax11.lines[-1].get_color())

ax11.errorbar(watanabe13_pcuo[1]-watanabe13_pcuo[0], watanabe13_rgs_beech,
              xerr=((np.sqrt(watanabe13_pcuo_std[0]**2+watanabe13_pcuo_std[1][0]**2),),(np.sqrt(watanabe13_pcuo_std[0]**2+watanabe13_pcuo_std[1][1]**2),)), yerr=watanabe13_rgs_beech_sigma,
              ls='None', marker='s')
ax11.errorbar(watanabe13_pcuo[3]-watanabe13_pcuo[2], watanabe13_rgs_oak,
              xerr=((np.sqrt(watanabe13_pcuo_std[2]**2+watanabe13_pcuo_std[3][0]**2),),(np.sqrt(watanabe13_pcuo_std[2]**2+watanabe13_pcuo_std[3][1]**2),)), yerr=watanabe13_rgs_oak_sigma,
              ls='None', marker='s', color=ax11.lines[-1].get_color())
ax11.errorbar(gao_pcuo[1::2]-gao_pcuo[0::2], gao_rgs,
              xerr=np.sqrt(gao_pcuo_std[0::2]**2+gao_pcuo_std[1::2]**2), yerr=gao_rgs_sigma,
              ls='None', marker='v')

ax11.errorbar(harmens_pcuo[0::2][1::2]-harmens_pcuo[0::2][0::2], harmens_rgs_1,
              xerr=np.sqrt(harmens_pcuo_std[0::2][1::2]**2+harmens_pcuo_std[0::2][0::2]**2), yerr=harmens_rgs_1_sigma,
              ls='None', marker='*')
ax11.errorbar(harmens_pcuo[1::2][1::2]-harmens_pcuo[1::2][0::2], harmens_rgs_2,
              xerr=np.sqrt(harmens_pcuo_std[1::2][1::2]**2+harmens_pcuo_std[1::2][0::2]**2), yerr=harmens_rgs_2_sigma,
              ls='None', marker='*', color=ax11.lines[-1].get_color())

ax12 = plt.subplot(232)
ax12.errorbar(xu_pcuo, xu_rA,
              xerr=xu_pcuo_std, yerr=xu_rA_sigma,
              ls='None', marker='o')
ax12.errorbar(pelle_pcuo[1::3], pelle_rA,
              xerr=pelle_pcuo_std[1::3], yerr=pelle_rA_sigma,
              ls='None', marker='d')
ax12.errorbar(watanabe_pcuo-watanabe_pcuo_cf, watanabe_rA,
              xerr=np.sqrt(watanabe_pcuo_std**2+watanabe_pcuo_std_cf**2), yerr=watanabe_rA_sigma,
              ls='None', marker='s')
ax12.errorbar(watanabe_pcuo-watanabe_pcuo_oc, watanabe_rA_oc,
              xerr=np.sqrt(watanabe_pcuo_std**2+watanabe_pcuo_std_oc**2), yerr=watanabe_rA_sigma_oc,
              ls='None', marker='s', color=ax12.lines[-1].get_color())
ax12.errorbar(watanabe_pcuo-watanabe_pcuo_co, watanabe_rA_co,
              xerr=np.sqrt(watanabe_pcuo_std**2+watanabe_pcuo_std_co**2), yerr=watanabe_rA_sigma_co,
              ls='None', marker='s', color=ax12.lines[-1].get_color())
ax12.errorbar(np.take(kinose_pcuo[1], (2,4,7))-np.take(kinose_pcuo[0], (2,4,7)), kinose_rA,
              xerr=np.sqrt(np.take(kinose_pcuo_std[0], (2,4,7))**2+np.take(kinose_pcuo_std[1], (2,4,7))**2), yerr=kinose_rA_sigma,
              ls='None', marker='^')
ax12.errorbar(pelle14_pcuo[1::3], pelle14_rA,
              xerr=pelle14_pcuo_std[1::3], yerr=pelle14_rA_sigma,
              ls='None', marker='d')
ax12.errorbar(np.take(kinose_pcuo[2], (2,4,7))-np.take(kinose_pcuo[0], (2,4,7)), kinose_rA_s15,
              xerr=np.sqrt(np.take(kinose_pcuo_std[0], (2,4,7))**2+np.take(kinose_pcuo_std[2], (2,4,7))**2), yerr=kinose_rA_s15_sigma,
              ls='None', marker='^', color=ax12.lines[-1].get_color())

ax12.errorbar(watanabe13_pcuo[1]-watanabe13_pcuo[0], watanabe13_rA_beech,
              xerr=((np.sqrt(watanabe13_pcuo_std[0]**2+watanabe13_pcuo_std[1][0]**2),),(np.sqrt(watanabe13_pcuo_std[0]**2+watanabe13_pcuo_std[1][1]**2),)), yerr=watanabe13_rA_beech_sigma,
              ls='None', marker='s')
ax12.errorbar(watanabe13_pcuo[3]-watanabe13_pcuo[2], watanabe13_rA_oak,
              xerr=((np.sqrt(watanabe13_pcuo_std[2]**2+watanabe13_pcuo_std[3][0]**2),),(np.sqrt(watanabe13_pcuo_std[2]**2+watanabe13_pcuo_std[3][1]**2),)), yerr=watanabe13_rA_oak_sigma,
              ls='None', marker='s', color=ax12.lines[-1].get_color())
ax12.errorbar(gao_pcuo[1::2]-gao_pcuo[0::2], gao_rA,
              xerr=np.sqrt(gao_pcuo_std[0::2]**2+gao_pcuo_std[1::2]**2), yerr=gao_rA_sigma,
              ls='None', marker='v')

ax12.errorbar(harmens_pcuo[0::2][1::2]-harmens_pcuo[0::2][0::2], harmens_rA_1,
              xerr=np.sqrt(harmens_pcuo_std[0::2][1::2]**2+harmens_pcuo_std[0::2][0::2]**2), yerr=harmens_rA_1_sigma,
              ls='None', marker='*')
ax12.errorbar(harmens_pcuo[1::2][1::2]-harmens_pcuo[1::2][0::2], harmens_rA_2,
              xerr=np.sqrt(harmens_pcuo_std[1::2][1::2]**2+harmens_pcuo_std[1::2][0::2]**2), yerr=harmens_rA_2_sigma,
              ls='None', marker='*', color=ax12.lines[-1].get_color())

# Jmax
ax14 = plt.subplot(234)
ax14.errorbar(xu_pcuo, xu_rJmax,
              xerr=xu_pcuo_std, yerr=xu_rJmax_sigma,
              ls='None', marker='o')
ax14.errorbar(pelle_pcuo[1::3], pelle_rJmax,
              xerr=pelle_pcuo_std[1::3], yerr=pelle_rJmax_sigma,
              ls='None', marker='d')
ax14.plot(-1,-1, ls='None', marker='s')
ax14.errorbar(np.take(kinose_pcuo[1], (2,4,7))-np.take(kinose_pcuo[0], (2,4,7)), kinose_rJmax,
              xerr=np.sqrt(np.take(kinose_pcuo_std[0], (2,4,7))**2+np.take(kinose_pcuo_std[1], (2,4,7))**2), yerr=kinose_rJmax_sigma,
              ls='None', marker='^')
ax14.errorbar(np.take(kinose_pcuo[2], (2,4,7))-np.take(kinose_pcuo[0], (2,4,7)), kinose_rJmax_s15,
              xerr=np.sqrt(np.take(kinose_pcuo_std[0], (2,4,7))**2+np.take(kinose_pcuo_std[2], (2,4,7))**2), yerr=kinose_rJmax_s15_sigma,
              ls='None', marker='^', color=ax14.lines[-1].get_color())
ax14.errorbar(pelle14_pcuo[1::3], pelle14_rJmax,
              xerr=pelle14_pcuo_std[1::3], yerr=pelle14_rJmax_sigma,
              ls='None', marker='d')
ax14.errorbar(watanabe13_pcuo[1]-watanabe13_pcuo[0], watanabe13_rJmax_beech,
              xerr=((np.sqrt(watanabe13_pcuo_std[0]**2+watanabe13_pcuo_std[1][0]**2),),(np.sqrt(watanabe13_pcuo_std[0]**2+watanabe13_pcuo_std[1][1]**2),)), yerr=watanabe13_rJmax_beech_sigma,
              ls='None', marker='s')
ax14.errorbar(watanabe13_pcuo[3]-watanabe13_pcuo[2], watanabe13_rJmax_oak,
              xerr=((np.sqrt(watanabe13_pcuo_std[2]**2+watanabe13_pcuo_std[3][0]**2),),(np.sqrt(watanabe13_pcuo_std[2]**2+watanabe13_pcuo_std[3][1]**2),)), yerr=watanabe13_rJmax_oak_sigma,
              ls='None', marker='s', color=ax14.lines[-1].get_color())
ax14.errorbar(gao_pcuo[1::2]-gao_pcuo[0::2], gao_rJmax,
              xerr=np.sqrt(gao_pcuo_std[0::2]**2+gao_pcuo_std[1::2]**2), yerr=gao_rJmax_sigma,
              ls='None', marker='v')
ax14.errorbar(harmens_pcuo[0::2][1::2]-harmens_pcuo[0::2][0::2], harmens_rJmax_1,
              xerr=np.sqrt(harmens_pcuo_std[0::2][1::2]**2+harmens_pcuo_std[0::2][0::2]**2), yerr=harmens_rJmax_1_sigma,
              ls='None', marker='*')
ax14.errorbar(harmens_pcuo[1::2][1::2]-harmens_pcuo[1::2][0::2], harmens_rJmax_2,
              xerr=np.sqrt(harmens_pcuo_std[1::2][1::2]**2+harmens_pcuo_std[1::2][0::2]**2), yerr=harmens_rJmax_2_sigma,
              ls='None', marker='*', color=ax11.lines[-1].get_color())

# Vcmax
ax15 = plt.subplot(235)
ax15.errorbar(xu_pcuo, xu_rVcmax, 
              xerr=xu_pcuo_std, yerr=xu_rVcmax_sigma, 
              ls='None', marker='o')
ax15.errorbar(pelle_pcuo[1::3], pelle_rVcmax, 
              xerr=pelle_pcuo_std[1::3], yerr=pelle_rVcmax_sigma, 
              ls='None', marker='d')
ax15.errorbar(watanabe_pcuo-watanabe_pcuo_cf, watanabe_rVcmax, 
              xerr=np.sqrt(watanabe_pcuo_std**2+watanabe_pcuo_std_cf**2), yerr=watanabe_rVcmax_sigma, 
              ls='None', marker='s')
ax15.errorbar(watanabe_pcuo-watanabe_pcuo_oc, watanabe_rVcmax_oc, 
              xerr=np.sqrt(watanabe_pcuo_std**2+watanabe_pcuo_std_oc**2), yerr=watanabe_rVcmax_sigma_oc, 
              ls='None', marker='s', color='black')
ax15.errorbar(watanabe_pcuo-watanabe_pcuo_co, watanabe_rVcmax_co, 
              xerr=np.sqrt(watanabe_pcuo_std**2+watanabe_pcuo_std_co**2), yerr=watanabe_rVcmax_sigma_co, 
              ls='None', marker='s', color='black')
ax15.errorbar(np.take(kinose_pcuo[1], (2,4,7))-np.take(kinose_pcuo[0], (2,4,7)), kinose_rVcmax, 
              xerr=np.sqrt(np.take(kinose_pcuo_std[0], (2,4,7))**2+np.take(kinose_pcuo_std[1], (2,4,7))**2), yerr=kinose_rVcmax_sigma, 
              ls='None', marker='^')
ax15.errorbar(np.take(kinose_pcuo[2], (2,4,7))-np.take(kinose_pcuo[0], (2,4,7)), kinose_rVcmax_s15, 
              xerr=np.sqrt(np.take(kinose_pcuo_std[0], (2,4,7))**2+np.take(kinose_pcuo_std[2], (2,4,7))**2), yerr=kinose_rVcmax_s15_sigma, 
              ls='None', marker='^', color=ax15.lines[-1].get_color())
ax15.errorbar(pelle14_pcuo[1::3], pelle14_rVcmax, 
              xerr=pelle14_pcuo_std[1::3], yerr=pelle14_rVcmax_sigma, 
              ls='None', marker='d')
ax15.errorbar(watanabe13_pcuo[1]-watanabe13_pcuo[0], watanabe13_rVcmax_beech, 
              xerr=((np.sqrt(watanabe13_pcuo_std[0]**2+watanabe13_pcuo_std[1][0]**2),),(np.sqrt(watanabe13_pcuo_std[0]**2+watanabe13_pcuo_std[1][1]**2),)), yerr=watanabe13_rVcmax_beech_sigma, 
              ls='None', marker='s')
ax15.errorbar(watanabe13_pcuo[3]-watanabe13_pcuo[2], watanabe13_rVcmax_oak, 
              xerr=((np.sqrt(watanabe13_pcuo_std[2]**2+watanabe13_pcuo_std[3][0]**2),),(np.sqrt(watanabe13_pcuo_std[2]**2+watanabe13_pcuo_std[3][1]**2),)), yerr=watanabe13_rVcmax_oak_sigma, 
              ls='None', marker='s', color=ax15.lines[-1].get_color())
ax15.errorbar(gao_pcuo[1::2]-gao_pcuo[0::2], gao_rVcmax,
              xerr=np.sqrt(gao_pcuo_std[0::2]**2+gao_pcuo_std[1::2]**2), yerr=gao_rVcmax_sigma,
              ls='None', marker='v')
ax15.errorbar(harmens_pcuo[0::2][1::2]-harmens_pcuo[0::2][0::2], harmens_rVcmax_1,
              xerr=np.sqrt(harmens_pcuo_std[0::2][1::2]**2+harmens_pcuo_std[0::2][0::2]**2), yerr=harmens_rVcmax_1_sigma,
              ls='None', marker='*')
ax15.errorbar(harmens_pcuo[1::2][1::2]-harmens_pcuo[1::2][0::2], harmens_rVcmax_2,
              xerr=np.sqrt(harmens_pcuo_std[1::2][1::2]**2+harmens_pcuo_std[1::2][0::2]**2), yerr=harmens_rVcmax_2_sigma,
              ls='None', marker='*', color=ax11.lines[-1].get_color())

# Dark respiration Rd
ax13 = plt.subplot(233)
ax13.errorbar(xu_pcuo, xu_rRd, 
              xerr=xu_pcuo_std, yerr=xu_rRd_sigma, 
              ls='None', marker='o')
ax13.errorbar(pelle_pcuo[1::3], pelle_rRd, 
              xerr=pelle_pcuo_std[1::3], yerr=pelle_rRd_sigma, 
              ls='None', marker='d')
ax13.plot(-1,-1, 
              ls='None', marker='s')
ax13.errorbar(np.take(kinose_pcuo[1], (2,4,7))-np.take(kinose_pcuo[0], (2,4,7)), kinose_rRd, 
              xerr=np.sqrt(np.take(kinose_pcuo_std[0], (2,4,7))**2+np.take(kinose_pcuo_std[1], (2,4,7))**2), yerr=kinose_rRd_sigma, 
              ls='None', marker='^')
ax13.errorbar(np.take(kinose_pcuo[2], (2,4,7))-np.take(kinose_pcuo[0], (2,4,7)), kinose_rRd_s15, 
              xerr=np.sqrt(np.take(kinose_pcuo_std[0], (2,4,7))**2+np.take(kinose_pcuo_std[2], (2,4,7))**2), yerr=kinose_rRd_s15_sigma, 
              ls='None', marker='^', color=ax13.lines[-1].get_color())
ax13.errorbar(pelle14_pcuo[1::3], pelle14_rRd, 
              xerr=pelle14_pcuo_std[1::3], yerr=pelle14_rRd_sigma, 
              ls='None', marker='d')
ax13.errorbar(watanabe13_pcuo[1]-watanabe13_pcuo[0], watanabe13_rRd_beech, 
              xerr=((np.sqrt(watanabe13_pcuo_std[0]**2+watanabe13_pcuo_std[1][0]**2),),(np.sqrt(watanabe13_pcuo_std[0]**2+watanabe13_pcuo_std[1][1]**2),)), yerr=watanabe13_rRd_beech_sigma, 
              ls='None', marker='s')
ax13.errorbar(watanabe13_pcuo[3]-watanabe13_pcuo[2], watanabe13_rRd_oak, 
              xerr=((np.sqrt(watanabe13_pcuo_std[2]**2+watanabe13_pcuo_std[3][0]**2),),(np.sqrt(watanabe13_pcuo_std[2]**2+watanabe13_pcuo_std[3][1]**2),)), yerr=watanabe13_rRd_oak_sigma, 
              ls='None', marker='s', color=ax13.lines[-1].get_color())
ax13.plot(-1,-1,
              ls='None', marker='v')
ax13.plot(-1,-1,
              ls='None', marker='*')

# Chlorophyll A+B content 
ax16 = plt.subplot(236)
ax16.errorbar(xu_pcuo, xu_rChl, 
              xerr=xu_pcuo_std, yerr=xu_rChl_sigma, 
              ls='None', marker='o', label='Xu2019')
ax16.errorbar([0,]+pelle_pcuo, pelle_rChl, 
              xerr=[0,]+pelle_pcuo_std, yerr=pelle_rChl_sigma, 
              ls='None', marker='d', label='Pellegrini2011')
ax16.errorbar(-1,-1,
              xerr=0.1, yerr=0.1,
              ls='None', marker='s', label='Watanabe2014')
ax16.errorbar(-1,-1,
              xerr=0.1, yerr=0.1,
              ls='None', marker='^', label='Kinose2019')
ax16.errorbar(-1,-1,
              xerr=0.1, yerr=0.1,
              ls='None', marker='d', label='Pellegrini2014')
ax16.errorbar(watanabe13_pcuo[1]-watanabe13_pcuo[0], watanabe13_rChl_beech, 
              xerr=((np.sqrt(watanabe13_pcuo_std[0]**2+watanabe13_pcuo_std[1][0]**2),),(np.sqrt(watanabe13_pcuo_std[0]**2+watanabe13_pcuo_std[1][1]**2),)), yerr=watanabe13_rChl_beech_sigma, 
              ls='None', marker='s', label='Watanabe2013')
ax16.errorbar(watanabe13_pcuo[3]-watanabe13_pcuo[2], watanabe13_rChl_oak, 
              xerr=((np.sqrt(watanabe13_pcuo_std[2]**2+watanabe13_pcuo_std[3][0]**2),),(np.sqrt(watanabe13_pcuo_std[2]**2+watanabe13_pcuo_std[3][1]**2),)), yerr=watanabe13_rChl_oak_sigma, 
              ls='None', marker='s', color=ax16.lines[-1].get_color())
ax16.errorbar(gao_pcuo[1::2]-gao_pcuo[0::2], gao_rChl,
              xerr=np.sqrt(gao_pcuo_std[0::2]**2+gao_pcuo_std[1::2]**2), yerr=gao_rChl_sigma,
              ls='None', marker='v', label='Gao2016')
ax16.errorbar(-1,-1,
              xerr=0.1, yerr=0.1,
              ls='None', marker='*', label='Harmens2016')

ax15.set_xlabel("CUO (mmol $m^{-2}$)")
ax11.set_ylabel("$g_s^{O_3}/g_s^{CF}$")
ax12.set_ylabel("$A_n^{O_3}/A_n^{CF}$")

ax13.set_ylabel("$R_{d}^{O_3}/R_{d}^{CF}$")

ax14.set_ylabel("$J_{max}^{O_3}/J_{max}^{CF}$")
ax15.set_ylabel("$V_{cmax}^{O_3}/V_{cmax}^{CF}$")

ax16.set_ylabel("$Chl_{a+b}^{O_3}/Chl_{a+b}^{CF}$")



# Fits
ax11.plot(np.arange(0,100), yfit_gs, ls='--', color='black', label="lin. fit: no uncert.")
ax11.plot(np.arange(0,100), yfit_gs2, ls=':', color='red', label="lin. fit: y uncert.")
ax11.plot(np.arange(0,100), yfit_gs_2, ls='-.', color='blue', label="lin. fit: x-y uncert.")
ax11.plot(np.arange(0,100), yfit_gs_free, ls='--', color='orange', label="lin. fit; 2 dof; non")
ax11.plot(np.arange(0,100), yfit_gs2_free, ls=':', color='orange', label="lin. fit: 2 dof; y")
ax11.plot(np.arange(0,100), yfit_gs_2_free, ls='-.', color='orange', label="lin. fit: 2 dof; x-y")


ax12.plot(np.arange(0,100), yfit_A, ls='--', color='black', label="lin. fit: no uncert.")
ax12.plot(np.arange(0,100), yfit_A2, ls=':', color='red', label="lin. fit: y uncert.")
ax12.plot(np.arange(0,100), yfit_A_2, ls='-.', color='blue', label="lin. fit: x-y uncert.")
ax12.plot(np.arange(0,100), yfit_A_free, ls='--', color='orange', label="lin. fit: 2 dof; non")
ax12.plot(np.arange(0,100), yfit_A2_free, ls=':', color='orange', label="lin. fit: 2 dof; y")
ax12.plot(np.arange(0,100), yfit_A_2_free, ls='-.', color='orange', label="lin. fit: 2 dof; x-y")


ax13.plot(np.arange(0,100), yfit_Rd, ls='--', color='black', label="lin. fit: no uncert.")
ax13.plot(np.arange(0,100), yfit_Rd2, ls=':', color='red', label="lin. fit: y uncert.")
ax13.plot(np.arange(0,100), yfit_Rd_2, ls='--', color='black', alpha=0.5, label="exp. fit: no uncert. ")
ax13.plot(np.arange(0,100), yfit_Rd2_2, ls=':', color='red', alpha=0.5, label="exp. fit: y uncert. exp")

ax14.plot(np.arange(0,100), yfit_Jmax, ls='--', color='black', label="lin. fit: no uncert.")
ax14.plot(np.arange(0,100), yfit_Jmax2, ls=':', color='red', label="lin. fit: y uncert.")
ax14.plot(np.arange(0,100), yfit_Jmax_2, ls='-.', color='blue', label="lin. fit: x-y uncert.")

ax15.plot(np.arange(0,100), yfit_Vcmax, ls='--', color='black', label="lin. fit: no uncert.")
ax15.plot(np.arange(0,100), yfit_Vcmax2, ls=':', color='red', label="lin. fit: y uncert.")
ax15.plot(np.arange(0,100), yfit_Vcmax_2, ls='-.', color='blue', label="lin. fit: x-y uncert.")


#ax16.plot(np.arange(0,100), yfit_Chl, ls='--', color='black')
#ax16.plot(np.arange(0,100), yfit_Chl2, ls=':', color='red')


for ax in fig1.axes:
    ax.set_ylim(0.,3)
    ax.set_xlim(0.,100)
    ax.legend(ncol=2)

ax16.legend(bbox_to_anchor=(-2.35, 2.225), loc='lower left', borderaxespad=0., ncol=8) # y=-0.25
'''
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
ax21.plot(np.arange(0,50), np.arange(0,50), 
              ls='--', color='black')

ax21.legend()
'''

fig3 = plt.figure(3, figsize=(10,10))
fig3.canvas.set_window_title("ozone_respons_ratios_jmax_vcmax")
ax31 = plt.subplot()

#Robs = Jmax/Vcmax

ax31.errorbar(flunder(Jmax), flunder(Vcmax), xerr=flunder(Jmax_std), yerr=flunder(Vcmax_std), ls='None')
ax31.plot(np.arange(2.1), np.arange(2.1), ls='--', color='grey')
ax31.plot(np.arange(0, 2.1, 0.1), yfit_VcmaxVSJmax, ls='--', color='black', label="lin. fit: no uncert.")
ax31.plot(np.arange(0, 2.1, 0.1), yfit_VcmaxVSJmax2, ls=':', color='red', label="lin. fit: y uncert.")
ax31.plot(np.arange(0, 2.1, 0.1), yfit_VcmaxVSJmax_2, ls='-.', color='blue', label="lin. fit: x-y uncert.")
ax31.set_xlabel("$J_{max}^{O_3}/J_{max}^{CF}$")
ax31.set_ylabel("$V_{cmax}^{O_3}/V_{cmax}^{CF}$")
ax31.set_xlim(0,1.2)
ax31.set_ylim(0,1.2)

b_con = True
# Figures for paper
if b_con:
    fig4, ((ax41, ax42, ax43), (ax51, ax52, ax53)) = plt.subplots(2,3, gridspec_kw={'width_ratios': [6, 6, 1]}, figsize=(18,9))
    fig4.canvas.set_window_title("ozone_response_fits")
else:
    fig4, (ax41, ax42, ax43) = plt.subplots(1,3, gridspec_kw={'width_ratios': [6, 6, 1]}, figsize=(18,6))
    fig4.canvas.set_window_title("ozone_response_jmax_vcmax_fits")
    
ax41.set_title("(a)", x=0.1, y=0.9, size='xx-large')
ax42.set_title("(b)", x=0.1, y=0.9, size='xx-large')
ax43.remove()

# Jmax
ax41.errorbar(xu_pcuo, xu_rJmax,
              xerr=xu_pcuo_std, yerr=xu_rJmax_sigma,
              ls='None', marker='o', color='black')
ax41.errorbar(pelle_pcuo[1::3], pelle_rJmax,
              xerr=pelle_pcuo_std[1::3], yerr=pelle_rJmax_sigma,
              ls='None', marker='d', color='black')
ax41.plot(-1,-1, ls='None', marker='s')
ax41.errorbar(np.take(kinose_pcuo[1], (2,4,7))-np.take(kinose_pcuo[0], (2,4,7)), kinose_rJmax,
              xerr=np.sqrt(np.take(kinose_pcuo_std[0], (2,4,7))**2+np.take(kinose_pcuo_std[1], (2,4,7))**2), yerr=kinose_rJmax_sigma,
              ls='None', marker='^', color='black')
ax41.errorbar(np.take(kinose_pcuo[2], (2,4,7))-np.take(kinose_pcuo[0], (2,4,7)), kinose_rJmax_s15,
              xerr=np.sqrt(np.take(kinose_pcuo_std[0], (2,4,7))**2+np.take(kinose_pcuo_std[2], (2,4,7))**2), yerr=kinose_rJmax_s15_sigma,
              ls='None', marker='^', color=ax41.lines[-1].get_color())
ax41.errorbar(harmens_pcuo[0::2][1::2]-harmens_pcuo[0::2][0::2], harmens_rJmax_1,
              xerr=np.sqrt(harmens_pcuo_std[0::2][1::2]**2+harmens_pcuo_std[0::2][0::2]**2), yerr=harmens_rJmax_1_sigma,
              ls='None', marker='*', color='black')
ax41.errorbar(harmens_pcuo[1::2][1::2]-harmens_pcuo[1::2][0::2], harmens_rJmax_2,
              xerr=np.sqrt(harmens_pcuo_std[1::2][1::2]**2+harmens_pcuo_std[1::2][0::2]**2), yerr=harmens_rJmax_2_sigma,
              ls='None', marker='*', color=ax41.lines[-1].get_color())
ax41.errorbar(gao_pcuo[1::2]-gao_pcuo[0::2], gao_rJmax,
              xerr=np.sqrt(gao_pcuo_std[0::2]**2+gao_pcuo_std[1::2]**2), yerr=gao_rJmax_sigma,
              ls='None', marker='v', color='black')

ax41.errorbar(pelle14_pcuo[1::3], pelle14_rJmax,
              xerr=pelle14_pcuo_std[1::3], yerr=pelle14_rJmax_sigma,
              ls='None', marker='d', fillstyle='none')
ax41.errorbar(watanabe13_pcuo[1]-watanabe13_pcuo[0], watanabe13_rJmax_beech,
              xerr=((np.sqrt(watanabe13_pcuo_std[0]**2+watanabe13_pcuo_std[1][0]**2),),(np.sqrt(watanabe13_pcuo_std[0]**2+watanabe13_pcuo_std[1][1]**2),)), yerr=watanabe13_rJmax_beech_sigma,
              ls='None', marker='s', fillstyle='none')
ax41.errorbar(watanabe13_pcuo[3]-watanabe13_pcuo[2], watanabe13_rJmax_oak,
              xerr=((np.sqrt(watanabe13_pcuo_std[2]**2+watanabe13_pcuo_std[3][0]**2),),(np.sqrt(watanabe13_pcuo_std[2]**2+watanabe13_pcuo_std[3][1]**2),)), yerr=watanabe13_rJmax_oak_sigma,
              ls='None', marker='s', fillstyle='none', color=ax41.lines[-1].get_color())

# Vcmax
ax42.errorbar(xu_pcuo, xu_rVcmax, 
              xerr=xu_pcuo_std, yerr=xu_rVcmax_sigma, 
              ls='None', marker='o', color='black', label='Xu et al. (2019)')
ax42.errorbar(pelle_pcuo[1::3], pelle_rVcmax, 
              xerr=pelle_pcuo_std[1::3], yerr=pelle_rVcmax_sigma, 
              ls='None', marker='d', color='black', label='Pellegrini et al. (2011)')
ax42.errorbar(watanabe_pcuo-watanabe_pcuo_cf, watanabe_rVcmax, 
              xerr=np.sqrt(watanabe_pcuo_std**2+watanabe_pcuo_std_cf**2), yerr=watanabe_rVcmax_sigma, 
              ls='None', marker='s', color='black', label='Watanabe et al. (2014)')
ax42.errorbar(watanabe_pcuo-watanabe_pcuo_oc, watanabe_rVcmax_oc, 
              xerr=np.sqrt(watanabe_pcuo_std**2+watanabe_pcuo_std_oc**2), yerr=watanabe_rVcmax_sigma_oc, 
              ls='None', marker='s', color='black', label='_')
ax42.errorbar(watanabe_pcuo-watanabe_pcuo_co, watanabe_rVcmax_co, 
              xerr=np.sqrt(watanabe_pcuo_std**2+watanabe_pcuo_std_co**2), yerr=watanabe_rVcmax_sigma_co, 
              ls='None', marker='s', color='black', label='_')
ax42.errorbar(np.take(kinose_pcuo[1], (2,4,7))-np.take(kinose_pcuo[0], (2,4,7)), kinose_rVcmax, 
              xerr=np.sqrt(np.take(kinose_pcuo_std[0], (2,4,7))**2+np.take(kinose_pcuo_std[1], (2,4,7))**2), yerr=kinose_rVcmax_sigma, 
              ls='None', marker='^', color='black', label='Kinose et al. (2019)')
ax42.errorbar(np.take(kinose_pcuo[2], (2,4,7))-np.take(kinose_pcuo[0], (2,4,7)), kinose_rVcmax_s15, 
              xerr=np.sqrt(np.take(kinose_pcuo_std[0], (2,4,7))**2+np.take(kinose_pcuo_std[2], (2,4,7))**2), yerr=kinose_rVcmax_s15_sigma, 
              ls='None', marker='^', color=ax42.lines[-1].get_color())
ax42.errorbar(harmens_pcuo[0::2][1::2]-harmens_pcuo[0::2][0::2], harmens_rVcmax_1,
              xerr=np.sqrt(harmens_pcuo_std[0::2][1::2]**2+harmens_pcuo_std[0::2][0::2]**2), yerr=harmens_rVcmax_1_sigma,
              ls='None', marker='*', color='black')
ax42.errorbar(harmens_pcuo[1::2][1::2]-harmens_pcuo[1::2][0::2], harmens_rVcmax_2,
              xerr=np.sqrt(harmens_pcuo_std[1::2][1::2]**2+harmens_pcuo_std[1::2][0::2]**2), yerr=harmens_rVcmax_2_sigma,
              ls='None', marker='*', color=ax42.lines[-1].get_color(), label='Harmens et al. (2016)')
ax42.errorbar(gao_pcuo[1::2]-gao_pcuo[0::2], gao_rVcmax,
              xerr=np.sqrt(gao_pcuo_std[0::2]**2+gao_pcuo_std[1::2]**2), yerr=gao_rVcmax_sigma,
              ls='None', marker='v', color='black', label='Gao et al. (2016)')

ax42.errorbar(pelle14_pcuo[1::3], pelle14_rVcmax, 
              xerr=pelle14_pcuo_std[1::3], yerr=pelle14_rVcmax_sigma, 
              ls='None', marker='d', fillstyle='none', label='Pellegrini (2014)')
ax42.errorbar(watanabe13_pcuo[1]-watanabe13_pcuo[0], watanabe13_rVcmax_beech, 
              xerr=((np.sqrt(watanabe13_pcuo_std[0]**2+watanabe13_pcuo_std[1][0]**2),),(np.sqrt(watanabe13_pcuo_std[0]**2+watanabe13_pcuo_std[1][1]**2),)), yerr=watanabe13_rVcmax_beech_sigma, 
              ls='None', marker='s', fillstyle='none', label='Watanabe et al. (2013)')
ax42.errorbar(watanabe13_pcuo[3]-watanabe13_pcuo[2], watanabe13_rVcmax_oak, 
              xerr=((np.sqrt(watanabe13_pcuo_std[2]**2+watanabe13_pcuo_std[3][0]**2),),(np.sqrt(watanabe13_pcuo_std[2]**2+watanabe13_pcuo_std[3][1]**2),)), yerr=watanabe13_rVcmax_oak_sigma, 
              ls='None', marker='s', fillstyle='none', color=ax42.lines[-1].get_color(), label='_')

ax41.plot(np.arange(0,100), yfit_Jmax, ls='--', color='black', label="1 dof;  no uncert.")
ax41.plot(np.arange(0,100), yfit_Jmax2, ls=':', color='red', label="1 dof;  y uncert.")
ax41.plot(np.arange(0,100), yfit_Jmax_2, ls='-.', color='blue', label="1 dof;  x-y uncert.")

ax42.plot(np.arange(0,100), yfit_Vcmax, ls='--', color='black')#, label="lin. fit: no uncert.")
ax42.plot(np.arange(0,100), yfit_Vcmax2, ls=':', color='red')#, label="lin. fit: y uncert.")
ax42.plot(np.arange(0,100), yfit_Vcmax_2, ls='-.', color='blue')#, label="lin. fit: x-y uncert.")


for ax in fig4.axes:
    if not b_con:
        ax.set_ylim(0.,1.4)
    else:
        ax.set_ylim(0.,2.4)
    ax.set_xlim(0.,100)

if not b_con:
    ax41.legend(bbox_to_anchor=(2.15, 0.55), loc='upper left', borderaxespad=0.)
    ax42.legend(bbox_to_anchor=(1, 1), loc='upper left', borderaxespad=0.)
    ax41.set_xlabel("CUO (mmol $m^{-2}$)", x=1)

else:
    ax42.legend(bbox_to_anchor=(1, 1), loc='upper left', borderaxespad=0.)
    
ax41.set_ylabel("$J_{max}^{O_3}/J_{max}^{CF}$")
ax42.set_ylabel("$V_{cmax}^{O_3}/V_{cmax}^{CF}$")

# Figures for paper
if  not b_con:
    fig5, (ax51, ax52, ax53) = plt.subplots(1,3, gridspec_kw={'width_ratios': [6, 6, 1]}, figsize=(18,6))
    fig5.canvas.set_window_title("ozone_response_gsto_anet_fits")
    ax51.set_title("(a)", x=0.1, y=0.9, size='xx-large')
    ax52.set_title("(b)", x=0.1, y=0.9, size='xx-large')

ax51.set_title("(c)", x=0.1, y=0.9, size='xx-large')
ax52.set_title("(d)", x=0.1, y=0.9, size='xx-large')

ax53.remove()

# Gsto
ax51.errorbar(xu_pcuo, xu_rgs,
              xerr=xu_pcuo_std, yerr=xu_rgs_sigma,
              ls='None', marker='o', color='black')
ax51.errorbar(pelle_pcuo[1::3], pelle_rgs,
              xerr=pelle_pcuo_std[1::3], yerr=pelle_rgs_sigma,
              ls='None', marker='d', color='black')
ax51.errorbar(watanabe_pcuo-watanabe_pcuo_cf, watanabe_rgs,
              xerr=np.sqrt(watanabe_pcuo_std**2+watanabe_pcuo_std_cf**2), yerr=watanabe_rgs_sigma,
              ls='None', marker='s', color='black')
ax51.errorbar(watanabe_pcuo-watanabe_pcuo_oc, watanabe_rgs_oc,
              xerr=np.sqrt(watanabe_pcuo_std**2+watanabe_pcuo_std_oc**2), yerr=watanabe_rgs_sigma_oc,
              ls='None', marker='s', color=ax51.lines[-1].get_color())
ax51.errorbar(watanabe_pcuo-watanabe_pcuo_co, watanabe_rgs_co,
              xerr=np.sqrt(watanabe_pcuo_std**2+watanabe_pcuo_std_co**2), yerr=watanabe_rgs_sigma_co,
              ls='None', marker='s', color=ax51.lines[-1].get_color())
ax51.errorbar(np.take(kinose_pcuo[1], (2,4,7))-np.take(kinose_pcuo[0], (2,4,7)), kinose_rgs,
              xerr=np.sqrt(np.take(kinose_pcuo_std[0], (2,4,7))**2+np.take(kinose_pcuo_std[1], (2,4,7))**2), yerr=kinose_rgs_sigma,
              ls='None', marker='^', color='black')
ax51.errorbar(harmens_pcuo[0::2][1::2]-harmens_pcuo[0::2][0::2], harmens_rgs_1,
              xerr=np.sqrt(harmens_pcuo_std[0::2][1::2]**2+harmens_pcuo_std[0::2][0::2]**2), yerr=harmens_rgs_1_sigma,
              ls='None', marker='*', color='black')
ax51.errorbar(harmens_pcuo[1::2][1::2]-harmens_pcuo[1::2][0::2], harmens_rgs_2,
              xerr=np.sqrt(harmens_pcuo_std[1::2][1::2]**2+harmens_pcuo_std[1::2][0::2]**2), yerr=harmens_rgs_2_sigma,
              ls='None', marker='*', color=ax51.lines[-1].get_color())
ax51.errorbar(gao_pcuo[1::2]-gao_pcuo[0::2], gao_rgs,
              xerr=np.sqrt(gao_pcuo_std[0::2]**2+gao_pcuo_std[1::2]**2), yerr=gao_rgs_sigma,
              ls='None', marker='v', color='black')

ax51.errorbar(pelle14_pcuo[1::3], pelle14_rgs,
              xerr=pelle14_pcuo_std[1::3], yerr=pelle14_rgs_sigma,
              ls='None', marker='d', fillstyle='none')
ax51.errorbar(np.take(kinose_pcuo[2], (2,4,7))-np.take(kinose_pcuo[0], (2,4,7)), kinose_rgs_s15,
              xerr=np.sqrt(np.take(kinose_pcuo_std[0], (2,4,7))**2+np.take(kinose_pcuo_std[2], (2,4,7))**2), yerr=kinose_rgs_s15_sigma,
              ls='None', marker='^', fillstyle='none', color=ax51.lines[-1].get_color())

ax51.errorbar(watanabe13_pcuo[1]-watanabe13_pcuo[0], watanabe13_rgs_beech,
              xerr=((np.sqrt(watanabe13_pcuo_std[0]**2+watanabe13_pcuo_std[1][0]**2),),(np.sqrt(watanabe13_pcuo_std[0]**2+watanabe13_pcuo_std[1][1]**2),)), yerr=watanabe13_rgs_beech_sigma,
              ls='None', marker='s', fillstyle='none')
ax51.errorbar(watanabe13_pcuo[3]-watanabe13_pcuo[2], watanabe13_rgs_oak,
              xerr=((np.sqrt(watanabe13_pcuo_std[2]**2+watanabe13_pcuo_std[3][0]**2),),(np.sqrt(watanabe13_pcuo_std[2]**2+watanabe13_pcuo_std[3][1]**2),)), yerr=watanabe13_rgs_oak_sigma,
              ls='None', marker='s', fillstyle='none', color=ax51.lines[-1].get_color())


# Anet
ax52.errorbar(xu_pcuo, xu_rA,
              xerr=xu_pcuo_std, yerr=xu_rA_sigma,
              ls='None', marker='o', color='black', label='Xu et al. (2019)')
ax52.errorbar(pelle_pcuo[1::3], pelle_rA,
              xerr=pelle_pcuo_std[1::3], yerr=pelle_rA_sigma,
              ls='None', marker='d', color='black', label='Pellegrini et al. (2011)')
ax52.errorbar(watanabe_pcuo-watanabe_pcuo_cf, watanabe_rA,
              xerr=np.sqrt(watanabe_pcuo_std**2+watanabe_pcuo_std_cf**2), yerr=watanabe_rA_sigma,
              ls='None', marker='s', color='black', label='Watanabe et al. (2014)')
ax52.errorbar(watanabe_pcuo-watanabe_pcuo_oc, watanabe_rA_oc,
              xerr=np.sqrt(watanabe_pcuo_std**2+watanabe_pcuo_std_oc**2), yerr=watanabe_rA_sigma_oc,
              ls='None', marker='s', color=ax52.lines[-1].get_color())
ax52.errorbar(watanabe_pcuo-watanabe_pcuo_co, watanabe_rA_co,
              xerr=np.sqrt(watanabe_pcuo_std**2+watanabe_pcuo_std_co**2), yerr=watanabe_rA_sigma_co,
              ls='None', marker='s', color=ax52.lines[-1].get_color())
ax52.errorbar(np.take(kinose_pcuo[1], (2,4,7))-np.take(kinose_pcuo[0], (2,4,7)), kinose_rA,
              xerr=np.sqrt(np.take(kinose_pcuo_std[0], (2,4,7))**2+np.take(kinose_pcuo_std[1], (2,4,7))**2), yerr=kinose_rA_sigma,
              ls='None', marker='^', color='black', label='Kinose et al. (2019)')
ax52.errorbar(np.take(kinose_pcuo[2], (2,4,7))-np.take(kinose_pcuo[0], (2,4,7)), kinose_rA_s15,
              xerr=np.sqrt(np.take(kinose_pcuo_std[0], (2,4,7))**2+np.take(kinose_pcuo_std[2], (2,4,7))**2), yerr=kinose_rA_s15_sigma,
              ls='None', marker='^', color=ax52.lines[-1].get_color())

ax52.errorbar(harmens_pcuo[0::2][1::2]-harmens_pcuo[0::2][0::2], harmens_rA_1,
              xerr=np.sqrt(harmens_pcuo_std[0::2][1::2]**2+harmens_pcuo_std[0::2][0::2]**2), yerr=harmens_rA_1_sigma,
              ls='None', marker='*', color='black')
ax52.errorbar(harmens_pcuo[1::2][1::2]-harmens_pcuo[1::2][0::2], harmens_rA_2,
              xerr=np.sqrt(harmens_pcuo_std[1::2][1::2]**2+harmens_pcuo_std[1::2][0::2]**2), yerr=harmens_rA_2_sigma,
              ls='None', marker='*', color=ax52.lines[-1].get_color(), label='Harmens et al. (2016)')
ax52.errorbar(gao_pcuo[1::2]-gao_pcuo[0::2], gao_rA,
              xerr=np.sqrt(gao_pcuo_std[0::2]**2+gao_pcuo_std[1::2]**2), yerr=gao_rA_sigma,
              ls='None', marker='v', color='black', label='Gao et al. (2016)')

ax52.errorbar(pelle14_pcuo[1::3], pelle14_rA,
              xerr=pelle14_pcuo_std[1::3], yerr=pelle14_rA_sigma,
              ls='None', marker='d', fillstyle='none', label='Pellegrini (2014)')

ax52.errorbar(watanabe13_pcuo[1]-watanabe13_pcuo[0], watanabe13_rA_beech,
              xerr=((np.sqrt(watanabe13_pcuo_std[0]**2+watanabe13_pcuo_std[1][0]**2),),(np.sqrt(watanabe13_pcuo_std[0]**2+watanabe13_pcuo_std[1][1]**2),)), yerr=watanabe13_rA_beech_sigma,
              ls='None', marker='s', fillstyle='none', label='Watanabe et al. (2013)')
ax52.errorbar(watanabe13_pcuo[3]-watanabe13_pcuo[2], watanabe13_rA_oak,
              xerr=((np.sqrt(watanabe13_pcuo_std[2]**2+watanabe13_pcuo_std[3][0]**2),),(np.sqrt(watanabe13_pcuo_std[2]**2+watanabe13_pcuo_std[3][1]**2),)), yerr=watanabe13_rA_oak_sigma,
              ls='None', marker='s', fillstyle='none', color=ax52.lines[-1].get_color())

# Fits
ax51.plot(np.arange(0,100), yfit_gs, ls='--', color='black', label="1 dof; no uncert.")
ax51.plot(np.arange(0,100), yfit_gs2, ls=':', color='red', label="1 dof; y uncert.")
ax51.plot(np.arange(0,100), yfit_gs_2, ls='-.', color='blue', label="1 dof; x-y uncert.")
ax51.plot(np.arange(0,100), yfit_gs_free, ls='--', color='orange', label="2 dof; no uncert.")
#ax51.plot(np.arange(0,100), yfit_gs2_free, ls=':', color='orange', label="2 dof; y uncert.")
#ax51.plot(np.arange(0,100), yfit_gs_2_free, ls='-.', color='orange', label="2 dof; x-y uncert.")


ax52.plot(np.arange(0,100), yfit_A, ls='--', color='black')#, label="lin. fit: no uncert.")
ax52.plot(np.arange(0,100), yfit_A2, ls=':', color='red') #, label="lin. fit: y uncert.")
ax52.plot(np.arange(0,100), yfit_A_2, ls='-.', color='blue') #, label="lin. fit: x-y uncert.")
ax52.plot(np.arange(0,100), yfit_A_free, ls='--', color='orange') #, label="lin. fit: 2 dof; non")
#ax52.plot(np.arange(0,100), yfit_A2_free, ls=':', color='orange') #, label="lin. fit: 2 dof; y")
#ax52.plot(np.arange(0,100), yfit_A_2_free, ls='-.', color='orange') #, label="lin. fit: 2 dof; x-y")


for ax in fig5.axes:
    ax.set_ylim(0.,2.4)
    ax.set_xlim(0.,100)

if not b_con:
    ax51.legend(bbox_to_anchor=(2.15, 0.55), loc='upper left', borderaxespad=0.)
    ax52.legend(bbox_to_anchor=(1, 1), loc='upper left', borderaxespad=0.)
else:
    ax51.legend(bbox_to_anchor=(2.15, 1), loc='upper left', borderaxespad=0.)

ax51.set_xlabel("CUO (mmol $m^{-2}$)", x=1)
ax51.set_ylabel("$g_s^{O_3}/g_s^{CF}$")
ax52.set_ylabel("$A_n^{O_3}/A_n^{CF}$")


#Show it
plt.show(block=False)

# Save
save_data = pd.DataFrame({'Vcmax_ratio':flunder(Vcmax), 'Vcmax_ratio_std': flunder(Vcmax_std), 'Jmax_ratio':flunder(Jmax), 'Jmax_ratio_std':flunder(Jmax_std)})
save_data.to_csv("Vcmax_Jmax_ratios_articles.csv")
