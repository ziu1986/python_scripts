import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from sample_from_norm import compute_cuo
from mytools.met_tools import print_all

def ratio(x1, x2, sigma_x1, sigma_x2):
    rg = x1.values/x2.values
    rg_sigma = np.sqrt((1/x2.values*sigma_x1.values)**2 + (x1.values/x2.values**2*sigma_x2.values)**2)
    return(rg, rg_sigma)

def flunder(x, **kwarg):
    '''
    Flatten any kind of list of lists, numpy.arrays, and numbers.
    Returns a flat numpy array.
    '''
    verbose = kwarg.pop('verbose', False)
    result = []
    for elem in x:
        try:
            for num in elem:
                if verbose:
                    print(num)
                result.append(num)
        except TypeError:
            if verbose:
                print(elem)
            result.append(elem)
    return(np.array(result))
try:
    harmens_o3_mu
except NameError:
    execfile("ozone_response_read_data.py")

# Compute accumulaed ozone
execfile("ozone_response_pcuo.py")

# Stomatal contuctance
execfile("ozone_response_rgs.py")

# Maximum electron transport rate Jmax
execfile("ozone_response_rJmax.py")

# Maximum carboxylation rate Vcmax
execfile("ozone_response_rVcmax.py")

# Dark respiration
execfile("ozone_response_rRd.py")

# Chlorophyll A+B content
execfile("ozone_response_rChl.py")

pcuo = np.array((xu_pcuo,
                pelle_pcuo[1::3],
                watanabe_pcuo-watanabe_pcuo_cf,
                watanabe_pcuo-watanabe_pcuo_oc,
                watanabe_pcuo-watanabe_pcuo_co,
                pelle14_pcuo[1::3],
                np.take(kinose_pcuo[1], (2,4,7))-np.take(kinose_pcuo[0], (2,4,7)),
                np.take(kinose_pcuo[2], (2,4,7))-np.take(kinose_pcuo[0], (2,4,7)),
                watanabe13_pcuo[1]-watanabe13_pcuo[0],
                watanabe13_pcuo[3]-watanabe13_pcuo[2],
                gao_pcuo[1::2]-gao_pcuo[0::2],
                harmens_pcuo[0::2][1::2]-harmens_pcuo[0::2][0::2],
                harmens_pcuo[1::2][1::2]-harmens_pcuo[1::2][0::2]))

pcuo_std = np.array((xu_pcuo_std,
                     pelle_pcuo_std[1::3],
                     np.sqrt(watanabe_pcuo_std**2+watanabe_pcuo_std_cf**2),
                     np.sqrt(watanabe_pcuo_std**2+watanabe_pcuo_std_oc**2),
                     np.sqrt(watanabe_pcuo_std**2+watanabe_pcuo_std_co**2),
                     pelle14_pcuo_std[1::3],
                     np.sqrt(np.take(kinose_pcuo_std[0], (2,4,7))**2+np.take(kinose_pcuo_std[1], (2,4,7))**2),
                     np.sqrt(np.take(kinose_pcuo_std[0], (2,4,7))**2+np.take(kinose_pcuo_std[2], (2,4,7))**2),
                     ((np.sqrt(watanabe13_pcuo_std[0]**2+watanabe13_pcuo_std[1][0]**2),),(np.sqrt(watanabe13_pcuo_std[0]**2+watanabe13_pcuo_std[1][1]**2),)),
                     ((np.sqrt(watanabe13_pcuo_std[2]**2+watanabe13_pcuo_std[3][0]**2),),(np.sqrt(watanabe13_pcuo_std[2]**2+watanabe13_pcuo_std[3][1]**2),)),
                     np.sqrt(gao_pcuo_std[0::2]**2+gao_pcuo_std[1::2]**2),
                     np.sqrt(harmens_pcuo_std[0::2][1::2]**2+harmens_pcuo_std[0::2][0::2]**2),
                     np.sqrt(harmens_pcuo_std[1::2][1::2]**2+harmens_pcuo_std[1::2][0::2]**2)
))
# Clean up
plt.close('all')
# Plot data
fig1 = plt.figure(1, figsize=(16,9))
fig1.canvas.set_window_title("ozone_response")
ax11 = plt.subplot(231)
ax11.errorbar(xu_pcuo, xu_rgs,
              xerr=xu_pcuo_std, yerr=xu_rgs_sigma,
              ls='None', marker='o', label='Xu2019')
ax11.errorbar(pelle_pcuo[1::3], pelle_rgs,
              xerr=pelle_pcuo_std[1::3], yerr=pelle_rgs_sigma,
              ls='None', marker='d', label='Pellegrini2011')
ax11.errorbar(watanabe_pcuo-watanabe_pcuo_cf, watanabe_rgs,
              xerr=np.sqrt(watanabe_pcuo_std**2+watanabe_pcuo_std_cf**2), yerr=watanabe_rgs_sigma,
              ls='None', marker='s', label='Watanabe2014')
ax11.errorbar(watanabe_pcuo-watanabe_pcuo_oc, watanabe_rgs_oc,
              xerr=np.sqrt(watanabe_pcuo_std**2+watanabe_pcuo_std_oc**2), yerr=watanabe_rgs_sigma_oc,
              ls='None', marker='s', color=ax11.lines[-1].get_color())
ax11.errorbar(watanabe_pcuo-watanabe_pcuo_co, watanabe_rgs_co,
              xerr=np.sqrt(watanabe_pcuo_std**2+watanabe_pcuo_std_co**2), yerr=watanabe_rgs_sigma_co,
              ls='None', marker='s', color=ax11.lines[-1].get_color())
ax11.errorbar(pelle14_pcuo[1::3], pelle14_rgs,
              xerr=pelle14_pcuo_std[1::3], yerr=pelle14_rgs_sigma,
              ls='None', marker='d', label='Pellegrini2014')
ax11.errorbar(np.take(kinose_pcuo[1], (2,4,7))-np.take(kinose_pcuo[0], (2,4,7)), kinose_rgs,
              xerr=np.sqrt(np.take(kinose_pcuo_std[0], (2,4,7))**2+np.take(kinose_pcuo_std[1], (2,4,7))**2), yerr=kinose_rgs_sigma,
              ls='None', marker='^', label='Kinose2019')
ax11.errorbar(np.take(kinose_pcuo[2], (2,4,7))-np.take(kinose_pcuo[0], (2,4,7)), kinose_rgs_s15,
              xerr=np.sqrt(np.take(kinose_pcuo_std[0], (2,4,7))**2+np.take(kinose_pcuo_std[2], (2,4,7))**2), yerr=kinose_rgs_s15_sigma,
              ls='None', marker='^', color=ax11.lines[-1].get_color())

ax11.errorbar(watanabe13_pcuo[1]-watanabe13_pcuo[0], watanabe13_rgs_beech,
              xerr=((np.sqrt(watanabe13_pcuo_std[0]**2+watanabe13_pcuo_std[1][0]**2),),(np.sqrt(watanabe13_pcuo_std[0]**2+watanabe13_pcuo_std[1][1]**2),)), yerr=watanabe13_rgs_beech_sigma,
              ls='None', marker='s', label='Watanabe2013')
ax11.errorbar(watanabe13_pcuo[3]-watanabe13_pcuo[2], watanabe13_rgs_oak,
              xerr=((np.sqrt(watanabe13_pcuo_std[2]**2+watanabe13_pcuo_std[3][0]**2),),(np.sqrt(watanabe13_pcuo_std[2]**2+watanabe13_pcuo_std[3][1]**2),)), yerr=watanabe13_rgs_oak_sigma,
              ls='None', marker='s', color=ax11.lines[-1].get_color())
ax11.errorbar(gao_pcuo[1::2]-gao_pcuo[0::2], gao_rgs,
              xerr=np.sqrt(gao_pcuo_std[0::2]**2+gao_pcuo_std[1::2]**2), yerr=gao_rgs_sigma,
              ls='None', marker='v', label='Gao2016')

ax11.errorbar(harmens_pcuo[0::2][1::2]-harmens_pcuo[0::2][0::2], harmens_rgs_1,
              xerr=np.sqrt(harmens_pcuo_std[0::2][1::2]**2+harmens_pcuo_std[0::2][0::2]**2), yerr=harmens_rgs_1_sigma,
              ls='None', marker='*', label='Harmens2016')
ax11.errorbar(harmens_pcuo[1::2][1::2]-harmens_pcuo[1::2][0::2], harmens_rgs_2,
              xerr=np.sqrt(harmens_pcuo_std[1::2][1::2]**2+harmens_pcuo_std[1::2][0::2]**2), yerr=harmens_rgs_2_sigma,
              ls='None', marker='*', color=ax11.lines[-1].get_color())


ax12 = plt.subplot(232)
ax12.errorbar(xu_pcuo, xu_rJmax,
              xerr=xu_pcuo_std, yerr=xu_rJmax_sigma,
              ls='None', marker='o')
ax12.errorbar(pelle_pcuo[1::3], pelle_rJmax,
              xerr=pelle_pcuo_std[1::3], yerr=pelle_rJmax_sigma,
              ls='None', marker='d')
ax12.plot(-1,-1, ls='None', marker='s')
ax12.errorbar(pelle14_pcuo[1::3], pelle14_rJmax,
              xerr=pelle14_pcuo_std[1::3], yerr=pelle14_rJmax_sigma,
              ls='None', marker='d')
ax12.errorbar(np.take(kinose_pcuo[1], (2,4,7))-np.take(kinose_pcuo[0], (2,4,7)), kinose_rJmax,
              xerr=np.sqrt(np.take(kinose_pcuo_std[0], (2,4,7))**2+np.take(kinose_pcuo_std[1], (2,4,7))**2), yerr=kinose_rJmax_sigma,
              ls='None', marker='^', label='Kinose2019')
ax12.errorbar(np.take(kinose_pcuo[2], (2,4,7))-np.take(kinose_pcuo[0], (2,4,7)), kinose_rJmax_s15,
              xerr=np.sqrt(np.take(kinose_pcuo_std[0], (2,4,7))**2+np.take(kinose_pcuo_std[2], (2,4,7))**2), yerr=kinose_rJmax_s15_sigma,
              ls='None', marker='^', color=ax12.lines[-1].get_color())
ax12.errorbar(watanabe13_pcuo[1]-watanabe13_pcuo[0], watanabe13_rJmax_beech,
              xerr=((np.sqrt(watanabe13_pcuo_std[0]**2+watanabe13_pcuo_std[1][0]**2),),(np.sqrt(watanabe13_pcuo_std[0]**2+watanabe13_pcuo_std[1][1]**2),)), yerr=watanabe13_rJmax_beech_sigma,
              ls='None', marker='s', label='Watanabe2013')
ax12.errorbar(watanabe13_pcuo[3]-watanabe13_pcuo[2], watanabe13_rJmax_oak,
              xerr=((np.sqrt(watanabe13_pcuo_std[2]**2+watanabe13_pcuo_std[3][0]**2),),(np.sqrt(watanabe13_pcuo_std[2]**2+watanabe13_pcuo_std[3][1]**2),)), yerr=watanabe13_rJmax_oak_sigma,
              ls='None', marker='s', color=ax12.lines[-1].get_color())
ax12.errorbar(gao_pcuo[1::2]-gao_pcuo[0::2], gao_rJmax,
              xerr=np.sqrt(gao_pcuo_std[0::2]**2+gao_pcuo_std[1::2]**2), yerr=gao_rJmax_sigma,
              ls='None', marker='v', label='Gao2016')
ax12.errorbar(harmens_pcuo[0::2][1::2]-harmens_pcuo[0::2][0::2], harmens_rJmax_1,
              xerr=np.sqrt(harmens_pcuo_std[0::2][1::2]**2+harmens_pcuo_std[0::2][0::2]**2), yerr=harmens_rJmax_1_sigma,
              ls='None', marker='*', label='Harmens2016')
ax12.errorbar(harmens_pcuo[1::2][1::2]-harmens_pcuo[1::2][0::2], harmens_rJmax_2,
              xerr=np.sqrt(harmens_pcuo_std[1::2][1::2]**2+harmens_pcuo_std[1::2][0::2]**2), yerr=harmens_rJmax_2_sigma,
              ls='None', marker='*', color=ax11.lines[-1].get_color())

ax13 = plt.subplot(233)
ax13.errorbar(xu_pcuo, xu_rVcmax, 
              xerr=xu_pcuo_std, yerr=xu_rVcmax_sigma, 
              ls='None', marker='o')
ax13.errorbar(pelle_pcuo[1::3], pelle_rVcmax, 
              xerr=pelle_pcuo_std[1::3], yerr=pelle_rVcmax_sigma, 
              ls='None', marker='d')
ax13.errorbar(watanabe_pcuo-watanabe_pcuo_cf, watanabe_rVcmax, 
              xerr=np.sqrt(watanabe_pcuo_std**2+watanabe_pcuo_std_cf**2), yerr=watanabe_rVcmax_sigma, 
              ls='None', marker='s')
ax13.errorbar(watanabe_pcuo-watanabe_pcuo_oc, watanabe_rVcmax, 
              xerr=np.sqrt(watanabe_pcuo_std**2+watanabe_pcuo_std_oc**2), yerr=watanabe_rVcmax_sigma_oc, 
              ls='None', marker='s', color='black')
ax13.errorbar(watanabe_pcuo-watanabe_pcuo_co, watanabe_rVcmax, 
              xerr=np.sqrt(watanabe_pcuo_std**2+watanabe_pcuo_std_co**2), yerr=watanabe_rVcmax_sigma_co, 
              ls='None', marker='s', color='black')
ax13.errorbar(pelle14_pcuo[1::3], pelle14_rVcmax, 
              xerr=pelle14_pcuo_std[1::3], yerr=pelle14_rVcmax_sigma, 
              ls='None', marker='d')
ax13.errorbar(np.take(kinose_pcuo[1], (2,4,7))-np.take(kinose_pcuo[0], (2,4,7)), kinose_rVcmax, 
              xerr=np.sqrt(np.take(kinose_pcuo_std[0], (2,4,7))**2+np.take(kinose_pcuo_std[1], (2,4,7))**2), yerr=kinose_rVcmax_sigma, 
              ls='None', marker='^', label='Kinose2019')
ax13.errorbar(np.take(kinose_pcuo[2], (2,4,7))-np.take(kinose_pcuo[0], (2,4,7)), kinose_rVcmax_s15, 
              xerr=np.sqrt(np.take(kinose_pcuo_std[0], (2,4,7))**2+np.take(kinose_pcuo_std[2], (2,4,7))**2), yerr=kinose_rVcmax_s15_sigma, 
              ls='None', marker='^', color=ax13.lines[-1].get_color())
ax13.errorbar(watanabe13_pcuo[1]-watanabe13_pcuo[0], watanabe13_rVcmax_beech, 
              xerr=((np.sqrt(watanabe13_pcuo_std[0]**2+watanabe13_pcuo_std[1][0]**2),),(np.sqrt(watanabe13_pcuo_std[0]**2+watanabe13_pcuo_std[1][1]**2),)), yerr=watanabe13_rVcmax_beech_sigma, 
              ls='None', marker='s', label='Watanabe2013')
ax13.errorbar(watanabe13_pcuo[3]-watanabe13_pcuo[2], watanabe13_rVcmax_oak, 
              xerr=((np.sqrt(watanabe13_pcuo_std[2]**2+watanabe13_pcuo_std[3][0]**2),),(np.sqrt(watanabe13_pcuo_std[2]**2+watanabe13_pcuo_std[3][1]**2),)), yerr=watanabe13_rVcmax_oak_sigma, 
              ls='None', marker='s', color=ax13.lines[-1].get_color())
ax13.errorbar(gao_pcuo[1::2]-gao_pcuo[0::2], gao_rVcmax,
              xerr=np.sqrt(gao_pcuo_std[0::2]**2+gao_pcuo_std[1::2]**2), yerr=gao_rVcmax_sigma,
              ls='None', marker='v', label='Gao2016')
ax13.errorbar(harmens_pcuo[0::2][1::2]-harmens_pcuo[0::2][0::2], harmens_rVcmax_1,
              xerr=np.sqrt(harmens_pcuo_std[0::2][1::2]**2+harmens_pcuo_std[0::2][0::2]**2), yerr=harmens_rVcmax_1_sigma,
              ls='None', marker='*', label='Harmens2016')
ax13.errorbar(harmens_pcuo[1::2][1::2]-harmens_pcuo[1::2][0::2], harmens_rVcmax_2,
              xerr=np.sqrt(harmens_pcuo_std[1::2][1::2]**2+harmens_pcuo_std[1::2][0::2]**2), yerr=harmens_rVcmax_2_sigma,
              ls='None', marker='*', color=ax11.lines[-1].get_color())

ax14 = plt.subplot(234)
ax14.errorbar(xu_pcuo, xu_rRd, 
              xerr=xu_pcuo_std, yerr=xu_rRd_sigma, 
              ls='None', marker='o')
ax14.errorbar(pelle_pcuo[1::3], pelle_rRd, 
              xerr=pelle_pcuo_std[1::3], yerr=pelle_rRd_sigma, 
              ls='None', marker='d')
ax14.plot(-1,-1, 
              ls='None', marker='s')
ax14.errorbar(pelle14_pcuo[1::3], pelle14_rRd, 
              xerr=pelle14_pcuo_std[1::3], yerr=pelle14_rRd_sigma, 
              ls='None', marker='d')
ax14.errorbar(np.take(kinose_pcuo[1], (2,4,7))-np.take(kinose_pcuo[0], (2,4,7)), kinose_rRd, 
              xerr=np.sqrt(np.take(kinose_pcuo_std[0], (2,4,7))**2+np.take(kinose_pcuo_std[1], (2,4,7))**2), yerr=kinose_rRd_sigma, 
              ls='None', marker='^', label='Kinose2019')
ax14.errorbar(np.take(kinose_pcuo[2], (2,4,7))-np.take(kinose_pcuo[0], (2,4,7)), kinose_rRd_s15, 
              xerr=np.sqrt(np.take(kinose_pcuo_std[0], (2,4,7))**2+np.take(kinose_pcuo_std[2], (2,4,7))**2), yerr=kinose_rRd_s15_sigma, 
              ls='None', marker='^', color=ax14.lines[-1].get_color())
ax14.errorbar(watanabe13_pcuo[1]-watanabe13_pcuo[0], watanabe13_rRd_beech, 
              xerr=((np.sqrt(watanabe13_pcuo_std[0]**2+watanabe13_pcuo_std[1][0]**2),),(np.sqrt(watanabe13_pcuo_std[0]**2+watanabe13_pcuo_std[1][1]**2),)), yerr=watanabe13_rRd_beech_sigma, 
              ls='None', marker='s', label='Watanabe2013')
ax14.errorbar(watanabe13_pcuo[3]-watanabe13_pcuo[2], watanabe13_rRd_oak, 
              xerr=((np.sqrt(watanabe13_pcuo_std[2]**2+watanabe13_pcuo_std[3][0]**2),),(np.sqrt(watanabe13_pcuo_std[2]**2+watanabe13_pcuo_std[3][1]**2),)), yerr=watanabe13_rRd_oak_sigma, 
              ls='None', marker='s', color=ax14.lines[-1].get_color())
ax14.plot(-1,-1,
              ls='None', marker='v', label='Gao2016')
ax14.plot(-1,-1,
              ls='None', marker='*', label='Harmens2017')

ax15 = plt.subplot(235)
ax15.errorbar(xu_pcuo, xu_rChl, 
              xerr=xu_pcuo_std, yerr=xu_rChl_sigma, 
              ls='None', marker='o')
ax15.errorbar([0,]+pelle_pcuo, pelle_rChl, 
              xerr=[0,]+pelle_pcuo_std, yerr=pelle_rChl_sigma, 
              ls='None', marker='d')
ax15.plot(-1,-1, 
              ls='None', marker='s')
ax15.plot(-1,-1, 
              ls='None', marker='d')
ax15.plot(-1,-1, 
              ls='None', marker='^')
ax15.errorbar(watanabe13_pcuo[1]-watanabe13_pcuo[0], watanabe13_rChl_beech, 
              xerr=((np.sqrt(watanabe13_pcuo_std[0]**2+watanabe13_pcuo_std[1][0]**2),),(np.sqrt(watanabe13_pcuo_std[0]**2+watanabe13_pcuo_std[1][1]**2),)), yerr=watanabe13_rChl_beech_sigma, 
              ls='None', marker='s', label='Watanabe2013')
ax15.errorbar(watanabe13_pcuo[3]-watanabe13_pcuo[2], watanabe13_rChl_oak, 
              xerr=((np.sqrt(watanabe13_pcuo_std[2]**2+watanabe13_pcuo_std[3][0]**2),),(np.sqrt(watanabe13_pcuo_std[2]**2+watanabe13_pcuo_std[3][1]**2),)), yerr=watanabe13_rChl_oak_sigma, 
              ls='None', marker='s', color=ax15.lines[-1].get_color())
ax15.errorbar(gao_pcuo[1::2]-gao_pcuo[0::2], gao_rChl,
              xerr=np.sqrt(gao_pcuo_std[0::2]**2+gao_pcuo_std[1::2]**2), yerr=gao_rChl_sigma,
              ls='None', marker='v', label='Gao2016')
ax15.plot(-1,-1,
              ls='None', marker='*', label='Harmens2017')

ax15.set_xlabel("CUO (mmol $m^{-2}$)")
ax11.set_ylabel("$g_s^{O_3}/g_s^{CF}$")
ax12.set_ylabel("$J_{max}^{O_3}/J_{max}^{CF}$")
ax13.set_ylabel("$V_{cmax}^{O_3}/V_{cmax}^{CF}$")
ax14.set_ylabel("$R_{d}^{O_3}/R_{d}^{CF}$")
ax15.set_ylabel("$Chl_{a+b}^{O_3}/Chl_{a+b}^{CF}$")

ax11.legend(bbox_to_anchor=(3.05, -1.), loc='lower right', borderaxespad=0.)

for ax in fig1.axes:
    ax.set_ylim(0.,3)
    ax.set_xlim(0.,100)

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

plt.show(block=False)

