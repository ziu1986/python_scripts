import os, glob
import numpy as np
import matplotlib.pyplot as plt
from mytools.met_tools import *

# Results using Prather fit or Sutherland's Law
src = os.environ['DATA'] + "/astra_data/ctm_results/viscosity_tests/Rb0_*"
label = ('linear', 'sutherland')
species = ("O3", "SO2", "NO", "HNO3", "H2O2", "Aceta", "HCHO", "PAN", "NO2", "CH3OOH", "PAA", "HCO2H", "HNO2", "NH3")
try:
    raw_data
except NameError:
    raw_data = []
    for each in sorted(glob.glob(src)):
        infile = open(each)
        print("Reading %s" % (each))
        indata = infile.readlines()
        data = np.array(indata).astype(np.float)
        raw_data.append(np.array(data))
        infile.close()


plt.close("all")
# Plot it
#---------------------------------------------------------------------------------------------------------------
hist_lin = []
hist_suth = []

fig1 = plt.figure(1, figsize=(20,9))
fig1.canvas.set_window_title("dynamic_viscosity_Rbi_1")
for i in range(0,7):
    ax = plt.subplot(2,4,i+1)
    ax.set_title("%s" % (species[i]))
    hist_lin.append(ax.hist(raw_data[i], bins=240, range=(-200,40), histtype='step', linewidth=1.5, color='black', label="linear"))
    hist_suth.append(ax.hist(raw_data[i+14], bins=240, range=(-200,40), histtype='step', linewidth=1.5, label="Sutherland"))
    ax.axvline(10, ls='--', color='red')
    ax.set_ylim(0,20000)
    ratio = np.sum(hist_lin[-1][0][0:210])/np.sum(hist_lin[-1][0])*100
    ax.text(-30, 19000, "%2.1f%s" % (ratio, '%'), color='black', fontsize="large")
    ax.text(14, 19000, "%2.1f%s" % (100-ratio, '%'), color='black', fontsize="large")
    ratio = np.sum(hist_suth[-1][0][0:210])/np.sum(hist_suth[-1][0])*100
    ax.text(-30, 18000, "%2.1f%s" % (ratio, '%'), color='blue', fontsize="large")
    ax.text(14, 18000, "%2.1f%s" % (100-ratio, '%'), color='blue', fontsize="large")
ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
fig1.axes[0].set_ylabel("counts", y=0)
fig1.axes[-2].set_xlabel("%s (%s)" % ("$R_b$", "$sm^{-1}$"), x=1)
fig2 = plt.figure(2, figsize=(20,9))
fig2.canvas.set_window_title("dynamic_viscosity_Rbi_2")
for i in range(0,7):
    ax = plt.subplot(2,4,i+1)
    ax.set_title("%s" % (species[i+7]))
    hist_lin.append(ax.hist(raw_data[i+7], bins=240, range=(-200,40), histtype='step', linewidth=1.5, color='black', label="linear"))
    hist_suth.append(ax.hist(raw_data[i+7+14], bins=240, range=(-200,40), histtype='step', linewidth=1.5, label="Sutherland"))
    ax.axvline(10, ls='--', color='red')
    ax.set_ylim(0,20000)
    ratio = np.sum(hist_lin[-1][0][0:210])/np.sum(hist_lin[-1][0])*100
    ax.text(-30, 19000, "%2.1f%s" % (ratio, '%'), color='black', fontsize="large")
    ax.text(14, 19000, "%2.1f%s" % (100-ratio, '%'), color='black', fontsize="large")
    ratio = np.sum(hist_suth[-1][0][0:210])/np.sum(hist_suth[-1][0])*100
    ax.text(-30, 18000, "%2.1f%s" % (ratio, '%'), color='blue', fontsize="large")
    ax.text(14, 18000, "%2.1f%s" % (100-ratio, '%'), color='blue', fontsize="large")

ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
fig2.axes[0].set_ylabel("counts", y=0)
fig2.axes[-2].set_xlabel("%s (%s)" % ("$R_b$", "$sm^{-1}$"), x=1)

# Show it
plt.show(block=False)

print("<Rb> (scm^-1)", "<v_b> (cms^-1)")
for i in range(0,14):
    Rb = (np.sum(hist_lin[i][0][:210]*10)+np.sum(hist_lin[i][0][210:]*(hist_lin[i][1][210:-1]+0.5)))/np.sum(hist_lin[i][0])
    print("%4.2f, %1.5f" % (Rb*100, 1/(Rb*100)))
    print("---")
    Rb = (np.sum(hist_suth[i][0][:210]*10)+np.sum(hist_suth[i][0][210:]*(hist_suth[i][1][210:-1]+0.5)))/np.sum(hist_suth[i][0])
    print("%4.2f, %1.5f" % (Rb*100, 1/(Rb*100)))
    print("")
