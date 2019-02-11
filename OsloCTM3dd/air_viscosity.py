import os, glob
import numpy as np
import matplotlib.pyplot as plt
from mytools.met_tools import *

def sutherland(T):
    '''
    Sutherland, W. (1893), "The viscosity of gases and molecular force", Philosophical Magazine, S. 5, 36, pp. 507-531 (1893).
    '''
    C1 = 1.458e-6
    S = 110.4
    return(C1*T**(3/2.)/(T+S))
    #T0 = 273.15
    #mu0 = 1.716e-5
    #return(mu0*(T/T0)**(3/2.)*(T0+S)/(T+S))

# Viscosity computation
temperature = np.array((-75, -50, -25, -15, -10, -5, 0, 5, 10, 15, 20, 25, 30, 40, 50, 60, 80, 100, 125, 150, 175, 200, 225, 300, 412, 500, 600, 700, 800, 900, 1000, 1100))+273.15

dynamic_viscosity = np.array((13.18, 14.56, 15.88, 16.40, 16.65, 16.90, 17.15, 17.40, 17.64, 17.89, 18.13, 18.37, 18.60, 19.07, 19.53, 19.99, 20.88, 21.74, 22.79, 23.80, 24.78, 25.73, 26.66, 29.28, 32.87, 35.47, 38.25, 40.85, 43.32, 45.66, 47.88, 50.01))

fit = np.polyfit(temperature[2:13], dynamic_viscosity[2:13],1)
fit_function = np.poly1d(fit)

fit2 = np.linalg.lstsq(temperature[2:13,np.newaxis],dynamic_viscosity[2:13])

eval_range = np.linspace(-30,40)+273.15

# Results using Prather fit or Sutherland's Law
src = os.environ['DATA'] + "/astra_data/ctm_results/viscosity_tests/"
srcfile = ('z0w_linear_calm.txt', 'z0w_linear_rough.txt','z0w_sutherland_calm.txt', 'z0w_sutherland_rough.txt',
           'Rb_linear_calm.txt', 'Rb_linear_rough.txt', 'Rb_sutherland_calm.txt', 'Rb_sutherland_rough.txt') #
label = ('linear', 'sutherland', 'linear', 'sutherland')
try:
    raw_data
except NameError:
    raw_data = []
    for each in srcfile:
        infile = open(src+each)
        print("Reading %s" % (each))
        indata = infile.readlines()
        data = np.array(indata).astype(np.float)
        raw_data.append(np.array(data))
        infile.close()
    
plt.close("all")
# Plot it
#---------------------------------------------------------------------------------------------------------------
fig1 = plt.figure(1)
fig1.canvas.set_window_title("dynamic_viscosity")
ax11 = plt.subplot(211)

ax11.scatter(temperature, dynamic_viscosity, label='tabulated viscosity')
ax11.plot(eval_range,fit_function(eval_range), color='red', label='linear fit %s' % ("$\\mu(T) = %1.3f \\cdot T + %1.3f$" % (fit[0],fit[1])))
ax11.plot(eval_range, eval_range*6.2e-2, color='black', label='OsloCTM3 %s'% ('$\\mu(T) = %s \\cdot T$' % (0.062)))
ax11.scatter(temperature[2:13], sutherland(temperature[2:13])*1e6, label="Sutherland's law %s" % ("$\\mu(T) = \\frac{C_1\\cdot T^{3/2}}{T+S}$"), marker='x', color='grey')

ax11.set_ylabel("$\\mu$ ($\\times 10^{-6}$ kg $m^{-1}s^{-1}$) ")
ax11.set_xlim(220,330)
ax11.set_ylim(14,20)

ax11.legend()

ax12 = plt.subplot(212)
ax12.plot(eval_range, (eval_range*6.2e-2-sutherland(eval_range)*1e6)/(sutherland(eval_range)*1e6)*100)
ax12.axhline(0,ls='--',color='grey')
ax12.set_xlim(220,330)
ax12.set_ylim(-4,2)
ax12.set_ylabel("$\\Delta\\mu/\\mu_{sut}$ (%)")
ax12.set_xlabel("Temperature (K)")

#----------------------------------------------------------------------------------------------------------------
fig2 = plt.figure(2)
fig2.canvas.set_window_title("dynamic_viscosity_zo_Rb")
hist = []
for i in range(0,4):
    ax = plt.subplot(2,2,i+1)
    ax.set_title("%s" % (label[i]),y=0.9, color='blue')
    
    if i < 2:
        hist.append(ax.hist(raw_data[i*2], bins=1010, range=(-1e-6,80e-5), label='calm'))
        hist.append(ax.hist(raw_data[i*2+1], bins=1010, range=(-1e-6,80e-5), histtype='step', color='black', linewidth=1.5, label='rough'))
        ax.set_xlabel("%s (%s)" % ("$z_0$", "m"))
        ax.legend()
        #ax.axvline(1.5e-5, color='red', ls='--')
        #ax.axvline(2e-3, color='black', ls='--')
    else:
        hist.append(ax.hist(raw_data[i*2], bins=240, range=(-200,40), label="calm"))
        hist.append(ax.hist(raw_data[i*2+1], bins=240, range=(-200,40), histtype='step', color='black', linewidth=1.5, label="rough"))
        ax.set_xlabel("%s (%s)" % ("$R_b$", "$sm^{-1}$"))
        ax.legend()
        
        ratio = np.sum(hist[i*2][0][0:210])/np.sum(hist[i*2][0])*100
        ax.text(-22, 20000, "%2.1f %s" % (ratio, '%'), color='red', fontsize="large")
        ax.text(12, 20000, "%2.1f %s" % (100-ratio, '%'), color='red', fontsize="large")
for ax in fig2.axes[::2]:
    ax.set_ylabel("count")
for ax in fig2.axes[:2]:
    ax.set_xlim(0,10e-5)
for ax in fig2.axes[2:]:
    ax.set_xlim(-200,40)
    ax.set_ylim(0,2500)
    ax.axvline(10, ls='--', color='red')
    #ax.legend()

#hist.append(fig2.axes[1].hist(raw_data[-1], bins=240, range=(-200,40), histtype='step', color='black', linewidth=1.5))
#fig2.axes[1].plot(hist[1][1][:-1]+0.5,hist[1][0]-hist[-1][0])

fig3 = plt.figure(3,figsize=(14,8))
fig3.canvas.set_window_title("dynamic_viscosity_diffs")
ax31 = plt.subplot(121)
ax31.plot(hist[0][1][:-1]+1e-07, hist[2][0]-hist[0][0])
ax31.set_xlabel("%s (%s)" % ("$z_0^{calm}$", "m"))
ax31.set_ylabel("$\\Delta$count")
#ax31.set_ylim(-300,300)
ax32 = plt.subplot(122)
ax32.plot(hist[4][1][:-1]+0.5, hist[6][0]-hist[4][0])
ax32.set_xlabel("%s (%s)" % ("$R_b$", "$sm^{-1}$"))
ax32.set_ylabel("$\\Delta$count")
#ax32.set_xlim(-200,40)
#ax32.set_ylim(-30000,30000)
#ax32.axvline(10, ls='--', color='red')

print("<Rb> (scm^-1)", "<v_b> (cms^-1)")
for i in range(4,8):
    Rb = (np.sum(hist[i][0][:210]*10)+np.sum(hist[i][0][210:]*(hist[i][1][210:-1]+0.5)))/np.sum(hist[i][0])
    print(Rb*100, 1/(Rb*100))
   
   
# Show it
plt.show(block=False)
