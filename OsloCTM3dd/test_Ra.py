import numpy as np
import matplotlib.pyplot as plt
from mytools.met_tools import *

plt.close("all")

L = (-1e3, -1e2, -1e1, -1, -0.1, 0.1, 1, 1e1, 1e2, 1e3)#np.linspace(1e-3,1e3)
z0 = np.arange(0.01, 2.6, 0.1)
d = 0.7*1
z = 8
# Garratt factor from Kansas experiment
beta = 5
gamma1 = gamma2 = 16

Lm, z0m = np.meshgrid(L, z0)

def compute_Ra(z, d, z0, L, beta, gamma1, gamma2):
    
    def Psi_m(ceta):
        if ceta >= 0:
            Psim = -beta*ceta
        elif ceta > -2:
            x = (1-gamma1*ceta)**0.25
            Psim = np.log(((1+x**2)*0.5)*((1+x)*0.5)**2)-2*np.arctan(x+np.pi*0.5)
        else:
            Psim = np.nan
        return(Psim)
    try:
        Ra = (np.log((z-d)/z0m)-Psi_m((z-d)/Lm)+Psi_m(z0m/Lm))
    except ZeroDivisionError:
        Ra = np.inf
    return(Ra)
print("Ra" ,"L.size", "z0.size")
# Plot it
fig1 = plt.figure(1)
fig1.canvas.set_window_title("test_Ra_wrong_form")
ax11 = plt.subplot(311)
Ra = []
for Lm in L:
    for z0m in z0:
        Ra.append(compute_Ra(z, d, z0m, Lm, beta, gamma1, gamma2))
        #print("%4.1f %0.2f %3.2f" % (Lm, z0m, Ra[-1]))

Ra = np.array(Ra).reshape(len(L),len(z0))
print(Ra.shape, len(L), len(z0))

pcm11 = ax11.pcolormesh(Ra.transpose(), vmin=-3, vmax=3, cmap=plt.cm.coolwarm)

ax11.set_xlabel("")
ax11.set_ylabel("$z_0$")
#ax11.set_xticklabels(L[::2])
#ax11.set_yticklabels(z0[::5])

fig1.colorbar(pcm11, ax=ax11)
fig1.axes[-1].set_ylabel("$R_a \cdot \kappa u_{*}$")

ax11.text(0.1, 24,"$z_{ref} = 8\,m$")
ax11.text(0.1, 22,"$d = 0.7 \cdot 1\,m$")


ax12 = plt.subplot(312)
L = np.arange(0,101,10)
Ra = []
for Lm in L:
    for z0m in z0:
        Ra.append(compute_Ra(z, d, z0m, Lm, beta, gamma1, gamma2))
        #print("%4.1f %0.2f %3.2f" % (Lm, z0m, Ra[-1]))

Ra = np.array(Ra).reshape(len(L),len(z0))
print(Ra.shape, len(L), len(z0))

pcm12 = ax12.pcolormesh(Ra.transpose(), vmin=-3, vmax=3, cmap=plt.cm.coolwarm)

ax12.set_xlabel("")
ax12.set_ylabel("$z_0$")
#ax12.set_xticklabels(L[::2])
#ax12.set_yticklabels(z0[::5])

fig1.colorbar(pcm12, ax=ax12)
fig1.axes[-1].set_ylabel("$R_a \cdot \kappa u_{*}$")

ax13 = plt.subplot(313)
L = np.arange(0,11,1)

Ra = []
for Lm in L:
    for z0m in z0:
        Ra.append(compute_Ra(z, d, z0m, Lm, beta, gamma1, gamma2))
        #print("%4.1f %0.2f %3.2f" % (Lm, z0m, Ra[-1]))

Ra = np.array(Ra).reshape(len(L),len(z0))
print(Ra.shape, len(L), len(z0))

pcm13 = ax13.pcolormesh(Ra.transpose(), vmin=-3, vmax=3, cmap=plt.cm.coolwarm)

ax13.set_xlabel("L")
ax13.set_ylabel("$z_0$")
#ax13.set_xticklabels(L[::2])
#ax13.set_yticklabels(z0[::5])


fig1.colorbar(pcm13, ax=ax13)
fig1.axes[-1].set_ylabel("$R_a \cdot \kappa u_{*}$")

# Show it
plt.show(block=False)

