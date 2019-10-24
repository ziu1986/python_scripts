import numpy as np
import matplotlib.pyplot as plt # Plotting
from mytools.met_tools import print_all

def LAI_h(height,**kwarg):
    '''
    LAI dependent on height above ground.
    '''
    
    plant_height = kwarg.pop('plant_height', 10)
    max_lai = kwarg.pop('max_lai', 10)
    species = kwarg.pop('species', 'CF')
    #print("Species: %s" % species)
    if (species == 'CF'):
        if(height<0.3*plant_height):
            LAI = 0
        elif height<0.5*plant_height:
            LAI = max_lai
        elif height<0.7*plant_height:
            LAI = 0.7*max_lai
        elif height<0.9*plant_height:
            LAI = 0.5*max_lai
        elif height<=plant_height:
            LAI = 0.2*max_lai
        else:
            LAI = 0
    elif (species == 'CF2'):
        s = max_lai/np.sqrt(0.5*plant_height)
        if height<0.3*plant_height:
            LAI = 0
        elif height <= plant_height:
            LAI = s*np.sqrt(-(height-plant_height))
        else:
            LAI = 0
    elif (species == 'GR'):
        if(height <= plant_height):
            LAI = max_lai
        else:
            LAI = 0
    elif (species == 'SH'):
        s = max_lai/np.sqrt(0.5*plant_height)
        if height <= plant_height:
            LAI = s*np.sqrt(-(height-plant_height))
        else:
            LAI = 0
    return LAI

def ozone_t(ozone_conc0, time, Rate):
    '''
    Deposition of ozone.
    '''
    return(ozone_conc0*np.exp(-Rate*time))

def diffuse(ozone_vertical_profile, **kwarg):
    lower_boundary = kwarg.pop('lower_boundary',0)
    upper_boundary = kwarg.pop('upper_boundary',0)
    ozone_vertical_profile_tmp = np.insert(ozone_vertical_profile,0,lower_boundary)
    ozone_vertical_profile_tmp = np.insert(ozone_vertical_profile_tmp,0,upper_boundary)
    ozone_tot_ns = ((ozone_vertical_profile_tmp+np.roll(ozone_vertical_profile_tmp,1)+np.roll(ozone_vertical_profile_tmp,-1))/3.)[1:-1]
    return(ozone_tot_ns)


plt.close('all')
height = np.linspace(0,50, num=100)
lai_low_tree = []
lai_high_tree = []
lai_mid_tree = []
lai_grass = []
lai_shrub = []

for each in height:
    lai_low_tree.append(LAI_h(each, species='CF2'))
    lai_high_tree.append(LAI_h(each, plant_height=40, species='CF2'))
    lai_mid_tree.append(LAI_h(each, plant_height=20, species='CF2'))
    lai_grass.append(LAI_h(each, plant_height=1, species='GR'))
    lai_shrub.append(LAI_h(each, plant_height=5, species='SH'))
    
lai_low_tree = np.array(lai_low_tree)
lai_high_tree = np.array(lai_high_tree)
lai_mid_tree = np.array(lai_mid_tree)
lai_grass = np.array(lai_grass)
lai_shrub = np.array(lai_shrub)

ozone_low_tree = []
ozone_high_tree = []
ozone_mid_tree = []
ozone_grass = []
ozone_shrub = []
ozone_tot = []

ozone_diffu = np.repeat(40.,height.size)
ozone_tot_diffu = []

k_leaf = 0.01
k_ground = np.zeros_like(lai_low_tree).astype(float)
k_ground[0] = 0.2


K1 = k_leaf*lai_low_tree
K2 = k_leaf*lai_high_tree
K3 = k_leaf*lai_mid_tree
K4 = k_leaf*lai_grass
K5 = k_leaf*lai_shrub
Ktot = 2*K1 + 2*K2 + 2*K3 + K4 + 2*K5 + k_ground

for time in np.arange(0,25):
    ozone_low_tree.append(ozone_t(40, time, 2*K1))
    ozone_high_tree.append(ozone_t(40, time, 2*K2))
    ozone_mid_tree.append(ozone_t(40, time, 2*K3))
    ozone_grass.append(ozone_t(40, time, K4))
    ozone_shrub.append(ozone_t(40, time, 2*K5))
    ozone_tot.append(ozone_t(40, time, Ktot))
    
    ozone_tot_diffu.append(diffuse(ozone_tot[-1]))
    try:
        ozone_diffu = np.append(ozone_diffu, diffuse(ozone_diffu[-1]), axis=0)
    except IndexError:
        ozone_diffu = np.append(ozone_diffu, diffuse(ozone_diffu), axis=0)

# Plot it
fig1 = plt.figure(1, figsize=(16,9))
fig1.canvas.set_window_title("ozone_uptake_static")

ax11 = plt.subplot(171)
ax11.set_title("Vertical Tree Profile")
plt.plot(np.array(lai_low_tree), height, label='low tree')
plt.plot(np.array(lai_mid_tree), height, label='mid tree', color='black')
plt.plot(np.array(lai_high_tree), height, label='high tree', color='red')
ax11.set_ylim(0,50)
ax11.set_xlabel("arbitrary LAI(h)")
ax11.set_ylabel("Height (m)")
ax11.legend()

ax12 = plt.subplot(172)
ax12.set_title("low tree")
plt.pcolormesh(np.array(ozone_low_tree).transpose())

ax13 = plt.subplot(173)
ax13.set_title("mid tree")
plt.pcolormesh(np.array(ozone_mid_tree).transpose())

ax14 = plt.subplot(174)
ax14.set_title("high tree")
plt.pcolormesh(np.array(ozone_high_tree).transpose())

ax15 = plt.subplot(175)
ax15.set_title("Grass")
plt.pcolormesh(np.array(ozone_grass).transpose())

ax16 = plt.subplot(176)
ax16.set_title("Shrub")
plt.pcolormesh(np.array(ozone_shrub).transpose())

ax14.set_xlabel("Time (atu)")

ax17 = plt.subplot(177)
ax17.set_title("all trees")
pcm4 = plt.pcolormesh(np.array(ozone_tot).transpose())
plt.colorbar(pcm4)
fig1.axes[-1].set_ylabel("[$O_3$] (ppt)")

ax11.text(0,42.5,'k_ground = %1.2f atu$^{-1}$' % k_ground[0], color='black')
ax11.text(0,43.5,'k_leaf = %1.2f atu$^{-1}$' % k_leaf, color='black')

for ax in fig1.axes[1:-1]:
    ax.set_yticklabels('')

fig2 = plt.figure(2, figsize=(16,9))
fig2.canvas.set_window_title("ozone_diffusion")

ax21 = plt.subplot()
ax21.pcolormesh(ozone_diffu)
# Show it
plt.show(block=False)
