import numpy as np
import matplotlib.pyplot as plt
from mytools.met_tools import *
#from mytools.met_tools import *
# clean up
plt.close('all')

def f1(start_height, rate):
    y = start_height
    f = [start_height,]
    lat = range(61,74)
    for x in lat:
        y = y*(1-rate)
        f.append(y)
    return f

def f2(start_height, rate):
    f = [start_height,]
    lat = range(61,74)
    for x in lat:
        f.append((1-(x-60)*rate)*start_height)
    return f

def f3(start_height, end_height, rate):
    f = [start_height,]
    lat = range(61,74)
    for x in lat:
        f.append((75-x)*start_height/15.+(x-60)*end_height/15.)
    return f

def f4(start_height, rate):
    f = []
    lat = np.arange(-90,90.5,0.5)
    for ilat in lat:
        if(np.fabs(ilat) <= 10):
            f.append(start_height*2)
        elif(np.fabs(ilat) <= 40):
            f.append((1-(np.fabs(ilat)-10)*rate/3)*start_height*2)
        elif(np.fabs(ilat) <= 60):
            f.append(start_height)
        elif(np.fabs(ilat) <= 74):
            f.append((1-(np.fabs(ilat)-60)*rate)*start_height)
        else:
            f.append(start_height*3/10.)
    return f

lat = range(61,74)
# plot it
fig1 = plt.figure(1, figsize=(16,9)) #figsize=(16,9)
fig1.canvas.set_window_title("tree_height")
ax11 = plt.subplot()
#ax11.scatter([60,]+lat+[74,], f1(20,0.05)+[6,], label='exp')
#ax11.scatter([60,]+lat+[74,], f3(20,6,0.05)+[6,], label='Amund')
ax11.plot(np.arange(-90,90.5,0.5), f4(20,0.05), color='green', linewidth=2.75, ls='--',label = 'CF/DF')
ax11.plot(np.arange(-90,90.5,0.5), f4(8,0.05), color='y',ls='-.', linewidth=2.75, label = 'NF')
ax11.plot(np.arange(-90,90.5,0.5), f4(15,0.05), color='darkgreen', linewidth=2.75, ls=':', label = 'BF')
ax11.plot(np.arange(-90,90.5,0.5), f4(2,0.05), color='olive', linewidth=2.75, label = 'MS')
ax11.plot(np.arange(-90,90.5,0.5), (np.array(f4(2,0.05))+np.array(f4(15,0.05))+np.array(f4(8,0.05))+np.array(f4(20,0.05)))/4, color='grey', linewidth=4, label = 'average')
ax11.scatter([60,]+lat+[74,], f2(20,0.05)+[6,], label='EMEP', marker='x')
ax11.set_xlabel("Latitude (deg)")
ax11.set_ylabel("Vegetation Height (m)")
ax11.legend()
# show it
plt.show(block=False)
