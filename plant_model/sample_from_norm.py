import numpy as np
import matplotlib.pyplot as plt

   
# Clean up
plt.close('all')
# Constants
kO3 = 1.67
kCO2 =1.52
kappa = kO3/kCO2

# Ozone distributions
xu_o3_mu = np.array((32.8,29.1,27.4,25.1))
xu_o3_sigma = np.array((1.2,1.0,0.8,0.5))
xu_o3_label = np.array(("S1", "S2", "S3", "S4"))
xu_o3_days = np.array((26,26,26,54))
m = 0.5/(435.59-302.41)
b = -302.41*m
xu_gs_cf = np.array((426.08,370.98,365.04,366.62))*m+b
xu_gs_cf_sigma = np.array((435.59,378.51,373.36,380.89))*m+b
xu_gs_o3 = np.array((416.56,373.76,356.71,358.30))*m+b
xu_gs_o3_sigma = np.array((427.66,384.86,365.03,370.98))*m+b

# Draw samples

# Draw 10 samples from dist(M10) (mean ozone and std) and sum these
# Repeat 1000 times to get the distribution of accumulated ozone per day
ceo3_day = []

for i in range(4):
    for j in range(1000):
        s = np.random.normal(xu_o3_mu[i],xu_o3_sigma[i],10)
        ceo3_day.append(np.sum(s))
# Reshape
ceo3_day = np.reshape(np.array(ceo3_day),(4,1000))

# Draw samples corresponding to leaf age from the new distribution
# Repeat 1000 times to get the distribution of accumulated ozone per leaf age
ceo3 = []
ceo3_day_mean = []
ceo3_day_sigma = []
for i in range(4):
    # Get new distribution
    ceo3_day_mean.append(ceo3_day[i].mean())
    ceo3_day_sigma.append(ceo3_day[i].std())
    for j in range(1000):
        s = np.random.normal(ceo3_day_mean[-1],ceo3_day_sigma[-1],xu_o3_days[i])
        ceo3.append(np.sum(s))
# Reshape      
ceo3 = np.reshape(np.array(ceo3),(4,1000))

cuo = []
ceo3_mean = []
ceo3_sigma = []
for i in range(4):
    ceo3_mean.append(ceo3_day[i].mean())
    ceo3_sigma.append(ceo3_day[i].std())
    for j in range(10000):
        s = np.random.normal(ceo3_mean[-1],ceo3_sigma[-1],1)
        gs = np.random.normal(xu_gs_o3[i],xu_gs_o3_sigma[i],1)
        tmp = s*kappa*gs*60**2*1e-6
        cuo.append(tmp)
# Reshape
cuo = np.reshape(np.array(cuo),(4,10000))


# Plot it
fig1 = plt.figure(1)

for i in range(4):
    ax = plt.subplot(2,2,i+1)
    ax.hist(ceo3_day[i])
    #ax.set_xlim(240,360)
    ax.set_ylim(0,40)
    ax.set_title(xu_o3_label[i])
    ax.text(0.05, 0.95, "(%3.0f $\pm$%1.0f) $ppb \cdot h$" % (ceo3_day_mean[i], ceo3_day_sigma[i]), transform=ax.transAxes)

fig1.axes[3].set_xlabel("$\Sigma_{10h} O_3 (ppb \cdot h)$ ", x=-0.15)


fig2 = plt.figure(2)

for i in range(4):
    ax = plt.subplot(2,2,i+1)
    ax.hist(ceo3[i])
    #ax.set_xlim(7000,14000)
    ax.set_ylim(0,40)
    ax.set_title(xu_o3_label[i])
    
    ax.text(0.05, 0.95, "(%5.0f $\pm$%2.0f) $ppb \cdot h \cdot d$" % (np.around(ceo3[i].mean(),-2), ceo3[i].std()), transform=ax.transAxes)

fig2.axes[3].set_xlabel("$\Sigma_{10h} O_3 (ppb \cdot h \cdot d)$ ", x=-0.15)

fig3 = plt.figure(3)

for i in range(4):
    ax = plt.subplot(2,2,i+1)
    ax.hist(cuo[i])
    #ax.set_xlim(7000,14000)
    ax.set_ylim(0,400)
    ax.set_title(xu_o3_label[i])
    ax.text(0.05, 0.95, "(%1.2f $\pm$%1.2f) $mmol m^{-2}$" % (np.around(cuo[i].mean(),2), cuo[i].std()), transform=ax.transAxes)

fig3.axes[3].set_xlabel("CUO $(mmol\,m^{-2}$) ", x=-0.15)

# Show it
plt.show(block=False)

    
