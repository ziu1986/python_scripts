import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from mytools.met_tools import print_all

def compute_cuo(o3_mu, o3_sigma, gs_o3, gs_o3_sigma, exp_hours, exp_days, **kwarg):
    '''
    Compute the ozone exposure using MC.
    Parameters:
    Daily ozone exposure mean, std (ppb).
    Hours of exposure per day.
    Total number of days.
    Stomatal conductance (mmol O3 m-2)
    
    '''
    verbose = kwarg.pop('verbose', False)
    if (verbose):
        print("hours: %d days: %d" % (exp_hours, exp_days))
    # Draw samples
    
    # Draw n samples from dist(Mn) (mean ozone and std) and sum these
    # Repeat 1000 times to get the distribution of accumulated ozone per day
    ceo3_day = []
    for j in range(1000):
        s = np.random.normal(o3_mu,o3_sigma,exp_hours)
        ceo3_day.append(np.sum(s))
    # Reshape
    ceo3_day = np.array(ceo3_day)

    # Draw samples corresponding to leaf age from the new distribution
    # Repeat 1000 times to get the distribution of accumulated ozone per leaf age
    ceo3 = []
    # Get new distribution
    ceo3_day_mean = ceo3_day.mean()
    ceo3_day_sigma = ceo3_day.std()
    for j in range(1000):
        s = np.random.normal(ceo3_day_mean,ceo3_day_sigma,exp_days)
        ceo3.append(np.sum(s))
    # Reshape      
    ceo3 = np.array(ceo3)

    cuo = []
    ceo3_mean = ceo3.mean()
    ceo3_sigma = ceo3.std()

    for j in range(10000):
        s = np.random.normal(ceo3_mean,ceo3_sigma,1)
        gs = np.random.normal(gs_o3,gs_o3_sigma,1)
        tmp = s*gs*60**2*1e-9
        cuo.append(tmp)
    # Reshape
    cuo = np.array(cuo)
    if(verbose):
        print("(%1.2f +- %1.2f) nmol mol^-1 h (%1.2f +- %1.2f) nmol mol^-1 h (%1.2f +- %1.2f) mmol m^-2"  % (np.around(ceo3_day_mean,2), ceo3_day_sigma, np.around(ceo3_mean,2), ceo3_sigma, np.around(cuo.mean(),2), cuo.std()))
    
    
    if (verbose):
        # Plot it
        fig1 = plt.figure(10)
        fig1.canvas.set_window_title("sample_from_norm-ceo3_day")
        ax = plt.subplot()
        ax.hist(ceo3_day)
        #ax.set_xlim(240,360)
        ax.set_ylim(0,40)
        ax.text(0.05, 0.95, "(%3.0f $\pm$%1.0f) $ppb \cdot h$" % (ceo3_day_mean, ceo3_day_sigma), transform=ax.transAxes)

        ax.set_xlabel("$\Sigma_{10h} O_3 (ppb \cdot h)$ ", x=-0.15)


        fig2 = plt.figure(12)
        fig2.canvas.set_window_title("sample_from_norm-ceo3")
        ax = plt.subplot()
        ax.hist(ceo3)
        #ax.set_xlim(7000,14000)
        ax.set_ylim(0,40)
           
        ax.text(0.05, 0.95, "(%5.0f $\pm$%2.0f) $ppb \cdot h \cdot d$" % (np.around(ceo3.mean(),-2), ceo3.std()), transform=ax.transAxes)

        ax.set_xlabel("$\Sigma_{10h} O_3 (ppb \cdot h \cdot d)$ ", x=-0.15)

        fig3 = plt.figure(13)
        fig3.canvas.set_window_title("sample_from_norm-cuo")
        ax = plt.subplot()
        ax.hist(cuo)
        #ax.set_xlim(7000,14000)
        ax.set_ylim(0,400)
        ax.text(0.05, 0.95, "(%1.2f $\pm$%1.2f) $mmol m^{-2}$" % (np.around(cuo.mean(),2), cuo.std()), transform=ax.transAxes)

        ax.set_xlabel("CUO $(mmol\,m^{-2}$) ", x=-0.15)

        # Show it
        plt.show(block=False)
        print_all()

    return(cuo.mean(), cuo.std())
# Read file
#data = pd.read_csv("xu_2019.dat",index_col=0, sep=' ')
            
# Clean up
#plt.close('all')

# Ozone distributions
#xu_o3_mu = data['o3_mean'][1::2]
#xu_o3_sigma = data['o3_sigma'][1::2]
#xu_o3_label = data.index[1::2]
#xu_o3_days = data['leaf_age'][1::2]

#xu_gs_cf = data['gs_co_mean'][0::2]
#xu_gs_cf_sigma = data['gs_co_sigma'][0::2]
#xu_gs_o3 = data['gs_co_mean'][1::2]
#xu_gs_o3_sigma = data['gs_co_sigma'][1::2]


#for i in range(4):
    #compute_cuo(xu_o3_mu[i], xu_o3_sigma[i], xu_gs_o3[i], xu_gs_o3_sigma[i], 10, xu_o3_days[i])
