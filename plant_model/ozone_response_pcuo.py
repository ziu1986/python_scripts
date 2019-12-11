import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from sample_from_norm import compute_cuo
from mytools.met_tools import print_all

def cuo(o3_mu, o3_sigma, gs_o3, gs_o3_sigma, o3_fumi, o3_days, **kwarg):
    article = kwarg.pop("article", "general")
    if (article=="watanabe14"):
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


# Compute accumulated ozone for each article
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

watanabe_pcuo_cf, watanabe_pcuo_std_cf = cuo(watanabe_o3_mu, watanabe_o3_sigma, watanabe_gs_o3, watanabe_gs_o3_sigma, watanabe_o3_fumi, watanabe_o3_days, exp='CC', article='watanabe14')

watanabe_pcuo, watanabe_pcuo_std = cuo(watanabe_o3_mu, watanabe_o3_sigma, watanabe_gs_o3, watanabe_gs_o3_sigma, watanabe_o3_fumi, watanabe_o3_days, exp='OO', article='watanabe14')
watanabe_pcuo_oc, watanabe_pcuo_std_oc = cuo(watanabe_o3_mu, watanabe_o3_sigma, watanabe_gs_o3, watanabe_gs_o3_sigma, watanabe_o3_fumi, watanabe_o3_days, exp='OC', article='watanabe14')
watanabe_pcuo_co, watanabe_pcuo_std_co = cuo(watanabe_o3_mu, watanabe_o3_sigma, watanabe_gs_o3, watanabe_gs_o3_sigma, watanabe_o3_fumi, watanabe_o3_days, exp='CO', article='watanabe14')

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

watanabe13_pcuo = []
watanabe13_pcuo_std = []

for i in (0,2):
    cuo_mean, cuo_std = compute_cuo(watanabe13_o3_mu[i], watanabe13_o3_sigma[i], watanabe13_gs_o3[i], watanabe13_gs_o3_sigma[i], int(watanabe13_o3_fumi[i]), watanabe13_o3_days[i])
    watanabe13_pcuo.append(cuo_mean)
    watanabe13_pcuo_std.append(cuo_std)
    cuo_mean_1, cuo_std_1 = compute_cuo(watanabe13_o3_mu[i], watanabe13_o3_sigma[i], watanabe13_gs_o3[i], watanabe13_gs_o3_sigma[i], int(watanabe13_o3_fumi[i]), watanabe13_o3_days[i]-watanabe13_o3_days[i+1])
    cuo_mean_max, cuo_std_max = compute_cuo(watanabe13_o3_mu[i+1], watanabe13_o3_sigma[i+1], watanabe13_gs_o3[i], watanabe13_gs_o3_sigma[i], int(watanabe13_o3_fumi[i+1]), watanabe13_o3_days[i+1])
    cuo_mean_min, cuo_std_min = compute_cuo(watanabe13_o3_mu[i+1], watanabe13_o3_sigma[i+1], watanabe13_gs_o3[i+1], watanabe13_gs_o3_sigma[i+1], int(watanabe13_o3_fumi[i+1]), watanabe13_o3_days[i+1])
    cuo_mean_mean, cuo_std_mean = compute_cuo(watanabe13_o3_mu[i+1], watanabe13_o3_sigma[i+1], (watanabe13_gs_o3[i]+watanabe13_gs_o3[i+1])*0.5, np.sqrt(watanabe13_gs_o3_sigma[i+1]**2+watanabe13_gs_o3_sigma[i+1]**2)*0.5, int(watanabe13_o3_fumi[i+1]), watanabe13_o3_days[i+1])

    cuo_mean = cuo_mean_1+cuo_mean_mean
    # Max uncertainty estimation (variation of start- and endpoint gsto)
    cuo_std = (cuo_mean_mean-cuo_mean_min, cuo_mean_max-cuo_mean_mean)
    watanabe13_pcuo.append(cuo_mean)
    watanabe13_pcuo_std.append(cuo_std)
    
    #print(watanabe13_o3_days[i], cuo_mean, cuo_std)
    #print(watanabe13_o3_days[i]-watanabe13_o3_days[i+1], cuo_mean_1, cuo_std_1)
    #print(watanabe13_o3_days[i+1], cuo_mean_max, cuo_std_max)
    #print(watanabe13_o3_days[i+1], cuo_mean_min, cuo_std_min)
    #print(watanabe13_o3_days[i+1], cuo_mean_mean, cuo_std_mean)                          

gao_pcuo = []
gao_pcuo_std = []

for i in range(8):
    cuo_mean, cuo_std = compute_cuo(gao_o3_mu[0::2][i], gao_o3_sigma[0::2][i], gao_gs_o3[0::2][i], gao_gs_o3_sigma[0::2][i], int(gao_o3_fumi[0::2][i]), gao_o3_days[0::2][i])

    gao_pcuo.append(cuo_mean)
    gao_pcuo_std.append(cuo_std)
    
    cuo_mean, cuo_std = compute_cuo(gao_o3_mu[1::2][i], gao_o3_sigma[1::2][i], gao_gs_o3[1::2][i], gao_gs_o3_sigma[1::2][i], int(gao_o3_fumi[1::2][i]), gao_o3_days[1::2][i])

    gao_pcuo.append(cuo_mean)
    gao_pcuo_std.append(cuo_std)
   
    cuo_mean, cuo_std = compute_cuo(gao_o3_mu[0::2][i], gao_o3_sigma[0::2][i], gao_gs_o3[1::2][i], gao_gs_o3_sigma[1::2][i], int(gao_o3_fumi[0::2][i]-gao_o3_fumi[1::2][i]), gao_o3_days[1::2][i])

    gao_pcuo[-1] = gao_pcuo[-1] + cuo_mean
    gao_pcuo_std[-1] = gao_pcuo_std[-1] + cuo_std

    #print(cuo_mean, cuo_std)

gao_pcuo = np.array(gao_pcuo)
gao_pcuo_std = np.array(gao_pcuo_std)


harmens_pcuo = []
harmens_pcuo_std = []
for j in range(3):
    for i in range(4):
        
        if j==0:
            cuo_mean, cuo_std = compute_cuo(harmens_o3_mu[4*j:4*j+4][i], harmens_o3_sigma[4*j:4*j+4][i], harmens_gs_o3[4*j:4*j+4][i],harmens_gs_o3_sigma[4*j:4*j+4][i], harmens_o3_fumi[4*j:4*j+4][i].astype(int), harmens_o3_days[4*j:4*j+4][i])
            #print(i, j, cuo_mean, cuo_std)
        if j>0:
            cuo_mean, cuo_std = compute_cuo(harmens_o3_mu[4*j:4*j+4][i], harmens_o3_sigma[4*j:4*j+4][i], (harmens_gs_o3[4*j:4*j+4][i]-harmens_gs_o3[4*(j-1):4*(j-1)+4][i])*0.5+harmens_gs_o3[4*(j-1):4*(j-1)+4][i], 0.5*np.sqrt(harmens_gs_o3_sigma[4*j:4*j+4][i]**2+harmens_gs_o3_sigma[4*(j-1):4*(j-1)+4][i]**2), harmens_o3_fumi[4*j:4*j+4][i].astype(int), harmens_o3_days[4*j:4*j+4][i]-harmens_o3_days[4*(j-1):4*(j-1)+4][i])
            #print(harmens_o3_mu[4*j:4*j+4][i], harmens_o3_sigma[4*j:4*j+4][i], (harmens_gs_o3[4*j:4*j+4][i]-harmens_gs_o3[4*(j-1):4*(j-1)+4][i])*0.5+harmens_gs_o3[4*(j-1):4*(j-1)+4][i], 0.5*np.sqrt(harmens_gs_o3_sigma[4*j:4*j+4][i]**2+harmens_gs_o3_sigma[4*(j-1):4*(j-1)+4][i]**2), harmens_o3_fumi[4*j:4*j+4][i].astype(int), harmens_o3_days[4*j:4*j+4][i]-harmens_o3_days[4*(j-1):4*(j-1)+4][i])
            #print(i, j, cuo_mean, cuo_std)
            cuo_mean = cuo_mean+harmens_pcuo[4*(j-1):4*(j-1)+4][i]
            cuo_std = np.sqrt(cuo_std**2+harmens_pcuo_std[4*(j-1):4*(j-1)+4][i]**2)
            #print(i, j, cuo_mean, cuo_std)
        #print(harmens_pcuo)
        harmens_pcuo.append(cuo_mean)
        harmens_pcuo_std.append(cuo_std)    

harmens_pcuo = np.array(harmens_pcuo).reshape(6,2)
harmens_pcuo_std = np.array(harmens_pcuo_std).reshape(6,2)
    
