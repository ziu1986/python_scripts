import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# Read file
data_xu = pd.read_csv("xu_2019.dat",index_col=0, sep=' ', keep_default_na=False, na_values=['na'])
data_pelle = pd.read_csv("pellegrini_2011.dat",index_col=0, sep=' ', keep_default_na=False, na_values=['na'])
data_watanabe = pd.read_csv("watanabe_2014.dat",index_col=0, sep=' ', keep_default_na=False, na_values=['na'])
data_pelle14 = pd.read_csv("pellegrini_2014.dat",index_col=0, sep=' ', keep_default_na=False, na_values=['na'])
data_kinose = pd.read_csv("kinose_2019.dat",index_col=0, sep=' ', keep_default_na=False, na_values=['na'])

# Ozone distributions
xu_o3_mu = data_xu['o3_mean'][1::2]
xu_o3_sigma = data_xu['o3_sigma'][1::2]
xu_o3_label = data_xu.index[1::2]
xu_o3_days = data_xu['leaf_age'][1::2]

xu_gs_cf = data_xu['gs_mean'][0::2]
xu_gs_cf_sigma = data_xu['gs_sigma'][0::2]
xu_gs_o3 = data_xu['gs_mean'][1::2]
xu_gs_o3_sigma = data_xu['gs_sigma'][1::2]

xu_Jmax_cf = data_xu['Jmax_mean'][0::2]
xu_Jmax_cf_sigma = data_xu['Jmax_sigma'][0::2]
xu_Jmax_o3 = data_xu['Jmax_mean'][1::2]
xu_Jmax_o3_sigma = data_xu['Jmax_sigma'][1::2]

xu_Vcmax_cf = data_xu['Vcmax_mean'][0::2]
xu_Vcmax_cf_sigma = data_xu['Vcmax_sigma'][0::2]
xu_Vcmax_o3 = data_xu['Vcmax_mean'][1::2]
xu_Vcmax_o3_sigma = data_xu['Vcmax_sigma'][1::2]

xu_Rd_cf = data_xu['Rd_mean'][0::2]
xu_Rd_cf_sigma = data_xu['Rd_sigma'][0::2]
xu_Rd_o3 = data_xu['Rd_mean'][1::2]
xu_Rd_o3_sigma = data_xu['Rd_sigma'][1::2]

xu_Chl_cf = data_xu['Chl_a+b_mean'][0::2]
xu_Chl_cf_sigma = data_xu['Chl_a+b_sigma'][0::2]
xu_Chl_o3 = data_xu['Chl_a+b_mean'][1::2]
xu_Chl_o3_sigma = data_xu['Chl_a+b_sigma'][1::2]


pelle_o3_mu = data_pelle['o3_mean'][1::2]
pelle_o3_sigma = data_pelle['o3_sigma'][1::2]
pelle_o3_label = data_pelle.index[1::2]
pelle_o3_days = data_pelle['leaf_age'][1::2]

pelle_gs_cf = data_pelle['gs_mean'][0::2]
pelle_gs_cf_sigma = data_pelle['gs_sigma'][0::2]
pelle_gs_o3 = data_pelle['gs_mean'][1::2]
pelle_gs_o3_sigma = data_pelle['gs_sigma'][1::2]

pelle_Jmax_cf = data_pelle['Jmax_mean'][4::6]
pelle_Jmax_cf_sigma = data_pelle['Jmax_sigma'][4::6]
pelle_Jmax_o3 = data_pelle['Jmax_mean'][5::6]
pelle_Jmax_o3_sigma = data_pelle['Jmax_sigma'][5::6]

pelle_Vcmax_cf = data_pelle['Vcmax_mean'][4::6]
pelle_Vcmax_cf_sigma = data_pelle['Vcmax_sigma'][4::6]
pelle_Vcmax_o3 = data_pelle['Vcmax_mean'][5::6]
pelle_Vcmax_o3_sigma = data_pelle['Vcmax_sigma'][5::6]

pelle_Rd_cf = data_pelle['Rd_mean'][4::6]
pelle_Rd_cf_sigma = data_pelle['Rd_sigma'][4::6]
pelle_Rd_o3 = data_pelle['Rd_mean'][5::6]
pelle_Rd_o3_sigma = data_pelle['Rd_sigma'][5::6]

pelle_Chl_cf = data_pelle['Chl_a+b_mean'][0::2]
pelle_Chl_cf_sigma = data_pelle['Chl_a+b_sigma'][0::2]
pelle_Chl_o3 = data_pelle['Chl_a+b_mean'][1::2]
pelle_Chl_o3_sigma = data_pelle['Chl_a+b_sigma'][1::2]


watanabe_o3_mu = data_watanabe['o3_mean']
watanabe_o3_sigma = data_watanabe['o3_sigma']
watanabe_o3_label = data_watanabe.index
watanabe_o3_days = data_watanabe['leaf_age']
watanabe_o3_fumi = data_watanabe['fumigation']


watanabe_gs_o3 = data_watanabe['gs_mean']
watanabe_gs_o3_sigma = data_watanabe['gs_sigma']

watanabe_Vcmax_o3 = data_watanabe['Vcmax_mean']
watanabe_Vcmax_o3_sigma = data_watanabe['Vcmax_sigma']


pelle14_o3_mu = data_pelle14['o3_mean'][1::2]
pelle14_o3_sigma = data_pelle14['o3_sigma'][1::2]
pelle14_o3_label = data_pelle14.index[1::2]
pelle14_o3_days = data_pelle14['leaf_age'][1::2]

pelle14_gs_cf = data_pelle14['gs_mean'][0::2]
pelle14_gs_cf_sigma = data_pelle14['gs_sigma'][0::2]
pelle14_gs_o3 = data_pelle14['gs_mean'][1::2]
pelle14_gs_o3_sigma = data_pelle14['gs_sigma'][1::2]

pelle14_Jmax_cf = data_pelle14['Jmax_mean'][4::6]
pelle14_Jmax_cf_sigma = data_pelle14['Jmax_sigma'][4::6]
pelle14_Jmax_o3 = data_pelle14['Jmax_mean'][5::6]
pelle14_Jmax_o3_sigma = data_pelle14['Jmax_sigma'][5::6]

pelle14_Vcmax_cf = data_pelle14['Vcmax_mean'][4::6]
pelle14_Vcmax_cf_sigma = data_pelle14['Vcmax_sigma'][4::6]
pelle14_Vcmax_o3 = data_pelle14['Vcmax_mean'][5::6]
pelle14_Vcmax_o3_sigma = data_pelle14['Vcmax_sigma'][5::6]

pelle14_Rd_cf = data_pelle14['Rd_mean'][4::6]
pelle14_Rd_cf_sigma = data_pelle14['Rd_sigma'][4::6]
pelle14_Rd_o3 = data_pelle14['Rd_mean'][5::6]
pelle14_Rd_o3_sigma = data_pelle14['Rd_sigma'][5::6]

kinose_o3_mu = data_kinose['o3_mean']
kinose_o3_sigma = data_kinose['o3_sigma']
kinose_o3_label = data_kinose.index
kinose_o3_days = data_kinose['leaf_age']
kinose_o3_fumi = data_kinose['fumigation']

kinose_gs_o3 = data_kinose['gs_mean']
kinose_gs_o3_sigma = data_kinose['gs_sigma']

kinose_Vcmax_o3 = data_kinose['Vcmax_mean']
kinose_Vcmax_o3_sigma = data_kinose['Vcmax_sigma']
kinose_Jmax_o3 = data_kinose['Jmax_mean']
kinose_Jmax_o3_sigma = data_kinose['Jmax_sigma']
kinose_Rd_o3 = data_kinose['Rd_mean']
kinose_Rd_o3_sigma = data_kinose['Rd_sigma']
