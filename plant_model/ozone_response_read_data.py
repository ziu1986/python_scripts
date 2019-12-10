import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import datetime as dt
from sun_hours import *

# Read file
data_xu = pd.read_csv("xu_2019.dat",index_col=0, sep=' ', keep_default_na=False, na_values=['na'])
data_pelle = pd.read_csv("pellegrini_2011.dat",index_col=0, sep=' ', keep_default_na=False, na_values=['na'])
data_watanabe = pd.read_csv("watanabe_2014.dat",index_col=0, sep=' ', keep_default_na=False, na_values=['na'])
data_pelle14 = pd.read_csv("pellegrini_2014.dat",index_col=0, sep=' ', keep_default_na=False, na_values=['na'])
data_kinose = pd.read_csv("kinose_2019.dat",index_col=0, sep=' ', keep_default_na=False, na_values=['na'])
data_watanabe13 = pd.read_csv("watanabe_2013.dat",index_col=0, sep=' ', keep_default_na=False, na_values=['na'])
data_gao = pd.read_csv("gao_2016.dat",index_col=0, sep=' ', keep_default_na=False, na_values=['na'])
data_harmens = pd.read_csv("harmens_2017.dat",index_col=0, sep=' ', keep_default_na=False, na_values=['na'])

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

watanabe13_o3_mu = data_watanabe13['o3_mean']
watanabe13_o3_sigma = data_watanabe13['o3_sigma']
watanabe13_o3_label = data_watanabe13.index
watanabe13_o3_days = data_watanabe13['leaf_age']
watanabe13_o3_fumi = data_watanabe13['fumigation']

watanabe13_gs_o3 = data_watanabe13['gs_mean']
watanabe13_gs_o3_sigma = data_watanabe13['gs_sigma']

watanabe13_Vcmax_o3 = data_watanabe13['Vcmax_mean']
watanabe13_Vcmax_o3_sigma = data_watanabe13['Vcmax_sigma']
watanabe13_Jmax_o3 = data_watanabe13['Jmax_mean']
watanabe13_Jmax_o3_sigma = data_watanabe13['Jmax_sigma']
watanabe13_Rd_o3 = data_watanabe13['Rd_mean']
watanabe13_Rd_o3_sigma = data_watanabe13['Rd_sigma']
watanabe13_Chl_o3 = data_watanabe13['Chl_a+b_mean']
watanabe13_Chl_o3_sigma = data_watanabe13['Chl_a+b_sigma']


gao_o3_mu = data_gao['o3_mean']
gao_o3_sigma = data_gao['o3_sigma']
gao_o3_label = data_gao.index
gao_o3_days = data_gao['leaf_age']
gao_o3_fumi = data_gao['fumigation']

gao_gs_o3 = data_gao['gs_mean']
gao_gs_o3_sigma = data_gao['gs_sigma']

gao_Vcmax_o3 = data_gao['Vcmax_mean']
gao_Vcmax_o3_sigma = data_gao['Vcmax_sigma']
gao_Jmax_o3 = data_gao['Jmax_mean']
gao_Jmax_o3_sigma = data_gao['Jmax_sigma']
gao_Chl_o3 = data_gao['Chl_a+b_mean']
gao_Chl_o3_sigma = data_gao['Chl_a+b_sigma']


harmens_o3_mu = data_harmens['o3_mean']
harmens_o3_sigma = data_harmens['o3_sigma']
harmens_o3_label = data_harmens.index
harmens_o3_days = data_harmens['leaf_age']
harmens_o3_fumi = data_harmens['fumigation']

harmens_gs_o3 = data_harmens['gs_mean']
harmens_gs_o3_sigma = data_harmens['gs_sigma']

harmens_Vcmax_o3 = data_harmens['Vcmax_mean']
harmens_Vcmax_o3_sigma = data_harmens['Vcmax_sigma']
harmens_Jmax_o3 = data_harmens['Jmax_mean']
harmens_Jmax_o3_sigma = data_harmens['Jmax_sigma']
