import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from sample_from_norm import compute_cuo
from mytools.met_tools import print_all

def ratio(x1, x2, sigma_x1, sigma_x2):
    rg = x1.values/x2.values
    rg_sigma = np.sqrt((1/x2.values*sigma_x1.values)**2 + (x1.values/x2.values**2*sigma_x2.values)**2)
    return(rg, rg_sigma)

def flunder(x, **kwarg):
    '''
    Flatten any kind of list of lists, numpy.arrays, and numbers.
    Returns a flat numpy array.
    '''
    verbose = kwarg.pop('verbose', False)
    result = []
    for elem in x:
        try:
            for num in elem:
                if verbose:
                    print(num)
                result.append(num)
        except TypeError:
            if verbose:
                print(elem)
            result.append(elem)
    return(np.array(result))
try:
    harmens_o3_mu
except NameError:
    execfile("ozone_response_read_data.py")

# Compute accumulaed ozone
execfile("ozone_response_pcuo.py")

# Stomatal contuctance
execfile("ozone_response_rgs.py")

# Maximum electron transport rate Jmax
execfile("ozone_response_rJmax.py")

# Maximum carboxylation rate Vcmax
execfile("ozone_response_rVcmax.py")

# Dark respiration
execfile("ozone_response_rRd.py")

# Chlorophyll A+B content
execfile("ozone_response_rChl.py")

pcuo = np.array((xu_pcuo,
                pelle_pcuo[1::3],
                watanabe_pcuo-watanabe_pcuo_cf,
                watanabe_pcuo-watanabe_pcuo_oc,
                watanabe_pcuo-watanabe_pcuo_co,
                pelle14_pcuo[1::3],
                np.take(kinose_pcuo[1], (2,4,7))-np.take(kinose_pcuo[0], (2,4,7)),
                np.take(kinose_pcuo[2], (2,4,7))-np.take(kinose_pcuo[0], (2,4,7)),
                watanabe13_pcuo[1]-watanabe13_pcuo[0],
                watanabe13_pcuo[3]-watanabe13_pcuo[2],
                gao_pcuo[1::2]-gao_pcuo[0::2],
                harmens_pcuo[0::2][1::2]-harmens_pcuo[0::2][0::2],
                harmens_pcuo[1::2][1::2]-harmens_pcuo[1::2][0::2]))

pcuo_std = np.array((xu_pcuo_std,
                     pelle_pcuo_std[1::3],
                     np.sqrt(watanabe_pcuo_std**2+watanabe_pcuo_std_cf**2),
                     np.sqrt(watanabe_pcuo_std**2+watanabe_pcuo_std_oc**2),
                     np.sqrt(watanabe_pcuo_std**2+watanabe_pcuo_std_co**2),
                     pelle14_pcuo_std[1::3],
                     np.sqrt(np.take(kinose_pcuo_std[0], (2,4,7))**2+np.take(kinose_pcuo_std[1], (2,4,7))**2),
                     np.sqrt(np.take(kinose_pcuo_std[0], (2,4,7))**2+np.take(kinose_pcuo_std[2], (2,4,7))**2),
                     ((np.sqrt(watanabe13_pcuo_std[0]**2+watanabe13_pcuo_std[1][0]**2))+(np.sqrt(watanabe13_pcuo_std[0]**2+watanabe13_pcuo_std[1][1]**2)))*0.5,
                     ((np.sqrt(watanabe13_pcuo_std[2]**2+watanabe13_pcuo_std[3][0]**2))+(np.sqrt(watanabe13_pcuo_std[2]**2+watanabe13_pcuo_std[3][1]**2)))*0.5,
                     np.sqrt(gao_pcuo_std[0::2]**2+gao_pcuo_std[1::2]**2),
                     np.sqrt(harmens_pcuo_std[0::2][1::2]**2+harmens_pcuo_std[0::2][0::2]**2),
                     np.sqrt(harmens_pcuo_std[1::2][1::2]**2+harmens_pcuo_std[1::2][0::2]**2)
))

Jmax = np.array((xu_rJmax,
                 pelle_rJmax,
                 np.full_like(watanabe_pcuo-watanabe_pcuo_cf,np.nan),
                 np.full_like(watanabe_pcuo-watanabe_pcuo_oc,np.nan),
                 np.full_like(watanabe_pcuo-watanabe_pcuo_co,np.nan),
                 pelle14_rJmax,
                 kinose_rJmax,
                 kinose_rJmax_s15,
                 watanabe13_rJmax_beech,
                 watanabe13_rJmax_oak,
                 gao_rJmax,
                 harmens_rJmax_1,
                 harmens_rJmax_2))

Jmax_std = np.array((xu_rJmax_sigma,
                 pelle_rJmax_sigma,
                 np.full_like(watanabe_pcuo-watanabe_pcuo_cf,np.nan),
                 np.full_like(watanabe_pcuo-watanabe_pcuo_oc,np.nan),
                 np.full_like(watanabe_pcuo-watanabe_pcuo_co,np.nan),
                 pelle14_rJmax_sigma,
                 kinose_rJmax_sigma,
                 kinose_rJmax_s15_sigma,
                 watanabe13_rJmax_beech_sigma,
                 watanabe13_rJmax_oak_sigma,
                 gao_rJmax_sigma,
                 harmens_rJmax_1_sigma,
                 harmens_rJmax_2_sigma))

Vcmax = np.array((xu_rVcmax,
                 pelle_rVcmax,
                 watanabe_rVcmax,
                 watanabe_rVcmax_oc,
                 watanabe_rVcmax_co,
                 pelle14_rVcmax,
                 kinose_rVcmax,
                 kinose_rVcmax_s15,
                 watanabe13_rVcmax_beech,
                 watanabe13_rVcmax_oak,
                 gao_rVcmax,
                 harmens_rVcmax_1,
                 harmens_rVcmax_2))

Vcmax_std = np.array((xu_rVcmax_sigma,
                      pelle_rVcmax_sigma,
                      watanabe_rVcmax,
                      watanabe_rVcmax_oc,
                      watanabe_rVcmax_co,
                      pelle14_rVcmax_sigma,
                      kinose_rVcmax_sigma,
                      kinose_rVcmax_s15_sigma,
                      watanabe13_rVcmax_beech_sigma,
                      watanabe13_rVcmax_oak_sigma,
                      gao_rVcmax_sigma,
                      harmens_rVcmax_1_sigma,
                      harmens_rVcmax_2_sigma))

Rd = np.array((xu_rRd,
               pelle_rRd,
               np.full_like(watanabe_pcuo,np.nan),
               np.full_like(watanabe_pcuo,np.nan),
               np.full_like(watanabe_pcuo,np.nan),
               pelle14_rRd,
               kinose_rRd,
               kinose_rRd_s15,
               watanabe13_rRd_beech,
               watanabe13_rRd_oak,
               np.full_like(gao_pcuo[1::2],np.nan),
               np.full_like(harmens_pcuo[0::2][1::2],np.nan),
               np.full_like(harmens_pcuo[1::2][1::2],np.nan)))

Rd_std = np.array((xu_rRd_sigma,
                   pelle_rRd_sigma,
                   np.full_like(watanabe_pcuo,np.nan),
                   np.full_like(watanabe_pcuo,np.nan),
                   np.full_like(watanabe_pcuo,np.nan),
                   pelle14_rRd_sigma,
                   kinose_rRd_sigma,
                   kinose_rRd_s15_sigma,
                   watanabe13_rRd_beech_sigma,
                   watanabe13_rRd_oak_sigma,
                   np.full_like(gao_pcuo[1::2],np.nan),
                   np.full_like(harmens_pcuo[0::2][1::2],np.nan),
                   np.full_like(harmens_pcuo[1::2][1::2],np.nan)))
Chl = np.array((xu_rChl,
                pelle_rChl,
                np.full_like(watanabe_pcuo,np.nan),
                np.full_like(watanabe_pcuo,np.nan),
                np.full_like(watanabe_pcuo,np.nan),
                np.full_like(pelle14_pcuo,np.nan),
                np.full_like(np.take(kinose_pcuo[1], (2,4,7)),np.nan),
                np.full_like(np.take(kinose_pcuo[2], (2,4,7)),np.nan),
                watanabe13_rChl_beech,
                watanabe13_rChl_oak,
                gao_rChl,
                np.full_like(harmens_pcuo[0::2][1::2],np.nan),
                np.full_like(harmens_pcuo[1::2][1::2],np.nan)))

Chl_std = np.array((xu_rChl_sigma,
                    pelle_rChl_sigma,
                    np.full_like(watanabe_pcuo,np.nan),
                    np.full_like(watanabe_pcuo,np.nan),
                    np.full_like(watanabe_pcuo,np.nan),
                    np.full_like(pelle14_pcuo,np.nan),
                    np.full_like(np.take(kinose_pcuo[1], (2,4,7)),np.nan),
                    np.full_like(np.take(kinose_pcuo[2], (2,4,7)),np.nan),
                    watanabe13_rChl_beech_sigma,
                    watanabe13_rChl_oak_sigma,
                    gao_rChl_sigma,
                    np.full_like(harmens_pcuo[0::2][1::2],np.nan),
                    np.full_like(harmens_pcuo[1::2][1::2],np.nan)))

execfile("ozone_response_fits.py")

execfile("plot_ozone_response.py")
