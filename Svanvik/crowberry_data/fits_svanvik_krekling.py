def f_light(M, alpha):
    '''
    Parameters:
    -----------
    M : array
    Irradiance [PPFD] = mumol_O3 m^-2 s^-1
    alpha: float
    Light activation, [alpha] = m^2 s mumol_O3^-1 
    
    Returns:
    --------
    relative stomatal conductance g_sto^rel e {0,1}
    '''
    gs = 1-np.exp(-alpha*M)

    return(gs)


stomatal_conductance_o3 = np.concatenate((((1e3*k_O3*data_krekling['June Cond'].where((np.log(data_krekling['PARo'])>3)).dropna())).values, ((1e3*k_O3*data_krekling['Aug Cond'].where((np.log(data_krekling['PARo.1'])>3)).dropna())).values, ((1e3*k_O3*data_krekling['Sept Cond'].where((np.log(data_krekling['PARo.2'])>3)).dropna())).values))

ppfd = np.concatenate((((data_krekling['PARo'].where((np.log(data_krekling['PARo'])>3)&(data_krekling['June Cond'].notna())).dropna())).values, ((data_krekling['PARo.1'].where((np.log(data_krekling['PARo.1'])>3)&(data_krekling['Aug Cond'].notna())).dropna())).values, ((data_krekling['PARo.2'].where((np.log(data_krekling['PARo.2'])>3)&(data_krekling['Sept Cond'].notna())).dropna())).values))
 

gmax = np.max(stomatal_conductance_o3)

print("gmax = %3.2f mu mol O_3 m^-2 s^-1 kg^-1" % gmax )

from scipy.optimize import curve_fit

fit_params, cov_mat = curve_fit(f_light, ppfd, stomatal_conductance_o3/gmax, p0=[1e-3,])
print("Fit results")
print("alpha", fit_params)
print("cov mat", cov_mat)
