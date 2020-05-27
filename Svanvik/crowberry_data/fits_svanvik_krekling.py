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

def f_temp(M, Topt, Tmin, Tmax):
    '''
    Parameters:
    -----------
    M : array
    Temperature [T] = deg C
    Topt, Twidth: float
    [Topt] = deg C 
    
    Returns:
    --------
    relative stomatal conductance g_sto^rel
    '''
    T = M[0]
    fmin = M[1]
    
    bt = float(Tmax-Topt)/(Topt-Tmin)
    gs = (T-Tmin)/float(Topt-Tmin)*((Tmax-T)/float(Tmax-Topt))**bt
    
    gs[np.where((T<Tmin) | (T>Tmax) | (gs<fmin))[0]] = fmin

    return(gs)

# gsto 

photosynth_loglimit = 4

stomatal_conductance_o3 = np.concatenate(((1e3*k_O3*data_krekling['June Cond']).where(np.log(data_krekling['PARo'])>photosynth_loglimit).dropna().values,(1e3*k_O3*data_krekling['Aug Cond']).where(np.log(data_krekling['PARo.1'])>photosynth_loglimit).dropna().values,(1e3*k_O3*data_krekling['Sept Cond']).where(np.log(data_krekling['PARo.2'])>photosynth_loglimit).dropna().values))
x_sample, pdf, fit, stat = fit_skew_normal(stomatal_conductance_o3)

stomatal_conductance_o3_aug = (1e3*k_O3*data_krekling['Aug Cond']).where(np.log(data_krekling['PARo.1'])>photosynth_loglimit).dropna().values
x_sample_aug, pdf_aug, fit_aug, stat_aug = fit_skew_normal(stomatal_conductance_o3_aug[np.where(stomatal_conductance_o3_aug>500)[0]])

quantiles = stats.skewnorm.interval(0.5, fit["shape"], fit["location"], fit["scale"])
gmax = stat_aug['mean']
fmin = stomatal_conductance_o3.min()/gmax

print("gmax = %3.2f mu mol O_3 m^-2 s^-1 kg^-1" % gmax )
print("fmin = %3.2f" % fmin )

from scipy.optimize import curve_fit

# f_light
ppfd = np.concatenate((((data_krekling['PARo'].where((np.log(data_krekling['PARo'])>photosynth_loglimit)&(data_krekling['June Cond'].notna())).dropna())).values, ((data_krekling['PARo.1'].where((np.log(data_krekling['PARo.1'])>photosynth_loglimit)&(data_krekling['Aug Cond'].notna())).dropna())).values, ((data_krekling['PARo.2'].where((np.log(data_krekling['PARo.2'])>photosynth_loglimit)&(data_krekling['Sept Cond'].notna())).dropna())).values))

fit_params, cov_mat = curve_fit(f_light, ppfd, stomatal_conductance_o3/gmax, p0=[1e-3,])
print("Fit results")
print("alpha", fit_params)
print("cov mat", cov_mat)

# f_temp
temp = np.concatenate((((data_krekling['Tleaf (air)'].where((np.log(data_krekling['PARo'])>photosynth_loglimit)&(data_krekling['June Cond'].notna())).dropna())).values, ((data_krekling['Tleaf (air).1'].where((np.log(data_krekling['PARo.1'])>photosynth_loglimit)&(data_krekling['Aug Cond'].notna())).dropna())).values, ((data_krekling['Tleaf (air).2'].where((np.log(data_krekling['PARo.2'])>photosynth_loglimit)&(data_krekling['Sept Cond'].notna())).dropna())).values))

fit_params_temp, cov_mat_temp = curve_fit(f_temp, (temp, fmin), stomatal_conductance_o3/gmax, p0=[12,0,25])
print("Fit results")
print("alpha", fit_params_temp)
print("cov mat", cov_mat_temp)

# f_vpd
relHum = np.concatenate((((data_krekling['RH_R'].where((np.log(data_krekling['PARo'])>photosynth_loglimit)&(data_krekling['June Cond'].notna())).dropna())).values, ((data_krekling['RH_R.1'].where((np.log(data_krekling['PARo.1'])>photosynth_loglimit)&(data_krekling['Aug Cond'].notna())).dropna())).values, ((data_krekling['RH_R.2'].where((np.log(data_krekling['PARo.2'])>photosynth_loglimit)&(data_krekling['Sept Cond'].notna())).dropna())).values))

vpd = VPD(relHum, temp, version='buck')/kilo
