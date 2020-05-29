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
    gs : float
        relative stomatal conductance g_sto^rel
    '''
    gs = 1-np.exp(-alpha*M)

    return(gs)

def f_temp(M, Topt, Tmin, Tmax):
    '''
    Parameters:
    -----------
    M : array
        0: Temperature: [T] = deg C
        1: fmin
    Topt: float
        [Topt] = deg C 
    Tmin: float
        [Tmin] = deg C 
    Tmax: float
        [Tmax] = deg C 
    Returns:
    --------
    gs : float
        relative stomatal conductance g_sto^rel
    '''
    T = M[0]
    fmin = M[1]
    
    bt = float(Tmax-Topt)/(Topt-Tmin)
    gs = (T-Tmin)/float(Topt-Tmin)*((Tmax-T)/float(Tmax-Topt))**bt
    
    gs[np.where((T<Tmin) | (T>Tmax) | (gs<fmin))[0]] = fmin

    return(gs)

def f_vpd(M, VPDmax, VPDmin):
    '''
    Parameters:
    -----------
    M : array
        0: VPD : [VPD] = kPa
        1: fmin
    VPDmax: float
        [VPDmax] = kPa 
    VPDmin: float
        [VPDmin] = kPa 
    Returns:
    --------
    gs : float
        relative stomatal conductance g_sto^rel
    '''
    vpd = M[0]
    fmin = M[1]

    gs = ((1-fmin)*(VPDmin-vpd)/float(VPDmin-VPDmax)) + fmin
    gs[np.where(gs < fmin)] = fmin
    gs[np.where(gs > 1)] = 1.

    return(gs)
    

# gsto 

photosynth_loglimit = 4

g_sto_o3 = np.concatenate(((1e3*k_O3*data_krekling['June Cond']).where(np.log(data_krekling['PARo'])>photosynth_loglimit).dropna().values,(1e3*k_O3*data_krekling['Aug Cond']).where(np.log(data_krekling['PARo.1'])>photosynth_loglimit).dropna().values,(1e3*k_O3*data_krekling['Sept Cond']).where(np.log(data_krekling['PARo.2'])>photosynth_loglimit).dropna().values))
x_sample, pdf, fit, stat = fit_skew_normal(g_sto_o3)

g_sto_o3_month = {"Jun":(1e3*k_O3*data_krekling['June Cond'].where(np.log(data_krekling['PARo'])>photosynth_loglimit)),
                  "Aug":(1e3*k_O3*data_krekling['Aug Cond'].where(np.log(data_krekling['PARo.1'])>photosynth_loglimit)),
                  "Sep":(1e3*k_O3*data_krekling['Sept Cond'].where(np.log(data_krekling['PARo.2'])>photosynth_loglimit))}
fit_results_month = {}
for imonth in ('Jun', 'Aug', 'Sep'):
    x_sample_month, pdf_month, fit_month, stat_month = fit_skew_normal(g_sto_o3_month[imonth].dropna().values)
    fit_results_month[imonth] = {"x_sample":x_sample_month, "pdf":pdf_month, "fit":fit_month, "stat":stat_month}

quantiles = stats.skewnorm.interval(0.5, fit["shape"], fit["location"], fit["scale"])
gmax = g_sto_o3_month["Aug"].dropna().values[np.where(g_sto_o3_month["Aug"].dropna().values>500)[0]].mean()
fmin = g_sto_o3.min()/gmax

print("gmax = %3.2f mu mol O_3 m^-2 s^-1 kg^-1" % gmax )
print("fmin = %3.2f" % fmin )

A_net_o3_month = {"Jun":(data_krekling['June Photo'].where(np.log(data_krekling['PARo'])>photosynth_loglimit)),
                  "Aug":(data_krekling['Aug Photo'].where(np.log(data_krekling['PARo.1'])>photosynth_loglimit)),
                  "Sep":(data_krekling['Sept Photo'].where(np.log(data_krekling['PARo.2'])>photosynth_loglimit))}


from scipy.optimize import curve_fit

# f_light
ppfd = np.concatenate((((data_krekling['PARo'].where((np.log(data_krekling['PARo'])>photosynth_loglimit)&(data_krekling['June Cond'].notna())).dropna())).values, ((data_krekling['PARo.1'].where((np.log(data_krekling['PARo.1'])>photosynth_loglimit)&(data_krekling['Aug Cond'].notna())).dropna())).values, ((data_krekling['PARo.2'].where((np.log(data_krekling['PARo.2'])>photosynth_loglimit)&(data_krekling['Sept Cond'].notna())).dropna())).values))

fit_params, cov_mat = curve_fit(f_light, ppfd, g_sto_o3/gmax, p0=[1e-3,])
print("Fit results")
print("alpha", fit_params)
print("cov mat", cov_mat)

# f_temp
temp = np.concatenate((((data_krekling['Tleaf (air)'].where((np.log(data_krekling['PARo'])>photosynth_loglimit)&(data_krekling['June Cond'].notna())).dropna())).values, ((data_krekling['Tleaf (air).1'].where((np.log(data_krekling['PARo.1'])>photosynth_loglimit)&(data_krekling['Aug Cond'].notna())).dropna())).values, ((data_krekling['Tleaf (air).2'].where((np.log(data_krekling['PARo.2'])>photosynth_loglimit)&(data_krekling['Sept Cond'].notna())).dropna())).values))

fit_params_temp, cov_mat_temp = curve_fit(f_temp, (temp, fmin), g_sto_o3/gmax, p0=[12,0,25])
print("Fit results")
print("Topt, Tmin, Tmax", fit_params_temp)
print("cov mat", cov_mat_temp)

# f_vpd
relHum = np.concatenate((((data_krekling['RH_R'].where((np.log(data_krekling['PARo'])>photosynth_loglimit)&(data_krekling['June Cond'].notna())).dropna())).values, ((data_krekling['RH_R.1'].where((np.log(data_krekling['PARo.1'])>photosynth_loglimit)&(data_krekling['Aug Cond'].notna())).dropna())).values, ((data_krekling['RH_R.2'].where((np.log(data_krekling['PARo.2'])>photosynth_loglimit)&(data_krekling['Sept Cond'].notna())).dropna())).values))

vpd = VPD(relHum, temp, version='buck')/kilo

fit_params_vpd, cov_mat_vpd = curve_fit(f_vpd, (vpd, fmin), g_sto_o3/gmax, p0=[0.5, 2.])
print("Fit results")
print("VPDmax, VPDmin", fit_params_vpd)
print("cov mat", cov_mat_vpd)
