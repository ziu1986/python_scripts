def stomatal_conductance(M, g0, g1):
    '''
    x = VPD^-0.5
    y = An/[CO2]
    '''
    # Unpack the 2d coordinates
    x, y = M
    # Compute the stomtal conductance
    gs = g0 + 1.6*(1+g1*x)*y

    # flatten the output
    return(np.ravel(gs))

x = []
y = []
z = []

for each in data_list:
    selection = each.where((each['"Cond"']>0) & (each['"Photo"']>0) & (each['"Photo"']/each['"Ci"']<0.2) & (each['"Ci"']>0)).dropna()
    y.append((selection['"Photo"']/selection['"Ci"']).values)
    x.append((1/np.sqrt(selection['"VpdL"'])).values)
    z.append((selection['"Cond"']).values)

X, Y = np.meshgrid(flunder(x), flunder(y))
xdata = np.vstack((X.ravel(), Y.ravel()))

Z = np.meshgrid(flunder(z), flunder(z))[0]

from scipy.optimize import curve_fit

fit_params, cov_mat = curve_fit(stomatal_conductance, xdata, Z.ravel(), p0=[0, 0.5])

