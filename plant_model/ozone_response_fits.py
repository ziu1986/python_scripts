from scipy.optimize import curve_fit

def f_flunder(x, x_std, y, y_std):
    '''
    Bring flundered arrays to same shape and size.
    
    Parameters
    ----------
    x : array
    1. dimension of the data

    x_std : array
    1. dimension uncertainty of the data

    y : array
    2. dimension of the data

    y_std : array
    2. dimension uncertainty of the data

    Returns
    -------
    flat_x, flat_x_std, flat_y, flat_y_std
    Flattened arrays of same sizes.
    
    '''
    if (flunder(y)[~np.isnan(flunder(y))].size <= flunder(x)[~np.isnan(flunder(x))].size):
        flat_y = flunder(y)[~np.isnan(flunder(y))]
        flat_y_std = flunder(y_std)[~np.isnan(flunder(y))]
        flat_x = flunder(x)[~np.isnan(flunder(y))]
        try:
            flat_x_std = flunder(x_std)[~np.isnan(flunder(y))]
        except TypeError:
            flat_x_std = 0
    else:
        flat_y = flunder(y)[~np.isnan(flunder(x))]
        flat_y_std = flunder(y_std)[~np.isnan(flunder(x))]
        flat_x = flunder(x)[~np.isnan(flunder(x))]
        try:
            flat_x_std = flunder(x_std)[~np.isnan(flunder(x))]
        except TypeError:
            flat_x_std = 0
    #print(flat_y, flat_y_std, flat_x_y, flat_x_std_y)
    return(flat_x, flat_x_std, flat_y, flat_y_std)


# Define flipping function arguments
flip = lambda f: lambda *a: f(*reversed(a))

def poly1(x, m):
    '''
    Line
    '''
    return(m*x+1)

def poly_origin(x, m):
    '''
    Line through origin
    '''
    return(m*x)

def poly_free(x, m, b=None):
    '''
    Line through with free ordinate
    '''
    if b is not None:
        return(m*x+b)
    else:
        return(m[0]*x+m[1])

def poly2(x, m):
    '''
    Parabol
    '''
    return(m*(x-0)**2+1)

def expo(x, m):
    '''
    Exponential
    '''
    return(np.exp(x*m))

def rms(y, yfit):
    return np.sqrt(np.sum((y-yfit)**2)/np.size(y))

def or_fit(x, y, x_std, y_std, **karg):

    deg = karg.pop("deg", 1)
    fit_range = karg.pop('range', np.arange(0,100))

    flat_x, flat_x_std, flat_y, flat_y_std = f_flunder(x, x_std, y, y_std)

    p0 = [(flat_y.max()-flat_y.min())/(flat_x.min()-flat_x.max()),]
        
    if deg==1:
        func = poly1
    elif deg==2:
        func = poly2
    elif deg=='exp':
        func = expo
    elif deg=='origin':
        func = poly_origin
    elif deg=='free':
        p0.append(1)
        func = poly_free
    print('first guess', p0)

    if type(flat_x_std) == int:
        
        # Unweighted fit
        popt, pcov = curve_fit(func, flat_x, flat_y, p0)
        yfit = func(fit_range, *popt)

        print('Unweighted fit parameters:', popt)
        print('Covariance matrix:'); print(pcov)
        print('rms error in fit:', rms(flat_y, func(flat_x, *popt)))

        # Weighted fit
        popt2, pcov2 = curve_fit(func, flat_x, flat_y, p0, sigma=flat_y_std, absolute_sigma=True)
        yfit2 = func(fit_range, *popt2)

        print('Weighted fit parameters:', popt2)
        print('Covariance matrix:'); print(pcov2)
        print('rms error in fit:', rms(flat_y, func(flat_x, *popt2)))

        return(yfit, yfit2)
    else:
        print('X-Y-uncertainty fit parameters:')
        # Need to flip function arguments, because ODR expects parameters (list!) to come first,
        # while the other fitters expect them to be last
        flip_func = flip(func)
        x_fit, y_fit = xy_uncertainty_regression(flat_x, flat_y, flat_x_std, flat_y_std, flip_func, p0, range=fit_range)
        
        
        return(y_fit)
        

    
def xy_uncertainty_regression(x, y, x_err, y_err, func, par, **karg):
    from scipy.odr import *

    fit_range = karg.pop('range', np.arange(0,100))
    # Create a model for fitting.
    model = Model(func)

    # Create a RealData object using our initiated data from above.
    data = RealData(x, y, sx=x_err, sy=y_err)

    # Set up ODR with the model and data.
    odr = ODR(data, model, beta0=par)

    # Run the regression.
    out = odr.run()

    # Use the in-built pprint method to give us results.
    out.pprint()
    print('rms error in fit:', rms(y, func(out.beta, x)))
    x_fit = fit_range
    y_fit = func(out.beta, x_fit)

    return(x_fit, y_fit)

print("+++ gs +++")
yfit_gs, yfit_gs2 = or_fit(pcuo, gs, 0, gs_std)
yfit_gs_2 = or_fit(pcuo, gs, pcuo_std, gs_std)
yfit_gs_free, yfit_gs2_free = or_fit(pcuo, gs, 0, gs_std, deg='free')
yfit_gs_2_free = or_fit(pcuo, gs, pcuo_std, gs_std, deg='free')
print("+++ An +++")
yfit_A, yfit_A2 = or_fit(pcuo, A, 0, A_std)
yfit_A_2 = or_fit(pcuo, A, pcuo_std, A_std)
yfit_A_free, yfit_A2_free = or_fit(pcuo, A, 0, A_std, deg='free')
yfit_A_2_free = or_fit(pcuo, A, pcuo_std, A_std, deg='free')
print("+++ Jmax +++")
yfit_Jmax, yfit_Jmax2 = or_fit(pcuo, Jmax, 0, Jmax_std)
yfit_Jmax_2 = or_fit(pcuo, Jmax, pcuo_std, Jmax_std)
print("+++ Vcmax +++")
yfit_Vcmax, yfit_Vcmax2 = or_fit(pcuo, Vcmax, 0, Vcmax_std)
yfit_Vcmax_2 = or_fit(pcuo, Vcmax, pcuo_std, Vcmax_std)
print("+++ Rd +++")
yfit_Rd, yfit_Rd2 = or_fit(pcuo, Rd, 0, Rd_std)
yfit_Rd_2, yfit_Rd2_2 = or_fit(pcuo, Rd, 0, Rd_std, deg='exp')

#yfit_Chl, yfit_Chl2 = or_fit(pcuo, Chl, 0, Chl_std)

# Fit Jmax-Vcmax ratios
print("+++ JmaxO3/Jmax0/VcmaxO3/Vcmax0 +++")
yfit_VcmaxVSJmax, yfit_VcmaxVSJmax2 = or_fit(Jmax, Vcmax, 0, Vcmax_std, deg='origin', range=np.arange(0, 2.1, 0.1))
yfit_VcmaxVSJmax_2 = or_fit(Jmax, Vcmax, Jmax_std, Vcmax_std, deg='origin', range=np.arange(0, 2.1, 0.1))
