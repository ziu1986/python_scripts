from scipy.optimize import curve_fit

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

def poly_free(x, m, b):
    '''
    Line through with free ordinate
    '''
    return(m*x+b)

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

    f_y, f_y_std, f_x_y, f_x_std_y = f_flunder(x, y, x_std, y_std)
            
    if deg==1:
        p0 = [(f_y.max()-f_y.min())/(f_x_y.min()-f_x_y.max()),]
        func = poly1
    elif deg==2:
        p0 = ((f_y.max()-f_y.min())/(f_x_y.min()-f_x_y.max()),)
        func = poly2
    elif deg=='exp':
        p0 = ((f_y.max()-f_y.min())/(f_x_y.min()-f_x_y.max()),)
        func = expo
    elif deg=='origin':
        p0 = [(f_y.max()-f_y.min())/(f_x_y.min()-f_x_y.max()),]
        func = poly_origin
    elif deg=='free':
        p0 = [(f_y.max()-f_y.min())/(f_x_y.min()-f_x_y.max()), 1]
        func = poly_free
    

    if type(f_x_std_y) == int:
        
        # Unweighted fit
        popt, pcov = curve_fit(func, f_x_y, f_y, p0)
        yfit = func(fit_range, *popt)

        print('Unweighted fit parameters:', popt)
        print('Covariance matrix:'); print(pcov)
        print('rms error in fit:', rms(f_y, func(f_x_y, *popt)))

        # Weighted fit
        popt2, pcov2 = curve_fit(func, f_x_y, f_y, p0, sigma=f_y_std, absolute_sigma=True)
        yfit2 = func(fit_range, *popt2)

        print('Weighted fit parameters:', popt2)
        print('Covariance matrix:'); print(pcov2)
        print('rms error in fit:', rms(f_y, func(f_x_y, *popt2)))

        return(yfit, yfit2)
    else:
        print('X-Y-uncertainty fit parameters:')
        flip_func = flip(func)
        x_fit, y_fit = xy_uncertainty_regression(f_x_y, f_y, f_x_std_y, f_y_std, flip_func, range=fit_range)
        
        
        return(y_fit)
        

def f_flunder(x, y, x_std, y_std):
    if (flunder(y)[~np.isnan(flunder(y))].size <= flunder(x)[~np.isnan(flunder(x))].size):
        f_y = flunder(y)[~np.isnan(flunder(y))]
        f_y_std = flunder(y_std)[~np.isnan(flunder(y))]
        f_x_y = flunder(x)[~np.isnan(flunder(y))]
        try:
            f_x_std_y = flunder(x_std)[~np.isnan(flunder(y))]
        except TypeError:
            f_x_std_y = 0
    else:
        f_y = flunder(y)[~np.isnan(flunder(x))]
        f_y_std = flunder(y_std)[~np.isnan(flunder(x))]
        f_x_y = flunder(x)[~np.isnan(flunder(x))]
        try:
            f_x_std_y = flunder(x_std)[~np.isnan(flunder(x))]
        except TypeError:
            f_x_std_y = 0
    #print(f_y, f_y_std, f_x_y, f_x_std_y)
    return(f_y, f_y_std, f_x_y, f_x_std_y)
    
def xy_uncertainty_regression(x, y, x_err, y_err, func, **karg):
    from scipy.odr import *

    fit_range = karg.pop('range', np.arange(0,100))
    # Create a model for fitting.
    model = Model(func)

    # Create a RealData object using our initiated data from above.
    data = RealData(x, y, sx=x_err, sy=y_err)

    # Set up ODR with the model and data.
    odr = ODR(data, model, beta0=[0.,])

    # Run the regression.
    out = odr.run()

    # Use the in-built pprint method to give us results.
    out.pprint()
    '''Beta: [ 1.01781493  0.48498006]
    Beta Std Error: [ 0.00390799  0.03660941]
    Beta Covariance: [[ 0.00241322 -0.01420883]
    [-0.01420883  0.21177597]]
    Residual Variance: 0.00632861634898189
    Inverse Condition #: 0.4195196193536024
    Reason(s) for Halting:
    Sum of squares convergence'''
    print('rms error in fit:', rms(y, func(out.beta, x)))
    x_fit = fit_range
    y_fit = func(out.beta, x_fit)

    return(x_fit, y_fit)

print("+++ gs +++")
yfit_gs, yfit_gs2 = or_fit(pcuo, gs, 0, gs_std)
yfit_gs_2 = or_fit(pcuo, gs, pcuo_std, gs_std)
yfit_gs_free, yfit_gs2_free = or_fit(pcuo, gs, 0, gs_std, deg='free')
print("+++ An +++")
yfit_A, yfit_A2 = or_fit(pcuo, A, 0, A_std)
yfit_A_2 = or_fit(pcuo, A, pcuo_std, A_std)
yfit_A_free, yfit_A2_free = or_fit(pcuo, A, 0, A_std, deg='free')
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
