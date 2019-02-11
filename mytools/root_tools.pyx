import numpy as np
from ROOT import gROOT, TCanvas, TH2D, TProfile, TLine

def atm_var_time_series(datetime, atm_var, atm_var_label):
    '''
    Get interpolation function, TH2D and TProfile from pressure level slice of atm_var
    '''
    gROOT.Reset()
    cdef int delta_time = (datetime[0][1]-datetime[0][0]).seconds/60**2
    cdef int n_max = np.max([len(each) for each in datetime[1:-1]]) # number of hours October 1 to April 30 per time step
    cdef int timeBinStart = 0 
    cdef int timeBinStop = timeBinStart+delta_time*n_max+10 # data point each 10 hours
    cdef int atmVarBinStart = round(np.min([np.min(each) for each in atm_var[1:-1]]))
    cdef int atmVarBinStop = round(np.max([np.max(each) for each in atm_var[1:-1]]))
       
    hWinter = TH2D("h%s" % (atm_var_label[0][:3]),"", n_max+10, timeBinStart, timeBinStop,
                   236, atmVarBinStart-abs(0.10*atmVarBinStart), atmVarBinStop+abs(0.10*atmVarBinStop))
    
    october_hours = [each[0].hour for each in datetime]
     
    cdef int i = 0
    cdef int j = 0
    for iUm1 in atm_var[1:-1]:
        for j in np.arange(len(iUm1)):
            date = october_hours[i]+j*delta_time
            hWinter.Fill(date, iUm1[j])
            
        i = i+1

    prWinter = hWinter.ProfileX("pr%s" % (atm_var_label[0][:3]), timeBinStart, timeBinStop)
    hWinter.GetYaxis().SetTitle('%s (%s)' % atm_var_label)
    hWinter.GetYaxis().CenterTitle()
    
    hWinter.GetXaxis().SetTitle('%s (hours since Oct, 1)' % ('Time'))
    hWinter.GetXaxis().CenterTitle()

    mean_x = [0]
    mean_y = [prWinter.GetBinContent(1), ]
    sem = [prWinter.GetBinError(1), ]
    sd = [sem[-1]*np.sqrt(prWinter.GetBinEntries(1)), ]
    
    # For all histogram types: nbins, xlow, xup
    #    bin = 0;       underflow bin
    #    bin = 1;       first bin with low-edge xlow INCLUDED
    #    bin = nbins;   last bin with upper-edge xup EXCLUDED
    #    bin = nbins+1; overflow bin
    cdef int iBin
    for iBin in np.arange(1, prWinter.GetNbinsX()+1):
        #if prWinter.GetBinEntries(iBin) > 0:
        mean_x.append(prWinter.GetBinCenter(iBin))
        mean_y.append(prWinter.GetBinContent(iBin))
        sem.append(prWinter.GetBinError(iBin))
        sd.append(sem[-1]*np.sqrt(prWinter.GetBinEntries(iBin)))

    # Now get a interpolation function
    from scipy.interpolate import interp1d
    # Define interpolation function
    f_hour_to_um1 = interp1d(mean_x, mean_y)

    return hWinter, prWinter, f_hour_to_um1, np.array(sem), np.array(sd)
