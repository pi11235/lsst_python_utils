# Import required libraries
import logging
from astropy.io import fits
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
print('This appears when you import fiberSpectrograph_utils')


def FiberSpectrograph_Fit_and_Analysis(filename,
                                       estimates=[None, None, None, None],
                                       return_data=False,
                                       plot=True,
                                       xlim=[None,None],
                                       ylim=[None,None]):

    # Function returns an object with fit data
    # return_data option will return the raw data as well
    # Plot parameter plots the data and fit

    # Define Gaussian Fitting Function
    def gaus(x, a, x0, sigma, offset):
        return a * np.exp(-(x - x0) ** 2 / (2 * sigma ** 2)) + offset

    # read in file
    logging.info('Opening file {}'.format(filename))
    hdul = fits.open(filename)
    xdata=(hdul[0].data)[:,0]
    ydata=(hdul[0].data)[:,1]
    #print(repr(hdul[0].header))

    # Fit a guassian to the profile
    # check if estimates are given
    if len(estimates) !=4 :
        raise TypeError('estimates parameter requires a length of 4, {} were given'.format(len(estimates)))
    
    # determine amplitude estimate
    if estimates[0] == None:
        Amp = np.max(ydata) # assume the peak we want is the maximum value (NOT ALWAYS THE CASE!)
        logging.info('No Amplitude estimate given, using {}'.format(Amp))
    else: 
        Amp = estimates[0]
    
    # determine center wavelength estimate
    if estimates[1] == None:
        mean = xdata[np.argmax(ydata)] # assume the peak we want is the maximum value (NOT ALWAYS THE CASE!)
        logging.info('No center wavelength estimate given, using max value of data of {}'.format(mean))
    else: 
        mean = estimates[1]
        
    # determine sigma estimate
    if estimates[2] == None:
        sigma = 2.0
        logging.info('No estimate for gaussian sigma, using {}'.format(sigma))
    else: 
        sigma = estimates[2]
        
    # determine zero point estimate
    if estimates[3] == None:
        offset=np.median(ydata)
        logging.info('No estimate for zero point estimate, using median of data which is {}'.format(offset))
    else: 
        offset = estimates[3]
    
    popt,pcov = curve_fit(gaus,xdata,ydata,p0=[Amp, mean, sigma, offset])
    
    if plot == True:
        # Plot data
        plt.plot(xdata,ydata,color='r', linestyle=':', label='data')
        plt.plot(xdata, gaus(xdata, *popt ), color='b', label='fit')
        if xlim != None:
            plt.xlim(xlim)
        if ylim == True:
            plt.xlim(xlim)
        plt.legend()
        plt.title('Fit to fiber spectrograph data')
        plt.ylabel('Intensity [ADU]')
        plt.xlabel('Wavelength [nm]')
    
    # create return object
    if return_data == True:
        return(popt, xdata, ydata)
    else:
        return(popt)

