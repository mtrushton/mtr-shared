import numpy as np

def dred(wavelength, flux, flux_error, EBV):
    """
    Dereddens flux for a given E(B-V) assuming galactic 
    extinction law of Howarth et al. (1983), MNRAS, 203, 301 
    
    Inputs:
    Wavelength: microns
    Flux
    Flux_error
    EBV: Colour excess E(B-V)
    
    Outputs:
    tuple: (flux_corr, flux_corr_error) arrays
    """
    
    R = 3.1
    x = []
    
    x = 1 / wavelength
    
    if 2.75 >= x >= 1.83:
        Xx = R + 2.56 * (x - 1.83) + 0.993 * (x - 1.83)**2
    elif 9.0 >= x >= 2.75:
        Xx = R - 0.236 + (0.462 * x) + (0.105 * x**2) + 0.454 / ((x - 4.557)**2 + 0.293)
    else:
        Xx = ((1.86 - 0.48 * x) * x - 0.1) * x

    flux_corr = flux * 10**(0.4 * Xx * EBV)
    flux_corr_error = flux_error / flux * flux_corr
    
    return flux_corr, flux_corr_error

