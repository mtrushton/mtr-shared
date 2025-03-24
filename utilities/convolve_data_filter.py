# -*- coding: utf-8 -*-
import pandas as pd
import numpy as np
from scipy import integrate
from scipy import interpolate
import matplotlib.pyplot as plt

def normalise_response(wavelength, response):
    """
    Takes a response curve for a filter and outputs
    a normalised curve.
    
    INPUTS:
     wavelength : array-like
        Wavelength values in Angstroms or nm
    response : array-like
        Filter response at each wavelength
    
    Returns:
    normalised_response : array-like 
        response curve
 
    """
    integrated_response = integrate.simpson(response, wavelength)
    normalised_response = response / integrated_response
    
    return wavelength, normalised_response
    
def get_effective_wavelength(wavelength, response):
    """
    Calculate the effective wavelength of a filter.
    
    INPUTS:
    wavelength : array-like
        Wavelength values in Angstroms or nm
    response : array-like
        Filter response at each wavelength
        
    Returns:
    float : The effective wavelength of the filter
    """
    # Normalise the response first (although for this calculation, it's not strictly necessary)
    wavelength, normalised_response = normalise_response(wavelength, response)
    
    # Calculate the effective wavelength
    numerator = integrate.simpson(wavelength * normalised_response, wavelength)
    denominator = integrate.simpson(normalised_response, wavelength)
    
    effective_wavelength = numerator / denominator
    
    return effective_wavelength
    
def conv_model(wave_model, flux_model, band):
    """
    Convolves FLUX_MODEL with Johnson filters and derives 
    convolved flux and magnitude
    Note file_filter needs to be accessible via dir
    
    INPUTS:
    wave_model : array-like
        wavelengths in FLUX_MODEL
    flux_model : array-like
        fluxes for convolution
   
    Returns:
    effective wavelength : array-lie
        the derived effective wavelength for each filter (microns)
    total_flux : array-like
        the integrated flux for each filter
    mag : array-like
        derived magnitude for each filter 
    
    """
    dir = 'filters/'
    
    band = band.upper()
    
    file_filter = {
    'band': ['U', 'B', 'V', 'R', 'I'],
    'file_response': ['Generic_Johnson.U.dat', 'Generic_Johnson.B.dat', 
                      'Generic_Johnson.V.dat', 'Generic_Johnson.R.dat', 'Generic_Johnson.I.dat'],
    'zero_point': [3.49719e-9, 6.72553e-9, 3.5833e-9, 1.87529e-9, 9.23651e-10]
    }
    
   
    if band not in file_filter['band']:
        print(f"Error: Band '{band}' not found in file_filter.")
        return None, None
    
    index = file_filter['band'].index(band)
    file_response = file_filter['file_response'][index]
    zero_point = file_filter['zero_point'][index]
    
    df = pd.read_csv(dir + file_response, sep='\s+', names=['wavelength', 'response'])
    wavelength = df['wavelength']
    response = df['response']
    
    wavelength, normalised_response = normalise_response(wavelength, response)
   
    f = interpolate.interp1d(wavelength, normalised_response, kind='linear', bounds_error=False, fill_value=np.nan)
    response_resampled = f(wave_model)
    
    mask = ~np.isnan(response_resampled)
    
    total_flux = integrate.simpson(response_resampled[mask] * flux_model[mask], wave_model[mask])
    
    mag = -2.5 * np.log10(total_flux / zero_point)
    
    return get_effective_wavelength(wavelength, response)/ 10000, total_flux, mag
    
