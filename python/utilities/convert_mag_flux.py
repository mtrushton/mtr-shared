import pandas as pd
import numpy as np

def convert_mag_flux(band_obs, mag_obs, mag_obs_error=None):
    """
    Calculates FLUX corresponding to mag_obs at band_obs, and
    returns FLUX_ERROR when mag_obs_error input.
    
    INPUTS
    BAND_OBS (str): Photomertic band
    MAG_OBS (float): BAND_OBS magnitude
    MAG_OBS_ERROR (float): ERROR in MAG_OBS (if known)
    
    OUTPUTS
    tuple: (wavelength, flux, flux_error)
    """
     
    # Johnson-Cousins-Glass System (Bessell et al. 1998)
    filter_name = ['U', 'B', 'V', 'R', 'I', 'J', 'H', 'K']
    lambda_eff = [0.36, 0.438, 0.545, 0.641, 0.798, 1.22, 1.63, 2.19] 
    f_zero = [1810, 4260, 3640, 3080, 2550, 1600, 1080, 670]  # Jy
    
    band_obs = band_obs.upper() # allows for lowercase input
    
    # Check band_obs is in filter list
    if band_obs not in filter_name:
        raise ValueError(f"Band '{band_obs}' not found in filter list: {filter_name}")
    
    band_idx = filter_name.index(band_obs)
    
    wavelength = lambda_eff[band_idx]
    
    zero_point = f_zero[band_idx]
    
    # Calculate flux and error (if provided) in Jy 
    flux = zero_point * 10**(-0.4 * mag_obs)
    
     if mag_obs_error is not None:
        flux_error = flux * 0.4 * np.log(10) * mag_obs_error
        return wavelength, flux, flux_error
    else:
        return wavelength, flux, None
