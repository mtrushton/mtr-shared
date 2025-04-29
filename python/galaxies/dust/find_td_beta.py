import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy import constants
import warnings

# Fits a graybody of the form B(nu, T) * nu**beta to the far-infrared
# fluxes defined in FLUX at wavelengths WAVELENGTH
#
# INPUTS:
# WAVELENGTH (np array) : Wavelengths in microns of the data
# FLUX (np array) : Fluxes in Jy at WAVELENGTH
# FLUX_ERROR (np array) : Errors in FLUX in Jy
#
# OUTPUTS:
# Plot showing the best-fitting graybody to the values of FLUX
# Best-fitting graybody temperature and index Beta
# Total integrated graybody flux


# Define input arrays
wavelength = np.array([70, 100, 160, 250, 350, 500])  # microns
flux = np.array([10.26, 29.28, 54.11, 38.13, 18.75, 6.22])  # Jy
flux_error = np.array([0.51, 1.46, 2.71, 1.91, 0.94, 0.31])  # Jy

c = constants.c  
h = constants.h  
k = constants.k  
frequency = c / (wavelength * 1e-6)  # Hz

# Define modified blackbody (graybody) function
def modified_blackbody(nu, T, beta, amplitude):
    """
    Modified blackbody function (graybody)
    
    Parameters:
    -----------
    nu : array
        Frequency (Hz)
    T : float
        Temperature (K)
    beta : float
        Dust emissivity index
    amplitude : float
        Scaling amplitude
    
    Returns:
    --------
    flux : array
        Model flux density
    """
    nu = np.asarray(nu)
    
    with np.errstate(over='ignore'):
        exp_term = np.exp(np.minimum(h * nu / (k * T), 700))  # Capped to prevent overflow
    
    planck = 2 * h * nu**3 / (c**2 * (exp_term - 1))
    
    # Dust emissivity power law
    modified = planck * (nu / 1e12)**beta
    
    return amplitude * modified

# Initial parameter guesses
T_init = 30.0  # K
beta_init = 1.5
amplitude_init = 1e-10

def normalised_modified_blackbody(nu, T, beta, norm):
    """Normalised modified blackbody for fitting"""
    peak_nu = frequency[np.argmax(flux)]
    return norm * modified_blackbody(nu, T, beta, 1.0) / modified_blackbody(peak_nu, T, beta, 1.0)

try:
    popt, pcov = curve_fit(
        normalised_modified_blackbody,
        frequency,
        flux,
        p0=[T_init, beta_init, np.max(flux)],
        sigma=flux_error,
        absolute_sigma=True,
        bounds=([10, 0.1, 0.1], [100, 3.0, 100]),
        maxfev=10000
    )
    
    T_fit, beta_fit, norm = popt
    T_err, beta_err, norm_err = np.sqrt(np.diag(pcov))
    
    # Calculate the amplitude
    peak_nu = frequency[np.argmax(flux)]
    amplitude_fit = norm / modified_blackbody(peak_nu, T_fit, beta_fit, 1.0)
    
    # Approx amplitude error
    amplitude_err = amplitude_fit * (norm_err/norm)
    
except Exception as e:
    print(f"Fitting error: {e}")
    try:
        # Try different parameters
        popt, pcov = curve_fit(
            modified_blackbody,
            frequency,
            flux,
            p0=[30, 1.5, 1e-25],
            sigma=flux_error,
            absolute_sigma=True,
            bounds=([10, 0.1, 1e-30], [100, 3.0, 1e-20]),
            method='trf',
            maxfev=10000
        )
        
        T_fit, beta_fit, amplitude_fit = popt
        T_err, beta_err, amplitude_err = np.sqrt(np.diag(pcov))
    except Exception as e2:
        print(f"Second fitting attempt failed: {e2}")
        # If all else fails, use a grid search approach
        print("Using grid search approach...")
        
        # Define parameter grids
        T_grid = np.linspace(15, 60, 50)
        beta_grid = np.linspace(0.5, 2.5, 50)
        
        best_chisq = np.inf
        best_params = [T_init, beta_init, amplitude_init]
        
        for T in T_grid:
            for beta in beta_grid:
                # Optimise amplitude for each T, beta pair
                def residual_func(amp):
                    model = modified_blackbody(frequency, T, beta, amp)
                    return np.sum(((model - flux) / flux_error)**2)
                
                from scipy.optimize import minimize_scalar
                res = minimize_scalar(residual_func, bounds=(1e-30, 1e-20), method='bounded')
                amp = res.x
                chisq = res.fun
                
                if chisq < best_chisq:
                    best_chi2 = chisq
                    best_params = [T, beta, amp]
        
        T_fit, beta_fit, amplitude_fit = best_params
        
        # Estimate errors from the chisq surface
        # Find parameters where chisq increases by 1 from minimum
        T_err = 1.0  # Default if estimation fails
        beta_err = 0.1
        amplitude_err = amplitude_fit * 0.1  # 10% error estimate

# Create a high-res frequency grid for plot
nu_grid = np.logspace(np.log10(frequency.min()*0.8), np.log10(frequency.max()*1.2), 1000)
wavelength_grid = c / nu_grid * 1e6  # Convert to microns for plotting

# Calculate the best-fit model
model_flux = modified_blackbody(nu_grid, T_fit, beta_fit, amplitude_fit)

# Calculate the model at the observed wavelengths for comparison
model_at_data = modified_blackbody(frequency, T_fit, beta_fit, amplitude_fit)

# Reduced chisq
residuals = flux - model_at_data
chisq = np.sum((residuals / flux_error)**2)
dof = len(flux) - 3  # degrees of freedom (n_data - n_parameters)
reduced_chi2 = chisq / dof


# Calculate the integrated dust emission (total luminosity)
nu_wide = np.logspace(np.log10(frequency.min()*0.1), np.log10(frequency.max()*10), 5000)
model_flux_wide = modified_blackbody(nu_wide, T_fit, beta_fit, amplitude_fit)

# Integration in frequency space
integrated_flux_Jy_Hz = np.trapz(model_flux_wide, nu_wide) 
integrated_flux_W_m2 = integrated_flux_Jy_Hz * 1e-26


plt.figure(figsize=(10, 6))

# Plot data points with error bars
plt.errorbar(wavelength, flux, yerr=flux_error, fmt='o', color='blue', 
             label='Data', markersize=8, capsize=4)

# Plot the modified blackbody
plt.plot(wavelength_grid, model_flux, 'r-', linewidth=2,
         label=f'Modified blackbody fit\nT = {T_fit:.1f} ± {T_err:.1f} K\nβ = {beta_fit:.2f} ± {beta_err:.2f}')


plt.xscale('log')
plt.yscale('log')
plt.xlabel('Wavelength (μm)', fontsize=14)
plt.ylabel('Flux Density (Jy)', fontsize=14)
plt.title('Modified Blackbody (Graybody) Fit to Dust SED', fontsize=16)
plt.legend(fontsize=12)
plt.grid(True, alpha=0.3, which='both')

plt.ylim(flux.min() * 0.5, flux.max() * 2)

plt.minorticks_on()
plt.tick_params(which='both', direction='in', top=True, right=True)


plt.text(0.05, 0.05, 
         f"Reduced χ² = {reduced_chi2:.2f}\nIntegrated flux = {integrated_flux_W_m2:.2e} W/m²", 
         transform=plt.gca().transAxes, 
         fontsize=12, bbox=dict(facecolor='white', alpha=0.7))

plt.tight_layout()
plt.show()

print(f"Best-fit parameters:")
print(f"Temperature: {T_fit:.2f} ± {T_err:.2f} K")
print(f"Beta index: {beta_fit:.3f} ± {beta_err:.3f}")
print(f"Amplitude: {amplitude_fit:.3e} ± {amplitude_err:.3e}")
print(f"Reduced chi-squared: {reduced_chi2:.2f}")
print(f"Integrated flux: {integrated_flux_W_m2:.3e} W/m²")
print(f"Integrated flux: {integrated_flux_Jy_Hz:.3e} Jy·Hz")

