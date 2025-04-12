import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter
from dred import dred

# Plots SED of V745 Sco from Andicam (SMARTS), GALEX, 2MASS, AKARI and 
# WISE. Additionally, An optical spectrum from SMARTS is also included. 
# An objective here is to find a suitable model for the quiescent SED which 
# serves as an input for a dusty model of Spitzer data (also included here). 
# This model only needs to be an approximate representation of the optical/
# near IR data, providing a reasonable fit rather than a perfect one. The 
# model should be a good representation of the data below about 0.5 micron, 
# as the RT code is  sensitive to the flux at these wavelengths. Note that 
# the flux at WISE/AKARI wavelengths shows an excess due to dust emission, 
# and should be excluded when evaluating the goodness of the fit to the 
# 'stellar' SED. 

# The model for the non-dust SED consists of a power-law in wavelength at 
# short wavelengths and a blackbody at lower wavelengths. The blackbody is a 
# suitable choice, as the model does not require a perfect representation of 
# the source, given that longer wavelengths are less crucial for RT modelling.
 

# INPUTS
#--------
# EBV (scalar): Colour excess (interstellar reddening)
# POWER_LAW_INDEX (scalar) : Index of the power law at short wavelengths
# SCALING_POWER_LAW (scalar) : value of the scaling needed for the power law 
#     to match data
# TEMPERATURE (scalar) : temperature of the blackbody approximating long 
#     wavelength data
# SCALING_BLACKBODY (scalar) : value of the scaling needed for the blackbody
#     to match data 

ebv = 1.0
# SED model parameters. Fine-tune as necessary
temperature=2700,      # Blackbody temperature
bb_scaling=1.3e-25,   # Blackbody scaling factor
power_law_index=-1.6,    # Spectral index
power_law_scaling=1.0e-14

# (assumed) effective wavelengths for Andicam bands (in microns)
ANDICAM_WAVELENGTHS = {
    'B': 0.44,   # Blue band
    'V': 0.55,   # Visual band
    'R': 0.64,   # Red band
    'I': 0.80,   # Infrared band
    'J': 1.22,   # J band
    'H': 1.63,   # H band
    'K': 2.19    # K band
}

dir_data = '/home/mtr/Downloads/V745 Sco/V745Sco/'
data_file = 'broad-fluxes_galex.dat'
spec_file = '170618_28.ascii'
andicam_file = 'V745Sco_SMARTS_BVRIJHK_Fluxes.dat'
merged_file = 'V745Sco_merged.dat'

# Manually read the data file
def read_data_file(filepath):
    wavelengths = []
    fluxes = []
    errors = []
    surveys = []
    
    with open(filepath, 'r') as f:
        for line in f:
            # Skip comment and header lines
            if line.startswith('#') or line.strip() == '':
                continue
            
            # Split the line and strip whitespace
            parts = line.split()
            
            # Check if we have at least 3 parts
            if len(parts) >= 3:
                wavelengths.append(float(parts[0]))
                fluxes.append(float(parts[1]))
                errors.append(float(parts[2]))
                
                # Add survey if available, otherwise empty string
                surveys.append(parts[4] if len(parts) > 4 else '')
    
    # Convert to numpy arrays
    return {
        'wavelength': np.array(wavelengths),
        'flux': np.array(fluxes),
        'error': np.array(errors),
        'survey': np.array(surveys, dtype='U20')
    }

# Read Andicam data
def read_andicam_file(filepath):
    wavelengths = []
    fluxes = []
    errors = []
    bands = []
    
    with open(filepath, 'r') as f:
        for line in f:
            # Skip header and empty lines
            if line.startswith('#') or line.strip() == '' or 'JD' in line:
                continue
            
            parts = line.split()
            
            # Ensure we have enough columns and flux is not empty
            if len(parts) >= 5 and parts[4] != '':
                band = parts[1]
                flux = float(parts[4])
                # Use flux error if available, otherwise assume 10% error
                error = float(parts[3]) * flux if len(parts) > 3 and parts[3] != '-' else flux * 0.1
                
                wavelengths.append(ANDICAM_WAVELENGTHS[band])
                fluxes.append(flux)
                errors.append(error)
                bands.append(band)
    
    return {
        'wavelength': np.array(wavelengths),
        'flux': np.array(fluxes),
        'error': np.array(errors),
        'band': np.array(bands, dtype='U10')
    }

# Read merged spectrum data (wavelength in microns, flux in Jy)
def read_merged_file(filepath):
    wavelengths = []
    fluxes = []
    
    with open(filepath, 'r') as f:
        for line in f:
            # Skip comment lines
            if line.startswith('#') or line.strip() == '':
                continue
      
            parts = line.split()
            
            # Check if we have at least 2 parts
            if len(parts) >= 2:
                wavelengths.append(float(parts[0]))
                fluxes.append(float(parts[1]))
    
    return {
        'wavelength': np.array(wavelengths),
        'flux': np.array(fluxes)
    }

# Read optical spectrum
def read_spectrum_file(filepath):
    wavelengths = []
    fluxes = []
    
    with open(filepath, 'r') as f:
        for line in f:
            # Skip comment lines
            if line.startswith('#'):
                continue
            
            # Split the line and strip whitespace
            parts = line.split()
            
            # Check if we have at least 2 parts
            if len(parts) >= 2:
                wavelengths.append(float(parts[0]))
                fluxes.append(float(parts[1]))
    
    return {
        'wavelength': np.array(wavelengths),
        'flux': np.array(fluxes)
    }

def apply_dred(wavelengths, fluxes, errors, err_fraction=0.05, ebv=1.0):
    """
    Vectorized dereddenning function
    
    Args:
        wavelengths (np.array): Wavelengths in microns
        fluxes (np.array): Flux values
        errors (np.array): Error values (or error fraction if None)
        err_fraction (float): Fractional error if errors not provided
    
    Returns:
        tuple: Dreddened fluxes and errors
    """
    # If errors not provided, use fractional error
    if errors is None:
        errors = fluxes * err_fraction
    
    # Vectorized dreddenning
    dred_results = np.array([dred(wave, flux, err, ebv) 
                              for wave, flux, err in zip(wavelengths, fluxes, errors)])
    
    return dred_results[:, 0], dred_results[:, 1]


def create_sed_model(wave_model, temperature=2700, bb_scaling=1.34e-25, 
                     power_law_index=-2, power_law_scaling=5.0e-15):
    """
    Create a Spectral Energy Distribution (SED) model combining 
    a blackbody radiation curve with a power-law component below a defined wavelength.
    
    INPUTS:
        wave_model (np.array): Wavelength grid in microns
        temperature (float): Blackbody temperature in Kelvin
        bb_scaling (float): Blackbody scaling factor
        power_law_index (float): Power-law spectral index
        power_law_scaling (float): Power-law scaling factor
    
    Returns:
        np.array: SED model flux values
    """
    # Physical constants
    h = 6.626e-34  # Planck's constant (J·s)
    c = 3e8        # Speed of light (m/s)
    k = 1.380e-23  # Boltzmann constant (J/K)
    
    # Convert wavelength to meters
    wavelength_m = wave_model * 1e-6
    
    # Blackbody radiation component
    numerator = 2 * h * c**2
    denominator = (wavelength_m**5) * (np.exp((h * c) / (wavelength_m * k * temperature)) - 1)
    blackbody_component = numerator / denominator * bb_scaling
    
    # Power-law component active only below 1 micron
    power_law_component = np.where(
        wave_model < 1.2, 
        power_law_scaling * wave_model**power_law_index, 
        0
    )
    
    # Combined model
    sed_model = blackbody_component + power_law_component
    
    return sed_model
    

# Read infrared data
ir_data = read_data_file(dir_data + data_file)

# Read Andicam data
andicam_data = read_andicam_file(dir_data + andicam_file)

# Read optical spectrum
optical_data = read_spectrum_file(dir_data + spec_file)

# Read merged spectrum (Spitzer)
merged_data = read_merged_file(dir_data + merged_file)

# Remove negative flux values and filter for wavelengths > 3500
wavelength_mask = (optical_data['wavelength'] > 3500) & (optical_data['flux'] > 0)
optical_wavelength = optical_data['wavelength'][wavelength_mask]
optical_flux = optical_data['flux'][wavelength_mask]

# The Spitzer spectrum is in Jy. Converting to mJy for consistency
merged_data['flux'] = merged_data['flux'] * 1000

# Apply Savitzky-Golay filter for smoothing
window_length = min(51, len(optical_flux) // 2 * 2 - 1)  # Ensure odd number
if window_length > 0:
    smoothed_flux = savgol_filter(optical_flux, window_length, 3)
else:
    smoothed_flux = optical_flux

# Convert wavelength to microns
ir_wavelength_microns = ir_data['wavelength']
optical_wavelength_microns = optical_wavelength * 1e-4

# Conversion function from mJy to ergs/cm²/Å
def mJy_to_ergs_per_cm2_per_angstrom(flux_mJy, wavelength_microns):
    return flux_mJy * 3e-8 / (wavelength_microns * 1e4)**2 
    

# Convert infrared flux to ergs/cm²/Å
ir_flux_ergs_per_cm2_per_angstrom = np.array([
    mJy_to_ergs_per_cm2_per_angstrom(flux, wavelength) 
    for flux, wavelength in zip(ir_data['flux'], ir_data['wavelength'])
])

# Convert infrared error to ergs/cm²/Å
ir_error_ergs_per_cm2_per_angstrom = np.array([
    mJy_to_ergs_per_cm2_per_angstrom(err, wavelength) 
    for err, wavelength in zip(ir_data['error'], ir_data['wavelength'])
])

# Convert Andicam flux to ergs/cm²/Å (assuming F_lambda input)
andicam_flux_ergs_per_cm2_per_angstrom = andicam_data['flux']
andicam_error_ergs_per_cm2_per_angstrom = andicam_data['error']

# Convert merged spectrum flux to ergs/cm²/Å
merged_flux_ergs_per_cm2_per_angstrom = np.array([
    mJy_to_ergs_per_cm2_per_angstrom(flux, wavelength) 
    for flux, wavelength in zip(merged_data['flux'], merged_data['wavelength'])
])

# Apply dereddenning to different datasets
smoothed_flux_dred, smoothed_flux_dred_err = apply_dred(
    optical_wavelength_microns, 
    smoothed_flux, 
    smoothed_flux * 0.05, ebv
)

ir_flux_dred, ir_flux_dred_err = apply_dred(
    ir_wavelength_microns, 
    ir_flux_ergs_per_cm2_per_angstrom, 
    ir_error_ergs_per_cm2_per_angstrom, ebv
)

andicam_flux_dred, andicam_flux_dred_err = apply_dred(
    andicam_data['wavelength'], 
    andicam_flux_ergs_per_cm2_per_angstrom, 
    andicam_error_ergs_per_cm2_per_angstrom, ebv
)

# Create SED model wavelength grid
wave_model = np.linspace(0.2, 50, 10000)

# Compute SED model
sed_model = create_sed_model(
    wave_model, 
    temperature,      # Blackbody temperature
    bb_scaling,   # Blackbody scaling factor
    power_law_index,    # Spectral index
    power_law_scaling  # Power-law scaling factor
)


plt.figure(figsize=(16, 10), dpi=300)

plt.rcParams['axes.linewidth'] = 2.5  # Bolder frame

# Plot the Spitzer spectrum as a black continuous line
plt.plot(merged_data['wavelength'], merged_flux_ergs_per_cm2_per_angstrom,  
         '-', 
         color='black', 
         alpha=0.7, 
         label='Spitzer/IRS',
         linewidth=3)

# Optical spectra
plt.plot(optical_wavelength_microns, smoothed_flux, 
         color='gray', 
         alpha=0.7, 
         label='CTIO',
         linewidth=3)

plt.plot(optical_wavelength_microns, smoothed_flux_dred, 
         ':', 
         color='gray', 
         alpha=0.7, 
         label='dreddened Spectra',
         linewidth=3)

# Define markers and colors for different surveys
survey_styles = {
    'GALEX': {'marker': 'x', 'markerfacecolor': 'none', 'color': 'blue', 'label': 'GALEX'},
    '2MASS': {'marker': 'o', 'markerfacecolor': 'none', 'color': 'red', 'label': '2MASS'},
    'WISE': {'marker': 'v', 'markerfacecolor': 'none', 'color': 'purple', 'label': 'WISE'},
    'AKARI': {'marker': 'X', 'markerfacecolor': 'none', 'color': 'cyan', 'label': 'AKARI'},
    'Andicam': {'marker': 's', 'markerfacecolor': 'none', 'color': 'green', 'label': 'Andicam'}
}

# Plot infrared surveys
for survey_type in ['GALEX', '2MASS', 'WISE', 'AKARI']:
    mask = ir_data['survey'] == survey_type
    if np.any(mask):
        plt.errorbar(
            ir_wavelength_microns[mask], 
            ir_flux_ergs_per_cm2_per_angstrom[mask], 
            yerr=ir_error_ergs_per_cm2_per_angstrom[mask], 
            fmt=survey_styles[survey_type]['marker'], 
            color=survey_styles[survey_type]['color'], 
            label=survey_styles[survey_type]['label'],
            capsize=5,
            markersize=12,  
            linewidth=2.5,   
            elinewidth=2.5,  
            markerfacecolor=survey_styles[survey_type]['markerfacecolor'],
            markeredgewidth=2.5
        )

# Plot Andicam data
plt.errorbar(
    andicam_data['wavelength'], 
    andicam_flux_ergs_per_cm2_per_angstrom, 
    yerr=andicam_error_ergs_per_cm2_per_angstrom, 
    fmt=survey_styles['Andicam']['marker'], 
    color=survey_styles['Andicam']['color'], 
    label=survey_styles['Andicam']['label'],
    capsize=5,
    markersize=12,  
    linewidth=2.5,   
    elinewidth=2.5,  
    markerfacecolor=survey_styles['Andicam']['markerfacecolor'],
    markeredgewidth=2.5
)

# Dreddened data plots with black empty markers
plt.errorbar(ir_wavelength_microns, ir_flux_dred, yerr=ir_flux_dred_err, 
             marker='o', linestyle='None', capsize=5, elinewidth=2.5, 
             markersize=12, markerfacecolor='None', markeredgecolor='k', 
             ecolor='k', label='dreddened fluxes')

plt.errorbar(andicam_data['wavelength'], andicam_flux_dred, 
             yerr=andicam_flux_dred_err, marker='o', linestyle='None', 
             capsize=5, elinewidth=2.5, markersize=12, 
             markerfacecolor='None', markeredgecolor='k', ecolor='k')
             
plt.plot(wave_model, sed_model, 
         color='red', 
         linestyle='--', 
         linewidth=3, 
         alpha=0.7, 
         label='SED Model')
plt.plot(wave_model, sed_model)


plt.xlabel('Wavelength (μm)', fontsize=24, fontweight='bold')
plt.ylabel('Flux (ergs/cm²/Å)', fontsize=24, fontweight='bold')
plt.title('V745 Sco', fontsize=30, fontweight='bold')

plt.xscale('log')
plt.yscale('log')

plt.tick_params(axis='both', which='major', labelsize=18, width=2.5, length=8)
plt.tick_params(axis='both', which='minor', labelsize=16, width=2, length=5)

# Place legend inside the plot with larger font
plt.legend(loc='best', fontsize=14, frameon=True, framealpha=0.7, edgecolor='black')

# Tight layout and save
plt.tight_layout()
plt.savefig('wavelength_flux_plot_with_surveys.png', dpi=300, bbox_inches='tight')
print("Plot saved as wavelength_flux_plot_with_surveys.png")

# Print some stats
print("\nData Summary:")
print("Infrared Surveys:")
for survey_type in ['GALEX', '2MASS', 'WISE', 'AKARI']:
    mask = ir_data['survey'] == survey_type
    if np.any(mask):
        print(f"{survey_type}: {np.sum(mask)} points")

print("\nAndicam Data:")
for band in np.unique(andicam_data['band']):
    mask = andicam_data['band'] == band
    print(f"{band} band: {np.sum(mask)} points")

print("\nMerged Spectrum:")
print(f"Points: {len(merged_data['wavelength'])}")
print(f"Wavelength range: {merged_data['wavelength'].min():.2f} - {merged_data['wavelength'].max():.2f} μm")

combined_wavelengths = np.concatenate([andicam_data['wavelength'], ir_wavelength_microns])
combined_dred_fluxes = np.concatenate([andicam_flux_dred, ir_flux_dred])

# Below is done for convenient export

df_export = pd.DataFrame({
    'Wavelength (μm)': combined_wavelengths,
    'Dreddened Flux (ergs/cm²/Å)': combined_dred_fluxes
})

export_filename = 'V745Sco_dreddened_fluxes.csv'
df_export.to_csv(export_filename, index=False)
print(f"Exported dreddened fluxes to {export_filename}")
print("\nData Export Preview:")
print(df_export.head())


