import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

# Prepares the V745 Sco Spitzer data for modelling with the Dusty code.

# In V745 Sco, the stellar continuum makes a significant contribution to the
# emission at the wavelengths covered by the Spitzer data. To isolate the dust 
# emission, which is prerequisite to fitting with Dusty, the model of the 
# V745 Sco SED determined in the V745_SCO_SED_PLOT program is used to predict 
# the stellar emission at the Spitzer wavelengths. This model is placed on the 
# same wavelength grid as the Spitzer data and converted to mJy. It is then 
# subtracted from the data to isolate the dust emission. 

# INPUTS:
# Spitzer data file (str) : V745Sco_merged.dat
# Temperature (float) : temperature of the BB component in the SED
# BB_scaling (float) : scaling factor for the BB component to match data in 
#     ergs/cm^2/A
# Power_law_index (float) : index of the power law describing the hot component 
#     in the SED
# power_law_scaling (float) : scaling factor for the got component to match data 
#     in ergs/cm^2/A

# OUTPUTS: 
# Output file (v745_sco_stellar_sub_data.dat) containing wavelength, dust 
# isolated data emission, data, and stellar model for the wavelengths covered by
#  Spitzer.

# Read Spitzer data
data = np.loadtxt('V745Sco_merged.dat')
wave_data = data[:, 0]  
flux_data = data[:, 1] 

# Parameters for the (approximate) model of the SED (see V745_SCO_SED_PLOT program) 
temperature = 2700
# I have reduced this scaling slightly to avoid negative values in the model subtracted
# data at 5 micron. Note we are assuming that the SED is non-variable which may be causing
# us to overestimate the stellar contribution at 5 micron
bb_scaling = 1.32e-25 
power_law_index = -2
power_law_scaling = 5.0e-15

# Create SED model wavelength grid 
wave_model_orig = np.linspace(0.2, 50, 10000)
   
# Physical constants
h = 6.626e-34  
c = 3e8       
k = 1.380e-23 

# Function to calculate SED of V745 Sco (see V745_SCO_SED_PLOT program)
def calculate_sed(wavelengths):

    wavelength_m = wavelengths * 1e-6
    
    # Blackbody (secondary) component
    numerator = 2 * h * c**2
    denominator = (wavelength_m**5) * (np.exp((h * c) / (wavelength_m * k * temperature)) - 1)
    blackbody_component = numerator / denominator * bb_scaling
    
    # Power-law component below 1.2 micron
    power_law_component = np.where(
        wavelengths < 1.2, 
        power_law_scaling * wavelengths**power_law_index, 
        0
    )
    
    # Combined model
    return blackbody_component + power_law_component

# Calculate model on original grid
sed_model_orig = calculate_sed(wave_model_orig)

# Calculate model directly on data wavelength grid
sed_model_data = calculate_sed(wave_data)

# Convert model to mJy
sed_model_jy = sed_model_data / (wave_data)**2  * 3e-8  * 1e23

# Subtract the 'stellar' continuum from the data to leave the dust continuum
flux_data_sub_stellar = flux_data - sed_model_jy

# Create plot
plt.figure(figsize=(12, 8))

# SED model
plt.subplot(2, 1, 1)
plt.plot(wave_model_orig, sed_model_orig, 'b-', label='Power law + BB')
plt.title('SED Model')
plt.xlabel('Wavelength (microns)')
plt.ylabel('Flux')
plt.legend()
plt.grid(True)

# Plot showing data, SED model, and 'stellar subtracted data' on the data 
# wavelength grid
plt.subplot(2, 1, 2)
plt.plot(wave_data, flux_data_sub_stellar, 'ro', label='Dust only (stellar subtracted)')
plt.plot(wave_data, sed_model_jy, 'g-', label='Model (on Data Grid)')
plt.plot(wave_data, flux_data, 'b-', label='Observed Data')
plt.title('Stellar Subtracted Data')
plt.xlabel('Wavelength (microns)')
plt.ylabel('Flux')
plt.legend()
plt.grid(True)

plt.tight_layout()
plt.savefig('v745_sco_stellar_sub_data.png')
plt.show()

output = np.column_stack((wave_data, flux_data_sub_stellar, flux_data, sed_model_jy))
np.savetxt('v745_sco_stellar_sub_data.dat', output, header='Wavelength(micron) Data_Flux_Model_Sub Data_Flux Model_Flux')
