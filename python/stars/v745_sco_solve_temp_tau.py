import os
import re
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from scipy import interpolate
from calc_chisq import calc_chisq  

# Finds the best-fitting Dusty model to the Spitzer data of V745 Sco from a 
# grid of models covering the dust temperature and opacity (tau) paramter space.
#
# The program takes the 'stellar subtracted' Spitzer spectra and searches through 
# a grid of models to find the one giving the lowest reduced chisq. The models file
# names contain the parameters temp and tau which are extracted and stored alongide
# the value of chisq. 
#
# The Dusty models contain model flux as flux*wavelength (flux in wavelength units)
# which are converted to flux in frequency units for comparison with the Spitzer 
# spectra. The CALC_CHISQ function takes this flux and scales it to the flux in the 
# Spitzer spectra (in Jy), returning chisq.
#
# INPUTS: 
# MODEL_DIRECTORY (str) : A directory containing all of the model files in the form
#     'model_TEMPK_TAU.s001' where TEMP and TAU are the values describing each model
# DATA_FILEPATH (str) : Data file name. Note that this is the 'stellar subtracted'
#     Spitzer data created by V745_SCO_DATA_PREP (i.e. dust emission only)
#
# OUTPUTS:
# CHISQ_RESULTS.TXT : Output file containing TEMP, TAU, REDUCED_CHISQ
# BEST_FIT_MODEL.PNG : Plot showing the data and best-fitting dusty model

def extract_model_params(filename):
    """Extract temp and tau from filename"""
    match = re.search(r'model_(\d+)k_(\d+\.\d+)', filename)
    if match:
        temp = int(match.group(1))
        tau = float(match.group(2))
        return temp, tau
    return None, None

def process_model_file(filepath):
    """Return wavelength and flux from model file.
       Note the Dusty outputs flux as flux*lambda,
       where flux is in wavelength units. Here we
       convert to flux in Jy"""
    # Skip first 5 lines in the model files (header)
    with open(filepath, 'r') as f:
        for _ in range(5):
            next(f)
        
        wavelengths = []
        fluxes_wl = []
        
        for line in f:
            values = line.strip().split()
            if len(values) >= 2:  # at least 2 columns
                wavelength = float(values[0])
                flux_wl = float(values[1])
                wavelengths.append(wavelength)
                fluxes_wl.append(flux_wl)
    
    wavelengths = np.array(wavelengths)
    fluxes_wl = np.array(fluxes_wl)
    flux = fluxes_wl * wavelengths  
    
    return wavelengths, fluxes_wl, flux

def read_data_file(data_filepath, error_col=None, min_wavelength=9.0):
    """
    Read data file with wavelength, flux, and error
    If error_col is None, errors are estimated as sqrt(flux)
    Only wavelengths greater than min_wavelength considered for
    fitting
    """
    wavelengths = []
    fluxes = []
    errors = []
    
    with open(data_filepath, 'r') as f:
        # Skip header
        next(f)
        
        for line in f:
            values = line.strip().split()
            if len(values) >= 2:
                wavelength = float(values[0])
                
                # Consider wavelengths > min_wavelength
                if wavelength > min_wavelength:
                    flux = float(values[1])
                    wavelengths.append(wavelength)
                    fluxes.append(flux)
                    
                    # Read error column if supplied; otherwise use sqrt(flux)
                    if error_col is not None and len(values) > error_col:
                        error = float(values[error_col])
                    else:
                        error = np.sqrt(abs(flux)) if flux > 0 else 1.0
                    errors.append(error)
    
    return np.array(wavelengths), np.array(fluxes), np.array(errors)

def interpolate_model_to_data(data_wavelength, model_wavelength, model_flux):
    """Regrid model flux onto data wavelength grid"""
    # Create interpolation function from model
    f = interpolate.interp1d(model_wavelength, model_flux, kind='linear', 
                            bounds_error=False, fill_value=np.nan)
    
    # Apply interpolation to get model flux on the data wavelength grid
    interpolated_flux = f(data_wavelength)
    
    # Filter NaN values
    valid_mask = ~np.isnan(interpolated_flux)
    
    return valid_mask, interpolated_flux

def process_all_models(model_directory, data_filepath, error_col=None, min_wavelength=9.0):
    """Process model  and calculate chisq against data"""
    results = {}
    
    # Read and filter (by wavelength) data
    data_wavelength, data_flux, data_error = read_data_file(data_filepath, error_col, min_wavelength)
    
    # Check for data points after filtering
    if len(data_wavelength) == 0:
        print(f"ERROR: No data points found with wavelength > {min_wavelength} microns")
        return {}, None, None, None, None
    
    
    print(f"Read data file with {len(data_wavelength)} data points (wavelength > {min_wavelength} microns)")
    
    # Get model files
    model_files = [f for f in os.listdir(model_directory) if f.startswith('model_') and '.s' in f]
    print(f"Found {len(model_files)} model files")
    
    # Sort model files
    model_files.sort()
    
    for model_file in model_files:
        filepath = os.path.join(model_directory, model_file)
        temp, tau = extract_model_params(model_file)
        
        if temp is not None and tau is not None:
            try:
                # Process model file
                model_wavelength, model_flux_wl, model_flux = process_model_file(filepath)
                
                # Filter model data to include only wavelengths greater than min_wavelength
                wl_mask = model_wavelength > min_wavelength
                filtered_model_wavelength = model_wavelength[wl_mask]
                filtered_model_flux_wl = model_flux_wl[wl_mask]
                filtered_model_flux = model_flux[wl_mask]
                
                # Check if any model points remain after filtering
                if len(filtered_model_wavelength) == 0:
                    print(f"Skipping {model_file}: No model points with wavelength > {min_wavelength} microns")
                    continue
                
                # Interpolate model flux to the data wavelength grid
                valid_mask, interp_model_flux = interpolate_model_to_data(
                    data_wavelength, filtered_model_wavelength, filtered_model_flux)
                
                # Use valid data for chisq calculation
                valid_data_wavelength = data_wavelength[valid_mask]
                valid_data_flux = data_flux[valid_mask]
                valid_data_error = data_error[valid_mask]
                valid_model_flux = interp_model_flux[valid_mask]
                
                # Calculate chisq for the model
                chisq, scaled_model_flux = calc_chisq(valid_data_flux, valid_data_error, valid_model_flux)
                
                # Degrees of freedom (data points - 2 free parameters)
                dof = len(valid_data_flux) - 2
                reduced_chisq = chisq / dof if dof > 0 else np.inf
                
                # Gather results
                key = (temp, tau)
                results[key] = {
                    'model_wavelength': filtered_model_wavelength,
                    'model_flux_wl': filtered_model_flux_wl,
                    'model_flux': filtered_model_flux,
                    'chisq': chisq,
                    'reduced_chisq': reduced_chisq,
                    'dof': dof,
                    'scaled_model_flux': scaled_model_flux,
                    'valid_data_wavelength': valid_data_wavelength,
                    'valid_data_flux': valid_data_flux,
                    'filename': model_file
                }
                
                print(f"Processed: {model_file} (T={temp}K, tau={tau}): χ²={chisq:.2f}, χ²/dof={reduced_chisq:.2f}")
                
            except Exception as e:
                print(f"Error processing {model_file}: {e}")
    
    # Read the full data file for plotting 
    full_data_wavelength, full_data_flux, _ = read_data_file(data_filepath, error_col, min_wavelength=0.0)
    
    return results, data_wavelength, data_flux, data_error, full_data_wavelength, full_data_flux

def find_best_fit_model(results):
    """Find the model (minimum reduced chisq)"""
    if not results:
        return None, None
    
    best_key = min(results.keys(), key=lambda k: results[k]['reduced_chisq'])
    return best_key, results[best_key]

def plot_best_fit(results, data_wavelength, data_flux, data_error, 
                 full_data_wavelength, full_data_flux, 
                 min_wavelength=9.0, output_dir='results'):
    """Plot the best fit model and data"""
    # Create output directory if it doesn't exist
    Path(output_dir).mkdir(exist_ok=True)
    
    # Find best fit model
    best_key, best_model = find_best_fit_model(results)
    
    if best_key is None:
        print("No valid models found for plotting")
        return None
    
    temp, tau = best_key
    
    # Set global font properties for matplotlib
    plt.rcParams.update({
        'font.size': 16,
        'font.weight': 'bold',
        'axes.titlesize': 20,
        'axes.labelsize': 18,
        'axes.titleweight': 'bold',
        'axes.labelweight': 'bold',
        'xtick.labelsize': 16,
        'ytick.labelsize': 16,
        'legend.fontsize': 16,
        'axes.linewidth': 2.0,
        'lines.linewidth': 3.0,
        'lines.markersize': 10,
        'xtick.major.width': 2.0,
        'ytick.major.width': 2.0,
        'xtick.major.size': 7,
        'ytick.major.size': 7,
        'xtick.minor.size': 4,
        'ytick.minor.size': 4,
        'xtick.minor.width': 1.5,
        'ytick.minor.width': 1.5,
    })
    
    plt.figure(figsize=(12, 10))
    
    # Plot all data points
    plt.plot(full_data_wavelength, full_data_flux, 'ko', alpha=0.3, label='All data', markersize=6)
    
    f = interpolate.interp1d(best_model['model_wavelength'], best_model['model_flux'], 
                           kind='linear', bounds_error=False, fill_value=np.nan)
    interp_model_flux = f(data_wavelength)
    
    plt.plot(best_model['valid_data_wavelength'], best_model['scaled_model_flux'], 
             'r-', linewidth=4, label=f'Dusty model: T={temp}K, τ={tau}')
    
    plt.title(f'V745 Sco')
    plt.xlabel('Wavelength (μm)')
    plt.ylabel('Flux (Jy)')
    plt.legend(frameon=True, fancybox=True, shadow=True)
    
    plt.grid(True, linestyle='--', alpha=0.7, linewidth=1.5)
    plt.tight_layout()
    
    plt.savefig(os.path.join(output_dir, 'best_fit_model.png'), dpi=300)
    plt.close()
    
    return best_key


def save_results_to_file(results, min_wavelength, output_dir='results'):
    """Save chisq results to file"""
    Path(output_dir).mkdir(exist_ok=True)
    
    output_file = os.path.join(output_dir, f'chisq_results_wl_gt_{min_wavelength}um.txt')
    
    with open(output_file, 'w') as f:
        f.write(f"# Chi-square results for wavelengths > {min_wavelength} microns\n")
        f.write("# Temperature(K)  Tau  Chi-square  Reduced_Chi-square  DOF  Filename\n")
        
        # Sort results by reduced chisq
        sorted_results = sorted(results.items(), key=lambda x: x[1]['reduced_chisq'])
        
        for (temp, tau), data in sorted_results:
            f.write(f"{temp:<12}  {tau:<6.3f}  {data['chisq']:<10.2f}  {data['reduced_chisq']:<18.2f}  {data['dof']:<4}  {data['filename']}\n")
    
    print(f"Results saved to {output_file}")


if __name__ == "__main__":
    # Directory of model files
    model_directory = "sed_v2_models/"
    
    # Path to data file
    data_filepath = "v745_sco_stellar_sub_data.dat"  
    
    # Set the error column (0-indexed) if available in your data file
    # Set to None if you want to use sqrt(flux) as error estimate
    error_col = None  # Change if you have an error column (e.g., 2 for third column)
    
    # Minimum wavelength for fitting (in microns). This constraint is applied because
    # the data at the shortest wavelengths are affected by SiO emission
    min_wavelength = 9.0
    
    # Process and calculate chisq
    results, data_wavelength, data_flux, data_error, full_data_wavelength, full_data_flux = process_all_models(
        model_directory, data_filepath, error_col, min_wavelength)
    
    save_results_to_file(results, min_wavelength)
    
    best_key = plot_best_fit(
        results, data_wavelength, data_flux, data_error,
        full_data_wavelength, full_data_flux, min_wavelength)
    
    print(f"Processed {len(results)} model files")
    if best_key:
        print(f"Best fit model: T={best_key[0]}K, τ={best_key[1]}")
    else:
        print("No valid models found")
