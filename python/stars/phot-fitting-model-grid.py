from dred import dred
import numpy as np
from read_nextgen_model import read_model
from convolve_data_filter import conv_model
from calc_chisq import calc_chisq
import matplotlib.pyplot as plt
import pandas as pd
import os

# Convolves next_gen models downloaded from
# http://svo2.cab.inta-csic.es/theory/newov2/
# with Johnson generic filters and compares to 
# observed dereddened fluxes from photometric
# observations.
# The model fluxes are scaled to the data and
# chisq is derived.
# The best-fitting model is found from the 
# minimum chisq and a 'chisq heatmap' is 
# generated for the parameter space

#
# ***At the moment this considers stellar temp
# and log(g) for solar metalicity***

# INPUTS
#--------
# FLUX_TARGET (ARRAY): Observed fluxes
# FILTER_TARGET (ARRAY): Filters used for FLUX_TARGET
# EBV (FLOAT): Colour excess of target due to interstellar reddening
#
# OUTPUTS
#---------
# CHISQ (FLOAT): The value of chisq for the scaled fluxes and 
#    dereddened FLUX_TARGET


def process_model_grid(temp_grid, logg_grid, filter_target, flux_corr_values, flux_corr_err_values, model_dir):
    """
    Process a grid of NextGen models and calculate chi-square values for each.
    
    INPUTS:
    -----------
    temp_grid : list
        List of temperature values (as strings)
    logg_grid : list
        List of log(g) values (as strings)
    filter_target : list
        List of filter names
    flux_corr_values : array
        Dereddened observed fluxes
    flux_corr_err_values : array
        Errors on dereddened fluxes
    model_dir : str
        Directory containing the NextGen models
        
    Returns:
    --------
    pd.DataFrame
        DataFrame containing temperature, log(g), and chisq
    """
    results = []
    
    # Loop through temperature and log(g) in our grid
    for temp in temp_grid:
        for logg in logg_grid:
            try:
 
                model_file = f'lte{temp}-{logg}-0.0.BT-NextGen.7.dat.txt'
                model_path = os.path.join(model_dir, model_file)
                
                wave_model, flux_model = read_model(model_path)
                
                # Calculate integrated fluxes for each filter
                model_fluxes = []
                for band in filter_target:
                    wav_eff, flux_int_model, mag_model = conv_model(wave_model, flux_model, band)
                    model_fluxes.append(flux_int_model)
                
                # Convert to numpy array
                model_fluxes = np.array(model_fluxes)
                
                # Calculate chisq
                chisq, flux_scaled = calc_chisq(flux_corr_values, flux_corr_err_values, model_fluxes)
                
                results.append({
                    'temperature': temp,
                    'logg': logg,
                    'chisq': chisq,
                    'scaled_fluxes': flux_scaled
                })
                
                print(f"Processed model: Temp={temp}, log(g)={logg}, Chi-Square={chisq:.4f}")
                
            except Exception as e:
                print(f"Error processing model: Temp={temp}, log(g)={logg}")
                print(f"Error message: {str(e)}")
    
    
    results_df = pd.DataFrame(results)
    return results_df

def plot_chisq_heatmap(results_df):
    """
    Create a heatmap of chisq for the model grid.
    
    INPUTS:
    -----------
    results_df : pd.DataFrame
        DataFrame with temperature, log(g), and chi-square values
    """
    # Pivot the data to create a 2D grid suitable for heatmap
    pivot_df = results_df.pivot(index='logg', columns='temperature', values='chisq')
    
    plt.figure(figsize=(12, 8))
    plt.imshow(pivot_df.values, cmap='viridis_r')
    plt.colorbar(label='Chi-Square')
    
    plt.xticks(np.arange(len(pivot_df.columns)), pivot_df.columns, rotation=45)
    plt.yticks(np.arange(len(pivot_df.index)), pivot_df.index)
    
    plt.xlabel('Temperature')
    plt.ylabel('log(g)')
    plt.title('Chi-Square Values for Model Grid')
    
    # Add chisq
    for i in range(len(pivot_df.index)):
        for j in range(len(pivot_df.columns)):
            plt.text(j, i, f"{pivot_df.iloc[i, j]:.2f}", 
                     ha="center", va="center", color="white", fontsize=8)
    
    plt.tight_layout()
    return plt.gcf()

def plot_best_fit(best_model, filter_target, flux_corr_values, wavelengths):
    """
    Plot the best-fitting model against the observed data.
    
    INPUTS:
    -----------
    best_model : dict
        Dictionary containing information about the best-fitting model
    filter_target: list
        List of filter names
    flux_corr_values : array
        Dereddened observed fluxes
    wavelengths : array
        Effective wavelengths of filters
    """
    plt.figure(figsize=(10, 6))
    
    # Plot observed data
    plt.scatter(wavelengths, flux_corr_values, label='Observed (dereddened)', color='blue', s=80, marker='o')
    
    # Plot best model
    plt.scatter(wavelengths, best_model['scaled_fluxes'], label=f"Best model (T={best_model['temperature']}, log(g)={best_model['logg']})", 
                color='red', s=80, marker='x')
    
    plt.xlabel('Wavelength (microns)')
    plt.ylabel('Flux')
    plt.title('Best Fitting Model vs. Observations')
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    # Add filter labels
    for i, filt in enumerate(filter_target):
        plt.annotate(filt, (wavelengths[i], flux_corr_values[i]), 
                    xytext=(5, 5), textcoords='offset points')
    
    plt.tight_layout()
    return plt.gcf()

# Main execution code
if __name__ == "__main__":
    # Define the model grid
    temp_grid = ['026', '027', '028', '029', '030', '031', '032', '033', '034', '035']
    logg_grid = ['0.0', '1.0', '1.5', '2.0']
    
    # Directory containing next gen models
    dir_model = '/home/mtr/Downloads/models/nextgen/'
    
    # Observed fluxes
    flux_target = np.array([3.614e-16, 1.50e-15, 5.984e-15, 2.993e-14])
    # Errors in FLUX_TARGET
    flux_target_err = 0.05 * flux_target  # Without access to errors I assume 5%
    # Filters of observation
    filter_target = ['B', 'V', 'R', 'I']
    # Color excess of the target (i.e. interstellar reddening)
    ebv = 0.32
    
    # Calculate effective wavelengths for each filter
    # We need to read one model to calculate the effective wavelengths
    ref_model_file = f'lte031-0.0-0.0.BT-NextGen.7.dat.txt'
    ref_model_path = os.path.join(dir_model, ref_model_file)
    wave_model, flux_model = read_model(ref_model_path)
    
    # Calculate effective wavelengths for each filter
    wavelengths = []
    for band in filter_target:
        wav_eff, _, _ = conv_model(wave_model, flux_model, band)
        wavelengths.append(wav_eff)
    wavelengths = np.array(wavelengths)
    
    # Create a dictionary to store data for each filter
    data = {}
    for filt, wav, flux, flux_err in zip(filter_target, wavelengths, flux_target, flux_target_err):
        data[filt] = {'wavelength': wav, 'flux': flux, 'flux_err': flux_err}
    
    # Deredden the observed fluxes
    for filt, values in data.items():
        flux_corr, flux_corr_err = dred(values['wavelength'], values['flux'], values['flux_err'], ebv)
        data[filt]['flux_corr'] = flux_corr
        data[filt]['flux_corr_err'] = flux_corr_err
        
        print(f"Filter: {filt}")
        print(f"Wavelength: {values['wavelength']:.3f} microns")
        print(f"Observed Flux: {values['flux']:.2e}")
        print(f"Dereddened Flux: {flux_corr:.2e}")
        print(f"Dereddened Flux Error: {flux_corr_err:.2e}")
    
    # Extract corrected fluxes and errors
    flux_corr_values = np.array([data[filt]['flux_corr'] for filt in filter_target])
    flux_corr_err_values = np.array([data[filt]['flux_corr_err'] for filt in filter_target])
    
    # Process the model grid
    results_df = process_model_grid(temp_grid, logg_grid, filter_target, 
                                   flux_corr_values, flux_corr_err_values, dir_model)
    
    # Find the best fitting model (minimum chi-square)
    best_model = results_df.loc[results_df['chisq'].idxmin()]
    print("\nBest fitting model:")
    print(f"Temperature: {best_model['temperature']}")
    print(f"log(g): {best_model['logg']}")
    print(f"Chi-Square: {best_model['chisq']:.4f}")
    
    # Save results to CSV
    results_df.to_csv('model_grid_results.csv', index=False)
    print("\nResults saved to 'model_grid_results.csv'")
    
    # Create and save plots
    chisq_heatmap = plot_chisq_heatmap(results_df)
    chisq_heatmap.savefig('chisq_heatmap.png', dpi=300, bbox_inches='tight')
    print("Chi-square heatmap saved to 'chisq_heatmap.png'")
    
    # Plot best fit
    best_fit_plot = plot_best_fit(best_model, filter_target, flux_corr_values, wavelengths)
    best_fit_plot.savefig('best_fit_plot.png', dpi=300, bbox_inches='tight')
    print("Best-fit plot saved to 'best_fit_plot.png'")
