# -*- coding: utf-8 -*-
#
# Determines CO rotational excitation temperature and 12C/13C ratio
# using the method described in Brittain et al. (2005, ApJ,626, 283).
# Note that this method applies for optically thick emission lines. 
# If the lines are optically thick, they will lie on a curve on the
# output plot as opposed to a straight line. The inverse of the line 
# gradient gives the excitation temperature. Here we assume the same
# excitation temperature for the 12CO and 13CO lines. This is reasonable
# as all lines originate from low-lying levels, close to the ground
# rotational state. The vertical offset between the 12CO and 13CO lines
# on the plot gives the 12C/13C ratio.
#
# INPUTS
# file_factors (dict): file names containing the data
# 
# The input contains two files, one for 12CO data and one for 13CO data. 
# They files contain line data as well as the fluxes and their errors
# derived from observational data. The CO line list is 
# Goorvitch et al. (1994, ApJS, 95, 535).
# 
# A plot is output showing the data and the best-fitting line. Reference
# lines are shown corresponding to 12C/13C = 4 (equilibrium) and
# 12C/13C = 90 (solar). Here I refer to the y axis quantity as 'reduced
# flux' (see Brittain et al. for details). 

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import linregress

# fundamental constants
h = 6.62e-34 
kb = 1.38e-23  
c = 2.998e8   

# file names and 'factors' (the factors are needed to convert fluxes to W/m^2)
file_factors = {'12co_lines.dat': 1e-16, '13co_lines.dat': 1e-17}

# Map file names to names used on the plot legend
file_display_names = {
    '12co_lines.dat': r'$^{12}$CO',
    '13co_lines.dat': r'$^{13}$CO'
}

def calc_para(filename, factor):
    """
    Calculate upper energy level and 'reduced flux'.
    
    Inputs:
    filename (str): Data file name
    factor (float): Multiplicative factor for flux values in W/m^2
    
    Outputs:
    dict: Dictionary containing energy_k, log_rd, filename, and rotational quantum number j
    """
    try:
        df = pd.read_csv(filename, delimiter=r'\s+', names=['wav', 'A', 'J', 'energy', 'flux', 'error'])
        
        jrot = df['J'].values
        # Calculate upper energy level in temperature units (K)
        energy_k = (df['energy'].values * 1e2 * h * c + h * c / (1e4 / df['wav'].values * 1e-6)) / kb
        
        # Calculate 'reduced flux'
        if filename == '12co_lines.dat':
            # main isotopologue (statistical weight for P branch lines)
            rd = df['flux'].values * factor * (1e4 / df['wav'].values)**2 / ((2 * (jrot + 1) + 1) * df['A'].values) 
        else:
            # 13C isotopologue (statistical weight for R branch lines)
            rd = df['flux'].values * factor * (1e4 / df['wav'].values)**2 / ((2 * (jrot - 1) + 1) * df['A'].values)
            
        log_rd = np.log(rd)
        
        return {'energy_k': energy_k, 'log_rd': log_rd,'filename': filename, 'j': jrot}
    except Exception as e:
        print(f"calc_para: error processing file {filename}: {e}")
        return None

def calc_intercept_error(energy_k, log_rd, slope):
    """
    Calculate plot intercept and its error for plotted line
    
    Inputs:
    energy_k (array): Energy levels in K
    log_rd (array): Log of reduced flux values
    slope (float): Fixed slope value
    
    Outputs:
    tuple: (intercept, error) pair
    """
    intercept = np.mean(log_rd - slope * energy_k)
    residuals = log_rd - (slope * energy_k + intercept)
    n = len(energy_k)
    mean_x = np.mean(energy_k)
    sum_x_squared = np.sum((energy_k - mean_x)**2)
    error = np.sqrt(np.sum(residuals**2)/(n-2)) * np.sqrt(1/n + mean_x**2/sum_x_squared)
    
    return intercept, error

def process_data():
    """Process data files and calculate plot parameters"""
    # Calculate parameters for each file
    para = []
    for filename, factor in file_factors.items():
        result = calc_para(filename, factor)
        if result:
            para.append(result)
    
    if not para:
        raise ValueError("process_data: no data files processed")
    
    slopes = []
    std_errs = []
    
    # For each dataset, calculate the slope
    for d in para:
        # For 12CO, use only data points with jrot > 2 for the fit
        if d['filename'] == '12co_lines.dat':
            mask = d['j'] > 2
            # Check if we have enough data points for the fit
            if np.sum(mask) < 2:
                print(f"Warning: Not enough data points with jrot > 2 for {d['filename']}")
                continue
            # Perform linear regression only on filtered data
            slope, intercept, r_value, p_value, std_err = linregress(
                d['energy_k'][mask], d['log_rd'][mask]
            )
        else:
            # For 13CO, use all data points
            slope, intercept, r_value, p_value, std_err = linregress(
                d['energy_k'], d['log_rd']
            )
        
        slopes.append(slope)
        std_errs.append(std_err)
    
    # Calculate weighted average slope
    weighted_avg_slope = np.average(slopes, weights=1/np.array(std_errs)**2)
    weighted_avg_std_err = 1/np.sqrt(np.sum(1/np.array(std_errs)**2))
    
    # Calculate temperature and its error
    temperature = 1/weighted_avg_slope
    temperature_error = weighted_avg_std_err/weighted_avg_slope**2
    
    # Calculate intercepts for each dataset
    results = []
    for d in para:
        # For 12CO, jrot < 2 lines are optically thick and should be excluded here
        if d['filename'] == '12co_lines.dat':
            mask = d['j'] > 2
  #          if np.sum(mask) < 2:
  #              continue
            intercept, error = calc_intercept_error(d['energy_k'][mask], d['log_rd'][mask], weighted_avg_slope)
        else:
            intercept, error = calc_intercept_error(d['energy_k'], d['log_rd'], weighted_avg_slope)
        
        results.append({'filename': d['filename'], 'intercept': intercept, 'error': error})
    
    # Calculate flux ratio between the two datasets (the 12C/13C ratio)
    if len(results) >= 2:
        intercept_diff = results[0]['intercept'] - results[1]['intercept']
        diff_error = np.sqrt(results[0]['error']**2 + results[1]['error']**2)
        flux_ratio = np.exp(intercept_diff)
        flux_ratio_error = flux_ratio * diff_error
    else:
        flux_ratio = flux_ratio_error = None
    
    # Calculate plot reference lines for key 12C/13C ratios
    ref_line = np.exp(results[0]['intercept']) if results else 0
    i1 = ref_line / 4  # Equilibrium value
    i2 = ref_line / 90  # Solar value
    cc_equl = np.log(i1)
    cc_solar = np.log(i2)
    
    return {
        'para': para,
        'weighted_avg_slope': weighted_avg_slope,
        'weighted_avg_std_err': weighted_avg_std_err,
        'temperature': temperature,
        'temperature_error': temperature_error,
        'results': results,
        'flux_ratio': flux_ratio,
        'flux_ratio_error': flux_ratio_error,
        'cc_equl': cc_equl,
        'cc_solar': cc_solar
    }

def create_plot(data):
    """
    Create and save the plot of upper energy level
    against reduced flux
    
    Inputs:
    data (dict): Dictionary containing processed data
    """
    para = data['para']
    weighted_avg_slope = data['weighted_avg_slope']
    cc_equl = data['cc_equl']
    cc_solar = data['cc_solar']
    results = data['results']
    flux_ratio = data['flux_ratio']
    
    plt.figure(figsize=(10, 6))
    
    for i, d in enumerate(para):
        display_name = file_display_names[d['filename']] 
        plt.scatter(d['energy_k'], d['log_rd'], label=display_name, s=30, 
                    marker='o' if i == 0 else 's')  # Different markers for different datasets
    
    # Create a combined label for both fitted lines
    combined_label = f"$T_{{rot}}$ = {abs(data['temperature']):.0f}K, $^{{12}}$C/$^{{13}}$C = {flux_ratio:.0f}"
    
    # Plot fitted lines for each dataset
    colors = ['red', 'red', 'green', 'purple']  # Different colors for different datasets
    for i, result in enumerate(results):
        intercept = result['intercept']
        
        # Find the corresponding data in para
        d = next((item for item in para if item['filename'] == result['filename']), None)
        if d:
            x_line = np.linspace(min(d['energy_k']), max(d['energy_k']), 100)
            y_line = weighted_avg_slope * x_line + intercept
            
            if i == 0:
                plt.plot(x_line, y_line, '--', linewidth=3, label=combined_label, color=colors[i])
            else:
                plt.plot(x_line, y_line, '--', linewidth=3, color=colors[i])
    
    # Plot reference lines for different 12C/13C ratios
    all_energy_k = np.concatenate([d['energy_k'] for d in para])
    x_line = np.linspace(min(all_energy_k), max(all_energy_k), 100)
    
    # Plot equilibrium 12C/13C ratio line
    y_line = weighted_avg_slope * x_line + cc_equl
    plt.plot(x_line, y_line, '-.', linewidth=3, label=f"$T_{{rot}}$ = {abs(data['temperature']):.0f}K, $^{{12}}$C/$^{{13}}$C = 4", color='green')
    
    # Plot solar 12C/13C ratio line
    y_line = weighted_avg_slope * x_line + cc_solar
    plt.plot(x_line, y_line, ':', linewidth=3, label=f"$T_{{ rot}}$ = {abs(data['temperature']):.0f}K, $^{{12}}$C/$^{{13}}$C = 90", color='purple')
    
    plt.xlabel('Energy (K)', fontsize=14)
    plt.ylabel(r'$\ln(F~\lambda^2 / g_J~A)$', fontsize=14)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.legend(fontsize=10, loc='best')
    plt.tight_layout() 
    
    plt.savefig('lines.pdf')

    plt.show()

def print_results(data):
    """Print derived parameters"""
    print("="*50)
    print("RESULTS:")
    print("="*50)
    print(f"Temperature: {abs(data['temperature']):.0f} K")
    print(f"Temperature error: {data['temperature_error']:.0f} K")
    
    if data['flux_ratio'] is not None:
        print(f"12C/13C ratio: {data['flux_ratio']:.4e} ± {data['flux_ratio_error']:.4e}")
    
    print("="*50)

def main():
    """Main function"""
    try:
        data = process_data()
        print_results(data)
        create_plot(data)
    except Exception as e:
        print(f"Error: {e}")

if __name__ == "__main__":
    main()
