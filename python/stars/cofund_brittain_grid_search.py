import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import linregress
from scipy.optimize import curve_fit 
from matplotlib import cm
from matplotlib.colors import LogNorm
from mpl_toolkits.mplot3d import Axes3D

# *** Note that this is an update to CONFUND_BRITTAIN with alterations
# to the fitting procedure outlined below***

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
# The files contain line data as well as the fluxes and their errors
# derived from observational data. The CO line list is 
# Goorvitch et al. (1994, ApJS, 95, 535).
# 
# A plot is output showing the data and the best-fitting line. Reference
# lines are shown corresponding to 12C/13C = 4 (equilibrium) and
# 12C/13C = 90 (solar). Here I refer to the y axis quantity as 'reduced
# flux' (see Brittain et al. for details). 

# In this version of the code the fitting is implemented as follows:

# 1) An initial guess of the excitation temperature is found from the slope
# of the 12CO lines (in the optically thin regime).

# 2) A 'coarse grid' in the parameter space is defined  and minimum chisq is 
# returned.

# 3) A finer grid is then defined around that minimum, and chisq
# is calculated at each position, as with the coarse grid. 

# 4) The parameter space around the chisq minimum in this finer grid is then 
# fitted with a quadratic. The best fit is found from the follwoing condition:
# dchisq/ dx = 0, 
# where x is a fitted parameter, and the error is given by 
# d^2chisq/dx^2. 

# Unlike the previous version of this code, the chisq topology is plotted 
# in the fine grid, with the best-fit marked by a star


# Fundamental constants
h = 6.62e-34 
kb = 1.38e-23  
c = 2.998e8   

# Calibration factors for different data files
file_factors = {'12co_lines.dat': 1e-16, '13co_lines.dat': 1e-17}

# Mapping of file names to display names for the legend
file_display_names = {
    '12co_lines.dat': r'$^{12}$CO',
    '13co_lines.dat': r'$^{13}$CO'
}

def calc_para(filename, factor):
    """
    Calculate parameters from CO line data.
    
    INPUTS:
    filename (str): Path to the data file
    factor (float): Calibration factor for flux values in W/m^2
    
    Returns:
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
            rd_err = df['error'].values / df['flux'].values * rd
        else:
            # 13C isotopologue (statistical weight for R branch lines)
            rd = df['flux'].values * factor * (1e4 / df['wav'].values)**2 / ((2 * (jrot - 1) + 1) * df['A'].values)
            rd_err = df['error'].values / df['flux'].values * rd
        log_rd = np.log(rd)
        log_rd_err = rd_err / rd
        
        return {'energy_k': energy_k, 'log_rd': log_rd, 'log_rd_err': log_rd_err, 'filename': filename, 'j': jrot}
    except Exception as e:
        print(f"Error processing file {filename}: {e}")
        return None
            
def  calc_para_error(fine_grid_results_df):
     """Takes a slice through the fine grid for quadratic fit"""
     fine_min_chisq_row = fine_grid_results_df.loc[fine_grid_results_df['reduced_chisq'].idxmin()]
        
     best_cc_ratio = fine_min_chisq_row['cc_ratio']
     temp_slice = fine_grid_results_df[fine_grid_results_df['cc_ratio'] == best_cc_ratio].sort_values('temp')
     best_temp = fine_min_chisq_row['temp']
     cc_ratio_slice = fine_grid_results_df[fine_grid_results_df['temp'] == best_temp].sort_values('cc_ratio')

     chisq_temp_values = temp_slice['reduced_chisq'].values
     chisq_cc_ratio_values = cc_ratio_slice['reduced_chisq'].values
     cc_ratio_values = cc_ratio_slice['cc_ratio'].values
        
     chisq_temp_values = temp_slice['reduced_chisq'].values
     temp_values = temp_slice['temp'].values
     
     best_temp, temp_err, best_temp_chisq = quad_fit(temp_values, chisq_temp_values)
     best_cc_ratio, cc_ratio_err, best_cc_ratio_chisq = quad_fit(cc_ratio_values, chisq_cc_ratio_values)
     best_chisq = (best_temp_chisq + best_cc_ratio_chisq) / 2
     
     quad_result_dict = {'min reduced chisq': best_chisq, 'best temp': best_temp, 'best cc ratio': best_cc_ratio, 'temp error': temp_err, 'cc ratio error': cc_ratio_err}

     return quad_result_dict
     
def quad_fit(para, chisq):
     """Fits a quadratic to the parameter space close to the minimum chisq returned by the fine grid search.
        The parameter error is given by the second order differential in the fitted function. For a function
        chisq = a + bx + cx^2, this is equal to 2c^-0.5. The minimum in chisq is returned from the condition
        dchisq/dx = b + 2x = 0. """
        
     best_index = np.argmin(chisq)
     min_index = best_index - 3
     max_index = best_index + 3
     para_sample = para[min_index:max_index]
     chisq_sample = chisq[min_index:max_index]
     
     quad_model = np.poly1d(np.polyfit(para_sample, chisq_sample, 2))
     best_para = (quad_model[1] / (2 * quad_model[2]))
     para_err = np.sqrt(1 / quad_model[2])
     reduced_chisq_min = quad_model[2] * best_para**2  - quad_model[1] * best_para + quad_model[0]
     
     return best_para, para_err, reduced_chisq_min
     
def create_plot(para, quad_result_dict, filtered_energy_k, filtered_log_rd, filtered_log_rd_err, lines_iso_data):
    """
    Creates plot of upper energy level against reduced flux. The best fitting model is shown, along
    with reference lines at solar and equilibrium 12C/13C ratios.
    
    Inputs:
    -----------
    para (list): List of dictionaries containing processed data
    quad_result_dict (dict): Dictionary with quadratic fit results
    filtered_energy_k (array): Filtered energy values for 12CO
    filtered_log_rd (array): Filtered log reduced flux values for 12CO
    filtered_log_rd_err (array): Filtered errors for 12CO
    lines_iso_data (dict): Data for 13CO
    """
    # Reference ratios
    cc_equl = 4    # Equilibrium 12C/13C ratio
    cc_solar = 90  # Solar 12C/13C ratio
    
    # Models for both isotopologues
    y_model_12co = filtered_energy_k / quad_result_dict['best temp']
    y_model_13co = lines_iso_data['energy_k'] / quad_result_dict['best temp'] - np.log(abs(quad_result_dict['best cc ratio']))
    
    # Determine vertical offset for the best-fitting model returned by the quad search
    comb_energy_k = np.concatenate([filtered_energy_k, lines_iso_data['energy_k']])
    comb_y_model = np.concatenate([y_model_12co, y_model_13co])
    comb_y_data = np.concatenate([filtered_log_rd, lines_iso_data['log_rd']])       
    cdiff = np.mean(comb_y_model) - np.mean(comb_y_data)
    
    plt.figure(figsize=(10, 6))
    
    best_fit_label = f"$T_{{rot}}$ = {abs(quad_result_dict['best temp']):.0f}K, $^{{12}}$C/$^{{13}}$C = {abs(quad_result_dict['best cc ratio']):.0f}"
    
    for i, item in enumerate(para):
    
        marker = 'o' if i == 0 else 's' 
        plt.scatter(item['energy_k'], item['log_rd'], 
                    label=file_display_names[item['filename']], 
                    s=30, marker=marker)
        
        x_line = np.linspace(min(comb_energy_k), max(comb_energy_k), 100)
        
    # Plot best-fit model lines for the whole energy range
    x_line = np.linspace(min(comb_energy_k), max(comb_energy_k), 100)
    # 12CO 
    y_12co_line = x_line / quad_result_dict['best temp'] - cdiff
    plt.plot(x_line, y_12co_line, '--', linewidth=3, label=best_fit_label, color='red')
    
    # 13CO 
    y_13co_line = x_line / quad_result_dict['best temp'] - cdiff - np.log(abs(quad_result_dict['best cc ratio']))
    plt.plot(x_line, y_13co_line, '--', linewidth=3, color='red')
    
    # Reference ratios
    y_equl_line = x_line / quad_result_dict['best temp'] - cdiff - np.log(cc_equl)
    plt.plot(x_line, y_equl_line, '-.', linewidth=3, 
             label=f"$T_{{rot}}$ = {abs(quad_result_dict['best temp']):.0f}K, $^{{12}}$C/$^{{13}}$C = {cc_equl}", 
             color='green')
    
    y_solar_line = x_line / quad_result_dict['best temp'] - cdiff - np.log(cc_solar)
    plt.plot(x_line, y_solar_line, '-.', linewidth=3, 
             label=f"$T_{{rot}}$ = {abs(quad_result_dict['best temp']):.0f}K, $^{{12}}$C/$^{{13}}$C = {cc_solar}", 
             color='purple')
    
    plt.xlabel('Energy Level (K)', fontsize=14)
    plt.ylabel(r'$\ln(F~\lambda^2 / g_J~A)$', fontsize=14)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.legend(fontsize=10, loc='best')
    plt.tight_layout() 
    
    plt.savefig('co_fund_fit.pdf')
    plt.show()
            
def process_data():
    """Process data files and calculate parameters"""
    # Calculate parameters for each file
    para = []
    for filename, factor in file_factors.items():
        result_unfiltered = calc_para(filename, factor)
        if result_unfiltered:
            para.append(result_unfiltered)
    
    if not para:
        raise ValueError("No valid data files processed")
	
    for d in para:
        if d['filename'] == '12co_lines.dat':
            lines_data = d
        else:
            lines_iso_data = d
            break
	    
    # Determine the slope from the 12^CO data for an initial guess of temperature
    # The mask here is necessary as low J 12CO lines are optically thick and should
    # be excluded from the model
    mask = lines_data['j'] > 2
    filtered_energy_k = lines_data['energy_k'][mask]
    filtered_log_rd = lines_data['log_rd'][mask]
    filtered_log_rd_err = lines_data['log_rd_err'][mask]
    
    # Calculate linear regression on filtered data
    slope, intercept, r_value, p_value, std_err = linregress(filtered_energy_k, filtered_log_rd)
        
    temp_init = round(abs(1 / slope), 0) # initial guess based on slope
    temp_init_err = abs(1 / (slope + std_err) - (1 / slope))
        
    # define temperature grid 
    temp_grid_factor = 1  # Sampling defined in multiples of error bar
    temp_grid_range = 3  # Sample range as a multiple of error bar
    temp_grid_step = temp_grid_factor * temp_init_err

    # Generate grid points with non-zero/negative check
    temp_min = max(1, temp_init - temp_grid_range * temp_init_err)
    temp_max = temp_init + temp_grid_range * temp_init_err
    temp_grid = np.arange(temp_min, temp_max + temp_grid_step/2, temp_grid_step)
    
    # define 12C/12C ratio grid
    cc_ratio_min = 4
    cc_ratio_max = 90
    cc_ratio_grid_step = 1
    cc_ratio_grid = np.arange(cc_ratio_min, cc_ratio_max + cc_ratio_grid_step/2, cc_ratio_grid_step)
    
    grid_results = []    
    
    chisq_grid = np.zeros((len(temp_grid), len(cc_ratio_grid)))
    
    # Degrees of freedom (number of data points minus number of fitted parameters)
    n_params = 2
    # Get total number of data points from both isotopologues
    n_data_points = len(filtered_energy_k) + len(lines_iso_data['energy_k'])
    # Degrees of freedom
    dof = n_data_points - n_params
    
    for i, temp in enumerate(temp_grid):
        grad = -1 / temp
        for j, cc_ratio in enumerate(cc_ratio_grid):
            y_model = {}
            y_model['12co_lines.dat'] = grad * filtered_energy_k
            y_model['13co_lines.dat'] = grad * lines_iso_data['energy_k'] - np.log(cc_ratio)
            comb_energy_k = np.concatenate([filtered_energy_k, lines_iso_data['energy_k']])
            comb_y_model = np.concatenate([y_model['12co_lines.dat'], y_model['13co_lines.dat']])
            comb_y_data = np.concatenate([filtered_log_rd, lines_iso_data['log_rd']])
            comb_y_err = np.concatenate([filtered_log_rd_err, lines_iso_data['log_rd_err']])
            cdiff = np.mean(comb_y_model) - np.mean(comb_y_data)

            scaled_model = comb_y_model - cdiff
            chisq = sum(((comb_y_data - scaled_model) / comb_y_err)**2)
            
            # Calculate reduced chi-square
            reduced_chisq = chisq / dof

            # Store results in a dictionary and append to results list
            grid_result_dict = {
                'chisq': chisq,
                'reduced_chisq': reduced_chisq,
                'temp': temp,
                'cc_ratio': cc_ratio
            }
            grid_results.append(grid_result_dict)
            
            # Store reduced chisq in grid for contour plot
            chisq_grid[i, j] = reduced_chisq
            
    # Convert results list to DataFrame for easier analysis
    grid_results_df = pd.DataFrame(grid_results)
    
    # Return results and grid data for plotting
    return grid_results_df, temp_grid, cc_ratio_grid, chisq_grid, filtered_energy_k, filtered_log_rd, filtered_log_rd_err, lines_iso_data, dof, para

def create_fine_grid(grid_results_df, temp_grid, cc_ratio_grid, filtered_energy_k, filtered_log_rd, filtered_log_rd_err, lines_iso_data, dof):
    """
    Create a fine grid around the minimum reduced chi-square value from the coarse grid.
    
    INPUTS:
    -----------
    grid_results_df (DataFrame): Results from coarse grid search
    temp_grid (array): Array of temperature values from coarse grid
    cc_ratio_grid (array): Array of 12C/13C ratio values from coarse grid
    filtered_energy_k, filtered_log_rd, filtered_log_rd_err: Data for model calculation
    lines_iso_data: Isotope data for model calculation
    dof (int): Degrees of freedom for reduced chi-square calculation
    
    Returns:
    --------
    fine_grid_results_df (DataFrame): Results from fine grid search
    fine_temp_grid (array): Array of temperature values for fine grid
    fine_cc_ratio_grid (array): Array of 12C/13C ratio values for fine grid
    fine_chisq_grid (array): 2D array of reduced chi-square values for fine grid
    """
    # Find the minimum reduced chi-square value and corresponding parameters
    min_chisq_row = grid_results_df.loc[grid_results_df['reduced_chisq'].idxmin()]
    best_temp = min_chisq_row['temp']
    best_cc_ratio = min_chisq_row['cc_ratio']
    
    # Get temperature step size from coarse grid
    temp_grid_step = temp_grid[1] - temp_grid[0] if len(temp_grid) > 1 else 1.0
    
    # Define fine grid around best values
    # Extend to two coarse grid steps in temperature
    fine_temp_min = best_temp - 1 * temp_grid_step
    fine_temp_max = best_temp + 1 * temp_grid_step
    fine_temp_step = temp_grid_step / 5  # 5x finer resolution
    
    # Extend ±10 in cc_ratio
    fine_cc_ratio_min = max(1, best_cc_ratio - 6)  # Ensure ratio > 0
    fine_cc_ratio_max = best_cc_ratio + 6
    fine_cc_ratio_step = 0.2  # 5x finer than coarse grid step of 1
    
    # Create fine grids
    fine_temp_grid = np.arange(fine_temp_min, fine_temp_max + fine_temp_step/2, fine_temp_step)
    fine_cc_ratio_grid = np.arange(fine_cc_ratio_min, fine_cc_ratio_max + fine_cc_ratio_step/2, fine_cc_ratio_step)
    
    print(f"Fine grid: {len(fine_temp_grid)} temperature points × {len(fine_cc_ratio_grid)} ratio points")
    
    # Create arrays to store fine grid reduced chisq values
    fine_chisq_grid = np.zeros((len(fine_temp_grid), len(fine_cc_ratio_grid)))
    fine_grid_results = []
    
    # Calculate reduced chi-square for fine grid
    for i, temp in enumerate(fine_temp_grid):
        grad = -1 / temp
        for j, cc_ratio in enumerate(fine_cc_ratio_grid):
            y_model = {}
            y_model['12co_lines.dat'] = grad * filtered_energy_k
            y_model['13co_lines.dat'] = grad * lines_iso_data['energy_k'] - np.log(cc_ratio)
            comb_energy_k = np.concatenate([filtered_energy_k, lines_iso_data['energy_k']])
            comb_y_model = np.concatenate([y_model['12co_lines.dat'], y_model['13co_lines.dat']])
            comb_y_data = np.concatenate([filtered_log_rd, lines_iso_data['log_rd']])
            comb_y_err = np.concatenate([filtered_log_rd_err, lines_iso_data['log_rd_err']])
            cdiff = np.mean(comb_y_model) - np.mean(comb_y_data)

            scaled_model = comb_y_model - cdiff
            chisq = sum(((comb_y_data - scaled_model) / comb_y_err)**2)
            
            # Calculate reduced chi-square
            reduced_chisq = chisq / dof
            
            # Store results
            fine_grid_result_dict = {
                'chisq': chisq,
                'reduced_chisq': reduced_chisq,
                'temp': temp,
                'cc_ratio': cc_ratio
            }
            fine_grid_results.append(fine_grid_result_dict)
            
            # Store reduced chisq in grid for plotting
            fine_chisq_grid[i, j] = reduced_chisq
    
    # Convert results list to DataFrame
    fine_grid_results_df = pd.DataFrame(fine_grid_results)
    
    return fine_grid_results_df, fine_temp_grid, fine_cc_ratio_grid, fine_chisq_grid

def plot_fine_chisq_topology(fine_temp_grid, fine_cc_ratio_grid, fine_chisq_grid, best_temp, best_cc_ratio):
    """
    Create contour and 3D plots of reduced chi-square values in the fine grid
    
    INPUTS:
    -----------
    fine_temp_grid (array): Array of temperature values for fine grid
    fine_cc_ratio_grid (array): Array of 12C/13C ratio values for fine grid
    fine_chisq_grid (array): 2D array of reduced chi-square values for fine grid
    best_temp (float): Best-fit temperature from fine grid
    best_cc_ratio (float): Best-fit 12C/13C ratio from fine grid
    """
    # Create a meshgrid for plotting
    X, Y = np.meshgrid(fine_temp_grid, fine_cc_ratio_grid, indexing='ij')
    
    # Set up a figure with two subplots
    fig = plt.figure(figsize=(18, 8))
    
    # 1. Create 2D contour plot
    ax1 = fig.add_subplot(121)
    
    # Create contour plot
    contour = ax1.contourf(X, Y, fine_chisq_grid, levels=20, cmap='viridis')
    
    cbar1 = plt.colorbar(contour, ax=ax1)
    cbar1.set_label(r'Reduced $\chi^2$', fontsize=14)
    
    # Add contour lines
    CS = ax1.contour(X, Y, fine_chisq_grid, colors='white', alpha=0.5, linewidths=0.5)
    ax1.clabel(CS, inline=True, fontsize=8, fmt='%.2f')
    
    # Mark the minimum reduced chi-square point
    ax1.plot(best_temp, best_cc_ratio, 'r*', markersize=15, label='Best fit')
    
    ax1.set_xlabel('Temperature (K)', fontsize=14)
    ax1.set_ylabel('$^{12}$C/$^{13}$C', fontsize=14)
    #ax1.set_title(r'Fine Grid Reduced $\chi^2$ Contour Plot', fontsize=16)
    ax1.legend(loc='upper right')
    ax1.grid(True, linestyle='--', alpha=0.7)
    
    # 2. Create 3D surface plot
    ax2 = fig.add_subplot(122, projection='3d')
    
    # Create the surface plot
    surf = ax2.plot_surface(X, Y, fine_chisq_grid, cmap='viridis', 
                           linewidth=0, antialiased=True, alpha=0.8)
    
    cbar2 = fig.colorbar(surf, ax=ax2, shrink=0.6, aspect=10)
    cbar2.set_label(r'Reduced $\chi^2$', fontsize=14)
    
    # Mark the minimum point in 3D
    best_idx = (np.abs(fine_temp_grid - best_temp)).argmin()
    best_jdx = (np.abs(fine_cc_ratio_grid - best_cc_ratio)).argmin()
    best_chisq = fine_chisq_grid[best_idx, best_jdx]
    
    # Plot the best fit point
    ax2.scatter([best_temp], [best_cc_ratio], [best_chisq], 
               color='red', s=200, marker='*', label='Best fit')
    
    ax2.set_xlabel('Temperature (K)', fontsize=14)
    ax2.set_ylabel('$^{12}$C/$^{13}$C', fontsize=14)
    ax2.set_zlabel(r'Reduced $\chi^2$', fontsize=14)
    #ax2.set_title(r'Fine Grid 3D Reduced $\chi^2$ Surface', fontsize=16)
    
    # Set view angle for better visualization
    ax2.view_init(elev=30, azim=45)
    
    plt.tight_layout()
    plt.savefig('fine_grid_reduced_chisq_topology.png', dpi=300)
    print("Fine grid reduced chi-square topology plots saved as 'fine_grid_reduced_chisq_topology.png'")
    
    # Display the plot
    plt.show()

def main():
    """Main function to run the analysis"""
    try:
        # Process data and get coarse grid results
        grid_results_df, temp_grid, cc_ratio_grid, chisq_grid, filtered_energy_k, filtered_log_rd, filtered_log_rd_err, lines_iso_data, dof, para = process_data()
        
        # Print out summary of coarse grid results
        print("\nCoarse Grid Summary:")
        print(f"Total number of coarse grid points: {len(grid_results_df)}")
        print(f"Degrees of freedom: {dof}")
        
        # Find the minimum reduced chi-square value and corresponding parameters in coarse grid
        coarse_min_chisq_row = grid_results_df.loc[grid_results_df['reduced_chisq'].idxmin()]
        print("\nCoarse grid best fit parameters:")

        print(f"Temperature: {coarse_min_chisq_row['temp']:.2f} K")
        print(f"12C/13C ratio: {coarse_min_chisq_row['cc_ratio']:.2f}")
        print(f"Chi-square: {coarse_min_chisq_row['chisq']:.4f}")
        print(f"Reduced chi-square: {coarse_min_chisq_row['reduced_chisq']:.4f}")
        
        # Create fine grid around coarse minimum
        fine_grid_results_df, fine_temp_grid, fine_cc_ratio_grid, fine_chisq_grid = create_fine_grid(
            grid_results_df, temp_grid, cc_ratio_grid, filtered_energy_k, filtered_log_rd, filtered_log_rd_err, lines_iso_data, dof)

        # Find the minimum reduced chi-square value and corresponding parameters in fine grid
        fine_min_chisq_row = fine_grid_results_df.loc[fine_grid_results_df['reduced_chisq'].idxmin()]
        print("\nFine Grid Summary:")
        print(f"Total number of fine grid points: {len(fine_grid_results_df)}")
        print("\nFine grid best fit parameters:")
        print(f"Temperature: {fine_min_chisq_row['temp']:.2f} K")
        print(f"12C/13C ratio: {fine_min_chisq_row['cc_ratio']:.2f}")
        print(f"Chi-square: {fine_min_chisq_row['chisq']:.4f}")
        print(f"Reduced chi-square: {fine_min_chisq_row['reduced_chisq']:.4f}")
        
        # Fit a quadratic to the parameter space around the chisq minimum returned above. Return updated values and errors
        quad_result_dict = calc_para_error(fine_grid_results_df)
        print("\nQuadratic Fit Summary:")
        print(f"Temperature: {round(abs(quad_result_dict['best temp']))} +/- {round(quad_result_dict['temp error'])} K")
        print(f"12C/13C ratio: : {round(abs(quad_result_dict['best cc ratio']))} +/- {round(quad_result_dict['cc ratio error'])}")
        print(f"Reduced chi-square: {round(quad_result_dict['min reduced chisq'])}")
        
        # Plot fine grid reduced chi-square topology as both contour and 3D plot
        plot_fine_chisq_topology(
            fine_temp_grid, fine_cc_ratio_grid, fine_chisq_grid, 
            fine_min_chisq_row['temp'], fine_min_chisq_row['cc_ratio'])
            
        create_plot(para, quad_result_dict, filtered_energy_k, filtered_log_rd, filtered_log_rd_err, lines_iso_data)
        
    except Exception as e:
        import traceback
        print(f"Error in analysis: {e}")
        print(traceback.format_exc())

if __name__ == "__main__":
    main()
