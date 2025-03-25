import numpy as np
import matplotlib.pyplot as plt

def read_model(filename):
     """
     Reads Nextgen stellar models downloaded from the following link:
    http://svo2.cab.inta-csic.es/theory/newov2/
     """
    wavelength = []
    flux = []
    
    with open(filename, 'r') as file:
        for line in file:
            # Skip comment lines that start with #
            if line.startswith('#'):
                continue
           
            columns = line.strip().split()
            
            # Check if the line has enough data
            if len(columns) >= 2:
                try:
                    # Convert string values to float
                    wavelength.append(float(columns[0]))
                    flux.append(float(columns[1]))
                except ValueError:
                    # Skip lines that can't be converted to float
                    continue
    
    # Convert lists to numpy arrays 
    wavelength = np.array(wavelength)
    flux = np.array(flux)
    
    return wavelength, flux


