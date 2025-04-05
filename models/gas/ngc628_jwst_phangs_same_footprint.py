from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
from astropy.wcs import WCS
from reproject import reproject_interp
from astropy.convolution import Gaussian2DKernel, convolve_fft
import matplotlib.font_manager as fm
from matplotlib.patches import Rectangle

# Takes the 21 micron JWST image of NGC628 and places it on the ALMA/PHANGS 
# CO(2-1) footprint. The JWST image is rotated and regridded before being
# convolved to the resolution of the ALMA data. 

# The field of view for the PHANGS dataset is then reduced to that of the JWST,
# allowing for a direct comparison of the emission in both images.

# INPUTS:
# DIR_JWST (str) : Directory containing the JWST data
# DIR_PHANGS (str) : Directory containing the PHANGS data
# IMAGE_FILE_JWST (str) : JWST file name
# IMAGE_FILE_PHANGS (str) : PHANGS file name
# GALAXY_DISTANCE_MPC (float) : Distance to NGC628 in Mpc
# JWST_FWHM (float) : Resolution of the JWST image (not present in header)

# OUTPUTS:
# NGC628_JWST_PHANGS.PNG : Plot showing the rotated/regridded/convolved JWST
#     image along side the cropped PHANGS image
# CROPPED_PHANGS.FITS and CONVOLVED_JWST.FITS : Fits files of the data on the
#     same footprint and with the same resolution

# Inputs
dir_jwst = '/home/mtr/Downloads/ngc628/jw02107/level_3/'
dir_phangs = '/home/mtr/Downloads/ngc628/'
image_file_jwst = 'jw02107-o039_t018_miri_f2100w_i2d.fits'
image_file_phangs = 'ngc0628_12m+7m+tp_co21_broad_mom0.fits'

# Galaxy parameters
galaxy_distance_mpc = 9.5  # Distance to NGC 628 in Mpc

# FWHM of PSFs
jwst_fwhm = 0.7  # JWST MIRI F2100W PSF FWHM (this is not found in the Header)

#--------------------------------------------------------------------------------

# Read JWST data and header
hdu_list = fits.open(dir_jwst + image_file_jwst)
image_data_jwst = hdu_list[1].data
image_header_jwst = hdu_list[1].header

# Read PHANGS data and header
hdu_list = fits.open(dir_phangs + image_file_phangs)
image_data_phangs = hdu_list[0].data
image_header_phangs = hdu_list[0].header

# Retrieve pixel scales and resolution from Headers (convert deg-to-arcsecs)
jwst_pixscale = np.abs(image_header_jwst['CDELT1']) * 3600
phangs_pixscale = np.abs(image_header_phangs['CDELT1']) * 3600
phangs_fwhm = np.abs(image_header_phangs['BMIN']) * 3600

# Create a new HDU with the JWST data and its original header
jwst_hdu = fits.PrimaryHDU(data=image_data_jwst, header=image_header_jwst)
jwst_hdul = fits.HDUList([jwst_hdu])

# Reproject JWST data to PHANGS coordinate system
reprojected_jwst, footprint = reproject_interp(jwst_hdul, image_header_phangs)

# Create a mask of valid (non-NaN and non-zero) pixels in the reprojected JWST data
jwst_valid_mask = (~np.isnan(reprojected_jwst)) & (reprojected_jwst != 0)

# Calculate the required sigma for the kernel (in PHANGS pixels)
target_fwhm = np.sqrt(phangs_fwhm**2 - jwst_fwhm**2) 
sigma_pix = target_fwhm / 2.355 / phangs_pixscale

# Create a Gaussian kernel for convolving JWST data
kernel = Gaussian2DKernel(x_stddev=sigma_pix)

# Find the bounding box of the valid JWST data
nonzero_indices = np.argwhere(jwst_valid_mask)
if len(nonzero_indices) > 0:  # Make sure there are valid pixels
    min_y, min_x = nonzero_indices.min(axis=0)
    max_y, max_x = nonzero_indices.max(axis=0)
    
    # Crop both datasets to this bounding box
    cropped_phangs = image_data_phangs[min_y:max_y+1, min_x:max_x+1].copy()
    cropped_jwst = reprojected_jwst[min_y:max_y+1, min_x:max_x+1].copy()
    
    # Create a combined mask for NaN values AND zeros in JWST data
    combined_mask = np.isnan(cropped_jwst) | (cropped_jwst == 0)
    
    # Set zeros in JWST data to NaN
    cropped_jwst[cropped_jwst == 0] = np.nan
    
    # Apply convolution to JWST data
    convolved_jwst = convolve_fft(cropped_jwst, kernel, fill_value=0, normalize_kernel=True, 
        allow_huge=True, nan_treatment='fill')
    
    # Apply the same mask to both datasets
    cropped_phangs[combined_mask] = np.nan
    convolved_jwst[combined_mask] = np.nan 
    
    # Update the PHANGS WCS information for the cropped image
    phangs_wcs = WCS(image_header_phangs)
    cropped_wcs = phangs_wcs[min_y:max_y+1, min_x:max_x+1]
    
    # Create a new header for the cropped PHANGS data
    cropped_header = cropped_wcs.to_header()
    
    # Calculate scale bar length in pixels
    # 1 kpc in arcseconds at galaxy distance
    kpc_in_arcsec = 1 * 1000 / galaxy_distance_mpc / 206265 * 1e6  # Convert Mpc to pc, then to arcsec
    scale_bar_pixels = kpc_in_arcsec / phangs_pixscale
    
    plt.figure(figsize=(14, 7))
    plt.gcf().set_facecolor("black")
    plt.subplots_adjust(wspace=0.05)  # Reduced space between subplots
    
    # Try to use a more visually pleasing font
    try:
        prop = fm.FontProperties(family='DejaVu Sans', weight='bold', size=14)
    except:
        prop = fm.FontProperties(weight='bold', size=14)
    
    # Plot JWST data
    ax1 = plt.subplot(1, 2, 1)
    plt.title('JWST', color='white', fontsize=18, fontweight='bold', fontproperties=prop)
    img1 = plt.imshow(convolved_jwst, cmap='gist_heat', vmin=273, vmax=290)
    plt.axis('off')
    
    # Plot PHANGS data
    ax2 = plt.subplot(1, 2, 2)
    plt.title('PHANGS', color='white', fontsize=18, fontweight='bold', fontproperties=prop)
    img2 = plt.imshow(cropped_phangs, cmap='gist_heat', vmin=0, vmax=30)
    plt.axis('off')
    
    # Get the WCS to determine the correct North direction
    # For PHANGS data - we use the second image for the arrow
    w = phangs_wcs
    
    # Extract the direction of North from the WCS
    # North direction should be determined by the CD/CDELT matrix
    # For simplicity, we'll determine it from pixel coordinates
    # This is a simplified approach and might need refinement
    pixel_y = np.array([0, 1])
    pixel_x = np.array([0, 0])
    
    # Convert pixel to world coordinates
    world = w.wcs_pix2world(np.column_stack([pixel_x, pixel_y]), 0)
    
    # Calculate the direction of North in pixel coordinates
    delta_ra = world[1, 0] - world[0, 0]
    delta_dec = world[1, 1] - world[0, 1]
    
    # Account for cos(dec) factor in RA
    delta_ra *= np.cos(np.radians(world[0, 1])) 
    
    # Calculate angle of North direction in display coordinates
    # We need to flip the dx because RA increases to the left
    north_angle = np.degrees(np.arctan2(-delta_ra, delta_dec))
    
    # Add North arrow (second plot only)
    arrow_props = dict(facecolor='white', edgecolor='white', width=2, headwidth=8, headlength=10)
    
    # Calculate arrow vector components based on the north_angle
    dx = np.sin(np.radians(north_angle)) * 0.15
    dy = np.cos(np.radians(north_angle)) * 0.15
    
    # North arrow 
    ax2.annotate('N', xy=(0.95, 0.15), xycoords='axes fraction', 
               xytext=(0.95 - dx, 0.15 - dy), textcoords='axes fraction',
               color='white', fontsize=16, fontweight='bold', ha='center',
               arrowprops=arrow_props)
    
    # Add scale bar 
    bar_x = 0.75
    bar_y = 0.05
    bar_length = 0.2  # proportion of axes width
    ax2.add_patch(Rectangle((bar_x, bar_y), bar_length, 0.01, 
                          transform=ax2.transAxes, color='white'))
    
    ax2.text(bar_x + bar_length/2, bar_y + 0.03, '1 kpc', 
           transform=ax2.transAxes, color='white', 
           fontsize=14, fontweight='bold', ha='center')
    
    plt.tight_layout()
    plt.savefig('ngc628_jwst_phangs.png', dpi=300)
    plt.show()
    
    fits.writeto('cropped_phangs.fits', cropped_phangs, cropped_header, overwrite=True)
    fits.writeto('convolved_jwst.fits', convolved_jwst, cropped_header, overwrite=True)
else:
    print("No valid JWST data found after reprojection")
