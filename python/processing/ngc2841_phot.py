from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
from cog import cog
import os

def process_galaxy_image(config):
    """
    Determines the 'sky' background for a galaxy fits image and calculates galaxy
    flux.
    
    The background is estimated by defining RMAX - the radius beyond which there 
    is no galaxy emission. RMAX is best estimated from the 'leveling' of the 
    surface brightness profile, where the galaxy counts drops below the background
    level. Some judgement maybe required to ensure the best estimate of the 
    background. Examine the *sky* images output to DIR_OUTPUT.
    
    The surface brightness profile is derived by the COG routine which divides
    the fits image into annuli (width WA) centred on XCEN and YCEN. NBACKGROUND
    is the number of radii over which the background is determined: the mean
    counts in the RMAX + NBACKGROUND annuli is the background method used here 
    and the std error in the RMAX + NBACKGROUND annuli is the background error.
    """
    # Extract configuration parameters
    dir_data = config['dir_data']
    dir_mask = config['dir_mask']
    file_data = config['file_data']
    file_mask = config['file_mask']
    pa = config['pa']
    wa = config['wa']
    xcen = config['xcen']
    ycen = config['ycen']
    ellipt = config['ellipt']
    gal_dist_mpc = config['gal_dist_mpc']
    pixel_scale = config['pixel_scale']
    fluxconv = config['fluxconv']
    calib_error = config['calib_error']
    rmax = config['rmax']
    nbackground = config['nbackground']
    xlab = config['xlab']
    ylab = config['ylab']
    title = config['title']
    xmin = config['xmin']
    xmax = config['xmax']
    ymin = config['ymin']
    ymax = config['ymax']
    gal_name = config['gal_name']
    dir_output = config['dir_output']
    filter_label = config['filter_label']
    file_tag_bcog = config['file_tag_bcog']
    file_tag_scog = config['file_tag_scog']
    
    # Calculate derived values
    area_of_pixel = (pixel_scale / (3600 * 180) * 3.14159)**2 
    pix_to_kpc = pixel_scale / 3600 * 3.14159 / 180
    
    # Ensure output directory exists
    os.makedirs(dir_output, exist_ok=True)
    
    # Load data and mask
    hdu_list = fits.open(os.path.join(dir_data, file_data))
    image_data = hdu_list[0].data
    hdu_list.close()
    
    hdu_list = fits.open(os.path.join(dir_mask, file_mask))
    image_data_mask = hdu_list[0].data
    hdu_list.close()
    
    # Run curve of growth analysis
    rcog_data, bcog_data, ncog, scog_data, r2pi = cog(
        image_data, image_data_mask, pa, ellipt, xcen, ycen, wa
    )
    
    # Convert to MJy/sr and r in kpc
    bcog_data_funits = fluxconv * bcog_data
    scog_data_funits = fluxconv * scog_data
    rcog_kpc = rcog_data * pix_to_kpc * gal_dist_mpc * 1e3
    
    # Find the annulus closest to Rmax (max radius at which there is galaxy emission)
    closest_idx = np.argmin(np.abs(rcog_kpc - rmax))
    total_galaxy_pixels = sum(ncog[0:closest_idx])  # total number of pixels within Rmax
    
    # Get nbackground annuli starting from the closest one
    if closest_idx + nbackground <= len(rcog_kpc):
        selected_indices = range(closest_idx, closest_idx + nbackground)
    else:
        # If not enough annuli at larger radii, adjust range
        selected_indices = range(closest_idx, len(rcog_kpc))
        print(f"Warning: Only {len(selected_indices)} annuli available from radius {rcog_kpc[closest_idx]:.2f} kpc onwards")
    
    # Calculate stats for the background annuli 
    mean_background = np.mean(bcog_data[selected_indices])
    std_background = np.std(bcog_data[selected_indices])
    background_error = std_background / np.sqrt(len(selected_indices))
    
    # Convert these stats to MJy/sr
    mean_background_funits = np.mean(bcog_data_funits[selected_indices])
    std_background_funits = np.std(bcog_data_funits[selected_indices])
    background_error_funits = std_background_funits / np.sqrt(len(selected_indices))
    
    # OUTPUTS
    print(f"\n===== Results for {filter_label} =====")
    print("\nSummary Statistics for All Selected Annuli (raw data units):")
    print(f"  Mean background: {mean_background:.6f}")
    print(f"  Std background:   {std_background:.6f}")
    print(f"  Background error:   {background_error:.6f}")
     
    # Convert scog and background to Jy / pixel
    scog_data_calib = scog_data_funits * 1e6 * area_of_pixel
    background_level_jy = mean_background_funits * 1e6 * area_of_pixel
    background_error_jy = background_error_funits * 1e6 * area_of_pixel
    
    # Total background error in galaxy flux 
    total_background_error = background_error_jy * total_galaxy_pixels
    scog_calib_error = max(scog_data_calib) * calib_error
    
    # Add in quadrature to derive total error
    total_error = np.sqrt(scog_calib_error**2 + total_background_error**2)
    
    # Find the galaxy flux and subtract the background
    flux_gal = scog_data_calib[closest_idx] - sum(ncog[0:closest_idx])*background_level_jy
    
    # The total Galaxy flux is given by the peak in sky-subtracted scog
    print(f"  Galaxy flux: {flux_gal:.6f} +/- {total_error:.6f} Jy")
    
    # Create plots
    # 1. Surface brightness profile
    plt.figure(figsize=(10, 6))
    plt.plot(rcog_kpc, bcog_data, '-', alpha=0.7, label='Surface brightness')
    plt.axhline(y=mean_background, color='g', linestyle='--', 
                label=f'Mean of background annuli ({mean_background:.6f})')
    plt.axvline(x=rmax, color='r', linestyle='--')
    
    plt.xlim(xmin, xmax)
    plt.ylim(ymin, ymax)
    plt.xlabel(xlab)
    plt.ylabel(ylab)
    plt.legend()
    plt.title(title)
    plt.grid(True, alpha=0.3)
    output_path = os.path.join(dir_output, f"{gal_name}_{filter_label}{file_tag_bcog}")
    plt.savefig(output_path)
    print(f"  Saved brightness profile to: {output_path}")
    plt.close()
    
    # 2. Curve of growth
    plt.figure()
    plt.plot(rcog_kpc, scog_data_calib, '-', alpha=0.7)
    plt.xlabel(xlab)
    plt.ylabel(ylab)
    plt.title(title)
    output_path = os.path.join(dir_output, f"{gal_name}_{filter_label}{file_tag_scog}")
    plt.savefig(output_path)
    print(f"  Saved curve of growth to: {output_path}")
    plt.close()
    
    return {
        'filter': filter_label,
        'mean_background': mean_background,
        'std_background': std_background,
        'background_error': background_error,
        'galaxy_flux': flux_gal,
        'flux_error': total_error
    }

def main():
    """
    This script processes multiple astronomical images to determine the sky background
    and galaxy flux for NGC2841 in different wavelength bands. Values are output to the
    screen.
    """
    # Define directories
    dir_data = '/home/mtr/Downloads/n2841/data/'
    dir_mask = '/home/mtr/Downloads/n2841/masks/'
    dir_output = '../figures/'
    gal_name = 'NGC2841'
    gal_dist_mpc = 14  # galaxy distance in Mpc
    
    # File tags
    file_tag_bcog = '_bcog_sky.png'
    file_tag_scog = '_scog.png'
    
    # Plot settings
    xlab = 'r (kpc)'
    ylab = 'Brightness (counts)'
    
    # Configuration for each image
    configs = [
        # GALEX FUV
        {
            'dir_data': dir_data,
            'dir_mask': dir_mask,
            'file_data': 'NGC_2841-I-FUV-g2006.fits',
            'file_mask': 'NGC2841_GALEX_FUV_mask.fits',
            'pa': -25.00,
            'wa': 1,
            'xcen': 243.50,
            'ycen': 242.50,
            'ellipt': 0.56,
            'gal_dist_mpc': gal_dist_mpc,
            'pixel_scale': 1.5,
            'fluxconv': 1.0724e-4 * 1.0e-6 / ((1.5 / 3600 * 3.14159 / 180.)**2),
            'calib_error': 0.1,
            'rmax': 25,
            'nbackground': 6,
            'xlab': xlab,
            'ylab': ylab,
            'title': 'GALEX FUV',
            'xmin': 5,
            'xmax': 30,
            'ymin': 1e-4,
            'ymax': 0.002,
            'gal_name': gal_name,
            'dir_output': dir_output,
            'filter_label': 'GALEX_FUV',
            'file_tag_bcog': file_tag_bcog,
            'file_tag_scog': file_tag_scog,
        },
	# GALEX NUV
        {
            'dir_data': dir_data,
            'dir_mask': dir_mask,
            'file_data': 'NGC_2841-I-NUV-g2006.fits',
            'file_mask': 'NGC2841_GALEX_NUV_mask.fits',
            'pa': -29.67,
            'wa': 1,
            'xcen': 248.17,
            'ycen': 241.40,
            'ellipt': 0.42,
            'gal_dist_mpc': gal_dist_mpc,
            'pixel_scale': 1.5,
            'fluxconv': 3.529e-5 * 1.0e-6 / ((1.5 / 3600 * 3.14159 / 180.)**2),
            'calib_error': 0.1,
            'rmax': 25,
            'nbackground': 6,
            'xlab': xlab,
            'ylab': ylab,
            'title': 'GALEX NUV',
            'xmin': 5,
            'xmax': 30,
            'ymin': 1e-4,
            'ymax': 0.0125,
            'gal_name': gal_name,
            'dir_output': dir_output,
            'filter_label': 'GALEX_NUV',
            'file_tag_bcog': file_tag_bcog,
            'file_tag_scog': file_tag_scog,
        },
        # SDSS u
        {
            'dir_data': dir_data,
            'dir_mask': dir_mask,
            'file_data': 'frame-u-002137-3-0267.fits',
            'file_mask': 'NGC2841_SDSS_U_mask.fits',
            'pa': -83.34,
            'wa': 1,
            'xcen': 1711.84,
            'ycen': 1374.37,
            'ellipt': 0.43,
            'gal_dist_mpc': gal_dist_mpc,
            'pixel_scale': 0.396,
            'fluxconv': 3.631e-6 * 1.0e-6 / ((0.396 / 3600 * 3.14159 / 180.)**2),
            'calib_error': 0.01,
            'rmax': 19,
            'nbackground': 6,
            'xlab': xlab,
            'ylab': ylab,
            'title': 'SDSS u',
            'xmin': 5,
            'xmax': 25,
            'ymin': -0.01,
            'ymax': 0.04,
            'gal_name': gal_name,
            'dir_output': dir_output,
            'filter_label': 'SDSS_u',
            'file_tag_bcog': file_tag_bcog,
            'file_tag_scog': file_tag_scog,
        },
	# B
        {
            'dir_data': dir_data,
            'dir_mask': dir_mask,
            'file_data': 'NGC_2841-I-B-kab2003.fits',
            'file_mask': 'NGC2841_B_mask.fits',
            'pa': -29.27,
            'wa': 1,
            'xcen': 1072.40,
            'ycen': 1013.87,
            'ellipt': 0.50,
            'gal_dist_mpc': gal_dist_mpc,
            'pixel_scale': 0.305,
            'fluxconv': 2.702938878e-06 * 1.0e-6 / ((0.305 / 3600 * 3.14159 / 180.)**2),
            'calib_error': 0.01,
            'rmax': 19,
            'nbackground': 6,
            'xlab': xlab,
            'ylab': ylab,
            'title': 'B',
            'xmin': 5,
            'xmax': 25,
            'ymin': 0.4,
            'ymax': 0.6,
            'gal_name': gal_name,
            'dir_output': dir_output,
            'filter_label': 'B',
            'file_tag_bcog': file_tag_bcog,
            'file_tag_scog': file_tag_scog,
        },
	# SDSS g
        {
            'dir_data': dir_data,
            'dir_mask': dir_mask,
            'file_data': 'frame-g-002137-3-0267.fits',
            'file_mask': 'NGC2841_SDSS_G_mask.fits',
            'pa': -82.92,
            'wa': 1,
            'xcen': 1709.40,
            'ycen': 1379.52,
            'ellipt': 0.44,
            'gal_dist_mpc': gal_dist_mpc,
            'pixel_scale': 0.396,
            'fluxconv': 3.631e-6 * 1.0e-6 / ((0.396 / 3600 * 3.14159 / 180.)**2),
            'calib_error': 0.01,
            'rmax': 19,
            'nbackground': 6,
            'xlab': xlab,
            'ylab': ylab,
            'title': 'SDSS g',
            'xmin': 5,
            'xmax': 25,
            'ymin': -0.01,
            'ymax': 0.02,
            'gal_name': gal_name,
            'dir_output': dir_output,
            'filter_label': 'SDSS_g',
            'file_tag_bcog': file_tag_bcog,
            'file_tag_scog': file_tag_scog,
        },
	# V
        {
            'dir_data': dir_data,
            'dir_mask': dir_mask,
            'file_data': 'NGC_2841-I-V-kab2003.fits',
            'file_mask': 'NGC2841_V_mask.fits',
            'pa': -30.08,
            'wa': 1,
            'xcen': 1065.58,
            'ycen': 1026.42,
            'ellipt': 0.51,
            'gal_dist_mpc': gal_dist_mpc,
            'pixel_scale': 0.305,
            'fluxconv': 6.301862642e-07 * 1.0e-6 / ((0.305 / 3600 * 3.14159 / 180.)**2),
            'calib_error': 0.01,
            'rmax': 25,
            'nbackground': 6,
            'xlab': xlab,
            'ylab': ylab,
            'title': 'V',
            'xmin': 5,
            'xmax': 30,
            'ymin': 1.0,
            'ymax': 1.7,
            'gal_name': gal_name,
            'dir_output': dir_output,
            'filter_label': 'V',
            'file_tag_bcog': file_tag_bcog,
            'file_tag_scog': file_tag_scog,
        },
	# SDSS r
        {
            'dir_data': dir_data,
            'dir_mask': dir_mask,
            'file_data': 'frame-r-002137-3-0267.fits',
            'file_mask': 'NGC2841_SDSS_R_mask.fits',
            'pa': -83.73,
            'wa': 1,
            'xcen': 1715.59,
            'ycen': 1371.09,
            'ellipt': 0.44,
            'gal_dist_mpc': gal_dist_mpc,
            'pixel_scale': 0.396,
            'fluxconv': 3.631e-6 * 1.0e-6 / ((0.396 / 3600 * 3.14159 / 180.)**2),
            'calib_error': 0.01,
            'rmax': 32,
            'nbackground': 6,
            'xlab': xlab,
            'ylab': ylab,
            'title': 'SDSS r',
            'xmin': 5,
            'xmax': 40,
            'ymin': -0.005,
            'ymax': 0.01,
            'gal_name': gal_name,
            'dir_output': dir_output,
            'filter_label': 'SDSS_r',
            'file_tag_bcog': file_tag_bcog,
            'file_tag_scog': file_tag_scog,
        },
	# I
        {
            'dir_data': dir_data,
            'dir_mask': dir_mask,
            'file_data': 'NGC_2841-I-I-kab2003.fits',
            'file_mask': 'NGC2841_I_mask.fits',
            'pa': -29.70,
            'wa': 1,
            'xcen': 1069.41,
            'ycen': 1019.08,
            'ellipt': 0.51,
            'gal_dist_mpc': gal_dist_mpc,
            'pixel_scale': 0.305,
            'fluxconv': 1.056433766e-06 * 1.0e-6 / ((0.305 / 3600 * 3.14159 / 180.)**2),
            'calib_error': 0.01,
            'rmax': 25,
            'nbackground': 6,
            'xlab': xlab,
            'ylab': ylab,
            'title': 'I',
            'xmin': 5,
            'xmax': 30,
            'ymin': 3.25,
            'ymax': 4.25,
            'gal_name': gal_name,
            'dir_output': dir_output,
            'filter_label': 'I',
            'file_tag_bcog': file_tag_bcog,
            'file_tag_scog': file_tag_scog,
        },
	# SDSS i
        {
            'dir_data': dir_data,
            'dir_mask': dir_mask,
            'file_data': 'frame-i-002137-3-0267.fits',
            'file_mask': 'NGC2841_SDSS_I_mask.fits',
            'pa': -83.67,
            'wa': 1,
            'xcen': 1713.00,
            'ycen': 1373.87,
            'ellipt': 0.45,
            'gal_dist_mpc': gal_dist_mpc,
            'pixel_scale': 0.396,
            'fluxconv': 3.631e-6 * 1.0e-6 / ((0.396 / 3600 * 3.14159 / 180.)**2),
            'calib_error': 0.01,
            'rmax': 29,
            'nbackground': 6,
            'xlab': xlab,
            'ylab': ylab,
            'title': 'SDSS i',
            'xmin': 5,
            'xmax': 50,
            'ymin': -0.005,
            'ymax': 0.01,
            'gal_name': gal_name,
            'dir_output': dir_output,
            'filter_label': 'SDSS_i',
            'file_tag_bcog': file_tag_bcog,
            'file_tag_scog': file_tag_scog,
	},
	# SDSS z
        {
            'dir_data': dir_data,
            'dir_mask': dir_mask,
            'file_data': 'frame-z-002137-3-0267.fits',
            'file_mask': 'NGC2841_SDSS_Z_mask.fits',
            'pa': -83.37,
            'wa': 1,
            'xcen': 1710.06,
            'ycen': 1378.14,
            'ellipt': 0.43,
            'gal_dist_mpc': gal_dist_mpc,
            'pixel_scale': 0.396,
            'fluxconv': 3.631e-6 * 1.0e-6 / ((0.396 / 3600 * 3.14159 / 180.)**2),
            'calib_error': 0.01,
            'rmax': 20,
            'nbackground': 6,
            'xlab': xlab,
            'ylab': ylab,
            'title': 'SDSS z',
            'xmin': 5,
            'xmax': 25,
            'ymin': -0.005,
            'ymax': 0.04,
            'gal_name': gal_name,
            'dir_output': dir_output,
            'filter_label': 'SDSS_z',
            'file_tag_bcog': file_tag_bcog,
            'file_tag_scog': file_tag_scog,
        },
        # 2MASS J
        {
            'dir_data': dir_data,
            'dir_mask': dir_mask,
            'file_data': '2MASS_NGC_2841_J.fits',
            'file_mask': 'NGC2841_2MASS_J_mask.fits',
            'pa': -29.71,
            'wa': 1,
            'xcen': 600.48,
            'ycen': 600.60,
            'ellipt': 0.44,
            'gal_dist_mpc': gal_dist_mpc,
            'pixel_scale': 1.0,
            'fluxconv': 7.8285e-6 * 1.0e-6 / ((1.0 / 3600 * 3.14159 / 180.)**2),
            'calib_error': 0.011,
            'rmax': 20,
            'nbackground': 6,
            'xlab': xlab,
            'ylab': ylab,
            'title': '2MASS J',
            'xmin': 5,
            'xmax': 25,
            'ymin': -0.05,
            'ymax': 0.5,
            'gal_name': gal_name,
            'dir_output': dir_output,
            'filter_label': '2MASS_J',
            'file_tag_bcog': file_tag_bcog,
            'file_tag_scog': file_tag_scog,
        },
	# 2MASS H
        {
            'dir_data': dir_data,
            'dir_mask': dir_mask,
            'file_data': '2MASS_NGC_2841_H.fits',
            'file_mask': 'NGC2841_2MASS_H_mask.fits',
            'pa': -29.77,
            'wa': 1,
            'xcen': 600.50,
            'ycen': 600.56,
            'ellipt': 0.45,
            'gal_dist_mpc': gal_dist_mpc,
            'pixel_scale': 1.0,
            'fluxconv': 5.7056e-6 * 1.0e-6 / ((1.0 / 3600 * 3.14159 / 180.)**2),
            'calib_error': 0.007,
            'rmax': 20.7,
            'nbackground': 6,
            'xlab': xlab,
            'ylab': ylab,
            'title': '2MASS H',
            'xmin': 5,
            'xmax': 25,
            'ymin': -0.1,
            'ymax': 0.5,
            'gal_name': gal_name,
            'dir_output': dir_output,
            'filter_label': '2MASS_H',
            'file_tag_bcog': file_tag_bcog,
            'file_tag_scog': file_tag_scog,
        },
	# 2MASS K
        {
            'dir_data': dir_data,
            'dir_mask': dir_mask,
            'file_data': '2MASS_NGC_2841_K.fits',
            'file_mask': 'NGC2841_2MASS_K_mask.fits',
            'pa': -29.71,
            'wa': 1,
            'xcen': 600.48,
            'ycen': 600.60,
            'ellipt': 0.44,
            'gal_dist_mpc': gal_dist_mpc,
            'pixel_scale': 1.0,
            'fluxconv': 6.5267e-6 * 1.0e-6 / ((1.0 / 3600 * 3.14159 / 180.)**2),
            'calib_error': 0.007,
            'rmax': 19,
            'nbackground': 6,
            'xlab': xlab,
            'ylab': ylab,
            'title': '2MASS K',
            'xmin': 5,
            'xmax': 25,
            'ymin': -0.1,
            'ymax': 0.5,
            'gal_name': gal_name,
            'dir_output': dir_output,
            'filter_label': '2MASS_K',
            'file_tag_bcog': file_tag_bcog,
            'file_tag_scog': file_tag_scog,
        },
        # IRAC I1
        {
            'dir_data': dir_data,
            'dir_mask': dir_mask,
            'file_data': 'NGC_2841-I-IRAC_3.6-srh2010.fits',
            'file_mask': 'NGC2841_IRAC_I1_mask.fits',
            'pa': -30.98,
            'wa': 1,
            'xcen': 507.15,
            'ycen': 1061.10,
            'ellipt': 0.46,
            'gal_dist_mpc': gal_dist_mpc,
            'pixel_scale': 0.75,
            'fluxconv': 1.0,
            'calib_error': 0.02,
            'rmax': 19,
            'nbackground': 6,
            'xlab': xlab,
            'ylab': ylab,
            'title': 'IRAC I1',
            'xmin': 5,
            'xmax': 25,
            'ymin': 0,
            'ymax': 0.4,
            'gal_name': gal_name,
            'dir_output': dir_output,
            'filter_label': 'IRAC_I1',
            'file_tag_bcog': file_tag_bcog,
            'file_tag_scog': file_tag_scog,
        },
	# IRAC I2
        {
            'dir_data': dir_data,
            'dir_mask': dir_mask,
            'file_data': 'NGC_2841-I-IRAC_4.5-srh2010.fits',
            'file_mask': 'NGC2841_IRAC_I2_mask.fits',
            'pa': -31.33,
            'wa': 1,
            'xcen': 704.20,
            'ycen': 551.78,
            'ellipt': 0.45,
            'gal_dist_mpc': gal_dist_mpc,
            'pixel_scale': 0.75,
            'fluxconv': 1.0,
            'calib_error': 0.02,
            'rmax': 19,
            'nbackground': 6,
            'xlab': xlab,
            'ylab': ylab,
            'title': 'IRAC I2',
            'xmin': 5,
            'xmax': 25,
            'ymin': 0.125,
            'ymax': 0.2,
            'gal_name': gal_name,
            'dir_output': dir_output,
            'filter_label': 'IRAC_I2',
            'file_tag_bcog': file_tag_bcog,
            'file_tag_scog': file_tag_scog,
        },
	# IRAC I3
        {
            'dir_data': dir_data,
            'dir_mask': dir_mask,
            'file_data': 'NGC_2841-I-IRAC_5.8-kab2003.fits',
            'file_mask': 'NGC2841_IRAC_I3_mask.fits',
            'pa': -31.34,
            'wa': 1,
            'xcen': 508.02,
            'ycen': 1047.50,
            'ellipt': 0.46,
            'gal_dist_mpc': gal_dist_mpc,
            'pixel_scale': 0.75,
            'fluxconv': 1.0,
            'calib_error': 0.02,
            'rmax': 20,
            'nbackground': 6,
            'xlab': xlab,
            'ylab': ylab,
            'title': 'IRAC I3',
            'xmin': 5,
            'xmax': 25,
            'ymin': -0.05,
            'ymax': 0.3,
            'gal_name': gal_name,
            'dir_output': dir_output,
            'filter_label': 'IRAC_I3',
            'file_tag_bcog': file_tag_bcog,
            'file_tag_scog': file_tag_scog,
        },
	# IRAC I4
        {
            'dir_data': dir_data,
            'dir_mask': dir_mask,
            'file_data': 'NGC_2841-I-IRAC_8.0-kab2003.fits',
            'file_mask': 'NGC2841_IRAC_I4_mask.fits',
            'pa': -27.11,
            'wa': 1,
            'xcen': 707.43,
            'ycen': 552.94,
            'ellipt': 0.45,
            'gal_dist_mpc': gal_dist_mpc,
            'pixel_scale': 0.75,
            'fluxconv': 1.0,
            'calib_error': 0.02,
            'rmax': 20,
            'nbackground': 6,
            'xlab': xlab,
            'ylab': ylab,
            'title': 'IRAC I4',
            'xmin': 5,
            'xmax': 25,
            'ymin': -0.05,
            'ymax': 0.3,
            'gal_name': gal_name,
            'dir_output': dir_output,
            'filter_label': 'IRAC_I4',
            'file_tag_bcog': file_tag_bcog,
            'file_tag_scog': file_tag_scog,
        },
    ]
    
    # Process each image
    results = []
    for config in configs:
        print(f"\nProcessing {config['filter_label']}...")
        result = process_galaxy_image(config)
        results.append(result)
    
    # Summary report
    print("\n===== SUMMARY =====")
    print("Filter     | Background  | Error     | Galaxy Flux     | Flux Error")
    print("-----------------------------------------------------------------------")
    for r in results:
        print(f"{r['filter']:<10} | {r['mean_background']:.6f} | {r['background_error']:.6f} | {r['galaxy_flux']:.6f} | {r['flux_error']:.6f}")

if __name__ == "__main__":
    main()
