import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.patches import Rectangle
import matplotlib.lines as mlines
from astropy.io import fits
import numpy as np

# Takes the Fits file outputs from NGC28_JWST_PHANGS_SAME_FOOTPRINT and
# creates an overlay plot to highlight differences in spatial emission.

# INPUTS:
# DIR (str) : Directory containing the Fits files
# JWST_FILE (str) : File name of the JWST data
# PHANGS_FILE (str) : File name of the PHANGS data
# OUTPUT_FILENAME (str) : File name of the overlay output plot
# PIXEL_SCALE (float) : Pixel scale in arcsecs of the two input images. This 
#    is used for creating a 1 kpc scale bar
# GALAXY_DISTANCE (float) : Galaxy distance in Mpc. This is also used for the
#     scale bar
# RESOLUTION (int) : Resolution of the output image
# SAVE_OUTPUT (bol) : True if the plot is to be output to a file. False for
#     display
#
# OUTPUTS:
# Overlay of the datasets in the form of a plot with file name OUTPUT_FILENAME

# INPUTS
dir = './'
jwst_file = dir + 'convolved_jwst.fits'
phangs_file = dir + 'cropped_phangs.fits'
extension = 0

save_output = True  # Set to False to display instead
output_filename = 'ngc628_jwst_phangs_overlay.png'
resolution = 300

# Scale bar parameters
pixel_scale = 0.1  # arcsec/pixel
galaxy_distance = 9.5  # Mpc
kpc_per_arcsec = galaxy_distance * 1000 / 206265  # Convert Mpc to kpc and arcsec to radians
pixels_per_kpc = 1 / (pixel_scale * kpc_per_arcsec)  # pixels per kpc

# Visualisation 
threshold = 4.0
colour_settings = {
    'jwst': {'cmap': 'Reds', 'vmin': 276, 'vmax': 280, 'alpha': 0.8, 'label': 'JWST'},
    'phangs': {'cmap': 'Blues', 'vmin': 0, 'vmax': 25, 'alpha': 0.8, 'label': 'PHANGS'}
}

def load_and_mask_image(filename, threshold=-5000):
    """Load FITS image and mask values below threshold"""
    with fits.open(filename) as hdul:
        data = hdul[extension].data.copy()  # Copy to avoid modifying original
        data[data < threshold] = -5000
        return data

def normalise_and_colourise(data, vmin, vmax, cmap_name):
    """Apply normalisation and colourmap to data"""
    cmap = plt.get_cmap(cmap_name)
    norm = Normalize(vmin=vmin, vmax=vmax)
    colored = cmap(norm(data))
    colored[data == -5000, :] = [0, 0, 0, 0]  # Make masked areas transparent
    return colored

def add_scale_bar(ax, pixels_per_kpc, length_kpc=1, position=(0.8, 0.1), color='black'):
    """Add a scale bar showing 1 kpc"""
    # Convert kpc to pixels
    length_pixels = length_kpc * pixels_per_kpc
    
    # Get figure dimensions
    transform = ax.transAxes
    
    # Create scale bar
    scalebar_x = position[0]
    scalebar_y = position[1]
    scalebar_width = length_pixels / ax.get_xlim()[1]
    scalebar_height = 0.01
    
    # Add scale bar as a thin rectangle
    scale_bar = Rectangle((scalebar_x, scalebar_y), scalebar_width, scalebar_height, 
                        transform=transform, color=color)
    ax.add_patch(scale_bar)
    
    ax.text(scalebar_x + scalebar_width/2, scalebar_y + 0.03, f'{length_kpc} kpc', 
            transform=transform, color=color, ha='center')

# Load and process JWST image
jwst_data = load_and_mask_image(jwst_file, threshold)
jwst_coloured = normalise_and_colourise(
    jwst_data, 
    colour_settings['jwst']['vmin'],
    colour_settings['jwst']['vmax'],
    colour_settings['jwst']['cmap']
)

# Load and process PHANGS image
phangs_data = load_and_mask_image(phangs_file, threshold)
phangs_coloured = normalise_and_colourise(
    phangs_data,
    colour_settings['phangs']['vmin'],
    colour_settings['phangs']['vmax'],
    colour_settings['phangs']['cmap']
)

phangs_coloured[phangs_data == -5000, :] = [1, 1, 1, 1]

fig, ax = plt.subplots(figsize=(10, 8), facecolor='white')
plt.margins(0, 0)
ax.imshow(jwst_coloured, alpha=colour_settings['jwst']['alpha'])
ax.imshow(phangs_coloured, alpha=colour_settings['phangs']['alpha'])

ax.set_xticks([])
ax.set_yticks([])
ax.set_xticklabels([])
ax.set_yticklabels([])
for spine in ax.spines.values():
    spine.set_visible(False)

# Legend
legend = ax.legend(
    handles=[
        mlines.Line2D([], [], color='red', marker='s', linestyle='None', markersize=15, label=colour_settings['jwst']['label']),
        mlines.Line2D([], [], color='blue', marker='s', linestyle='None', markersize=15, label=colour_settings['phangs']['label'])
    ],
    loc='upper right',
    frameon=True,
    labelcolor='black',
    edgecolor='black'
)
legend.get_frame().set_facecolor('white')

# Add scale bar (1 kpc) 
add_scale_bar(ax, pixels_per_kpc, color='black')

plt.tight_layout()
if save_output:
    plt.savefig(output_filename, dpi=resolution, bbox_inches='tight', facecolor='white')
else:
    plt.show()
