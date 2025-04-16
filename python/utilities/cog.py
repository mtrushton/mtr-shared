import numpy as np

def dist_ellipse(nx, ny, xcen, ycen, axis_ratio, pa):
    """
    Calculate elliptical radii for pixels in an image.
    
    Parameters:
    -----------
    nx, ny : int
        Image dimensions
    xcen, ycen : float
        Center coordinates
    axis_ratio : float
        Ratio of minor to major axis (b/a)
    pa : float
        Position angle (radians, positive clockwise from vertical)
    
    Returns:
    --------
    r : ndarray
        2D array of elliptical radii
    """
    y, x = np.indices((ny, nx))
    x = x.astype(float) - xcen
    y = y.astype(float) - ycen
    
    # Rotate coordinates
    xp = x * np.cos(pa) + y * np.sin(pa)
    yp = -x * np.sin(pa) + y * np.cos(pa)
    
    # Calculate elliptical radius
    r = np.sqrt((xp**2) + (yp/axis_ratio)**2)
    
    return r

def cog(x_in, mask_in, pa, ellipt, xcen, ycen, wa):
    """
    Computes surface brightness profiles and curve of growth in elliptical annuli
    for a 2D image array
    
    Parameters:
    -----------
    x_in : ndarray
        2D image array (flux/pixel)
    mask_in : ndarray
        2D mask array (0=bad pixels, 1=good pixels)
    pa : float
        Position angle of major axis (radians, positive clockwise from vertical)
    ellipt : float
        Ratio of minor to major axis (consistent with GALFIT output)
    xcen, ycen : float
        Coordinates of galaxy center
    wa : float
        Desired width of annulus in pixels
    
    Returns:
    --------
    rcog : ndarray
        Outer positions of elliptical annuli along major axis
    bcog : ndarray
        Average brightness (flux/pixel) in each annulus
    ncog : ndarray
        Number of good pixels in each annulus
    scog : ndarray
        Cumulative flux within radius rcog
    r2pi : float
        Largest element of RCOG with full azimuthal coverage
    """
    # Image dimensions
    ny, nx = x_in.shape
    
    # Calculate elliptical radii for all pixels
    r = dist_ellipse(nx, ny, xcen, ycen, 1.0/ellipt, pa)
    
    # Determine n of annuli and boundaries
    xnr = np.max(r) / wa
    nr = int(np.ceil(xnr))
    waout = np.max(r) / float(nr)
    rcog = (1.0 + np.arange(nr, dtype=float)) * waout
    
    # Radius with full azimuthal coverage
    a1 = np.min(r[0, :])
    a2 = np.min(r[nx-1, :])
    a3 = np.min(r[:, 0])
    a4 = np.min(r[:, ny-1])
    r2pi = min(a1, a2, a3, a4)
    
    bcog = np.zeros(nr, dtype=float)
    ncog = np.zeros(nr, dtype=float)
    scog = np.zeros(nr, dtype=float)
    
    # Sort pixels by radius
    r_flat = r.flatten()
    x_flat = x_in.flatten()
    mask_flat = mask_in.flatten()
    
    sort_indices = np.argsort(r_flat)
    r_sorted = r_flat[sort_indices]
    x_sorted = x_flat[sort_indices]
    mask_sorted = mask_flat[sort_indices]
    
    npix = nx * ny
    
 
    ir = 0  
    ipix = -1  
    r1 = 0.0
    r2 = rcog[ir]
    n_inbin = 0  # Good pixels in bin
    npix_inbin = 0  # Total pixels (good and bad) in bin
    flux_inbin = 0.0
    nr = nr - 1
    
    while True:
        ipix += 1
        
        # Handle special case for central pixel
        if ipix == 0:
            r_sorted[ipix] = r1 + 0.001 * (r2 - r1)
        
        if ipix >= len(r_sorted):
            break
            
        if r_sorted[ipix] > r1 and r_sorted[ipix] <= r2:
            npix_inbin += 1
            # Use mask == 1 to identify good pixels
            if mask_sorted[ipix] == 1:
                n_inbin += 1
                flux_inbin += x_sorted[ipix]
            continue
        
        if r_sorted[ipix] > r2:
            ncog[ir] = n_inbin
            
            if n_inbin == 0:
                bcog[ir] = 0.0
                scog[ir] = 0.0
            else:
                bcog[ir] = flux_inbin / float(n_inbin)
                s_inbin = bcog[ir] * npix_inbin
                if ir == 0:
                    scog[ir] = s_inbin
                else:
                    scog[ir] = scog[ir-1] + s_inbin
            
            ir += 1
            if ir >= nr:
                break
                
            r1 = rcog[ir-1]
            r2 = rcog[ir]
            n_inbin = 0
            npix_inbin = 0
            flux_inbin = 0.0
    
    return rcog, bcog, ncog, scog, r2pi
