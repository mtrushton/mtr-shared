PRO cog, x_in, mask_in, pa, ellipt, xcen, ycen, wa, rcog, bcog, ncog, scog, r2pi
  ;
  ;   Computes surface brightness profiles and curve of growth in elliptical annuli
  ;   for a 2D brightness map.
  ;
  ; INPUTS:
  ;   x_in    - 2D image array (flux/pixel)
  ;   mask_in - 2D mask array (0=good pixels, non-zero=bad pixels)
  ;   pa      - Position angle of major axis (radians, positive clockwise from vertical)
  ;   ellipt  - Ratio of minor to major axis (consistent with GALFIT output)
  ;   xcen    - X-coordinate of galaxy center
  ;   ycen    - Y-coordinate of galaxy center
  ;   wa      - Desired width of annulus in pixels
  ;
  ; OUTPUTS:
  ;   rcog    - Outer positions of elliptical annuli along major axis
  ;   bcog    - Average brightness (flux/pixel) in each annulus
  ;   ncog    - Number of good pixels in each annulus
  ;   scog    - Cumulative flux within radius rcog
  ;   r2pi    - Largest element of RCOG with full azimuthal coverage
  ;

  ; Get image dimensions
  nx = DOUBLE(N_ELEMENTS(x_in[*,0]))
  ny = DOUBLE(N_ELEMENTS(x_in[0,*]))
  
  ; Calculate elliptical radii for all pixels
  dist_ellipse, r, [nx, ny], xcen, ycen, 1.0D/ellipt, pa
  
  ; Determine number of annuli and their boundaries
  xnr = MAX(r)/wa
  nr = CEIL(xnr)
  waout = MAX(r)/FLOAT(nr)
  rcog = (1.0D + DINDGEN(nr)) * waout
  
  ; Find r2pi - radius with full azimuthal coverage
  a1 = MIN(r[0,0:ny-1L])
  a2 = MIN(r[nx-1L,0:ny-1L])
  a3 = MIN(r[0:nx-1L,0])
  a4 = MIN(r[0:nx-1L,ny-1L])
  r2pi = MIN([a1,a2,a3,a4])
  
  ; Initialise output arrays
  bcog = FLTARR(nr)
  ncog = FLTARR(nr)
  scog = bcog
  
  ; Sort pixels by radius
  irsort = SORT(r)
  r = r[irsort]
  x = x_in[irsort]
  mask = mask_in[irsort]
  npix = nx * ny
  
  ; Process each annulus 
  ir = 0L          ; Counter for the annulus
  ipix = -1L       ; Counter for the (radius sorted) map pixels
  r1 = 0.0D
  r2 = rcog[ir]
  n_inbin = 0L     ; Good pixels in bin
  npix_inbin = 0L  ; Total pixels (good and bad) in bin
  flux_inbin = 0.0D
  nr = nr - 1L     
  
  c555:
    ipix = ipix + 1L
    
    ; Handle special case for central pixel
    IF (ipix EQ 0) THEN r[ipix] = r1 + 0.001*(r2-r1)
    
    IF (r[ipix] GT r1 AND r[ipix] LE r2) THEN BEGIN
      npix_inbin = npix_inbin + 1L
      ; Use mask == 1 to identify good pixels
      IF (mask[ipix] EQ 1) THEN BEGIN
        n_inbin = n_inbin + 1L
        flux_inbin = flux_inbin + x[ipix]
      ENDIF ; mask[ipix] EQ 1
      GOTO, c555
    ENDIF ; r[ipix] GT r1 AND r[ipix] LE r2
    
    IF (r[ipix] GT r2) THEN BEGIN
      ncog[ir] = n_inbin
      
      IF (n_inbin EQ 0) THEN BEGIN
        bcog[ir] = 0.0D
        scog[ir] = 0.0D
        ir = ir + 1L
        IF (ir GE nr) THEN GOTO, c666
        r1 = rcog[ir-1L]
        r2 = rcog[ir]
        n_inbin = 0L
        npix_inbin = 0L
        flux_inbin = 0.0D
        GOTO, c555
      ENDIF ; n_inbin EQ 0
      
      IF (n_inbin GT 0) THEN BEGIN
        bcog[ir] = flux_inbin/FLOAT(n_inbin)
        s_inbin = bcog[ir] * npix_inbin
        IF (ir EQ 0) THEN scog[ir] = s_inbin ELSE scog[ir] = scog[ir-1L] + s_inbin
        ir = ir + 1L
        IF (ir GE nr) THEN GOTO, c666
        r1 = rcog[ir-1L]
        r2 = rcog[ir]
        n_inbin = 0L
        npix_inbin = 0L
        flux_inbin = 0.0D
        GOTO, c555
      ENDIF ; n_inbin GT 0
    ENDIF ; r[ipix] GT r2
  
  c666:
  
  RETURN
END
