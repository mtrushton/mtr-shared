pro inclination_corr, pix_size, q_obs, Rs_pix, Re_bulge, dist_gal, band, $
                    q_corr, incl_corr, Rs_arcsec, Rs_kpc, Reff_b_arcsec, Rb_kpc, status
;
;     Performs inclination corrections for galaxy disk parameters and calculates
;     physical sizes from pixel measurements.
;
; INPUTS:
;     PIX_SIZE    - Pixel size in arcseconds
;     Q_OBS       - Observed axis ratio (b/a)
;     RS_PIX      - Disk scale length in pixels
;     RE_BULGE    - Bulge effective radius in pixels
;     DIST_GAL    - Galaxy distance in Mpc
;     BAND        - Photometric band ('B','V','I','J','K')
;
; OUTPUTS:
;     Q_CORR      - Intrinsic (corrected) axis ratio
;     INCL_CORR   - Corrected inclination in degrees
;     RS_ARCSEC   - Disk scale length in arcseconds
;     RS_KPC      - Disk scale length in kpc
;     REFF_B_ARCSEC - Bulge effective radius in arcseconds
;     RB_KPC      - Bulge effective radius in kpc
;     STATUS      - 0 for success, negative for error
;-
  status = 0
  
  ; Define photometric bands and their associated intrinsic axis ratios
  band_index = ['B', 'V', 'I', 'J', 'K']
  q_model_wave = [0.074, 0.076, 0.085, 0.095, 0.108]
  
  ; Input validation
  if pix_size le 0 then begin
    print, 'Error: Pixel size must be positive'
    status = -1
    return
  endif
  
  if q_obs le 0 || q_obs gt 1 then begin
    print, 'Error: Observed axis ratio must be between 0 and 1'
    status = -2
    return
  endif
  
  if Rs_pix le 0 then begin
    print, 'Error: Disk scale length must be positive'
    status = -3
    return
  endif
  
  if Re_bulge le 0 then begin
    print, 'Error: Bulge effective radius must be positive'
    status = -4
    return
  endif
  
  if dist_gal le 0 then begin
    print, 'Error: Galaxy distance must be positive'
    status = -5
    return
  endif
  
  ; Find model axis ratio for the given band
  band_index = strupcase(band_index)
  band_upper = strupcase(band)
  band_idx = where(band_index eq band_upper, count)
  
  if count eq 0 then begin
    print, 'Error: Unrecognized photometric band: ', band
    print, 'Valid bands are: ', band_index
    status = -6
    return
  endif
  
  q_model = q_model_wave[band_idx[0]]
  
  ; Calculate observed inclination
  incl = acos(q_obs) * 180.0 / !PI
  print, 'Disk inclination (degrees) = ', incl
  
  ; Define constants for coordinate conversion
  arcsec_to_kpc = !PI * dist_gal * 1000.0 / (3600.0 * 180.0)
  
  ; Calculate physical sizes
  Rs_arcsec = Rs_pix * pix_size
  Reff_b_arcsec = Re_bulge * pix_size
  Rs_kpc = Rs_arcsec * arcsec_to_kpc
  Rb_kpc = Reff_b_arcsec * arcsec_to_kpc
  
  ; Print physical sizes
  print, 'Disk scale-length (arcsecs) = ', Rs_arcsec
  print, 'Disk scale-length (kpc) = ', Rs_kpc
  print, 'Bulge eff radius (arcsecs) = ', Reff_b_arcsec
  print, 'Bulge eff radius (kpc) = ', Rb_kpc
  
  ; Calculate intrinsic parameters
  if q_obs le q_model then begin
    print, 'Warning: Observed axis ratio is less than model intrinsic ratio'
    print, 'Setting corrected axis ratio to model value: ', q_model
    q_corr = q_model
  endif else begin
    q_corr = sqrt((q_obs^2 - q_model^2) / (1.0 - q_model^2))
  endelse
  
  incl_corr = acos(q_corr) * 180.0 / !PI
  
  ; Print corrected values
  print, 'Corrected axis-ratio = ', q_corr
  print, 'Corrected disk inclination = ', incl_corr
  
  return
end
