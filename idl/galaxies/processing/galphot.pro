pro galphot
;
; Plots the surface brightness profile and curve of growth (COG) for a galaxy
; from GALFIT output. 
;
; Simply run the procedure with no parameters: IDL> gal_phot
;
; ROUTINES CALLED:
; COG_INPUT    - Input galaxy names, GALFIT output file, mask file, and define range and labels for the plots
; RHEADER      - Reads galaxy parameters and data from the GALFIT output file
; COG          - Calculates the surface brightness profile and cog
; INCLINATION_CORR - Calculates the galaxy inclination and morphological parameters in arcsecs and kpc
; COG_PLOT     - Plots the surface brightness profile and cog
; ERR_COG      - Calculates errors in the surface brightness profile and cog
;
; OUTPUTS:
; Generates plots of the surface brightness profile and curve of growth
; Files are saved in the current directory according to the naming convention:
; [gal_name]_[band]_cog.ps - For the curve of growth
; [gal_name]_[band]_prof.ps - For the surface brightness profile
;
;
 compile_opt idl2  ; Enforce strict IDL compilation rules
 
 ; Define inputs
 cog_input, gal_name, fname, mname, dist_gal, band, pix_size, wa, sky_sub, r_max, ps, comps, sky_lev, $
            plot_reff, plot_rmask, rmask, plot_xrange_prof, plot_yrange_prof, plot_xrange_cog, plot_yrange_cog, $
            xlab, prof_ylab, cog_ylab, sky, plot_log, x_lab_prof, x_lab_cog, y_lab_prof, y_lab_cog, ierr
 
 ; Check for errors in input parameters
 if(ierr ne 0) then begin
   print, 'TEST: Error in input parameters. Please correct and try again.'
   return
 endif
 
 ; Print info about the analysis
 print, '======================================================='
 print, 'Processing galaxy: ' + gal_name + ' in band: ' + band
 print, 'Using GALFIT output file: ' + fname
 print, 'Using mask file: ' + mname
 print, '======================================================='
 
 ; Read header from the GALFIT output fits file
 catch, error_status  ; Enable error handling
 if (error_status ne 0) then begin
   print, 'TEST: Error reading header from file ' + fname
   print, 'Error message: ', !ERROR_STATE.MSG
   return
 endif
 
 rheader, fname, mname, pa, ellipt, xcen, ycen, x_data, x_bulge, x_disk, mask_in, Rs_pix, Re_bulge, S_index
 
 catch, /cancel  ; Disable error handling
 
 ; Check if the profile and cog plot are to be generated for sky subtracted data
 sky_sub = strlowcase(sky_sub)  ; Convert to lowercase for consistent comparison
 
 if(sky_sub eq 'y') then begin
   x_in = x_data - sky
   print, 'Using sky-subtracted data (sky value: ', strtrim(sky,2), ')'
 endif else begin
   x_in = x_data
   print, 'Using original data (no sky subtraction)'
 endelse
 
 ; Calculate surface brightness profile and cog
 print, 'Calculating surface brightness profile and curve of growth...'
 cog, x_in, mask_in, pa, ellipt, xcen, ycen, wa, rcog, bcog, ncog, scog, r2pi
 
 ; Check if bulge and disk components are to be plotted
 comps = strlowcase(comps)  ; Convert to lowercase for consistent comparison
 
 if(comps eq 'y') then begin
   print, 'Calculating bulge and disk components...'
   cog, x_bulge, mask_in, pa, ellipt, xcen, ycen, wa, rcog, bcog_bulge, ncog, scog_bulge, r2pi
   cog, x_disk, mask_in, pa, ellipt, xcen, ycen, wa, rcog, bcog_disk, ncog, scog_disk, r2pi
 endif else begin
   ; Initialize variables to avoid undefined variable errors
   bcog_bulge = dblarr(n_elements(bcog))
   scog_bulge = dblarr(n_elements(scog))
   bcog_disk = dblarr(n_elements(bcog))
   scog_disk = dblarr(n_elements(scog))
 endelse
 
 ; Calculate galaxy inclination and morphology parameters in arcsecs and kpc
 print, 'Calculating galaxy inclination and morphological parameters...'
 inclination_corr, pix_size, ellipt, Rs_pix, Re_bulge, dist_gal, band, Q_corr, incl, Rs_arcsec, Rs_kpc, $
                   Reff_b_arcsec, Rb_kpc  
 
 ; Convert radius to arcseconds
 r_arcsecs = rcog * pix_size
 
 ; Plot the surface brightness profile and cog
 print, 'Generating plots...'
 cog_plot, gal_name, pix_size, r_max, rcog, bcog, scog, bcog_bulge, scog_bulge, bcog_disk, scog_disk, $
           ps, comps, sky_lev, plot_reff, plot_rmask, rmask, rs_arcsec, r_arcsecs, band, ellipt, incl, Reff_b_arcsec, S_index, r2pi, $
           plot_xrange_prof, plot_yrange_prof, plot_xrange_cog, plot_yrange_cog, xlab, prof_ylab, cog_ylab, $
           sky, plot_log, x_lab_prof, x_lab_cog, y_lab_prof, y_lab_cog
 
 ; Calculate errors in the surface brightness profile and cog
 print, 'Calculating errors...'
 err_cog, sky, rcog, bcog, scog, sky_sub, r_arcsecs, ellipt, wa
 
 print, 'Analysis complete!'
 if(ps eq 'y') then begin
   print, 'Output files:'
   print, '  ' + gal_name + '_' + band + '_prof.ps'
   print, '  ' + gal_name + '_' + band + '_cog.ps'
 endif
 
 print, '======================================================='
 
end
