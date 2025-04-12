pro cog_plot, gal_name, pix_size, r_max, rcog, bcog, scog, bcog_bulge, scog_bulge, bcog_disk, scog_disk, $
              ps, comps, sky_lev, plot_reff, plot_rmask, rmask, rs_arcsec, r_arcsecs, band, ellipt, incl, $
              Reff_b_arcsec, S_index, r2pi, plot_xrange_prof, plot_yrange_prof, plot_xrange_cog, $
              plot_yrange_cog, xlab, prof_ylab, cog_ylab, sky, plot_log, x_lab_prof, x_lab_cog, $
              y_lab_prof, y_lab_cog, status=status
      
;     Plots surface brightness profile and curve of growth (COG) of a galaxy
;     with optional component decomposition.
;
; INPUTS:
;     GAL_NAME      - String identifier of the galaxy
;     PIX_SIZE      - Scalar pixel size in arcsecs
;     R_MAX         - Integer. Maximum radius in pixels
;     RCOG          - Radii of the annuli in pixels
;     BCOG          - Surface brightness of each annuli
;     SCOG          - Total flux encompassed by annuli
;     BCOG_BULGE    - Surface brightness of annuli for bulge component
;     SCOG_BULGE    - Total flux in annuli for bulge component
;     BCOG_DISK     - Surface brightness of annuli for disk component
;     SCOG_DISK     - Total flux in annuli for disk component
;     PS            - String ('y' or 'n') for postscript output
;     COMPS         - String ('y' or 'n') to plot bulge/disk components
;     SKY_LEV       - String ('y' or 'n') to plot horizontal sky level line
;     PLOT_REFF     - String ('y' or 'n') to plot half-light radius line
;     PLOT_RMASK    - String ('y' or 'n') to plot mask radius line
;     RMASK         - Integer. Limiting radius in arcsecs if used in fitting
;     RS_ARCSEC     - Scalar scale length of the disk in arcsecs
;     BAND          - String containing the observed band
;     ELLIPT        - Scalar axis ratio b/a of the disk
;     INCL          - Scalar inclination angle in degrees
;     REFF_B_ARCSEC - Scalar effective radius of the bulge
;     S_INDEX       - Scalar sersic index of the bulge
;     R2PI          - Scalar largest element for RCOG with 2pi azimuthal coverage
;     PLOT_XRANGE_PROF - 2D vector for x range in brightness profile
;     PLOT_YRANGE_PROF - 2D vector for y range in brightness profile
;     PLOT_XRANGE_COG  - 2D vector for x range in COG plot
;     PLOT_YRANGE_COG  - 2D vector for y range in COG plot
;     XLAB          - String label of the x axis
;     PROF_YLAB     - String y axis label for brightness profile
;     COG_YLAB      - String y axis label for COG plot
;     SKY           - Scalar value of the sky level
;     PLOT_LOG      - Integer. 1 for log plots in y, 0 for linear
;     X_LAB_PROF    - Vector of values for x positions of profile labels
;     Y_LAB_PROF    - Vector of values for y positions of profile labels
;     X_LAB_COG     - Vector of values for x positions of COG labels
;     Y_LAB_COG     - Vector of values for y positions of COG labels
;
; OPTIONAL OUTPUTS:
;     STATUS        - 0 for success, negative for error
;-
  
  ; Initialise status
  status = 0
  
  ; Input validation
  if n_elements(rcog) eq 0 || n_elements(bcog) eq 0 || n_elements(scog) eq 0 then begin
    message, 'Required input arrays are missing or empty', /CONTINUE
    status = -1
    return
  endif
  
  ; Make sure r_arcsecs is available
  if n_elements(r_arcsecs) eq 0 then begin
    r_arcsecs = rcog * pix_size
  endif
  
  ; Calculate r_max in arcsecs
  r_max_arcsecs = r_max * pix_size * 3L
  
  ; Format display strings with better precision control
  str_incl = string(incl, FORMAT='(F7.2)')
  str_sky = string(sky, FORMAT='(G12.5)')
  str_rs = string(rs_arcsec, FORMAT='(F7.2)')
  str_reff = string(Reff_b_arcsec, FORMAT='(F7.2)')
  str_ellipt = string(ellipt, FORMAT='(F6.3)')
  str_sersic = string(S_index, FORMAT='(F5.2)')
  str_max_scog = string(max(scog), FORMAT='(E11.4)')
  
  ; Calculate bulge-to-disk ratio if components are requested
  if strlowcase(comps) eq 'y' then begin
    if max(scog_disk) gt 0 then begin
      bd = max(scog_bulge) / max(scog_disk)
      str_bd = string(bd, FORMAT='(F6.3)')
    endif else begin
      message, 'Warning: Disk component has zero or negative maximum flux', /CONTINUE
      str_bd = 'N/A'
    endelse
  endif
  
  ; Set up plot device
  orig_device = !D.NAME
  
  if strlowcase(ps) eq 'y' then begin
    set_plot, 'PS'
    device, file=gal_name+'_brprof.ps', /COLOR, BITS_PER_PIXEL=8
  endif else begin
    window, 0, XSIZE=600, YSIZE=500, TITLE='Surface Brightness Profile: ' + gal_name
  endelse
  
  ; Define colors for better visibility
  data_color = 0       ; Black for data
  bulge_color = 200    ; Red for bulge
  disk_color = 80      ; Green for disk
  model_color = 250    ; Blue for total model
  
  ; Plot surface brightness profile
  if plot_log eq 1L then begin
    ; Log plot
    plot, r_arcsecs, bcog, YSTYLE=1, XSTYLE=1, $
          XRANGE=plot_xrange_prof, YRANGE=plot_yrange_prof, $
          XTITLE=xlab, YTITLE=prof_ylab, YLOG=1, $
          THICK=2, CHARSIZE=1.2, CHARTHICK=1.5
          
    ; Add components if requested
    if strlowcase(comps) eq 'y' then begin
      oplot, r_arcsecs, bcog_bulge, COLOR=bulge_color, LINESTYLE=0, THICK=2
      oplot, r_arcsecs, bcog_disk, COLOR=disk_color, LINESTYLE=0, THICK=2
      oplot, r_arcsecs, bcog_bulge + bcog_disk, COLOR=model_color, LINESTYLE=2, THICK=2
      
      ; Add legend
      plots, [x_lab_prof-10, x_lab_prof-5], [y_lab_prof(7)*1.2, y_lab_prof(7)*1.2], COLOR=bulge_color, THICK=2
      xyouts, x_lab_prof-4, y_lab_prof(7)*1.2, 'Bulge', CHARSIZE=0.9
      plots, [x_lab_prof-10, x_lab_prof-5], [y_lab_prof(7)*1.1, y_lab_prof(7)*1.1], COLOR=disk_color, THICK=2
      xyouts, x_lab_prof-4, y_lab_prof(7)*1.1, 'Disk', CHARSIZE=0.9
      plots, [x_lab_prof-10, x_lab_prof-5], [y_lab_prof(7)*1.0, y_lab_prof(7)*1.0], COLOR=model_color, THICK=2, LINESTYLE=2
      xyouts, x_lab_prof-4, y_lab_prof(7)*1.0, 'Model', CHARSIZE=0.9
    endif
    
    ; Add reference lines
    if strlowcase(sky_lev) eq 'y' then begin
      oplot, plot_xrange_prof, [sky, sky], LINESTYLE=2, THICK=1.5
    endif
    
    oplot, [r2pi*pix_size, r2pi*pix_size], plot_yrange_prof, LINESTYLE=1, THICK=1.5
    oplot, [r_max_arcsecs, r_max_arcsecs], plot_yrange_prof, LINESTYLE=0, THICK=1.5
    
    ; Add labels for reference lines
    xyouts, 1.10*r2pi*pix_size, plot_yrange_prof(0)*1.6, 'R(2!4p!3)', ORIENTATION=90, CHARSIZE=1.0, CHARTHICK=1.2
    xyouts, 1.10*r_max_arcsecs, plot_yrange_prof(0)*1.7, 'R!Dmax!N', ORIENTATION=90, CHARSIZE=1.0, CHARTHICK=1.2
    
    ; Add mask radius line if requested
    if strlowcase(plot_rmask) eq 'y' then begin
      oplot, [rmask, rmask], plot_yrange_prof, LINESTYLE=2, THICK=1.5
      xyouts, 1.08*rmask, plot_yrange_prof(0)*2.0, 'R!Dmask!N', ORIENTATION=90, CHARSIZE=1.0, CHARTHICK=1.2
    endif
    
  endif else begin
    ; Linear plot
    plot, r_arcsecs, bcog, YSTYLE=1, XSTYLE=1, $
          XRANGE=plot_xrange_prof, YRANGE=plot_yrange_prof, $
          XTITLE=xlab, YTITLE=prof_ylab, $
          THICK=2, CHARSIZE=1.2, CHARTHICK=1.5
          
    ; Add components if requested
    if strlowcase(comps) eq 'y' then begin
      oplot, r_arcsecs, bcog_bulge, COLOR=bulge_color, LINESTYLE=0, THICK=2
      oplot, r_arcsecs, bcog_disk, COLOR=disk_color, LINESTYLE=0, THICK=2
      oplot, r_arcsecs, bcog_bulge + bcog_disk, COLOR=model_color, LINESTYLE=2, THICK=2
      
      ; Add legend
      plots, [x_lab_prof-10, x_lab_prof-5], [y_lab_prof(7)*1.2, y_lab_prof(7)*1.2], COLOR=bulge_color, THICK=2
      xyouts, x_lab_prof-4, y_lab_prof(7)*1.2, 'Bulge', CHARSIZE=0.9
      plots, [x_lab_prof-10, x_lab_prof-5], [y_lab_prof(7)*1.1, y_lab_prof(7)*1.1], COLOR=disk_color, THICK=2
      xyouts, x_lab_prof-4, y_lab_prof(7)*1.1, 'Disk', CHARSIZE=0.9
      plots, [x_lab_prof-10, x_lab_prof-5], [y_lab_prof(7)*1.0, y_lab_prof(7)*1.0], COLOR=model_color, THICK=2, LINESTYLE=2
      xyouts, x_lab_prof-4, y_lab_prof(7)*1.0, 'Model', CHARSIZE=0.9
    endif
    
    ; Add reference lines
    if strlowcase(sky_lev) eq 'y' then begin
      oplot, plot_xrange_prof, [sky, sky], LINESTYLE=2, THICK=1.5
    endif
    
    oplot, [r2pi*pix_size, r2pi*pix_size], plot_yrange_prof, LINESTYLE=1, THICK=1.5
    oplot, [r_max_arcsecs, r_max_arcsecs], plot_yrange_prof, LINESTYLE=0, THICK=1.5
    
    ; Add labels for reference lines
    xyouts, 1.05*r2pi*pix_size, plot_yrange_prof(0)*0.7, 'R(2!4p!3)', ORIENTATION=90, CHARSIZE=1.0, CHARTHICK=1.2
    xyouts, 1.10*r_max_arcsecs, plot_yrange_prof(0)*0.7, 'R!Dmax!N', ORIENTATION=90, CHARSIZE=1.0, CHARTHICK=1.2
    
    ; Add mask radius line if requested
    if strlowcase(plot_rmask) eq 'y' then begin
      oplot, [rmask, rmask], plot_yrange_prof, LINESTYLE=2, THICK=1.5
      xyouts, 1.08*rmask, plot_yrange_prof(0)*2.0, 'R!Dmask!N', ORIENTATION=90, CHARSIZE=1.0, CHARTHICK=1.2
    endif
  endelse
  
  ; Add galaxy properties as labels
  xyouts, x_lab_prof, y_lab_prof(0), gal_name, CHARSIZE=1.2, CHARTHICK=1.5
  xyouts, x_lab_prof, y_lab_prof(1), band + ' band', CHARSIZE=1.0
  xyouts, x_lab_prof, y_lab_prof(2), 'Sky = ' + strtrim(str_sky, 2), CHARSIZE=1.0
  xyouts, x_lab_prof, y_lab_prof(3), 'Rs = ' + strtrim(str_rs, 2), CHARSIZE=1.0
  xyouts, x_lab_prof, y_lab_prof(4), 'Rb = ' + strtrim(str_reff, 2), CHARSIZE=1.0
  xyouts, x_lab_prof, y_lab_prof(5), 'b/a = ' + strtrim(str_ellipt, 2), CHARSIZE=1.0
  xyouts, x_lab_prof, y_lab_prof(6), 'n = ' + strtrim(str_sersic, 2), CHARSIZE=1.0
  
  if strlowcase(comps) eq 'y' then begin
    xyouts, x_lab_prof, y_lab_prof(7), 'b/d = ' + strtrim(str_bd, 2), CHARSIZE=1.0
  endif
  
  ; Close first plot if in PS mode
  if strlowcase(ps) eq 'y' then begin
    device, /CLOSE
    
    ; Open new file for COG plot
    device, FILE=gal_name+'_cog.ps', /COLOR, BITS_PER_PIXEL=8
  endif else begin
    ; Create new window for COG plot
    window, 1, XSIZE=600, YSIZE=500, TITLE='Curve of Growth: ' + gal_name
  endelse
  
  ; Plot COG
  plot, r_arcsecs, scog, YSTYLE=1, XSTYLE=1, $
        XRANGE=plot_xrange_cog, YRANGE=plot_yrange_cog, $
        XTITLE=xlab, YTITLE=cog_ylab, YLOG=1, $
        THICK=2, CHARSIZE=1.2, CHARTHICK=1.5
  
  ; Add components if requested
  if strlowcase(comps) eq 'y' then begin
    oplot, r_arcsecs, scog_bulge, COLOR=bulge_color, LINESTYLE=0, THICK=2
    oplot, r_arcsecs, scog_disk, COLOR=disk_color, LINESTYLE=0, THICK=2
    oplot, r_arcsecs, scog_bulge + scog_disk, COLOR=model_color, LINESTYLE=2, THICK=2
    
    ; Add legend
    plots, [x_lab_cog-10, x_lab_cog-5], [y_lab_cog(7)*1.2, y_lab_cog(7)*1.2], COLOR=bulge_color, THICK=2
    xyouts, x_lab_cog-4, y_lab_cog(7)*1.2, 'Bulge', CHARSIZE=0.9
    plots, [x_lab_cog-10, x_lab_cog-5], [y_lab_cog(7)*1.1, y_lab_cog(7)*1.1], COLOR=disk_color, THICK=2
    xyouts, x_lab_cog-4, y_lab_cog(7)*1.1, 'Disk', CHARSIZE=0.9
    plots, [x_lab_cog-10, x_lab_cog-5], [y_lab_cog(7)*1.0, y_lab_cog(7)*1.0], COLOR=model_color, THICK=2, LINESTYLE=2
    xyouts, x_lab_cog-4, y_lab_cog(7)*1.0, 'Model', CHARSIZE=0.9
  endif
  
  ; Add galaxy properties as labels
  xyouts, x_lab_cog, y_lab_cog(0), gal_name, CHARSIZE=1.2, CHARTHICK=1.5
  xyouts, x_lab_cog, y_lab_cog(1), band + ' band', CHARSIZE=1.0
  xyouts, x_lab_cog, y_lab_cog(2), 'Sky = ' + strtrim(str_sky, 2), CHARSIZE=1.0
  xyouts, x_lab_cog, y_lab_cog(3), 'Rs = ' + strtrim(str_rs, 2), CHARSIZE=1.0
  xyouts, x_lab_cog, y_lab_cog(4), 'Rb = ' + strtrim(str_reff, 2), CHARSIZE=1.0
  xyouts, x_lab_cog, y_lab_cog(5), 'b/a = ' + strtrim(str_ellipt, 2), CHARSIZE=1.0
  xyouts, x_lab_cog, y_lab_cog(6), 'n = ' + strtrim(str_sersic, 2), CHARSIZE=1.0
  
  if strlowcase(comps) eq 'y' then begin
    xyouts, x_lab_cog, y_lab_cog(7), 'b/d = ' + strtrim(str_bd, 2), CHARSIZE=1.0
  endif
  
  ; Add total flux value
  xyouts, 0.75*plot_xrange_cog(1), plot_yrange_cog(1)*0.12, $
          'Total = ' + strtrim(str_max_scog, 2), CHARSIZE=1.0, CHARTHICK=1.2
  
  ; Add reference lines
  oplot, [r2pi*pix_size, r2pi*pix_size], plot_yrange_cog, LINESTYLE=1, THICK=1.5
  oplot, [r_max_arcsecs, r_max_arcsecs], plot_yrange_cog, LINESTYLE=0, THICK=1.5
  
  ; Add labels for reference lines
  xyouts, 1.05*r2pi*pix_size, plot_yrange_cog(0)*2.0, 'R(2!4p!3)', ORIENTATION=90, CHARSIZE=1.0, CHARTHICK=1.2
  xyouts, 1.10*r_max_arcsecs, plot_yrange_cog(0)*2.0, 'R!Dmax!N', ORIENTATION=90, CHARSIZE=1.0, CHARTHICK=1.2
  
  ; Add half-light radius line if requested
  if strlowcase(plot_reff) eq 'y' then begin
    ; Find half-light radius (radius containing half the total flux)
    half_flux = 0.5 * max(scog)
    ds = abs(scog - half_flux)
    ids = where(ds eq min(ds))
    
    ; Check if valid index was found
    if ids[0] ge 0 then begin
      r_half = r_arcsecs[ids[0]]
      oplot, [r_half, r_half], plot_yrange_cog, LINESTYLE=0, THICK=1.5
      xyouts, 1.10*r_half, plot_yrange_cog(0)*2.0, 'R!Deff!N', ORIENTATION=90, CHARSIZE=1.0, CHARTHICK=1.2
    endif
  endif
  
  ; Add mask radius line if requested
  if strlowcase(plot_rmask) eq 'y' then begin
    oplot, [rmask, rmask], plot_yrange_cog, LINESTYLE=2, THICK=1.5
    xyouts, 1.08*rmask, plot_yrange_cog(0)*2.0, 'R!Dmask!N', ORIENTATION=90, CHARSIZE=1.0, CHARTHICK=1.2
  endif
  
  ; Close the plot and restore original device if in PS mode
  if strlowcase(ps) eq 'y' then begin
    device, /CLOSE
    set_plot, orig_device
  endif
  
  return
end
