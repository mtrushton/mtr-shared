pro plot_sed, lam_uvo, lam_ir, s_uvo, sig_uvo, s_ir, sig_ir, gal_name, sedplot_xrange, $
             sedplot_yrange, q_uv, q_opt, bfit_para, bd, thet_gal, sfr_ran, old_ran, $
             error, para_axl, parafind
;
; Plots the SED of the galaxy GAL_NAME, showing observed fluxes and best fitting model,
; with diffuse dust emission and star formation template components.
;
; INPUTS:
;   LAM_UVO    - Observed wavelengths in UV and optical ranges
;   LAM_IR     - Observed wavelengths in IR/Sub-mm ranges
;   S_UVO      - Observed fluxes in UV and optical ranges
;   SIG_UVO    - Uncertainties in UV and optical fluxes
;   S_IR       - Observed fluxes in IR/Sub-mm ranges
;   SIG_IR     - Uncertainties in IR fluxes
;   GAL_NAME   - String identifier of the galaxy
;   SEDPLOT_XRANGE - X-axis range for plot [min, max]
;   SEDPLOT_YRANGE - Y-axis range for plot [min, max]
;   Q_UV       - Attenuation correction factors for UV data
;   Q_OPT      - Attenuation correction factors for optical data
;   BFIT_PARA  - Vector of best fitting parameters
;   BD         - Intrinsic bulge-disk ratio
;   THET_GAL   - Intrinsic B-band scalelength of old stellar disk
;   SFR_RAN    - SFR parameter range
;   OLD_RAN    - OLD parameter range
;   ERROR      - Vector of errors in free parameters
;   PARA_AXL   - Vector of parameter labels
;   PARAFIND   - Integer vector (0=fixed, 1=free parameter)
;
; OUTPUT:
;   SED plot saved as GAL_NAME_SED.ps
;
  @common_model
  @common_kern_ir
  
  ; Constants
  HREF = 5.67*3.09d19
  WAVE_EXT = 4.5d0  ; Extrapolation wavelength (microns)
  
  ; Set up color table
  tvlct, 0, 0, 0, 0     ; black
  tvlct, 0, 0, 255, 1   ; blue
  tvlct, 255, 0, 0, 2   ; red
  
  ; Convert observed UV/optical fluxes from mJy to Jy
  s_uvo = s_uvo/1000.0
  sig_uvo = sig_uvo/1000.0
  
  ; Combine observed data
  lambda = [lam_uvo, lam_ir]       ; wavelengths
  flux = [s_uvo, s_ir]             ; fluxes
  sigma = [sig_uvo, sig_ir]        ; uncertainties
  q = [q_uv, q_opt]/fpimpc2        ; deattenuations
  
  ; Check for data outside plotting range
  if (min(lambda) lt sedplot_xrange[0] || max(lambda) gt sedplot_xrange[1] || $
      min(flux) lt sedplot_yrange[0] || max(flux) gt sedplot_yrange[1]) then $
    print, 'SEDPLOT: warning - data exists outside the plotting range'
  
  ; Deattenuate the UV/optical fluxes
  s_uvo_q = s_uvo*q
  
  ; Scale the star formation template
  template = template_ir * (bfit_para[1]*clocal1 + clocal2) * bfit_para[2]
  
  ; Prepare for the stellar light extrapolation
  wave_direct = wave_sedmodel[where(wave_sedmodel ge WAVE_EXT)]
  template = template[where(wave_sedmodel ge WAVE_EXT)]
  n = where(lam_uvo eq WAVE_EXT)
  
  ; λ^(-2) extrapolation of stellar light
  s_model_direct = s_uvo[n] * lam_uvo[n]^2
  s_model_direct = s_model_direct # (1.0d0 / wave_direct^2)
  
  ; Calculate diffuse dust emission
  test, bfit_para[0], bd
  ierr = 0
  test2, d_gal_ran, sfr_ran, old_ran, thet_gal, bfit_para[1], bfit_para[2], bfit_para[3], ierr
  diff = model_diffuse * thet_gal^2 * 1.0d46 / (4*!pi*HREF^2)
  diff = diff[where(wave_sedmodel ge WAVE_EXT)]
  
  ; Create the plot
  set_plot, 'PS'
  device, file=gal_name+'_SED.ps'
  
  ; Plot observed fluxes with error bars
  plot, lambda, flux, xst=1, yst=1, xrange=sedplot_xrange, yrange=sedplot_yrange, $
        /xlog, /ylog, psym=4, xtitle='!4k!3 [!4l!3m]', ytitle='F!l!4m!3!n [Jy]', /noerase
  oploterror, lambda, flux, sigma, psym=4
  
  ; Plot stellar components
  oplot, lam_uvo, s_uvo, color=1                 ; attenuated stellar light
  oplot, lam_uvo, s_uvo_q, color=1, linestyle=1  ; deattenuated stellar light
  oplot, wave_direct, s_model_direct, color=1    ; extrapolated stellar light
  
  ; Plot dust emission components
  oplot, wave_direct, template, color=2, linestyle=1  ; star formation regions
  oplot, wave_direct, diff, color=2, linestyle=1      ; diffuse dust emission
  
  ; Plot total model
  model = s_model_direct + template + diff
  oplot, wave_direct, model, color=0
  oplot, lam_uvo, s_uvo, color=0
  
  ; Annotations
  x_text_pos = exp(alog(sedplot_xrange[0]) + (alog(sedplot_xrange[1])-alog(sedplot_xrange[0]))*0.1)
  y_text_pos = exp(alog(sedplot_yrange[1]) - (alog(sedplot_yrange[1])-alog(sedplot_yrange[0]))*0.1)
  
  xyouts, x_text_pos, y_text_pos, gal_name
  
  y_step = 0.05d0
  y_pos = 0.1d0
  
  for i=0, n_elements(bfit_para)-1 do begin
    y_pos += y_step
    y_text_pos = exp(alog(sedplot_yrange[1]) - (alog(sedplot_yrange[1])-alog(sedplot_yrange[0]))*y_pos)
    
    if (parafind[i] eq 0) then begin
      ; Fixed parameter
      xyouts, x_text_pos, y_text_pos, para_axl[i]+'='+string(bfit_para[i], format='(F4.2)')
    endif else begin
      ; Free parameter with error
      xyouts, x_text_pos, y_text_pos, para_axl[i]+'='+string(bfit_para[i], format='(F4.2)')+ $
              string('!9+!3')+string(error[i], format='(F4.2)')
    endelse
  endfor ; i=0, n_elements(bfit_para)-1
  
  device, /close
  return
end
