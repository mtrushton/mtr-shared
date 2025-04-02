forward_function calculate_totals
pro create_extc
;
; Creates extinction curves for a galaxy modelled with the RT analysis. The
; maps for the observed galaxy emission are ratioed with the corresonding 
; intrinsic maps to determine the flux ratio of the attenuated to unattenuated
; light. The magnitude of attenuation for each specific band is then calculated
; from the Pogson's relation.
;
; It is customary to present an extinction curve as a plot of extinction  
; magnitude against wavelength. This is normalised with respect to a reference
; band, usually B or V. 
;
; The extinction curve is determined for the whole galaxy from COG photometry on 
; the models for the observed and intrinsic emission (COG program used here).
; COG photometry is also performed out to certain wavelengths (RIN and REF1 below)
; and extinction curves plotted for those wavelengths. This is done to investigate
; the dependence of extinction on galactocentric radius.
;
; This can be adapted to any modelled galaxy by changing the inputs below and the 
; file suffixes. Extra components can be added or removed, as necessary. 
;
  ; INPUTS:
  wa = 1 ; width in pixels of annuli used for averaging over galactocentric radius
  i = 14 ; inclination angle of the galaxy
  ellipt = cos(i * !pi / 180) ; ellipticity of the disk
  pa = 14.33 ; postion angle
  dist = 17.9 ; galaxy distance in MPC
  pix_pc = 2.4 / 3600. * !pi / 180. * dist * 1.d6  
  rin = 2.1d3 ; the inner radius (an attenutation curve is plotted at this radius)
  ref1 = 14.0d3 ; a reference radius in kpc (an attenuation is also plotted at this radius)  
  dir = '/nfs/d58/vlk/sedmodel/RTOMA/NUrad/out/' ; directory with the models
  norm_wavelength = 4400. ; reference wavelength for normalising extinction curve
  
  ; Sampled wavelengths in the RT analysis (stellar emission)
  lambda = [912., 1350., 1500., 1650., 2000., 2200., 2500., 2800., 3650., 4400., 5640.]
  
  ; Define filters for the non-UV stellar emission included in the RT analysis
  filter_opt = ['uSDSS', 'uv36', 'b', 'g', 'v', 'rSDSS', 'iSDSS', 'i', 'z', 'j', 'h', 'k', $
   'ir34', 'ir36', 'ir45', 'ir46', 'ir57']
  
  ; Short wavelength filters used for the young stellar components in addition to FILTER_OPT
   filter_uv = ['uv09', 'uv13', 'uv15', 'uv16', 'uv20', 'uv22', 'uv25', 'uv28']

    ; Define file suffixes
  suffix_obs = '_wd01_q06_i14_t20_hd12000_zd160_hd1_3200_zd1_90'
  suffix_intr = '_wd01_q06_i14_t0_hd12000_zd160_hd1_3200_zd1_90'
  
  ; Since we are only interested in the extinction curve out to B (or V), we restrict the calculation to the
  ; bands defined below
  bv_files = ['uv36', 'b', 'v']
  uv_files = ['uv09', 'uv13', 'uv15', 'uv16', 'uv20', 'uv22', 'uv25', 'uv28', 'uv36', 'b', 'v']
  ; -----------------------------------------------------------------------------------------------------------
  
  nopt = n_elements(bv_files)
  nuv = n_elements(uv_files)
  
  ; Initialise arrays
  tot_b_obs = dblarr(nopt)
  tot_d_obs = dblarr(nopt)
  tot_td_obs = dblarr(nuv)
  tot_b_intr = dblarr(nopt)
  tot_d_intr = dblarr(nopt)
  tot_td_intr = dblarr(nuv)
  
  tot_b_obs_rin = dblarr(nuv)
  tot_d_obs_rin = dblarr(nuv)
  tot_td_obs_rin = dblarr(nuv)
  tot_b_intr_rin = dblarr(nuv)
  tot_d_intr_rin = dblarr(nuv)
  tot_td_intr_rin = dblarr(nuv)
  
  tot_b_obs_ref1 = dblarr(nuv)
  tot_d_obs_ref1 = dblarr(nuv)
  tot_td_obs_ref1 = dblarr(nuv)
  tot_b_intr_ref1 = dblarr(nuv)
  tot_d_intr_ref1 = dblarr(nuv)
  tot_td_intr_ref1 = dblarr(nuv)
  
  ; Create mask
  mask = dblarr(1, 1) 
  
  ; Process young stellar components
 
  ; Process observed td files
  for n = 0, nuv - 1 do begin
    file = dir + 'map_mtd_' + uv_files[n] + suffix_obs + '_hs1_3200_zs1_90_sca.fits'
    process_file, file, mask, pa, ellipt, wa, pix_pc, rin, ref1, tot_td_obs, tot_td_obs_rin, tot_td_obs_ref1, n
  endfor ;  n = 0, nuv - 1
  
  ; Process intrinsic td files
  for n = 0, nuv - 1 do begin
    file = dir + 'map_mtd_' + uv_files[n] + suffix_intr + '_hs1_3200_zs1_90_abs.fits'
    process_file, file, mask, pa, ellipt, wa, pix_pc, rin, ref1, tot_td_intr, tot_td_intr_rin, tot_td_intr_ref1, n
  endfor ;  n = 0, nuv - 1
  
 ; Process old stellar components
  
  ; Process observed d files
  for n = 0, nopt - 1 do begin
    file = dir + 'map_md_' + bv_files[n] + suffix_obs + '_hs3000_zs190_sca.fits'
    process_file, file, mask, pa, ellipt, wa, pix_pc, rin, ref1, tot_d_obs, tot_d_obs_rin, tot_d_obs_ref1, n
  endfor ;  n = 0, nopt - 1
  
  ; Process intrinsic d files
  for n = 0, nopt - 1 do begin
    file = dir + 'map_md_' + bv_files[n] + suffix_intr + '_hs3000_zs190_abs.fits'
    process_file, file, mask, pa, ellipt, wa, pix_pc, rin, ref1, tot_d_intr, tot_d_intr_rin, tot_d_intr_ref1, n
  endfor ;  n = 0, nopt - 1
  
  ; Process observed b files
  for n = 0, nopt - 1 do begin
    file = dir + 'map_mb_' + bv_files[n] + suffix_obs + '_reff528_ell94_n1_sca.fits'
    process_file, file, mask, pa, ellipt, wa, pix_pc, rin, ref1, tot_b_obs, tot_b_obs_rin, tot_b_obs_ref1, n
  endfor ;  n = 0, nopt - 1
  
  ; Process intrinsic b files
  for n = 0, nopt - 1 do begin
    file = dir + 'map_mb_' + bv_files[n] + suffix_intr + '_reff528_ell94_n1_abs.fits'
    process_file, file, mask, pa, ellipt, wa, pix_pc, rin, ref1, tot_b_intr, tot_b_intr_rin, tot_b_intr_ref1, n
  endfor ;  n = 0, nopt - 1
  
  ; Calculate totals and ratios
  tot_obs = calculate_totals(tot_td_obs, tot_b_obs, tot_d_obs, nuv)
  tot_intr = calculate_totals(tot_td_intr, tot_b_intr, tot_d_intr, nuv)
  tot_obs_rin = calculate_totals(tot_td_obs_rin, tot_b_obs_rin, tot_d_obs_rin, nuv)
  tot_intr_rin = calculate_totals(tot_td_intr_rin, tot_b_intr_rin, tot_d_intr_rin, nuv)
  tot_obs_ref1 = calculate_totals(tot_td_obs_ref1, tot_b_obs_ref1, tot_d_obs_ref1, nuv)
  tot_intr_ref1 = calculate_totals(tot_td_intr_ref1, tot_b_intr_ref1, tot_d_intr_ref1, nuv)
  
  ; Calculate attenuation and normalised attenuation
  ratio = tot_intr / tot_obs
  atten = -2.5 * alog10(tot_obs / tot_intr)
  atten_rin = -2.5 * alog10(tot_obs_rin / tot_intr_rin)
  atten_ref1 = -2.5 * alog10(tot_obs_ref1 / tot_intr_ref1)
  
  ; Normalise the data to a specific band/wavelength
  ; data point closest to reference wavelength: B Band (4400 Angstrom)
  ; some studies use V Band (5500 Angstrom)
  diff = abs(lambda - norm_wavelength)
  ind = where(diff eq min(diff))
  
  atten_norm = atten # (1. / atten[ind])
  atten_norm_rin = atten_rin # (1. / atten_rin[ind]) 
  atten_norm_ref1 = atten_ref1 # (1. / atten_ref1[ind])
  
  ; Plot results
  plot, lambda, atten_norm
  oplot, lambda, atten_norm_rin, linestyle=1
  plot, lambda, atten_norm_ref1, linestyle=1
  
  ; Save results
  save, lambda, atten_norm, atten_norm_rin, atten_norm_ref1, ratio, file='atten.sav'
  
end

; Helper function to process files
pro process_file, file, mask, pa, ellipt, wa, pix_pc, rin, ref1, tot_arr, tot_arr_rin, tot_arr_ref1, index
  x = mrdfits(file, 0, nh)
  
  ; Initialise mask if needed
  if n_elements(mask) eq 1 then begin
    nx = n_elements(x[*, 0])
    ny = n_elements(x[0, *])
    mask = dblarr(nx, ny)
    mask[*, *] = 1.0
  endif ; n_elements(mask) eq 1 
  
  ; Process data
  x_cen = 0
  y_cen = n_elements(x[0, *]) / 2
  cog, x, mask, pa, ellipt, x_cen, y_cen, wa, rcog, bcog, ncog, scog, r2pi
  
  ; Store results
  tot_arr[index] = max(scog)
  
  ; Calculate index rin
  ii = where(rcog * pix_pc lt rin)
  tot_arr_rin[index] = scog[max(ii)]
  
  ; Calculate index ref1
  ii = where(rcog * pix_pc gt rin and rcog * pix_pc le ref1)
  tot_arr_ref1[index] = scog[max(ii)] - scog[min(ii)]
end ; process_file

; Calculate totals
function calculate_totals, tot_td, tot_b, tot_d, nuv
  tot = dblarr(nuv)
  
  for n = 0, nuv - 1 do begin
    tot[n] = tot_td[n]
    if n gt 7 then tot[n] += tot_b[n-8] + tot_d[n-8]
  endfor
  
  return, tot
end ; Calculate totals
