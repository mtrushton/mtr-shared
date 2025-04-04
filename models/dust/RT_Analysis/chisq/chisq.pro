pro chisq

; Main procedure to calculate chisq for multiple RT modelled wavelength bands
; Two procedures are called as follows:
; CHISQ_INPUT to input parameters and file names
; CHISQ_CALC to perform calculations
; Outputs are chisq for each waveband and a 'global' chisq for all bands combined  

; Declare variables to hold input parameters
wa = 0.0
EBV = 0.0
dir_data = ''
dir_mask = ''
dir_model = ''
dir_psf = ''
lambda = dblarr(1)
n_obs = 0L
data_file = strarr(1)
mask_file = strarr(1)
sky = dblarr(1)
model_file = strarr(1)
psf_file = strarr(1)
pixel_scale = fltarr(1)
psf_pixel_scale = fltarr(1)
xcen = fltarr(1)
ycen = fltarr(1)
ellipt = 0.0
pa = 0.0
map_scale_arcsecs_uvopt = 0.0
map_scale_arcsecs_ir = 0.0
flux_conv = fltarr(1)
col_corr = fltarr(1)
rmax = fltarr(1)
model_scale_factor = fltarr(1)

; Get input parameters
chisq_input, wa, EBV, dir_data, dir_mask, dir_model, dir_psf, lambda, n_obs, $
             data_file, mask_file, sky, model_file, psf_file, pixel_scale, $
             psf_pixel_scale, xcen, ycen, ellipt, pa, map_scale_arcsecs_uvopt, $
             map_scale_arcsecs_ir, flux_conv, col_corr, rmax, model_scale_factor

; Calculate chisq values
chisq_values = dblarr(n_obs)
chisq_calc, n_obs, wa, lambda, EBV, sky, dir_data, dir_mask, dir_model, dir_psf, $
            data_file, mask_file, model_file, psf_file, pixel_scale, psf_pixel_scale, $
            xcen, ycen, ellipt, pa, map_scale_arcsecs_uvopt, map_scale_arcsecs_ir, $
            flux_conv, col_corr, rmax, model_scale_factor, chisq_values

; Print a summary of chi-square values for each waveband
print, '========== CHI-SQUARE SUMMARY =========='
print, 'WAVEBAND     LAMBDA     CHI-SQUARE'
print, '-----------------------------------'
for i = 0, n_obs - 1 do begin
    print, string(i+1, format='(I2)'), '         ', $
           string(lambda[i], format='(F7.3)'), '      ', $
           string(chisq_values[i], format='(F12.4)')
endfor
print, '-----------------------------------'
print, 'TOTAL CHI-SQUARE: ', total(chisq_values)
print, '========================================'

end
