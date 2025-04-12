pro chisq_input, wa, EBV, dir_data, dir_mask, dir_model, dir_psf, lambda, n_obs, $
                data_file, mask_file, sky, model_file, psf_file, pixel_scale, psf_pixel_scale, $
                xcen, ycen, ellipt, pa, map_scale_arcsecs_uvopt, map_scale_arcsecs_ir, $
                flux_conv, col_corr, rmax, model_scale_factor
;--------------------------------------------------------------------------------------------------
; INPUTS for the chisq program

; Galaxy parameters
d_gal = 9.5  ; Galaxy distance in Mpc
EBV = 0.254 - 0.192  ; Color excess E(B-V)
wa = 1.0  ; Radial profile width parameter

; Directory paths
dir_data = '/nfs/d58/vlk/sedmodel/MTR/N0628/data/'
dir_mask = '/nfs/d58/vlk/sedmodel/MTR/N0628/mask/'
dir_model = './std_model/'
dir_psf = 'PSF/'

; Wavelength bands to process
lambda = [0.1516, 0.2267, 0.3543, 0.4770, 0.6231, 0.7625, 0.9134, 1.2590, 1.6620, $
          2.200, 24.0, 70.0, 100.0, 160.0, 250.0, 350.0, 500.0]
;lambda = 2.2  ; If you want to run just one band, use this line instead

; Number of observations
n_obs = n_elements(lambda)

; Initialise arrays for each observation
data_file = strarr(n_obs)
mask_file = strarr(n_obs)
sky = dblarr(n_obs)
model_file = strarr(n_obs)
psf_file = strarr(n_obs)
pixel_scale = fltarr(n_obs)
psf_pixel_scale = fltarr(n_obs)
xcen = fltarr(n_obs)
ycen = fltarr(n_obs)
flux_conv = fltarr(n_obs)
col_corr = fltarr(n_obs)
rmax = fltarr(n_obs)
model_scale_factor = fltarr(n_obs)

;--------------------------------------------------------------------------------------------------
; Define observation-specific parameters using arrays instead of individual assignments
; This allows processing all bands at once

; Data file names
data_files = ['NGC628_GALEX_FUV.fits', 'NGC628_GALEX_NUV.fits', 'NGC628_SDSS_U.fits', $
              'NGC628_SDSS_G.fits', 'NGC628_SDSS_R.fits', 'NGC628_SDSS_I.fits', $
              'NGC628_SDSS_Z.fits', 'NGC628_2MASS_J.fits', 'NGC628_2MASS_H.fits', $
              'NGC628_2MASS_K.fits', 'NGC628_MIPS_24um.fits', 'NGC628_PACS_70um.fits', $
              'NGC628_PACS_100um.fits', 'NGC628_PACS_160um.fits', 'NGC628_SPIRE_250um.fits', $
              'NGC628_SPIRE_350um.fits', 'NGC628_SPIRE_500um.fits']

; Mask file names
mask_files = ['NGC628_MASK_GALEX_FUV.fits', 'NGC628_MASK_GALEX_NUV.fits', 'NGC628_MASK_SDSS_U.fits', $
              'NGC628_MASK_SDSS_G.fits', 'NGC628_MASK_SDSS_R.fits', 'NGC628_MASK_SDSS_I.fits', $
              'NGC628_MASK_SDSS_Z.fits', 'NGC628_MASK_2MASS_J.fits', 'NGC628_MASK_2MASS_H.fits', $
              'NGC628_MASK_2MASS_K.fits', 'NGC628_MASK_MIPS_24.fits', 'NGC628_MASK_PACS_70.fits', $
              'NGC628_MASK_PACS_100.fits', 'NGC628_MASK_PACS_160.fits', 'NGC628_MASK_SPIRE_250.fits', $
              'NGC628_MASK_SPIRE_350.fits', 'NGC628_MASK_SPIRE_500.fits']

; Model file names
model_files = ['map_total_FUV_wd01_q06_i0_t128_sca.fits', 'map_total_NUV_wd01_q06_i0_t128_sca.fits', $
               'map_total_uSDSS_wd01_q06_i0_t128_sca.fits', 'map_total_gSDSS_wd01_q06_i0_t128_sca.fits', $
               'map_total_rSDSS_wd01_q06_i0_t128_sca.fits', 'map_total_iSDSS_wd01_q06_i0_t128_sca.fits', $
               'map_total_zSDSS_wd01_q06_i0_t128_sca.fits', 'map_total_j_wd01_q06_i0_t128_sca.fits', $
               'map_total_h_wd01_q06_i0_t128_sca.fits', 'map_total_k_wd01_q06_i0_t128_sca.fits', $
	       'map_wd01_q06_t128_s230_no18_bd9_hd7800_zd140_hd1_2800_zd1_90_hs3000_zs215_hs1_2800_zs1_90_reff800_ell60_sca_l24um_tot.fits', $
               'map_wd01_q06_t128_s230_no18_bd9_hd7800_zd140_hd1_2800_zd1_90_hs3000_zs215_hs1_2800_zs1_90_reff800_ell60_sca_l70um_tot.fits', $
               'map_wd01_q06_t128_s230_no18_bd9_hd7800_zd140_hd1_2800_zd1_90_hs3000_zs215_hs1_2800_zs1_90_reff800_ell60_sca_l100um_tot.fits', $
               'map_wd01_q06_t128_s230_no18_bd9_hd7800_zd140_hd1_2800_zd1_90_hs3000_zs215_hs1_2800_zs1_90_reff800_ell60_sca_l160um_tot.fits', $
               'map_wd01_q06_t128_s230_no18_bd9_hd7800_zd140_hd1_2800_zd1_90_hs3000_zs215_hs1_2800_zs1_90_reff800_ell60_sca_l250um_tot.fits', $
               'map_wd01_q06_t128_s230_no18_bd9_hd7800_zd140_hd1_2800_zd1_90_hs3000_zs215_hs1_2800_zs1_90_reff800_ell60_sca_l350um_tot.fits', $
               'map_wd01_q06_t128_s230_no18_bd9_hd7800_zd140_hd1_2800_zd1_90_hs3000_zs215_hs1_2800_zs1_90_reff800_ell60_sca_l500um_tot.fits']

; PSF file names
psf_files = ['PSF_GALEX_FUV.fits', 'PSF_GALEX_NUV.fits', 'PSF_SDSS_U.fits', $
             'PSF_SDSS_G.fits', 'PSF_SDSS_R.fits', 'PSF_SDSS_I.fits', $
             'PSF_SDSS_Z.fits', 'PSF_2MASS_J.fits', 'PSF_2MASS_H.fits', $
             'PSF_2MASS_K.fits', 'PSF_MIPS_24.fits', 'PSF_PACS_70.fits', $
             'PSF_PACS_100.fits', 'PSF_PACS_160.fits', 'PSF_SPIRE_250.fits', $
             'PSF_SPIRE_350.fits', 'PSF_SPIRE_500.fits']

; Pixel scales in arcseconds
pixel_scales = [1.5, 1.5, 0.396, 0.396, 0.396, 0.396, 0.396, 1.0, 1.0, $
                1.0, 1.5, 1.4, 1.7, 2.85, 6.0, 10.0, 14.0]

; PSF pixel scales in arcseconds
psf_pixel_scales = [1.5, 1.5, 0.396, 0.396, 0.396, 0.396, 0.396, 1.0, 1.0, $
                   1.0, 0.498, 0.3996, 0.3996, 0.3996, 0.599, 0.599, 0.599]

; X centre coordinates
xcens = [840.0, 840.0, 470.0, 470.0, 470.0, 470.0, 470.0, 700.0, 700.0, $
         700.0, 302.710, 562.0, 468.0, 238.0, 203.0, 121.27, 87.73]

; Y centre coordinates
ycens = [840.0, 840.0, 783.0, 783.0, 783.0, 783.0, 783.0, 700.0, 700.0, $
         700.0, 288.60, 592.0, 488.0, 300.0, 215.0, 133.30, 95.87]

; Sky background values
sky_values = [4.01d-4, 3.45d-3, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, $
              0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

; Colour correction factors
col_corrs = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, $
             1.0, 0.93, 0.933, 0.972, 0.948, 1.0, 1.0, 1.0]

; Maximum radius values (radius in pixels beyond which there is no galaxy emission)
rmaxs = [270.0, 270.0, 270.0, 270.0, 270.0, 270.0, 290.0, 290.0, 240.0, $
         240.0, 240.0, 240.0, 240.0, 270.0, 270.0, 270.0, 300.0]

; Assign parameters to the output arrays
for i = 0, n_obs - 1 do begin
    data_file[i] = data_files[i]
    mask_file[i] = mask_files[i]
    model_file[i] = model_files[i]
    psf_file[i] = psf_files[i]
    pixel_scale[i] = pixel_scales[i]
    psf_pixel_scale[i] = psf_pixel_scales[i]
    xcen[i] = xcens[i]
    ycen[i] = ycens[i]
    sky[i] = sky_values[i]
    col_corr[i] = col_corrs[i]
    rmax[i] = rmaxs[i]
endfor ; i = 0, n_obs - 1

; Set galaxy geometry parameters
ellipt = 0.94
pa = -66.08

; Set map scale parameters
map_scale_arcsecs_uvopt = 1.0
map_scale_arcsecs_ir = 68544 / (1026 * d_gal * 1.d6) * 3600 * 180 / !pi

; Calculate flux conversion factors
for i = 0, n_obs - 1 do begin
    if i le 6 then begin
        ; SDSS and GALEX bands
        flux_conv[i] = 3.631d-6 * 1.d-6 / (pixel_scale[i] * !pi / (3600 * 180.))^2
        ; Special cases for GALEX
        if i eq 0 then flux_conv[i] = 1.0724d-4 * 1.d-6 / (pixel_scale[i] * !pi / (3600 * 180.))^2
        if i eq 1 then flux_conv[i] = 3.529d-5 * 1.d-6 / (pixel_scale[i] * !pi / (3600 * 180.))^2
    endif else if i ge 7 and i le 9 then begin
        ; 2MASS bands
        if i eq 7 then flux_conv[i] = 7.8285d-12 / (pixel_scale[i] * !pi / (3600 * 180.))^2
        if i eq 8 then flux_conv[i] = 5.7056e-12 / (pixel_scale[i] * !pi / (3600 * 180.))^2
        if i eq 9 then flux_conv[i] = 6.5267e-12 / (pixel_scale[i] * !pi / (3600 * 180.))^2
    endif else if i eq 10 then begin
        ; MIPS 24um
        flux_conv[i] = 1.0
    endif else if i ge 11 and i le 13 then begin
        ; PACS bands
        flux_conv[i] = 1.d-6 / (pixel_scale[i] * !pi/ (3600 * 180.))^2
    endif else begin
        ; SPIRE bands
        flux_conv[i] = 1.0
    endelse
endfor ; i = 0, n_obs - 1

; Calculate model scale factors
for i = 0, n_obs - 1 do begin
    if lambda[i] lt 3.0 then begin
        ; UV/optical/near-IR bands
        model_scale_factor[i] = 1. / (map_scale_arcsecs_uvopt^2 * 2.3d-11 * 1.d6)
    endif else begin
        ; IR bands
        model_scale_factor[i] = 1.d20 / 3.09d16^2
    endelse
endfor ; i = 0, n_obs - 1 

return
end
