pro chisq_calc, n_obs, wa, lambda, EBV, sky, dir_data, dir_mask, dir_model, dir_psf, $
                data_file, mask_file, model_file, psf_file, pixel_scale, psf_pixel_scale, $
                xcen, ycen, ellipt, pa, map_scale_arcsecs_uvopt, map_scale_arcsecs_ir, $
                flux_conv, col_corr, rmax, model_scale_factor, chisq_values
;---------------------------------------------------------------------------------------------------------------------

; Initialise arrays for file paths
data_f = strarr(n_obs)
mask_f = strarr(n_obs)
model_f = strarr(n_obs)
psf_f = strarr(n_obs)

; Initialise chisq array to store values for each observation
chisq_values = dblarr(n_obs)

; Process each observation
for n = 0, n_obs - 1L do begin
    ; Construct file paths
    data_f[n] = dir_data + data_file[n]
    mask_f[n] = dir_mask + mask_file[n]  
    model_f[n] = dir_model + model_file[n] 
    psf_f[n] = dir_psf + psf_file[n]
    
    print, 'Processing waveband:', lambda[n], ' μm'
    print, 'Model file:', model_f[n]

    ; Read data files
    data = mrdfits(data_f[n], 0, xm)
    mask = mrdfits(mask_f[n], 0, xm)
    model = mrdfits(model_f[n], 0, xm)
    psf = mrdfits(psf_f[n], 0, xm)
    
    ; Subtract sky background
    data_sky_sub = data - sky[n]
    
    ; Apply dereddening correction for UV and optical bands
    if(lambda[n] lt 3.0) then dered, lambda[n], EBV, data_sky_sub, data

    ; Get dimensions
    nx = n_elements(model[*, 0])
    ny = n_elements(model[0, *])
    
    nx_psf = n_elements(psf[*, 0])
    ny_psf = n_elements(psf[0, *])

    ; Set regridding factor based on wavelength
    if(lambda[n] lt 3.0) then $
        regrid_factor = map_scale_arcsecs_uvopt / pixel_scale[n] $
    else $
        regrid_factor = map_scale_arcsecs_ir / pixel_scale[n]
    
    ; Calculate regridded dimensions
    nx_reg = nx * regrid_factor
    ny_reg = ny * regrid_factor

    ; Set model centre coordinates
    xcen_model = nx_reg / 2.
    ycen_model = ny_reg / 2.

    ; Account for the 'half maps' in the uv/optical
    if(lambda[n] lt 3.0) then xcen_model = 0.

    ; Calculate PSF regridding dimensions
    nx_psf_reg = nx_psf * psf_pixel_scale[n] / pixel_scale[n]
    ny_psf_reg = ny_psf * psf_pixel_scale[n] / pixel_scale[n]

    ; Regrid model
    model_reg = congrid(model, nx_reg, ny_reg)

    ; Regrid PSF if needed
    psf_reg = psf
    if(pixel_scale[n] ne psf_pixel_scale[n]) then $ 
        psf_reg = congrid(psf, nx_psf_reg, ny_psf_reg)
    
    ; Normalise PSF
    psf_reg_tot = total(psf_reg)
    psf_reg_norm = psf_reg / psf_reg_tot

    ; Convolve model with PSF (using Gaussian smoothing as a substitute)
    model_reg_conv = gauss_smooth(model_reg, 1.5, /edge_truncate)
    model_reg_conv = model_scale_factor[n] * model_reg_conv

    ; Create model mask
    mask_model = intarr(nx_reg, ny_reg)
    mask_model[*, *] = 1L

    ; Apply flux conversion
    flux = data * flux_conv[n] / col_corr[n]

    ; Calculate radial profiles
    cog, flux, mask, pa, ellipt, xcen[n], ycen[n], wa, rcog_data, bcog_data, ncog, scog_data, r2pi
    cogsig, flux, mask, pa, ellipt, xcen[n], ycen[n], wa, rcog_data, bcog_data, ncog, sig
    cog, model_reg_conv, mask_model, pa, ellipt, xcen_model, ycen_model, wa, rcog_model, bcog_model, ncog, scog_model, r2pi

    ; Create diagnostic plot if desired (comment out for batch processing)
    ;set_plot, 'PS'
    ;device, file='lum_plot_'+string(lambda[n], format='(F6.3)')+'.ps', /color
    ;plot, rcog_model, bcog_model, /ylog, xrange=[0, max(rmax)], yrange=[0.00002, 10], $
    ;      xtitle='Radius', ytitle='Surface Brightness', title='λ = '+string(lambda[n], format='(F6.3)')+' μm'
    ;oplot, rcog_data, bcog_data, color=2
    ;device, /close

    ; Get array sizes for comparison
    nr_obs = n_elements(rcog_data)
    nr_mod = n_elements(rcog_model)
    nr = nr_obs < nr_mod  ; Use the smaller of the two

    ; Calculate chisq for this waveband
    chisq = 0.0
    nn = 0L
    for i = 0, nr - 1L do begin
        if(i gt 10 and rcog_data[i] lt rmax[n] and sig[i] gt 0.) then begin
            chisq = chisq + (bcog_data[i] - bcog_model[i])^2 / sig[i]^2
            nn = nn + 1
        endif
    endfor ; i = 0, nr - 1L
    
    ; Store chisq value for this waveband
    chisq_values[n] = chisq
    
    ; Print results for this waveband
    print, 'λ =', lambda[n], 'μm, Chi-square =', $
           string(chisq, format='(F12.4)'), ' (', nn, ' points)'
    print, '-----------------------------------'
    
endfor ; n = 0, n_obs - 1L 

return
end
