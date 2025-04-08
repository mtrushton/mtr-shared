pro plot_sb_profiles_dust
l
  ; Creates infrared surface brightness (SB) radial profile plots 
  ; for galaxy GAL_NAME and comparison models with irr1, irr2, 
  ; irr_HII components.
  
  ; Galaxy parameters
  gal_name = 'NGC628'
  d_gal = 9.5  ; galaxy distance in Mpc
  
  ; Model parameters
  tau = 12.4
  f = 0.1
  sfr = 2.30
  old = 0.18
  bd = 0.09
  
  ; Data file and plot settings
  file_data = 'N0628_data_par_v3.dat'
  plot_xrange = [0, 500]
  
  ; Define bands
  band_str = ['MIPS', 'PACS', 'PACS', 'PACS', 'SPIRE', 'SPIRE', 'SPIRE']
  nobs = n_elements(band_str)
  
  ; Initialise arrays
  lambda = intarr(nobs)
  sky = fltarr(nobs)
  xcen_dat = fltarr(nobs)
  ycen_dat = fltarr(nobs)
  ellipt = fltarr(nobs)
  pa = fltarr(nobs)
  pixel_scale = fltarr(nobs)
  psf_scale_arcsec = fltarr(nobs)
  
  ; Read data parameters
  gal_data_format = '(i3,f11.2,f10.3,f8.2,f7.2,f9.2,f9.3,f16.4)'
  openr, lun, file_data, /get_lun
  skip_lun, lun, 1, /lines
  readf, lun, ''
  
  ; Read data for each observation
  for n = 0, nobs-1 do begin
    readf, lun, dum1, dum2, dum3, dum4, dum5, dum6, dum7, dum8, format=gal_data_format
    lambda[n] = fix(dum1)
    sky[n] = dum2
    xcen_dat[n] = dum3
    ycen_dat[n] = dum4
    ellipt[n] = dum5
    pa[n] = dum6
    pixel_scale[n] = dum7
    psf_scale_arcsec[n] = dum8
    
    ; Process each band
    filter = band_str[n]
    wav = strcompress(lambda[n], /remove_all)
    band = filter + ' ' + wav + '!4l!3m'
    
    ; Dust map size
    map_size_pc = 68544.0
    map_size_pix = 1026L
    
    ; Data unit conversion - default is MJy/sr
    fluxconv = 1.0
    if(filter eq 'PACS') then fluxconv = 1.d-6 / (pixel_scale[n]/3600.*!pi/180.)^2
    
    ; Set up profiles and labels
    wa = 1  ; width of annuli in pixels
    xlab = 'arcsecs'
    ylab = 'MJy/sr'
    
    ; Model settings
    model = 'wd01'
    qyear = '06'
    scaabs = 'abs'
    
    ; Get model filename
    fname, model, qyear, tau, sfr, old, bd, scaabs, wav, namep
    
    ; Special case for MIPS
    if(filter eq 'MIPS' and lambda[n] eq 160) then namep = namep + '_mips'
    if(filter eq 'MIPS' and lambda[n] eq 70) then namep = namep + '_mips'
    
    ; Directories
    dir_base = '/nfs/d58/vlk/sedmodel/MTR/N0628/'
    dir_mod = dir_base + 'gn/PROFILE_RV/MAPS/'
    dir_dat = dir_base + 'N0628_' + filter + '_' + wav + '/'
    dir_mask = dir_dat + 'star_remove/'
    dir_psf = '../PROFILE/PSF/' + filter + '/'
    dir_plot = 'PLOTS/'
    
    ; File paths
    data_file = dir_dat + gal_name + '_' + filter + '_' + wav + 'um.fits'
    model_tot = dir_mod + 'map_' + namep + '_tot.fits'
    model_HII = dir_mod + 'map_' + namep + '_irr_HII.fits'
    model_irr1 = dir_mod + 'map_' + namep + '_irr1.fits'
    model_irr2 = dir_mod + 'map_' + namep + '_irr2.fits'
    model_irr3 = dir_mod + 'map_' + namep + '_irr3.fits'
    fmask = dir_mask + 'mask.fits'
    psf_file = dir_psf + 'PSF_' + filter + '_' + wav + '.fits'
    outname = 'lum_' + namep
    
    ; Calculate map scales
    map_scale_pc = map_size_pc / map_size_pix
    map_scale_rads = map_scale_pc / (d_gal * 1d6)
    map_scale_arcsecs = map_scale_rads * 180/!pi * 3600
    map_regrid_f = map_scale_arcsecs / pixel_scale[n]
    
    ; Read data and model files
    xdat = mrdfits(data_file, 0, xd)
    
    ; Handle NaN values in specific datasets
    if (filter eq 'PACS' or filter eq 'SPIRE' or $
        (filter eq 'MIPS' and lambda[n] eq 160)) then $
      xdat[where(finite(xdat) eq 0)] = 0.0d0
    
    ; Subtract sky
    xdat = xdat - sky[n]
    
    ; Read model components
    xmod_tot = mrdfits(model_tot, 0, xm)
    xmod_HII = mrdfits(model_HII, 0, xm)
    xmod_irr1 = mrdfits(model_irr1, 0, xm)
    xmod_irr2 = mrdfits(model_irr2, 0, xm)
    xmod_irr3 = mrdfits(model_irr3, 0, xm)
    
    ; Read mask
    mask = mrdfits(fmask, 0, xm)
    if(lambda[n] eq 23) then mask[*,*] = 1.
    
    ; Get dimensions
    nx = n_elements(xmod_tot[*,0])
    ny = n_elements(xmod_tot[0,*])
    
    nnx = nx * map_regrid_f
    nny = ny * map_regrid_f
    
    ; Define center pixels for model
    xcen_mod = nnx / 2L
    ycen_mod = nny / 2L
    
    ; PSF
    p = mrdfits(psf_file, 0, xp)
    nx_psf = n_elements(p[*,0])
    ny_psf = n_elements(p[0,*])
    
    ; Regrid to data pixel scale and normalise PSF
    nx_psf = nx_psf * psf_scale_arcsec[n] / pixel_scale[n]
    ny_psf = ny_psf * psf_scale_arcsec[n] / pixel_scale[n]
    
    p_regrid = congrid(p, nx_psf, ny_psf)
    p_regrid = p_regrid / total(p_regrid)
    
    ; Regrid models 
    rmod_tot = congrid(xmod_tot, nnx, nny)
    rmod_HII = congrid(xmod_HII, nnx, nny)
    rmod_irr1 = congrid(xmod_irr1, nnx, nny)
    rmod_irr2 = congrid(xmod_irr2, nnx, nny)
    rmod_irr3 = congrid(xmod_irr3, nnx, nny)
    
    ; Convolve models with PSF
    mod_conv_tot = convolve(rmod_tot, p_regrid, ft_psf=psf_ft)
    mod_conv_HII = convolve(rmod_HII, p_regrid, ft_psf=psf_ft)
    mod_conv_irr1 = convolve(rmod_irr1, p_regrid, ft_psf=psf_ft)
    mod_conv_irr2 = convolve(rmod_irr2, p_regrid, ft_psf=psf_ft)
    mod_conv_irr3 = convolve(rmod_irr3, p_regrid, ft_psf=psf_ft)
    
    ; Convert from W/Hz/pc^2/sr (dust maps) to MJy/sr
    f_conv = 1.d20 / 3.09d16^2
    mod_conv_tot *= f_conv
    mod_conv_HII *= f_conv
    mod_conv_irr1 *= f_conv
    mod_conv_irr2 *= f_conv
    mod_conv_irr3 *= f_conv
    
    ; Calculate sb profiles for data
    cog, xdat, mask, pa[n], ellipt[n], xcen_dat[n], ycen_dat[n], wa, rcog, bcog_dat, scog_dat, ncog
    
    ; Create uniform mask for model
    mask = fltarr(nnx, nny) + 1.
    
    ; Calculate sb profiles for models
    cog, mod_conv_tot, mask, pa[n], ellipt[n], xcen_mod, ycen_mod, wa, rcog_mod, bcog_mod_tot, scog_mod_tot, ncog
    cog, mod_conv_HII, mask, pa[n], ellipt[n], xcen_mod, ycen_mod, wa, rcog_mod, bcog_mod_HII, scog_mod_HII, ncog
    cog, mod_conv_irr1, mask, pa[n], ellipt[n], xcen_mod, ycen_mod, wa, rcog_mod, bcog_mod_irr1, scog_mod_irr1, ncog
    cog, mod_conv_irr2, mask, pa[n], ellipt[n], xcen_mod, ycen_mod, wa, rcog_mod, bcog_mod_irr2, scog_mod_irr2, ncog
    cog, mod_conv_irr3, mask, pa[n], ellipt[n], xcen_mod, ycen_mod, wa, rcog_mod, bcog_mod_irr3, scog_mod_irr3, ncog
    
    ; Apply flux conversion
    bcog_dat *= fluxconv
    scog_dat *= fluxconv
    
    ; Create plots
    mb = max(bcog_dat)
    plot_yrange = [mb*(-0.10), mb*1.3]
    
    plot_profile_residuals, rcog, rcog_mod, bcog_dat, bcog_mod_tot, bcog_mod_HII, bcog_mod_irr1, $
                  bcog_mod_irr2, bcog_mod_irr3, scog_dat, scog_mod_tot, scog_mod_HII, $
                  scog_mod_irr1, scog_mod_irr2, scog_mod_irr3, pixel_scale[n], $
                  map_scale_arcsecs, plot_xrange, plot_yrange, xlab, ylab, $
                  gal_name, band, outname, dir_plot
  endfor
  
  close, lun
  free_lun, lun
end
