pro plot_sb_profiles_stellar
;
; Plots the UV/OPT/NIR surface brightness (SB) profiles of galaxy 
; GAL_NAME, along with the best fitting model and model components 
; for a galaxy with SB profiles fitted by the PT11 RT code. 

; Includes residuals in a panel below the main plot. These are best 
; plotted as far as the noise level (see RESIDUAL_RESTRICT vector).
;
; Note that the UV only includes young disks, whereas the models for
; longer wavelength data include a bulge component, as well as the
; old stellar disk. At wavelengths (>3 micron) it is necessay to also
; include a dust component, as the dust continuum emission and PAHs
; can make a significant contribution in that region.
; 
; The procedure can be summarised as follows:
; The data is read and corrected for interstellar extinction. A 
; conversion factor is applied ensuring that the pixels are in MJy/sr.
;
; The models are then read and regridded onto the data pixel scale.
; The SB profiles for the data and model are derived using the COG
; program. 
;
; PSFs for each filter are read and regridded onto the data scale. The
; models are then convolved with the regridded PSFs.
;
; INPUTS:
; GALNAME (str): Galaxy name
; D_GAL (float) : Distance of the galaxy in MPC
; EBV (float) : Colour excess (interstellar reddening)
; TAU (float) : Best-fitting value of TAU at R=0
; F (float) : Best-fitting 'F factor'
; BD (float) : Bulge-to-disk ratio
; NSERSIC (INT) : Sersic index of the bulge
; FILTER (str) : vector of filters for the observational data
; LAMBDA (float) : vector of observed wavelengths
; CORR (float) : vector of colour corrections at each wavelength
; REF_NAME (str) : vector of labels used in the input file namenames
; PLOT_LAB (str) : vector of labels for the filter on the plots
; MAP_SCALE_ARCSECS (float) : pixel scale in arcsecs of the models
; NLAB (INT) : number of labels to be included on the plot legend

; OUTPUTS:
; PS files of the SB profiles at each wavelength

  ; Add Coyote directory to path
  !PATH = Expand_Path('+/nfs/d58/vlk/sedmodel/MTR/N0628/gn/PROFILE_RV/sca/coyote/') + ':' + !PATH 
  loadct, 13
  
  ; Set common plotting parameters
  !p.thick = 4
  !x.thick = 5
  !y.thick = 5
  !p.charthick = 5
  !p.charsize = 1.5
  
  ; Galaxy parameters
  gal_name = 'NGC628'  ; identifier of the galaxy
  d_gal = 9.5          ; galaxy distance in Mpc
  EBV = 0.254 - 0.192  ; interstellar reddening
  
  ; Initialise residuals for quality of fit
  average_residuals = 0.
  nresiduals = 0L
  
  ; Model parameters
  tau = 12.8
  f = 0.1
  sfr = 2.3
  old = 0.18
  bd = 0.09
  nsersic = 2
  
  ; Filter and wavelength setup
  filter = ['FUV', 'NUV', 'uSDSS', 'g', 'rSDSS', 'iSDSS', 'z', 'j', 'h', 'k', 'ir34', 'ir36', 'ir45', 'ir46', 'ir57']
  lambda = [1516, 2267, 3543, 4770, 6231, 7625, 9134, 12590, 16620, 22000, 33680, 35500, 44930, 46180, 57310] * 1.d ; in AA
  corr = [1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.1888, 1.0129, 1., 1.4049]
  nfilter = n_elements(filter)
  
  refname = ['GALEX_FUV', 'GALEX_NUV', 'SDSS_U', 'SDSS_G', 'SDSS_R', 'SDSS_I', 'SDSS_Z', '2MASS_J', '2MASS_H', '2MASS_K', 'WISE_34', 'IRAC_36', 'IRAC_45', 'WISE_46', 'IRAC_57']
  plot_lab = ['GALEX FUV', 'GALEX NUV', 'SDSS u', 'SDSS g', 'SDSS r', 'SDSS i', 'SDSS z', '2MASS J', '2MASS H', '2MASS K', 'WISE 3.4', 'IRAC 3.6', 'IRAC 4.5', 'WISE 4.6', 'IRAC 5.7']
  
  map_scale_arcsecs = 1.0
  
  ; Model setup
  model = 'wd01'
  qyear = '06'
  scaabs = 'sca'
  stau = strcompress(string(fix(tau*10)), /remove_all)
  snsersic = strcompress(string(fix(nsersic)), /remove_all)
  
  ; Initialise data arrays
  sky = dblarr(nfilter)
  xcen_data = dblarr(nfilter)
  ycen_data = dblarr(nfilter)
  ellipt = dblarr(nfilter)
  pa = dblarr(nfilter)
  pixel_scale = dblarr(nfilter)
  psf_scale_arcsec = dblarr(nfilter)
  
  ; Plot parameters
  wa = 1                    ; width of the annuli in pixels
  xlab = 'arcsecs'
  ylab = 'MJy/sr'
  nlab = 7L ; number of labels added to plot legend
  y_lab = findgen(nlab)
  
  ; Directory paths
  model_label_disk_2d = 'map_'
  dir_mod_2d = '/nfs/d58/vlk/sedmodel/urad_ngc0628_sca/out/'
  dir_plot = 'PLOTS/'       ; directory for output
  
  ; Read data image and galfit parameters
  file_data = 'N0628_uvo_par_v2.dat'
  gal_data_format = '(f6.4,f10.5,f8.2,f8.2,f6.2,f8.3,f12.3,f18.3)'
  
  openr, lun, file_data, /get_lun
  skip_lun, lun, 1, /lines
  readf, lun, ''
  
  ; Read parameters for each filter
  for n = 0, nfilter-1L do begin
    readf, lun, dum1, dum2, dum3, dum4, dum5, dum6, dum7, dum8, format=gal_data_format
    sky[n] = dum2
    xcen_data[n] = dum3
    ycen_data[n] = dum4
    ellipt[n] = dum5
    pa[n] = dum6
    pixel_scale[n] = dum7
    psf_scale_arcsec[n] = dum8
  endfor
  free_lun, lun
  
  ; Scale factors for different components
  scletd = [0.94, 0.76, 0.34, 0.27, 0.14, 0.14, 0.13, 1, 1, 1., 1., 1., 1., 1., 1.] ; young disk
  sclet2d = [3.1, 4.6, 1.8, 2.4, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0] * 1.2 * 0.65 ; second young disk
  scle = [1, 1, 1.70, 1.10, 1.10, 0.9, 0.70, 0.60, 0.60, 0.6, 0.25, 0.30, 0.35, 0.30, 0.3] * 1.1 ; bulge
  scled = [1, 1, 7.7, 5.0, 4.5, 3.8, 2.5, 1.2, 1.60, 1.25, 1.7, 1.9, 2.3, 2.0, 2.3] * 0.7 ; old disk
  
  ; Y-range parameters for log plots
  high_y = [0.1, 0.2, 1.0, 4.0, 10., 10., 10., 10., 10., 10., 10., 10., 10., 10., 10.]
  low_y = [1e-3, 1e-3, 1e-3, 2e-3, 3e-3, 1e-3, 1e-3, 1e-3, 1e-3, 1e-3, 1e-3, 1e-3, 1e-3, 1e-3, 1e-3]
  
  ; Process each filter
  for i = 0, nfilter-1L do begin
    ; Set up directories
    dir_data = '/nfs/d58/vlk/sedmodel/MTR/N0628/N0628_' + refname[i] + '/'
    dir_mask = '/nfs/d58/vlk/sedmodel/MTR/N0628/N0628_' + refname[i] + '/star_remove/' ; mask directory
    
    ; Wavelength and filenames
    slam = string(lambda[i]/1e4, format='(f4.2)')
    fname, model, qyear, tau, sfr, old, bd, scaabs, slam, namep
    
    ; File paths
    data = dir_data + gal_name + '_' + refname[i] + '.fits' ; data file
    fmask = dir_mask + 'mask.fits' ; mask file
    mod_tdisk_2d = dir_mod_2d + model_label_disk_2d + 'mtd_' + filter[i] + '_' + model + '_q' + qyear + '_i0_t' + stau + '_' + scaabs + '.fits'
    mod_t2disk_2d = dir_mod_2d + model_label_disk_2d + 'mtd2_' + filter[i] + '_' + model + '_q' + qyear + '_i0_t' + stau + '_' + scaabs + '.fits'
    
    ; Set PSF file path
    psf = dir_data + 'PSF_' + refname[i] + '.fits' ; default psf file
    if refname[i] eq '2MASS_J' then psf = '/nfs/d58/vlk/sedmodel/MTR/psf/psf_2mass_j.fits'
    if refname[i] eq '2MASS_H' then psf = '/nfs/d58/vlk/sedmodel/MTR/psf/psf_2mass_h.fits'
    if refname[i] eq '2MASS_K' then psf = '/nfs/d58/vlk/sedmodel/MTR/psf/psf_2mass_k.fits'
    
    ; Read data files
    xdata = mrdfits(data, 0, xm) ; read data
    mask = mrdfits(fmask, 0, xm) ; read mask for data
    
    ; Optical wavelength handling
    if lambda[i] gt 2800. then begin
      mod_disk_2d = dir_mod_2d + model_label_disk_2d + 'md_' + filter[i] + '_' + model + '_q' + qyear + '_i0_t' + stau + '_' + scaabs + '.fits'
      mod_bulge_2d = dir_mod_2d + model_label_disk_2d + 'mb_' + filter[i] + '_' + model + '_q' + qyear + '_i0_t' + stau + '_n' + snsersic + '_' + scaabs + '.fits'
      xmod_disk_2d = mrdfits(mod_disk_2d, 0, xm)
      xmod_bulge_2d = mrdfits(mod_bulge_2d, 0, xm)
      
      nx_bulge_2d = n_elements(xmod_bulge_2d[*,0])
      ny_bulge_2d = n_elements(xmod_bulge_2d[0,*])
    endif
    
    ; Read model files
    xmod_tdisk_2d = mrdfits(mod_tdisk_2d, 0, xm)
    xmod_t2disk_2d = mrdfits(mod_t2disk_2d, 0, xm)
    p = mrdfits(psf, 0, xp) ; read psf
    
    ; Model dimensions
    nx_2d = n_elements(xmod_tdisk_2d[*,0])
    ny_2d = n_elements(xmod_tdisk_2d[0,*])
    
    ; Below the models and PSFs are regridded onto the data scale
    ; Regridding factors
    map_regrid_f = map_scale_arcsecs/pixel_scale[i]
    nnx = nx_2d * map_regrid_f
    nny = ny_2d * map_regrid_f
    xcen_mod_2d = 0L
    ycen_mod_2d = nny/2L
    
    if(lambda[i] gt 2800.) then begin
      nnx_bulge = nx_bulge_2d * map_regrid_f
      nny_bulge = ny_bulge_2d * map_regrid_f
      xcen_mod_bulge_2d = 0L
      ycen_mod_bulge_2d = nny_bulge/2L
    endif
    
    ; PSF dimensions
    nx_psf = n_elements(p[*,0]) * psf_scale_arcsec[i]/pixel_scale[i]
    ny_psf = n_elements(p[0,*]) * psf_scale_arcsec[i]/pixel_scale[i]
    
    ; Normalise PSF
    p_regrid = p
    tot_p = total(p_regrid)
    p_regrid = p_regrid/tot_p
    
    if(lambda[i] gt 2800. and lambda[i] lt 10000.) then begin
      p_regrid_bulge = p
      p_regrid_bulge = p_regrid_bulge/total(p_regrid_bulge)
    endif
    
    ; Regrid models
    rmod_tdisk_2d = congrid(xmod_tdisk_2d, nnx, nny)
    rmod_t2disk_2d = congrid(xmod_t2disk_2d, nnx, nny)
    
    if(lambda[i] gt 2800.) then begin
      rmod_disk_2d = congrid(xmod_disk_2d, nnx, nny)
      rmod_bulge_2d = congrid(xmod_bulge_2d, nnx_bulge, nny_bulge)
      mod_conv_bulge = gauss_smooth(rmod_bulge_2d, 1.5, /edge_truncate)
      xcen_mod_bulge_2d = 0L
      ycen_mod_bulge_2d = nny_bulge/2L
    endif
    
    ; Convolve models with PSF
    mod_conv_tdisk = convolve(rmod_tdisk_2d, p_regrid, ft_psf=psf_ft)
    mod_conv_t2disk = convolve(rmod_t2disk_2d, p_regrid, ft_psf=psf_ft)
    
    ; Assuming here gaussian convolution
    if(lambda[i] eq 33680. or lambda[i] eq 46180.) then begin
      mod_conv_tdisk = gauss_smooth(rmod_tdisk_2d, 1.5, /edge_truncate)
      mod_conv_t2disk = gauss_smooth(rmod_t2disk_2d, 1.5, /edge_truncate)
    endif
    
    if(lambda[i] gt 2800.) then begin
      mod_conv_disk = convolve(rmod_disk_2d, p_regrid, ft_psf=psf_ft)
      if(lambda[i] eq 33680. or lambda[i] eq 46180.) then mod_conv_disk = gauss_smooth(rmod_disk_2d, 1.5, /edge_truncate)
    endif
    
    ; Set flux conversion factors based on filter
    CASE refname[i] OF
      'GALEX_FUV': fluxconv = 1.0724d-4*1.0d-6 /(pixel_scale[i]/3600*!pi/180.)^2
      'GALEX_NUV': fluxconv = 3.529d-5*1.0d-6 /(pixel_scale[i]/3600*!pi/180.)^2
      'SDSS_U': fluxconv = 3.631d-6*1.0d-6 /(pixel_scale[i]/3600*!pi/180.)^2
      'SDSS_G': fluxconv = 3.631d-6*1.0d-6 /(pixel_scale[i]/3600*!pi/180.)^2
      'SDSS_R': fluxconv = 3.631d-6*1.0d-6 /(pixel_scale[i]/3600*!pi/180.)^2
      'SDSS_I': fluxconv = 3.631d-6*1.0d-6 /(pixel_scale[i]/3600*!pi/180.)^2
      'SDSS_Z': fluxconv = 3.631d-6*1.0d-6 /(pixel_scale[i]/3600*!pi/180.)^2
      '2MASS_J': fluxconv = 7.8285e-12/(pixel_scale[i]/3600*!pi/180.)^2
      '2MASS_H': fluxconv = 5.7056e-12/(pixel_scale[i]/3600*!pi/180.)^2
      '2MASS_K': fluxconv = 6.5267e-12/(pixel_scale[i]/3600*!pi/180.)^2
      'IRAC_36': fluxconv = 1.0d0
      'IRAC_45': fluxconv = 1.0d0
      'IRAC_57': fluxconv = 1.0d0
      'WISE_34': fluxconv = 1.57e-6*1d-6/(pixel_scale[i]/3600*!pi/180.)^2
      'WISE_46': fluxconv = 2.54e-6*1d-6/(pixel_scale[i]/3600*!pi/180.)^2
      ELSE: fluxconv = 1.0d0
    ENDCASE
    
    ; Correct data for interstellar extinction
    x = 1/(lambda[i]/1d4)
    f = 3.1 - 0.236 + 0.462*x + 0.105*x^2 + 0.454/((x-4.557)^2 + 0.293)
    if(x gt 1.83 and x lt 2.75) then f = 3.1 + 2.56*(x-1.83) - 0.993*(x-1.83)^2
    if(x lt 1.83) then f = ((1.86-0.48*x)*x-0.1)*x
    A = f*EBV
    f0 = 10^(0.4*A)
    flux = (xdata-sky[i])*fluxconv*f0/corr[i]
    
    ; Scale models
    rmod_tdisk_2d = mod_conv_tdisk/(map_scale_arcsecs^2*2.3d-11*1d6)
    rmod_t2disk_2d = mod_conv_t2disk*0.008/(map_scale_arcsecs^2*2.3d-11*1d6)
    
    ; Calculate sb profiles
    cog, flux, mask, pa[i], ellipt[i], xcen_data[i], ycen_data[i], wa, rcog_data, bcog_data, ncog, scog_data, r2pi
    rarcsdat = rcog_data * pixel_scale[i]
    plot,rcog_data, bcog_data
    
    ; Create mask for models
    mask = fltarr(nnx, nny) + 1.

    ; Calculate circular profiles for models
    cog, rmod_tdisk_2d, mask, pa[i], ellipt[i], xcen_mod_2d, ycen_mod_2d, wa, rcog_mod, bcog_mod_tdisk_2d, ncog, scog_mod_tdisk, r2pi
    cog, rmod_t2disk_2d, mask, pa[i], ellipt[i], xcen_mod_2d, ycen_mod_2d, wa, rcog_mod, bcog_mod_t2disk_2d, ncog, scog_mod_t2disk, r2pi
    
    ; Apply scale factors
    bcog_mod_tdisk_2d = scletd[i] * bcog_mod_tdisk_2d
    bcog_mod_t2disk_2d = sclet2d[i] * bcog_mod_t2disk_2d
    scog_mod_tdisk = scog_mod_tdisk * scletd[i] * 1e6 * (pixel_scale[i]/3600*!pi/180.)^2
    scog_mod_t2disk = scog_mod_t2disk * sclet2d[i] * 1e6 * (pixel_scale[i]/3600*!pi/180.)^2
    scog_mod_tot = scog_mod_t2disk + scog_mod_tdisk
    
    ; Handle optical wavelengths
    if(lambda[i] gt 2800.) then begin
      mask = fltarr(nnx_bulge, nny_bulge) + 1.
      
      rmod_bulge_2d = mod_conv_bulge/(map_scale_arcsecs^2*2.3d-11*1d6)
      cog, rmod_bulge_2d, mask, pa[i], ellipt[i], xcen_mod_bulge_2d, ycen_mod_bulge_2d, wa, rcog_mod, bcog_mod_bulge_2d, ncog, scog_mod_bulge, r2pi
      
      ; Extend bulge profile if needed
      if n_elements(bcog_mod_bulge_2d) lt 300 then begin
        bcog_mod_bulge_2d = [bcog_mod_bulge_2d, replicate(0., 300-n_elements(bcog_mod_bulge_2d))]
      endif
      
      rmod_disk_2d = mod_conv_disk/(map_scale_arcsecs^2*2.3d-11*1d6)
      cog, rmod_disk_2d, mask, pa[i], ellipt[i], xcen_mod_2d, ycen_mod_2d, wa, rcog_mod, bcog_mod_disk_2d, ncog, scog_mod_disk, r2pi
      
      ; Apply scale factors
      bcog_mod_disk_2d = bcog_mod_disk_2d * scled[i]
      bcog_mod_bulge_2d = bcog_mod_bulge_2d * scle[i]
      scog_mod_disk = scog_mod_tdisk * scled[i] * 1e6 * (pixel_scale[i]/3600*!pi/180.)^2
      scog_mod_bulge = scog_mod_bulge * scle[i] * 1e6 * (pixel_scale[i]/3600*!pi/180.)^2
      
      ; Combine components
      bcog_mod_tot = bcog_mod_tdisk_2d + bcog_mod_disk_2d + bcog_mod_t2disk_2d + bcog_mod_bulge_2d * 1.0
      scog_mod_tot = scog_mod_t2disk + scog_mod_tdisk + scog_mod_disk + scog_mod_bulge
    endif else begin
      ; For UV wavelengths
      bcog_mod_tot = bcog_mod_tdisk_2d + bcog_mod_t2disk_2d
    endelse
    
    ; Dust emission handling for IR bands
    if(lambda[i] gt 30000.) then begin
      dir_ir = '/nfs/d58/vlk/sedmodel/MTR/N0628/gn/PROFILE_RV/sca/MAPS/old/2018_12_14/'
      
      ; Handle different IR bands
      case filter[i] of
        'ir34': begin
          dir_psf = 'PSF/WISE/'
          psf = dir_psf + 'PSF_WISE_34.fits'
          wav_ir = 3.55
          fname, model, qyear, tau, sfr, old, bd, scaabs, wav_ir, namet
          model_dust_tot = dir_ir + 'map_' + namet + '_tot.fits'
          xmod_dust_tot = mrdfits(model_dust_tot, 0, xm)
        end
        'ir46': begin
          dir_psf = 'PSF/WISE/'
          psf = dir_psf + 'PSF_WISE_45.fits'
          wav_ir = 4.60
          fname, model, qyear, tau, sfr, old, bd, scaabs, wav_ir, namet
          model_tot = dir_ir + 'map_' + namet + '_tot.fits'
          xmod_dust_tot = mrdfits(model_dust_tot, 0, xm)
        end
        'ir36': begin
          dir_psf = 'PSF/IRAC/'
          psf = dir_psf + 'PSF_IRAC_36.fits'
          wav_ir = 3.55
          fname, model, qyear, tau, sfr, old, bd, scaabs, wav_ir, namet
          model_tot = dir_ir + 'map_' + namet + '_tot.fits'
          xmod_dust_tot = mrdfits(model_dust_tot, 0, xm)
        end
        'ir45': begin
          dir_psf = 'PSF/IRAC/'
          psf = dir_psf + 'PSF_IRAC_45.fits'
          wav_ir = 4.49
          fname, model, qyear, tau, sfr, old, bd, scaabs, wav_ir, namet
          model_tot = dir_ir + 'map_' + namet + '_tot.fits'
          xmod_dust_tot = mrdfits(model_dust_tot, 0, xm)
        end
        'ir57': begin
          dir_psf = 'PSF/IRAC/'
          psf = dir_psf + 'PSF_IRAC_57.fits'
          wav_ir = 5.73
          fname, model, qyear, tau, sfr, old, bd, scaabs, wav_ir, namet
          model_tot = dir_ir + 'map_' + namet + '_tot.fits'
          xmod_dust_tot = mrdfits(model_dust_tot, 0, xm)
        end
      endcase
      
      ; Process dust model
      map_size_pc = 68544.0        ; printed model output
      map_size_pix = 2500L         ; n
      map_dust_scale_pc = map_size_pc/map_size_pix         ; map scale in parsecs
      map_dust_scale_arcsecs = map_dust_scale_pc/(d_gal*1d6)*180/!pi*3600  ; map scale in arcsecs
      map_regrid_f = map_dust_scale_arcsecs/pixel_scale[i]
      
      nx = n_elements(xmod_dust_tot[*,0])
      ny = n_elements(xmod_dust_tot[0,*])
      nnx = nx * map_regrid_f
      nny = ny * map_regrid_f
      xcen_mod = nnx/2L
      ycen_mod = nny/2L
      
      ; Read PSF
      p = mrdfits(psf, 0, xp)
      nx_psf = n_elements(p[*,0]) * psf_scale_arcsec[i]/pixel_scale[i]
      ny_psf = n_elements(p[0,*]) * psf_scale_arcsec[i]/pixel_scale[i]
      
      p_regrid = congrid(p, nx_psf, ny_psf)
      p_regrid = p_regrid/total(p_regrid)
      
      ; Process dust model
      rmod_dust_tot = congrid(xmod_dust_tot, nnx, nny)
      mod_conv_dust_tot = convolve(rmod_dust_tot, p_regrid, ft_psf=psf_ft)
      f = 1.d20/3.09d16^2     ; model scale factor
      mod_conv_dust_tot = mod_conv_dust_tot * f
      
      ; Calculate profile
      cog, mod_conv_dust_tot, mask, pa[i], ellipt[i], xcen_mod, ycen_mod, wa, rcog_dust_mod, bcog_mod_dust_tot, scog_mod_dust_tot, ncog
    endif
    
    ; Set plot parameters
    rarcs = map_scale_arcsecs * rcog_mod
    plot_xrange = [0, 400]
    plot_yrange = [low_y[i], high_y[i]]
    
    ; Calculate error bands (20% either side)
    high_error = bcog_data + (bcog_data / 100*20)
    low_error = bcog_data - (bcog_data / 100*20)
    
    ; Set label positions for log plot
    oft_x = 0.6d0 ; offset in x from right-side of the plot
    x_lab = oft_x * plot_xrange[1]
    oft_y = (alog10(plot_yrange(1)) - alog10(plot_yrange(0)))*0.93 ; vertical offset between labels
    y_lab_max_log = alog10(plot_yrange(0)) + oft_y ; position of the first ('highest') label
    y_lab_max_pos = 10^y_lab_max_log
    y_lab_pos_log = (alog10(plot_yrange(1))-alog10(y_lab_max_pos)) * findgen(nlab) 
    y_lab_pos = 10^y_lab_pos_log
    y_lab_pos(0) = 1 ; ensure first label is at y_lab_max_pos
    y_lab = y_lab_max_pos / y_lab_pos ; y label positions in linear space
    
    ; Create standard plot
    set_plot, 'PS'
    device, file='lum_'+namep+'_log.ps', /color
    
    ; Plot data and models
    plot, rarcsdat, bcog_data, ystyle=1, xstyle=1, xrange=plot_xrange, yrange=[low_y[i], high_y[i]], $
          xtickformat="(A1)", ytit=ylab, position=[.16, 0.37, .97, .97], /ylog
    
    ; Add shaded 20 percent error margin 
    n_points = n_elements(rarcsdat)
    x_poly = [rarcsdat, reverse(rarcsdat)]
    y_poly = [high_error, reverse(low_error)]
    polyfill, x_poly, y_poly, color=210, /line_fill, orientation=45, spacing=0.5, thick=1
    
    ; Plot data and model components
    oplot, rarcsdat, bcog_data, thick=3
    oplot, rarcsdat, bcog_mod_tdisk_2d, color=100
    oplot, rarcsdat, bcog_mod_t2disk_2d, color=80
    
    ; Add bulge for non-UV wavelengths
    if(lambda[i] gt 2800.) then begin
      oplot, rarcsdat, bcog_mod_disk_2d, color=220
      oplot, rarcsdat, bcog_mod_bulge_2d*1.0, color=250
    endif
    
    ; Plot total model
    if(lambda[i] lt 30000.) then oplot, rarcsdat, bcog_mod_tot, color=150
    if(lambda[i] gt 30000.) then begin
      oplot, rarcsdat, bcog_mod_tot+bcog_mod_dust_tot, color=150
      oplot, rarcsdat, bcog_mod_dust_tot, color=50
    endif
    
    ; Labels
    xyouts, x_lab, y_lab[0], plot_lab[i]+' band'
    xyouts, x_lab, y_lab[1], 'Total Emission', color=150
    xyouts,x_lab,y_lab(2),'Young Stellar Disk',color=100,charsize=1.4
    xyouts,x_lab,y_lab(3),'Young Stellar Ring',color=80,charsize=1.4
    if(lambda(i) gt 2800.) then begin
     xyouts,x_lab,y_lab(4),'Old Stellar Disk',color=220,charsize=1.4
     xyouts,x_lab,y_lab(5),'Bulge',color=250,charsize=1.4
    endif ; lambda(i) gt 2800.

   ; For some bands, it is necessary to restrict the residuals to radii below
   ; the noise level (there is no point showing residuals in the noise)
   residual_restrict = where(rarcsdat lt 290)
   if(lambda[i] eq 12590) then residual_restrict = where(rarcsdat lt 200)
   if(lambda[i] eq 16620) then residual_restrict = where(rarcsdat lt 180)
   if(lambda[i] eq 22000) then residual_restrict = where(rarcsdat lt 150)
   
   residuals=(bcog_data[residual_restrict]-bcog_mod_tot[residual_restrict]) / $ 
	   bcog_data[residual_restrict]*100 

   for n=0,n_elements(residuals)-1L do begin
       if(finite(residuals[n], /INFINITY) eq 0L) then begin
           average_residuals = average_residuals + abs(residuals[n])
           nresiduals = nresiduals + 1L
       endif ; (finite(residuals[n], /INFINITY) eq 0L)
   endfor ; n=0,n_elements(residuals)-1L 

  if(lambda[i] gt 30000.) then xyouts, x_lab, y_lab[6], 'Dust emission', color=50

  Plot, rarcsdat, residuals, yrange=[-100,100], position=[.16,0.12,.97,.36], xrange=plot_xrange, xtit='arcsecs', ytit='R (%)', ytickinterval=100, $
      /NOERASE, thick=5

  Plots, [0,400], [0,0], linestyle=1, thick=5
  Plots, [0,400], [20,20], linestyle=2, color='sky blue', thick=5
  Plots, [0,400], [-20,-20], linestyle=2, color='sky blue', thick=5

  endfor ;i=0,nfilter-1L 
  
  average_residuals = average_residuals/nresiduals
  print, 'Average residual is = ', average_residuals
end
