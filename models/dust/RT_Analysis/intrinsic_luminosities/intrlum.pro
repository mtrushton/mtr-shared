forward_function Calculate_Integrated_Luminosity
pro intrlum
  ; Calculates and plots the intrinsic luminosity density at each sampled wavelengths
  ; for a galaxy modelled with the PT11 RT code. 
  ;
  ; Note that the young stellar component should also include the clumpy component. 
  ; For this reason, CLUMPY.PRO should be run to get the intrinsic luminosity of the
  ; local component at each sampled wavelength.
  ;
  ; INPUTS:
  ; SFR (float) : best fitting sfr
  ; OLD (float) : best fitting value of old
  ; BD (float) : bulge-to-disk ratio
  ; LAMBA (float) : vector of sampled wavelengths (Angstroms)
  ; 
  ; OUTPUTS:
  ; Plot of the SED of intrinsic luminosities
  ; IDL save file containing the intrinsic luminosity density at
  ;    sampled wavelength for each component
  ; Ascii files containing the intrinsic luminosity density at
  ;    sampled wavelength for each component
  ; Integrated intrinsic luminosites printed to the screen
  ;
  ; Best fitting parameters, needed for scaling the template
  
  sfr = 4.2
  old = 0.3
  bd = 0.1

  ; Input files containing the scaling factors of the template at each sampled
  ;     wavelength
  file_young = '../scale_factors_young.dat'
  file_old = '../scale_factors_old.dat'
  file_clumpy = 'clumpy.dat'
  
  ; Sampled wavelengths
  lambda = [912.,1350.,1500.,1650.,2000.,2200.,2500.,2800.,3650.,4430.,5640., $
           8090.,12590.,22000.,36000.,45000.,58000.] * 1.d  ; in AA
  ; young template (W/m^2; see PT11)
  ll1 = [0.344,0.905,0.844,0.863,0.908,0.926,0.843,0.910,1.842,2.271,3.837, $
         5.734,0.931,0.728,0.141,0.141,0.141] * 1.d+21
  ; old template
  ll = [0.,0.,0.,0.,0.,0.,0.,0.,4.771,4.771,9.382,19.54,72.20,64.97,12.58, $
        12.58,12.58] * 1.d+21
  
  np = n_elements(lambda)
  freq_uvo = 2.998d8 / (lambda * 1d-10)

  ; Initialise arrays
  lum_young_disk = dblarr(np)
  lum_old_disk = dblarr(np)
  lum_bulge = dblarr(np)
  lum_clumpy = dblarr(np)

  ; Get scaling factors
  openr, unit, file_young, /get_lun
  openr, unit1, file_old, /get_lun
  openr, unit5, file_clumpy, /get_lun

  ; Header
  readf, unit, ss
  readf, unit1, ss

  for n = 0, np - 1 do begin
    readf, unit, dum1, dum2
    readf, unit1, dum3, dum4, dum5
    readf, unit5, dum6, dum7
    ; Scale the disks to get the intrinsic luminosities
    lum_young_disk[n] = ll1[n] * sfr * dum2
    lum_bulge[n] = ll[n] * old * bd * dum4
    lum_old_disk[n] = ll[n] * old * dum5
  endfor ; n = 0, np - 1 
  
  ; Add clumpy component to young disk to get the total young component
  lum_young_disk += lum_clumpy
  
  save, lambda, lum_young_disk, lum_old_disk, lum_bulge, file = 'NGC3938_intrinsic_SED.sav'
  
  ; Write intrinsic luminosities to files
  c = 2.998d8
  openw, unit2, 'young_stellar_disk.dat', /get_lun
  openw, unit3, 'old_stellar_disk.dat', /get_lun
  openw, unit4, 'old_stellar_bulge.dat', /get_lun
  
  for n = 0, np - 1 do begin
    freq = c / (lambda[n] * 1d-6)
    printf, unit2, freq, lum_young_disk[n]
    printf, unit3, freq, lum_old_disk[n]
    printf, unit4, freq, lum_bulge[n]
  endfor ; n = 0, np - 1
  
  ; Calculate plot settings
  lambda_micron = lambda / 1e4
  
  ymin = 1e19
  ymax = 2e24
  xmin = 0.15
  xmax = 5.73
  
  ; Plot parameters
  loadct, 13
  !p.thick = 4
  !x.thick = 5
  !y.thick = 5
  !p.charthick = 5
  !p.charsize = 1.5
  
  ; Create plot
  set_plot, 'PS'
  device, file='intrinsic_SED.ps', /color
  
  plot, lambda_micron, lum_young_disk, /xlog, /ylog, yrange=[ymin,ymax], xstyle=1, $
        ytit='F!D!7t!3!N (W/Hz)', xtit='!7k!3!N (!7l!3m)', ystyle=1, /nodata
  
  oplot, lambda_micron, lum_young_disk, color=100, psym=4
  oplot, lambda_micron, lum_young_disk, color=100
  oplot, lambda_micron, lum_old_disk, color=250, psym=4
  oplot, lambda_micron, lum_old_disk, color=250
  oplot, lambda_micron, lum_bulge, color=50, psym=4
  oplot, lambda_micron, lum_bulge, color=50
  
  ; Legend
  plot_label = ['stellar disk', 'bulge', 'thin stellar disk']
  plot_label_colour = [250, 50, 100]
  
  nlab = n_elements(plot_label)
  plot_label_xfrac = 0.02
  plot_label_xpos = plot_label_xfrac * xmax
  plot_xrange = xmax - xmin
  bar_length = 0.004 * plot_xrange 
  lab_oft_x = 0.003 * plot_xrange
  
  oft_ymax = (alog10(ymax) - alog10(ymin)) * 0.93
  delta_y = alog10(ymin) + oft_ymax
  oft_y = 10^delta_y
  bs = (alog10(ymax) - alog10(oft_y)) * findgen(nlab)
  y_lab = oft_y / 10^bs
  
  for n = 0, nlab - 1 do begin
    plots, [plot_label_xpos, plot_label_xpos + bar_length], $
           y_lab[n], color = plot_label_colour[n], thick = 6 
    xyouts, plot_label_xpos + bar_length + lab_oft_x, y_lab[n], plot_label[n]
  endfor ; n = 0, nlab - 1 
  
  device, /close
  
  ; Calculate integrated luminosities
  Calculate_Integrated_Luminosity, lum_young_disk, freq_uvo, sum1
  Calculate_Integrated_Luminosity, lum_old_disk, freq_uvo, sum2
  Calculate_Integrated_Luminosity, lum_bulge, freq_uvo, sum3
  
  print, 'Young disk luminosity: ', sum1, ' W'
  print, 'Old disk luminosity: ', sum2, ' W'
  print, 'Bulge luminosity: ', sum3, ' W'
  print, 'Total luminosity: ', sum1 + sum2 + sum3, ' W'
  
  ; Close open files
  free_lun, unit, unit1, unit2, unit3, unit4, unit5
end

; Calculate integrated luminosity
pro Calculate_Integrated_Luminosity, luminosity, frequency, integrated_lum
  integrated_lum = 0.0d
  np = n_elements(luminosity)
  
  for i = np-1, 1, -1 do begin
    s1 = luminosity[i]
    s2 = luminosity[i-1]
    nu1 = frequency[i]
    nu2 = frequency[i-1]
    integrated_lum += 0.5d0 * abs(nu2-nu1) * (s2+s1)
  endfor ; i = np-1, 1, -1 
end
