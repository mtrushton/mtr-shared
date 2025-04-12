pro lum_tex
  ; PRINTS intrinsic luminosities to aid creation of latex table
  ; Fitted parameters (used for scaling template)
  sfr = 2.0
  old = 3.0
  bd = 0.1

 ; The following files should contain the factors for scaling the template at 
 ; each sampled wavelength
  file_young = '../scale_factors_young.dat'
  file_old = '../scale_factors_old.dat'
  file_clumpy = 'clumpy.dat'

 ; List of sampled wavelengths and template luminosities in W/m^2 (see PT11)
  lambda = [912.,1350.,1500.,1650.,2000.,2200.,2500.,2800.,3650.,4430.,5640.,8090.,12590.,22000.,36000.,45000.,58000.] * 1.d ;in AA
  ll1 = [0.344,0.905,0.844,0.863,0.908,0.926,0.843,0.910,1.842,2.271,3.837,5.734,0.931,0.728,0.141,0.141,0.141] * 1.d+21
  ll = [0.,0.,0.,0.,0.,0.,0.,0.,4.771,4.771,9.382,19.54,72.20,64.97,12.58,12.58,12.58] * 1.d+21

  np = n_elements(lambda)
  
  ; Initialise arrays
  lum_young_disk = dblarr(np)
  lum_old_disk = dblarr(np)
  lum_bulge = dblarr(np)
  lum_clumpy = dblarr(np)

  ; Open files 
  openr, unit_young, file_young, /get_lun
  openr, unit_old, file_old, /get_lun
  openr, unit_clumpy, file_clumpy, /get_lun

  ; Headers
  readf, unit_young, ss
  readf, unit_old, ss

  for n = 0, np - 1 do begin
    readf, unit_young, dum1, dum2
    readf, unit_old, dum3, dum4, dum5
    readf, unit_clumpy, dum6, dum7
    
    lum_young_disk[n] = ll1[n] * sfr * dum2
    lum_bulge[n] = ll[n] * old * bd * dum4
    lum_old_disk[n] = ll[n] * old * dum5
    lum_clumpy[n] = dum7
    
    print, lambda[n] / 1e4, ' & ', lum_bulge[n], ' & ', lum_old_disk[n], ' & ', $
           lum_young_disk[n] + lum_clumpy[n], ' //'
  endfor ; n = 0, np - 1

  free_lun, unit_young
  free_lun, unit_old
  free_lun, unit_clumpy

end
