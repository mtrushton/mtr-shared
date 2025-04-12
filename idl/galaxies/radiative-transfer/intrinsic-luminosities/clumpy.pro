pro clumpy
  ; Determines the luminosity density of the local 'clumpy' dust component in 
  ; a PT11 RT fitted galaxy. 
  ;
  ; INPUTS:
  ; L_local (float) : bolometric luminosity of the clumpy component
  ; Lambda (float): vector of sampled wavelengths
  ; F (float) : values of the 'F factor' at F=0.35
  ; ll1 (float) : vector of Lunit for the young disks. This is also known as
  ;     the 'template'
  
  ; OUTPUTS:
  ; CLUMPY.DAT (file): two column file containing the luminosity density of the
  ;     clumpy component at each frequency  

  ; Luminosity of the local 'clumpy' component
  ; found in the PT11 RT fitting
  L_local = 1.7327888d35

  ; Sampled wavelengths and scaling factors (the F factor and Lunit; see PT11)
  lambda = [912.,1350.,1500.,1650.,2000.,2200.,2500.,2800.,3650.,4430.,5640.,8090.,12590.,22000.,36000.,45000.,58000.] * 1.d ;in AA
  f = [0.573,0.516,0.473,0.455,0.372,0.324,0.261,0.206,0.108,0.068,0.038,0.016,0.009,0.006,0.001,0.001,0.001] ;Fcal=0.35
  ll1 = [0.344,0.905,0.844,0.863,0.908,0.926,0.843,0.910,1.842,2.271,3.837,5.734,0.931,0.728,0.141,0.141,0.141] * 1.d+21

  freq = 2.998d8 / (lambda * 1e-10)
  nfreq = n_elements(freq)
  
  df = dblarr(nfreq)
  df[0:nfreq-2] = freq[0:nfreq-2] - freq[1:nfreq-1]
  
  ; Normalise fractions
  f_norm = f / total(f)
  f_freq = f_norm / df

  ; Sort frequencies and scale appropriately
  ii = sort(freq)
  sorted_freq = freq[ii]
  fdiff = f
  cnst = L_local / int_tabulated(sorted_freq, fdiff[ii])
  scale = fdiff[ii] * cnst
  
  ; integrate 
  print, int_tabulated(sorted_freq, scale)
  print, sorted_freq
  
  openw, unit, 'clumpy.dat', /get_lun
  reverse_freq = reverse(sorted_freq)
  reverse_scale = reverse(scale)
  
  for i = 0, nfreq - 1 do $
    printf, unit, reverse_freq[i], reverse_scale[i]
    
  free_lun, unit

end
