pro int, nposr, nposz, npos, scale, ncomp, file
  ; Calculate scaling factors for dust density components
  
  ; Files containing the opacities of each dust component
  file = ['opacity0.dat', 'opacity1.dat', 'opacity2.dat']
  dust_comp_mass = [65616992., 3999528.8, 0.0000] ; masses of each dust component (GLOB_PARAM.DAT)
  
  ncomp = n_elements(file)
  scale = dblarr(ncomp)
  npos = nposr * nposz
  
  for i = 0, ncomp - 1 do begin
    rr = dblarr(npos)
    zz = dblarr(npos)
    tau = dblarr(npos)
    
    openr, unit, file[i], /get_lun
    readf, unit, rr, zz, tau, format='(3F)'
    free_lun, unit
    
    ; Reshape arrays to 2D
    taudouble = reform(tau, nposr, nposz)
    rrr = reform(rr, nposr, nposz)
    zzz = reform(zz, nposr, nposz)
    
    ; Calculate integrated dust mass
    pass_sum = 0.0
    for Cdo begin
      for kr = 0, nposr-2 do begin
        pass = !pi * (rrr[kr+1,kz]^2 - rrr[kr,kz]^2) * (zzz[kr,kz+1] - zzz[kr,kz]) * $
               (taudouble[kr,kz] + taudouble[kr+1,kz+1])/2.
        pass_sum += pass
      endfor ; kr = 0, nposr-2 
    endfor ;  Cdo 
    
    pass_sum *= 2.  ; Account for symmetry
    scale[i] = (dust_comp_mass[i] > 0) / (pass_sum > 1e-10)  ; Avoid division by zero
  endfor ;  i = 0, ncomp - 1
end
