pro int_z, npos, nposr, nposz, ncomp, file, r, dust_dens_z, scale
; Integrate dust density over z and calculate scaling factors
  
; Files containing the opacities of each dust component (unscaled)
 file = ['opacity0.dat', 'opacity1.dat', 'opacity2.dat', 'opacity3.dat', 'opacity4.dat', 'opacity5.dat']
  
; Corresponding dust mass for each component (M_sol) used for scaling
 dust_comp_mass = [97950360., 3107165.5, 1567412.4, 26159.156, 27386412., 9599040.0]
  
; Get number of components
 ncomp = n_elements(file)
 scale = dblarr(ncomp)
  
; Calculate total number of grid points
 npos = nposr * nposz
  
; Initialise arrays
 dust_dens_z = dblarr(nposr, ncomp)
  
; Process each component
 for i = 0, ncomp - 1 do begin
   ; Read opacity data
   openr, unit, file[i], /get_lun
    
   ; Arrays to store grid data
   rr = dblarr(npos)
   zz = dblarr(npos)
   tau = dblarr(npos)
    
   ; Read data
   for n = 0L, npos - 1L do begin
     readf, unit, dum1, dum2, dum3
     rr[n] = dum1
     zz[n] = dum2
     tau[n] = dum3
   endfor ; n = 0L, npos - 1L 
    
    free_lun, unit
    
    ; Reshape data to 2D grid
    taudouble = dblarr(nposr, nposz)
    rrr = dblarr(nposr, nposz)
    zzz = dblarr(nposr, nposz)
    
    k = 0L
    for kr = 0L, nposr - 1L do begin
      for kz = 0L, nposz - 1L do begin
        taudouble[kr, kz] = tau[k]
        rrr[kr, kz] = rr[k]
        zzz[kr, kz] = zz[k]
        k = k + 1
      endfor ; kz = 0L, nposz - 1L 
    endfor ; kr = 0L, nposr - 1L 
    
    ; Integrate over z for each radial position
    for kr = 0L, nposr - 2L do begin
      pass_z = 0.0
      for kz = 0L, nposz - 2L do begin
        pass_z = pass_z + (zzz[kr, kz + 1] - zzz[kr, kz]) * $
                 (taudouble[kr, kz] + taudouble[kr + 1, kz + 1]) / 2.
      endfor ; kz = 0L, nposz - 2L 
      dust_dens_z[kr, i] = pass_z
    endfor ; kr = 0L, nposr - 2L 
    
    ; Account for symmetry about z=0
    dust_dens_z[*, i] = dust_dens_z[*, i] * 2.
    
    ; Integrate over positions to find scaling factor
    pass_sum = 0.0
    for kz = 0L, nposz - 2 do begin
      for kr = 0L, nposr - 2 do begin
        pass = !pi * (rrr[kr+1, kz]^2 - rrr[kr, kz]^2) * (zzz[kr, kz+1] - zzz[kr, kz]) * $
               (taudouble[kr, kz] + taudouble[kr+1, kz+1]) / 2.
        pass_sum = pass_sum + pass
      endfor ; kr = 0L, nposr - 2
    endfor ; kz = 0L, nposz - 2
    
    ; Account for symmetry and calculate scaling factor
    pass_sum = pass_sum * 2.
    scale[i] = dust_comp_mass[i] / pass_sum
  endfor ; n = 0L, npos - 1L 
  
  ; Extract unique radial positions
  ind = rem_dup(rr)
  r = rr[ind]
end
