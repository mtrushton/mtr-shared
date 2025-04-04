pro int_z, nposr, nposz, file, r, pass_energy_z
  ; Get dimensions and initialise arrays
  npos = nposr * nposz
  
  ; Read data from file
  openr, unit, file, /get_lun
  rr = dblarr(npos)
  zz = dblarr(npos)
  lum = dblarr(npos)
  readf, unit, rr, zz, lum, format='(3F)'
  free_lun, unit
  
  ; Reshape 1D arrays to 2D grid
  rrr = reform(rr, nposr, nposz)
  zzz = reform(zz, nposr, nposz)
  lumdouble = reform(lum, nposr, nposz)
  
  ; Integrate over z
  pass_energy_z = dblarr(nposr)
  
  for kr = 0L, nposr - 2L do begin
    ; Calculate z-differences and average luminosities
    z_diffs = zzz[kr, 1:nposz-1] - zzz[kr, 0:nposz-2]
    avg_lum = (lumdouble[kr, 0:nposz-2] + lumdouble[kr+1, 1:nposz-1]) / 2.0
    
    ; Sum over z dimension
    pass_energy_z[kr] = total(z_diffs * avg_lum)
  endfor ; kr = 0L, nposr - 2L 
  
  ; Apply multiplication factor
  pass_energy_z = pass_energy_z * 2.0
  
  ; Extract unique r values
  r = rr[uniq(rr, sort(rr))]
end
