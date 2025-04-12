function integrate_model, nposr, nposz, rr, zz, lum
; Reshape 1D arrays to 2D grid
 rrr = reform(rr, nposr, nposz)
 zzz = reform(zz, nposr, nposz)
 lumdouble = reform(lum, nposr, nposz)
  
; Calculate cell volumes and luminosities
 r_diffs = rrr[1:*,0:nposz-2]^2 - rrr[0:nposr-2,0:nposz-2]^2
 z_diffs = zzz[0:nposr-2,1:*] - zzz[0:nposr-2,0:nposz-2]
 avg_lum = (lumdouble[0:nposr-2,0:nposz-2] + lumdouble[1:*,1:*]) / 2.0
  
 total_energy = 2.0 * !pi * total(r_diffs * z_diffs * avg_lum)
  
 return, total_energy
end ; integrate_model

pro int_model, nposr, nposz, total_model_tdisk, total_model_tdisk2, total_model_disk, total_model_bulge
; Define arrays for model components
 component_files = ['emissivity_tdisk.dat', 'emissivity_disk.dat', 'emissivity_bulge.dat', 'emissivity_tdisk2.dat']
 model_totals = dblarr(n_elements(component_files)
 npos = nposr * nposz
  
  ; Process each component file
for i = 0, n_elements(component_files) do begin
; Read data from file
    rr = dblarr(npos)
    zz = dblarr(npos)
    lum = dblarr(npos)
    
    openr, unit, component_files[i], /get_lun
    readf, unit, rr, zz, lum, format='(3F)'
    free_lun, unit
    
    ; Model integration
    model_totals[i] = integrate_model(nposr, nposz, rr, zz, lum)
endfor ; i = 0, n_elements(component_files)
  
; Assign results to output variables
 total_model_tdisk = model_totals[0]
 total_model_disk = model_totals[1]
 total_model_bulge = model_totals[2]
 total_model_tdisk2 = model_totals[3]
end
