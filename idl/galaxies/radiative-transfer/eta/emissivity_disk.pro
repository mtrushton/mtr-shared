pro emissivity_disk, nposr, nposz, dr, dz, hstin, hs, hsin, hssolar, zs, zsin, zssolar, xis0, zs01, sha, rtrun
;
; Main old stellar disk emissivity
; 
; Constants
 bfactor = 1.
 bin = 1
 bin2 = bin * bin
 pixsize = 1.
 pixsize2 = pixsize * pixsize
 refeta = 20.
 eta0 = (10.^(-refeta / 2.5)) * pixsize2 * bin2 / (2.d0 * zs * bfactor)
  
 flares02 = alog10(hssolar / hsin)
 flares01 = alog10((zssolar - zs)/(zsin - zs))
 flares = flares01 / flares02
  
 sigtruncate = sha * rtrun
  
; Open output file
 openw, unit, 'emissivity_disk.dat', /get_lun
  
 for i = 0L, nposr - 1L do begin
    rho = i * dr  ; Calculate rho directly from i
    
    for n = 0L, nposz - 1L do begin
       z = n * dz  ; Calculate z directly from n
      
       truncate = 0.5 * (1.0 - erf((rho - rtrun) / sigtruncate))
      
      ; Calculate disk emissivity based on radial position
       if (rho lt hstin) then begin
         etadisk0 = 0
       endif else if (rho lt hsin) then begin
         etadisk0 = eta0 * (zs/zs01) * ((rho/hsin) * (1. - xis0) + xis0) * $
                  exp(-hsin/hs - abs(z)/zs01)
       endif else begin
         etadisk0 = eta0 * (zs/zs01) * exp(-rho/hs - abs(z)/zs01)
       endelse
      
       etadisk0 *= truncate
    
       printf, unit, rho, z, etadisk0
    endfor ; n = 0L, nposz - 1L 
 endfor ; i = 0L, nposr - 1L
  
 free_lun, unit
end
