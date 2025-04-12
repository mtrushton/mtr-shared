pro emissivity_tdisk, nposr, nposz, dr, dz, hs1tin, hs1, hs1in, hs1solar, zs1, zs1in, zs1solar, xis1, zs11, sha1, rtrun1
;
; Thin (young stellar) disk emissivity
;
; Constants
 bfactor = 1.
 bin = 1
 bin2 = bin * bin
 pixsize = 1.
 pixsize2 = pixsize * pixsize
 refeta = 20.
 eta1 = (10.^(-refeta / 2.5)) * pixsize2 * bin2 / (2.d0 * zs1 * bfactor)
  
 flares02 = alog10(hs1solar / hs1in)
 flares01 = alog10((zs1solar - zs1)/(zs1in - zs1))
 flares = flares01 / flares02
  
 sigtruncate = sha1 * rtrun1
  
; Open output file
 openw, unit, 'emissivity_tdisk.dat', /get_lun
  
 for i = 0L, nposr - 1L do begin
    rho = i * dr  ; Calculate rho directly from i
    
    for n = 0L, nposz - 1L do begin
       z = n * dz  ; Calculate z directly from n
      
       truncate = 0.5 * (1.0 - erf((rho - rtrun1) / sigtruncate))
      
     ; Calculate thick disk emissivity based on radial position
       if (rho lt hs1tin) then begin
         etadisk1 = 0
       endif else if (rho lt hs1in) then begin
        etadisk1 = eta1 * (zs1/zs11) * ((rho/hs1in) * (1. - xis1) + xis1) * $
                   exp(-hs1in/hs1 - abs(z)/zs11)
       endif else begin
        etadisk1 = eta1 * (zs1/zs11) * exp(-rho/hs1 - abs(z)/zs11)
       endelse
      
       etadisk1 *= truncate
      
       printf, unit, rho, z, etadisk1
    endfor ; n = 0L, nposz - 1L
 endfor ; i = 0L, nposr - 1L 
  
 free_lun, unit
end
