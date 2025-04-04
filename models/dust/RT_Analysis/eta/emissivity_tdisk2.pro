pro emissivity_tdisk2, nposr, nposz, dr, dz, hs2tin, hs2, hs2in, hs2solar, zs2, zs2in, zs2solar, xis2, zs21, sha2, rtrun2
;
; Second thin (young stellar) disk emissivity
;
; Constants
 bfactor = 1.
 bin = 1
 bin2 = bin * bin
 pixsize = 1.
 pixsize2 = pixsize * pixsize
 refeta = 20.
 eta2 = (10.^(-refeta / 2.5)) * pixsize2 * bin2 / (2.d0 * zs2 * bfactor)
  
 flares02 = alog10(hs2solar / hs2in)
 flares01 = alog10((zs2solar - zs2) / (zs2in - zs2))
 flares = flares01 / flares02
  
 sigtruncate = sha2 * rtrun2
  
; Open output file
 openw, unit, 'emissivity_tdisk2.dat', /get_lun
  
; Main calculation loops
 for i = 0L, nposr - 1L do begin
    rho = i * dr  ; Calculate rho directly from i
    
    for n = 0L, nposz - 1L do begin
       z = n * dz  ; Calculate z directly from n
      
       truncate = 0.5 * (1.0 - erf((rho - rtrun2) / sigtruncate))
      
       ; Calculate second thick disk emissivity based on radial position
       if (rho lt hs2tin) then begin
         etadisk2 = 0
       endif else if (rho lt hs2in) then begin
        etadisk2 = eta2 * (zs2/zs21) * ((rho/hs2in) * (1. - xis2) + xis2) * $
                   exp(-hs2in/hs2 - abs(z)/zs21)
       endif else begin
        etadisk2 = eta2 * (zs2/zs21) * exp(-rho/hs2 - abs(z)/zs21)
       endelse
      
       etadisk2 *= truncate
      
       printf, unit, rho, z, etadisk2
    endfor ; n = 0L, nposz - 1L
 endfor ; i = 0L, nposr - 1L 
  
 free_lun, unit
end
