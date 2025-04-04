pro emissivity_bulge, nposr, nposz, dr, dz, reff, ellipt, nsersic, pixsize, refeta, shab, rtrunb, bsersic
;
; Old Stellar Bulge emissivity
; 
; Define constants
 ellipt2 = ellipt * ellipt
 pixsize2 = pixsize * pixsize
 etab = 10.^(-refeta / 2.5) * pixsize2
 ref = reff / 3.
 sigtruncate = shab * rtrunb
 rin = 0.4
 esersic = 1.0d0 / float(nsersic)
 dsersic = 2.0d0 * float(nsersic)
 dsersic = (dsersic - 1.0d0) / dsersic
  
; Open output file
 openw, unit, 'emissivity_bulge.dat', /get_lun
  
 for i = 0, nposr-1 do begin
    rho = i * dr  ; Calculate rho directly from i
    
    for n = 0, nposz-1 do begin
       z = n * dz  ; Calculate z directly from n
      
       ; Determine truncation factor and acap based on position
       if(rho lt ref) then begin
         truncate = 0.5 * (1.0 - erf((ref - rtrunb) / sigtruncate))
         acap = sqrt(ref^2 + (z^2 / ellipt2)) / reff
       endif else begin
         truncate = 0.5 * (1.0 - erf((rho - rtrunb) / sigtruncate))
         acap = sqrt(rho^2 + (z^2 / ellipt2)) / reff
       endelse
      
       ; Apply minimum threshold to acap
       acapdum = acap > rin
       if(acapdum eq 0.0d0) then acapdum = 0.4 / reff
      
      ; Calculate bulge emissivity
       dum = bsersic * (acapdum^esersic)
       etabulge = exp(-dum) / (acapdum^dsersic) * etab * truncate
      
       printf, unit, rho, z, etabulge
    endfor ; n = 0, nposz-1 
 endfor ; i = 0, nposr-1 
  
 free_lun, unit
end
