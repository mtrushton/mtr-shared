pro calc_stellar_mass
;
; Calculates the total stellar mass for an RT galaxy using the ouputs of the 
; ETA program. 
;
restore, 'stellar_mass.xdr'

nposr = n_elements(r)
lumdouble = stellar_mass

pass_energy = 0.0
;
; Integrate the suface density of the stellar mass over radius
;
for kr = 0L,nposr - 2L do begin
   pass = !pi * (r(kr+ 1)^2 - r(kr)^2) $
    * (lumdouble(kr) + lumdouble(kr + 1)) / 2.
   pass_energy = pass_energy + pass  ;in W/Hz
endfor ; kr = 0L,nposr - 2L
print,'total stellar mass = ', pass_energy, 'solar masses'

end
