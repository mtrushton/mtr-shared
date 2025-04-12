pro calc_sfr
;
; Calculate the SFR for an RT modelled galaxy using outputs
; of the ETA program
;
restore, 'sfr.xdr'

nposr = n_elements(r)
lumdouble = FLUX_TOT_SFR

pass_energy = 0.0
;
; Integrate the surface density of sfr over radius
;
for kr = 0L,nposr - 2L do begin
   pass = !pi * (r(kr+ 1)^2 - r(kr)^2) $
    * (lumdouble(kr) + lumdouble(kr + 1)) / 2.
   pass_energy = pass_energy + pass  ;in W/Hz
endfor ; kr = 0L,nposr - 2L
print,pass_energy

end
