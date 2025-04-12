pro scale_factor,ncomp,total_energy,scale

; Determine scaling factors for converting
; emissivities to physical units. Results
; stored in SCALE vector
;-----------------------------------------
; file containing component luminosities
; (W/Hz) at relevant wavelength, with 
; component elements corresponding to those
; in file vector

restore,'../intrinsic_luminosities/lum_comp.sav'

scale = lum[0:ncomp-1] / total_energy

end
