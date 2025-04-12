pro ssfr
;
; Plots the radial profile for the surface density of the specfic star formation
; rate using the outputs from the ETA program, and saves the result to SSFR.SAV
;
restore,'sfr.xdr'
restore,'stellar_mass.xdr'

!p.thick=5
!x.thick=5
!y.thick=5
!p.charthick=5
!p.charsize=2.0
fname = 'ssfr.ps'

Device = !D.NAME
SET_PLOT, 'PS'
DEVICE, FILENAME = fname, set_font = 'HELVETICA'

ssfr_r = FLUX_TOT_SFR/stellar_mass
ind=where(finite(ssfr_r) eq 1)

plot,r(ind)/1e3,ssfr_r(ind),/ylog,xtit='R (kpc)',ytit = 'sSFR (yr!E-1!N)',xrange=[0,14], xstyle=1
save,r,ssfr_r,filename='ssfr.sav'

DEVICE, /CLOSE_FILE

end
