pro make_fig_stellar_mass, r, flux_bulge, flux_disk, flux_tdisk, flux_tdisk2, fname, stellar_mass, scale_factor_mass, Mval, scale_mass_prof
;
; Creates radial profile of the surface density of stellar mass. 
;
; radius in pc to kpc
 rkpc = r / 1d3
  
; Constants
 Lsolar = 3.83d26
  
; Combine flux components and convert to stellar mass using L band relation
 flux_tot = flux_bulge + flux_disk + flux_tdisk
 stellar_mass = flux_tot * 0.6 / Lsolar
  
; Calculate total stellar mass using vectorized operations
 nposr = n_elements(r)
 r_diffs = r[1:nposr-1]^2 - r[0:nposr-2]^2
  
; Total luminosity
 avg_flux = (flux_tot[0:nposr-2] + flux_tot[1:nposr-1]) / 2.0
 total_lum = !pi * total(r_diffs * avg_flux)
  
; Calculate initial total stellar mass
 avg_mass = (stellar_mass[0:nposr-2] + stellar_mass[1:nposr-1]) / 2.0
 total_mass = !pi * total(r_diffs * avg_mass)

;
; Scale total stellar mass to desired value if required
; 
 if (scale_mass_prof eq 1) then begin 
     target_mass = Mval  ; Target stellar mass
  
    ; Calculate scaling factor to reach target mass
     scale_factor = target_mass / total_mass
  
    ; Apply scaling factor
     stellar_mass = stellar_mass * scale_factor
 endif ; if scale_mass_prof
  
; Verify final total mass
 avg_mass = (stellar_mass[0:nposr-2] + stellar_mass[1:nposr-1]) / 2.0
 total_mass = !pi * total(r_diffs * avg_mass)
 print, 'stellar_mass = ', total_mass, ' solar masses'
  
; Find index for 17 kpc for y-range calculation
 diff = abs(17 - rkpc)
 ind = where(diff eq min(diff))
 ymin = stellar_mass[ind] * 0.8
 ymax = max(stellar_mass) * 1.2
  
; Plot settings
 !p.thick = 8
 !x.thick = 8
 !y.thick = 8
 !p.charthick = 8
 !p.charsize = 2.0
 chsize = 2
  
; Create plot
 SET_PLOT, 'PS'
 DEVICE, FILENAME=fname, SET_FONT='HELVETICA'
  
; Convert to Msolar per kpc
 stellar_mass = stellar_mass * 1e6
  
 plot, rkpc, stellar_mass, $
     XTITLE='R (kpc)', $
     YTITLE='!7R!X!BM!B*!N (M'+sunsymbol()+' kpc!E-2!N)', $
     CHARSIZE=chsize, $
     /YLOG, $
     YRANGE=[ymin*1e6, ymax*1e6], $
     YSTYLE=1
  
 DEVICE, /CLOSE_FILE
  
 save, r, stellar_mass, FILENAME='stellar_mass.xdr'
end
