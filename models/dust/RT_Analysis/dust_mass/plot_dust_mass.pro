pro plot_dust_mass, nposr, ncomp, file, nr, nz, dust_dens_z, scale, r
; Calculate and plot the surface density dust mass distribution.
; Results saved to 'DUST_MASS.SAV'
  
; Calculate total dust mass by summing relevant components
 dust_mass_total = dblarr(nposr)
 for i = 0, ncomp - 4 do begin
    dust_mass_total = dust_mass_total + dust_dens_z[*, i] * scale[i]
 endfor ; i = 0, ncomp - 4
  
; Calculate total dust mass by integrating over radius
 pass_energy = 0.0
 for kr = 0L, nposr - 2L do begin
    pass = !pi * (r[kr+1]^2 - r[kr]^2) * (dust_mass_total[kr] + dust_mass_total[kr+1]) / 2.
    pass_energy = pass_energy + pass
 endfor ;  kr = 0L, nposr - 2L
 print, 'Total dust mass: ', pass_energy, ' solar masses'
  
; Replace zeros with small values for log plotting
 ind = where(dust_mass_total eq 0.)
 if ind[0] ne -1 then dust_mass_total[ind] = 1e-5
  
; Set up plotting parameters
 !p.thick = 4
 !x.thick = 5
 !y.thick = 5
 !p.charthick = 5
 !p.charsize = 1.5
  
; Define plot range (kpc)
 xmin = 0.
 xmax = 30.
  
; Find value at 14 kpc for scaling reference
 diff = abs(14e3 - r)
 ind = where(diff eq min(diff))
 ymin = dust_mass_total[ind] * 0.03 * 1e6
 ymax = max(dust_mass_total) * 1.3 * 1e6
  
; Create plot
 set_plot, 'PS'
 device, file = 'dust_mass_r.ps', /color
  
; Plot radial profile of dust mass surface density
 plot, r/1e3, dust_mass_total*1e6, yrange=[ymin, ymax], ystyle=1, xrange=[xmin, xmax], $
       xtit='R [kpc]', ytit='!4R!X!ldust!N (M'+sunsymbol() +' kpc!U-2!N)', /ylog
  
; Galaxy identifier
 ya = (alog10(ymax) - alog10(ymin)) * 0.9
 yb = alog10(ymin) + ya
 yc = 10^yb
 xyouts, 25, yc, 'M101'
  
 device, /close

 save, r, dust_mass_total, file = 'dust_mass.sav'
end
