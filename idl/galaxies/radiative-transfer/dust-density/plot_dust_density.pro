pro plot_dust_density, npos, ncomp, file, nr, nz, scale
  ; Read and plot dust density at z=0
  
  ; Initialise arrays
  r = dblarr(nr)
  dust_density_comp = dblarr(nr, ncomp)
  
  ; Read data for each component
  for i = 0, ncomp - 1 do begin
    rr = dblarr(npos)
    zz = dblarr(npos)
    em = dblarr(npos)
    
    openr, unit, file[i], /get_lun
    readf, unit, rr, zz, em, format='(3F)'
    free_lun, unit
    
    ; Extract radial profile at z=0
    for n = 0, nr - 1 do begin
      r[n] = rr[n*nz]
      dust_density_comp[n, i] = em[n*nz] * scale[i]
    endfor ; n = 0, nr - 1 
  endfor ; i = 0, ncomp - 1 
  
  ; Extract individual components
  dust_density_irr1 = dust_density_comp[*, 0]
  dust_density_irr2 = dust_density_comp[*, 1]
  dust_density_irr3 = dust_density_comp[*, 2]
  
  ; Calculate total dust density
  dust_density_total = dust_density_irr1 + dust_density_irr2 ; + dust_density_irr3
  
  ; Set plot parameters
  !p.thick = 4
  !x.thick = 5
  !y.thick = 5
  !p.charthick = 5
  !p.charsize = 1.5
  
  xmin = 0
  xmax = 20
  ymin = 100
  ymax = 1000
  
  ; Create plot
  set_plot, 'PS'
  device, file = 'dust_density_r.ps', /color
  
  plot, r/1e3, dust_density_total*1e6, yrange=[ymin,ymax], ystyle=1, xrange=[xmin,xmax], $
        xtit='R [kpc]', ytit='!4q!3!Ddust!N(M'+sunsymbol()+' kpc!U-3!N)', /ylog
  
  ya = (alog10(ymax) - alog10(ymin)) * 0.9
  yb = alog10(ymin) + ya
  yc = 10^yb
  xyouts, 15, yc, 'NGC3938'
  
  ; plot individual components (uncomment if required)
  ; oplot, r/1e3, dust_density_irr1*1e6, linestyle=1
  ; oplot, r/1e3, dust_density_irr2*1e6, linestyle=2
  
  device, /close
  
  save, r, dust_density_irr1, dust_density_irr2, dust_density_irr3, dust_density_total, $
        filename='dust_density.sav'
end
