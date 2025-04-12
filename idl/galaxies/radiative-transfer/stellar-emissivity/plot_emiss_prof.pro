pro plot_emiss_prof, npos, nr, nz, ncomp, file, scale

  ; Scales component emissivities to physical units and plots radial profile for z = 0
  
  rr = dblarr(npos)
  zz = dblarr(npos)
  em = dblarr(npos)
  total_em = dblarr(nr)
  emiss_comp = dblarr(nr, ncomp)
  r = dblarr(nr)
  
  ; Read data and scale emissivities
  for i = 0, ncomp - 1 do begin
    openr, unit, file[i], /get_lun
    
    for n = 0L, npos - 1L do begin
      readf, unit, dum, dum1, dum2
      rr[n] = dum
      zz[n] = dum1
      em[n] = dum2
    endfor ; n = 0L, npos - 1L
    
    ; Extract radial profile at z = 0
    for n = 0L, nr - 1L do begin
      r[n] = rr[n*nz]
      emiss_comp[n, i] = em[n*nz] * scale[i]
    endfor ; n = 0L, nr - 1
    
    free_lun, unit
  endfor ; i = 0, ncomp - 1
  
  ; Calculate total emission
  for n = 0, nr - 1 do begin
    total_em[n] = total(emiss_comp[n, *])
  endfor ; n = 0, nr - 1
  
  ; Plot settings
  plot_label = ['g SDSS', 'bulge', 'stellar disk', 'thin stellar disk', 'outer stellar disk', $
               'inner thin stellar disk', 'inner stellar', 'outer thin', 'nuclear']
  plot_comp_colour = [50, 130, 90, 130, 90, 130, 90, 90]
  plot_label_colour = [0, 50, 130, 90, 130, 90, 90, 130, 90]
  plot_line_style = [0, 2, 2, 2, 2, 2, 2, 2, 2]
  plot_comp_line_style = [2, 2, 2, 2, 2, 2, 2, 2]
  
  ; Set plotting parameters
  !p.thick = 4
  !x.thick = 5
  !y.thick = 5
  !p.charthick = 5
  !p.charsize = 1.5
  loadct, 13
  
  ; Generate full range plot
  set_plot, 'PS'
  device, file='total_b_emissivity_log.ps', /color
  
  xmin = 0
  xmax = 30
  ymin = 1e8
  ymax = 2e14
  
  plot, r/1e3, total_em, yrange=[ymin, ymax], /ylog, ystyle=1, xrange=[xmin, xmax], $
        xtit='R [kpc]', ytit='W/Hz/pc!U3', title='M101'
  
  for i = 0, ncomp - 1 do begin
    oplot, r/1e3, emiss_comp[*, i], color=plot_comp_colour[i], linestyle=plot_comp_line_style[i]
  endfor ;  i = 0, ncomp - 1 
  
  oplot, r/1e3, total_em
  
  ; Add legend
  nlab = n_elements(plot_label)
  plot_xrange = xmax - xmin
  plot_label_xfrac = 0.50
  plot_label_xpos = plot_label_xfrac * xmax
  bar_length = 0.07 * plot_xrange
  lab_oft_x = 0.03 * plot_xrange
  plot_label_xpos2 = plot_label_xpos + bar_length
  
  ; Calculate y positions for legend (logarithmic spacing)
  oft_ymax = (alog10(ymax) - alog10(ymin)) * 0.93
  delta_y = alog10(ymin) + oft_ymax
  oft_y = 10^delta_y
  bs = (alog10(ymax) - alog10(oft_y)) * findgen(nlab)
  y_lab = oft_y / 10^bs
  
  ; Draw legend
  for n = 0, nlab - 1 do begin
    if(n ne 0 and n le 3) then begin
      plots, [plot_label_xpos, plot_label_xpos2], y_lab[n], $
             color=plot_label_colour[n], linestyle=plot_line_style[n], thick=6
    endif ; n ne 0 and n le 3
    if(n le 3) then xyouts, plot_label_xpos2 + lab_oft_x, y_lab[n], plot_label[n]
  endfor ; n = 0, nlab - 1 
  
  device, /close
  
  ; Generate zoomed plot
  set_plot, 'PS'
  device, file='total_b_emissivity_zoom_log.ps', /color
  
  xmin = 0
  xmax = 2.5
  
  plot, r/1e3, total_em, yrange=[ymin, ymax], /ylog, ystyle=1, xrange=[xmin, xmax], $
        xtit='R [kpc]', title='M101', ytickformat='(A1)'
  
  for i = 0, ncomp - 1 do begin
    oplot, r/1e3, emiss_comp[*, i], color=plot_comp_colour[i], linestyle=plot_comp_line_style[i]
  endfor ; i = 0, ncomp - 1 
  
  device, /close
end
