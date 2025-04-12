pro plot_profile_residuals, rcog, rcog_mod, bcog_dat, bcog_mod_tot, bcog_mod_HII, bcog_mod_irr1, $
                  bcog_mod_irr2, bcog_mod_irr3, scog_dat, scog_mod_tot, scog_mod_HII, $
                  scog_mod_irr1, scog_mod_irr2, scog_mod_irr3, pixel_scale, map_scale_arcsecs, $
                  plot_xrange, plot_yrange, xlab, ylab, gal_name, band, namep, dir_plot
  
  ; Plots SB radial profile of galaxy GAL_NAME
  
  ; Set up display
  loadct, 13
  !p.thick = 3
  !x.thick = 4
  !y.thick = 4
  !p.charthick = 3
  !p.charsize = 1.5
  
  ; Convert radii to arcseconds
  r_arcsecs = pixel_scale * rcog
  r_map_arcsecs = pixel_scale * rcog_mod
  
  ; Effective radius
  ds_dat = abs(scog_dat - 0.5*max(scog_dat))
  ds_mod = abs(scog_mod_tot - 0.5*max(scog_mod_tot))
  ids_dat = where(ds_dat eq min(ds_dat))
  ids_mod = where(ds_mod eq min(ds_mod))
  
  reff_dat = r_arcsecs[ids_dat]
  reff_mod = r_map_arcsecs[ids_mod]
  
  ; Position settings for plot labels
  nlab = 9L
  oft_x = 0.62d0
  x_lab = oft_x * plot_xrange[1]
  oft_y = 0.9d0
  d_lab = 0.10d0
  y_lab = (oft_y - d_lab*findgen(nlab)) * plot_yrange[1]
  
  ; Create main plot
  set_plot, 'PS'
  device, file=dir_plot + namep + '.ps', /color
  
  ; Linear plot
  plot, r_arcsecs, bcog_dat, ystyle=1, xstyle=1, xrange=plot_xrange, yrange=plot_yrange, $
        ytit=ylab, tit=gal_name, position=[0.15, 0.35, 0.97, 0.93], xtickformat="(A1)"
  
  ; Plot effective radius
  oplot, [reff_dat, reff_dat], plot_yrange, linestyle=1
  oplot, [reff_mod, reff_mod], plot_yrange, linestyle=1, color=250
  
  ; Labels
  xyouts, x_lab, y_lab[0], band
  xyouts, x_lab, y_lab[2], 'Total', color=250
  xyouts, x_lab, y_lab[3], 'Old Dust Disk', color=150
  xyouts, x_lab, y_lab[4], 'Young Dust Disk', color=100
  xyouts, x_lab, y_lab[5], 'HII', color=200
  xyouts, x_lab, y_lab[6], 'Young Dust Ring', color=50
  
  ; Plot model components
  oplot, r_map_arcsecs, bcog_mod_tot, color=250
  oplot, r_map_arcsecs, bcog_mod_HII, color=200
  oplot, r_map_arcsecs, bcog_mod_irr1, color=150
  oplot, r_map_arcsecs, bcog_mod_irr2, color=100
  oplot, r_map_arcsecs, bcog_mod_irr3, color=50
  
  ; Plot residuals
  plot_yrange[0] = 0.1
  plot, r_map_arcsecs, (bcog_dat-bcog_mod_tot)/bcog_dat*100, yrange=[-100, 100], $
        position=[0.15, 0.15, 0.97, 0.35], xrange=plot_xrange, xtit='arcsecs', $
        ytit='R (%)', ytickinterval=100, /NOERASE
  plots, [0, 500], [0, 0], linestyle=1
  plots, [0, 500], [20, 20], linestyle=1
  plots, [0, 500], [-20, -20], linestyle=1
  
  device, /close
  
  ; Log plot
  device, file=dir_plot + namep + '_log.ps', /color
  plot, r_arcsecs, bcog_dat, ystyle=1, xstyle=1, xrange=plot_xrange, yrange=plot_yrange, $
        ytit=ylab, tit=gal_name, /ylog, xtickformat="(A1)", position=[0.15, 0.35, 0.97, 0.93]
  
  ; Plot model components
  oplot, r_map_arcsecs, bcog_mod_tot, color=250
  oplot, r_map_arcsecs, bcog_mod_HII, color=200
  oplot, r_map_arcsecs, bcog_mod_irr1, color=150
  oplot, r_map_arcsecs, bcog_mod_irr2, color=100
  oplot, r_map_arcsecs, bcog_mod_irr3, color=50
  
  ; Labels for log plot
  oft_x = 0.60d0
  x_lab = oft_x * plot_xrange[1]
  oft_y = 0.9
  d_lab = 1.5
  y_lab = oft_y / (d_lab^findgen(nlab)) * plot_yrange[1]
  
  xyouts, x_lab, y_lab[1], band
  xyouts, x_lab, y_lab[2], 'Total', color=250
  xyouts, x_lab, y_lab[3], 'Old Dust Disk', color=150
  xyouts, x_lab, y_lab[4], 'Young Dust Disk', color=100
  xyouts, x_lab, y_lab[5], 'HII', color=200
  xyouts, x_lab, y_lab[6], 'Young Dust Ring', color=50
  
  ; Plot residuals
  plot, r_arcsecs, (bcog_dat-bcog_mod_tot)/bcog_dat*100, yrange=[-100, 100], $
        position=[0.15, 0.15, 0.97, 0.35], xrange=plot_xrange, xtit='arcsecs', $
        ytit='R (%)', ytickinterval=100, /NOERASE
  plots, [0, 500], [0, 0], linestyle=1
  plots, [0, 500], [20, 20], linestyle=1
  plots, [0, 500], [-20, -20], linestyle=1
  
  device, /close
  
  return
end
