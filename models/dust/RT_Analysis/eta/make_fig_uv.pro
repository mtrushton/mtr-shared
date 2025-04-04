pro make_fig_uv, r, flux_tdisk, flux_tdisk2, flux_tot, fname1, fname2
;
; Creates radial profile of the surface density of SFR from the intrinsic
; galaxy emission integrated over UV-to-B band. 
;
; Radius in pc to kpc
 rkpc = r / 1d3
  
; Constants
 sfr_unit = 2.241e36
 nposr = 2000
  
; Calculate total model flux using array operations
 rrr = rkpc
 lumdouble = flux_tdisk
  
 r_diffs = rrr[1:nposr-1]^2 - rrr[0:nposr-2]^2
 avg_flux = (lumdouble[0:nposr-2] + lumdouble[1:nposr-1]) / 2.0
 total_model = !pi * total(r_diffs * avg_flux)
  
; Scale factor
 scale = 1.10 / total_model
  
; Set total flux
 flux_tot = flux_tdisk
  
; Common plot settings
 common_settings = {$
 thick: 8, $
 charthick: 8, $
 charsize: 2.0, $
 lineThick: 8 $
 }
  
 ; Find index for 14 kpc
 diff = abs(14 - rkpc)
 ind = where(diff eq min(diff))
  
; First plot - raw flux
 SET_PLOT, 'PS'
 DEVICE, FILENAME=fname1, SET_FONT='HELVETICA'
  
 !p.thick = common_settings.thick
 !x.thick = common_settings.thick
 !y.thick = common_settings.thick
 !p.charthick = common_settings.charthick
 !p.charsize = common_settings.charsize
  
; Y-axis limits
 ymin = flux_tot[ind] * 0.8 
 ymax = max(flux_tot) * 1.4
  
; Create first plot
 plot, rkpc, flux_tot, $
     XTITLE='R (kpc)', $
     YTITLE='W/pc!E2', $
     CHARSIZE=common_settings.charsize, $
     /YLOG, $
     YRANGE=[ymin, ymax], $
     YSTYLE=1, $
     XRANGE=[0,15]
  
 oplot, rkpc, flux_tdisk, LINESTYLE=1, THICK=common_settings.lineThick, COLOR=cgcolor('blue')
 oplot, rkpc, flux_tot
  
; legend
 cglegend, TITLES=['Thin Disk', 'Thin Ring'], $
     COLORS=[cgcolor('blue'), cgcolor('green')], $
     LINESTYLES=[1, 1], $
     BOX=0, $
     CHARSIZE=common_settings.charsize, $
     THICK=common_settings.lineThick, $
     SYMTHICK=common_settings.lineThick, $
     LENGTH=0.04, $
     LOCATION=[0.60, 0.85], $
     VSPACE=1.6
  
 DEVICE, /CLOSE_FILE
  
; Convert to SFR units and apply scaling
 flux_tot_sfr = flux_tot * scale
 flux_tdisk_sfr = flux_tdisk * scale
  
; Y-axis limits for second plot
 ymin = flux_tot_sfr[ind] * 0.8
 ymax = max(flux_tot_sfr) * 1.4
  
; Second plot - SFR
 SET_PLOT, 'PS'
 DEVICE, FILENAME=fname2, SET_FONT='HELVETICA'
  
; Create second plot
 plot, rkpc, flux_tot_sfr, $
     XTITLE='R (kpc)', $
     YTITLE='!7R!X!BSFR!N (M'+sunsymbol()+' yr!E-1!Nkpc!E-2!N)', $
     CHARSIZE=common_settings.charsize, $
     /YLOG, $
     YSTYLE=1, $
     YRANGE=[ymin, ymax], $
     XRANGE=[0,15]
 oplot, rkpc, flux_tot_sfr
  
 DEVICE, /CLOSE_FILE
  
 save, r, flux_tot_sfr, flux_tdisk_sfr, FILENAME='sfr.xdr'
end
