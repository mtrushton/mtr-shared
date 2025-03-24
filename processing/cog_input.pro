pro cog_input,gal_name,fname,mname,dist_gal,band,pix_size,wa,sky_sub,r_max,ps,comps,sky_lev, $
     plot_reff,plot_rmask,rmask,plot_xrange_prof, plot_yrange_prof, plot_xrange_cog,plot_yrange_cog, $
     xlab,prof_ylab,cog_ylab,sky,plot_log,x_lab_prof,x_lab_cog,y_lab_prof,y_lab_cog,ierr
;
; Input the parameters for the surface brightness profile and curve of growth (cog) analysis for galaxy gal_name.
;
; INPUTS:
; None - All parameters are set within this function
;
; OUTPUTS:
;  GAL_NAME: Galaxy identifier string to be included in the plot file names.
;  FNAME: Galfit fits file output string. This should include the data and fitted image.
;  MNAME: Mask file name for stars (0 for bad pixels, 1 for good pixels)
;  DIST_GAL: Galaxy distance in Mpc
;  BAND: String containing the observed band (U,B,V,R,I,JHK)
;  PIX_SIZE: Scalar data pixel size in arcsecs
;  WA: Scalar width of the annuli in arcsecs for the brightness profiles and cog.
;  SKY_SUB: String. Equal to 'y' for sky subtracted surface brightness profile and cog plots
;  SKY: Scalar value of the sky
;  PS: String. Equal to 'y' for postscript plots
;  COMPS: String. Equal to 'y' for plots with bulge and disk components
;  R_MAX: Integer. Radius in pixels beyond which it is judged there to be no galaxy emission.
;  SKY_LEV: String. Equal to 'y' for plotting a horizontal line at the sky level on the surface brightness profile plot
;  PLOT_REFF: String. Equal to 'y' for plotting a vertical line on the cog at the radius where counts=0.5*maxCounts
;  PLOT_RMASK: String. Equal to 'y' if a limiting radius has been used in the fitting.
;  RMASK: Integer. The limiting radius in arcsecs if used in the fitting.
;  PLOT_XRANGE_PROF: 2D vector of the x range for the surface brightness profile plot [xmin,xmax].
;  PLOT_YRANGE_PROF: 2D vector of the y range for the surface brightness profile plot [ymin,ymax].
;  PLOT_XRANGE_COG: 2D vector of the x range for the cog plot [xmin,xmax].
;  PLOT_YRANGE_COG: 2D vector of the y range for the cog plot [ymin,ymax].
;  XLAB: String label for the x axis of the surface brightness profile and cog plots.
;  PROF_YLAB: String label for the y axis of the surface brightness profile plot. 
;  COG_YLAB: String label for the y axis of the cog plot. 
;  X_LAB_PROF: Vector of float values denoting the position of labels on the x axis of the brightness profile plot.
;  X_LAB_COG: Vector of float values denoting the position of labels on the x axis of the cog plot.
;  Y_LAB_PROF: Vector of float values denoting the position of labels on the y axis of the brightness profile plot.
;  Y_LAB_COG: Vector of float values denoting the position of labels on the y axis of the cog plot.
;  PLOT_LOG: Integer set to 0 for linear plots in y axis and 1 for a logarithmic y axis.
;  IERR: Error flag (0 = no error, 1 = error)
;
; HISTORY:
; 03/18/2025 - Updated with improved error handling and documentation
;
 ierr=0L
 
 ; Set default values where appropriate
 def_pix_size = 1.0
 def_wa = 5.0
 def_sky = 0.0
 def_rmask = 250L
 def_plot_log = 0L
 
 ; Galaxy parameters
 gal_name='NGC628' ; galaxy identifier
 fname='N0628_J_subcomps_fit.fits' ; data image file. Includes data and fitted image generated using subcomps in GALFIT
 mname='star_remove/mask.fits' ; name of the data image masked for stars (opposite GALFIT convention: 0 for bad, 1 for good)
 dist_gal=9.5d0 ; galaxy distance in Mpc
 band='J' ; band of observation
 
 ; Image parameters
 pix_size=def_pix_size ; size of the data pixels in arcsecs
 wa=def_wa ; width of the annuli for the brightness profiles and cog in arcsecs
 sky=def_sky ; sky value shown as a horizontal line on the brightness profile
 r_max=75L ; radius beyond which there is no galaxy emission (pixels)
 rmask=def_rmask ; radius in pixels beyond which no data is included in the fit
 
 ; Plot control parameters
 sky_sub='n' ; equal to 'y' for sky subtracted surface brightness profile and cog plots
 ps='y' ; equal to 'y' for postscript plots
 comps='n' ; equal to 'y' for plots with bulge and disk components
 sky_lev='n' ; equal to 'y' for plotting sky level on the surface brightness profile
 plot_reff='n' ; equal to 'y' for plotting a vertical line on the cog at half max counts
 plot_rmask='n' ; equal to 'y' if a limiting radius has been used in the fitting
 plot_log=def_plot_log ; 0 for linear plots and 1 for a logarithmic y axis
 
 ; Plot ranges and labels
 plot_xrange_prof=[0,1000] ; x range for the surface brightness profile plot
 plot_yrange_prof=[-0.3,0.74] ; y range for the surface brightness profile plot
 plot_xrange_cog=[0,1000] ; x range for the cog plot
 plot_yrange_cog=[0.4,1d7] ; y range for the cog plot
 
 ; Axis labels
 xlab='arcsecs' ; x label for the surface brightess profile and cog plots
 prof_ylab='counts/sec/pixel' ; y label for the surface brightness profile plot
 cog_ylab='counts/sec' ; y label for the cog plot
 
 ; Plot labels (x and y positions)
 nlab=8L ; number of labels on the plot
 
 ; Surface profile plot parameters
 oft_prof_x=0.84d0 ; determines the position of the label along the x axis
 x_lab_prof=oft_prof_x*plot_xrange_prof(1) ; position of labels on x axis
 oft_prof_y=0.9d0 ; determines the first label position along the y axis
 d_lab_prof=0.06d0 ; label offset on the y axis
 y_lab_prof=(oft_prof_y-d_lab_prof*findgen(nlab))*plot_yrange_prof(1) ; positions of the label on the y axis
 
 ; Adjust parameters for logarithmic plots
 if(plot_log eq 1) then begin
   oft_prof_x=0.81d0; determines the position of the label along the x axis
   x_lab_prof=oft_prof_x*plot_xrange_prof(1)
   oft_prof_y=0.5 ; determines the position of the first label along the y axis
   d_lab_prof=1.6d0 ; label offset on the y axis
   y_lab_prof=oft_prof_y/(d_lab_prof^findgen(nlab))*plot_yrange_prof(1) ; positions of the label on the y axis
 endif
 
 ; COG plot parameters
 oft_cog_x=0.78d0; determines the position of the label along the x axis
 x_lab_cog=oft_cog_x*plot_xrange_cog(1)
 oft_cog_y=5d-5 ; determines the position of the first label along the y axis
 d_lab_cog=2.2d0 ; label offset on the y axis
 y_lab_cog=oft_cog_y/(d_lab_cog^findgen(nlab))*plot_yrange_cog(1) ; positions of the label on the y axis
 
 ; Input validation
 if(sky lt plot_yrange_prof(0) or sky gt plot_yrange_prof(1)) then begin
   print, 'COG_INPUT: Warning - sky value outside plotted y range'
 endif
 
 ; Convert string parameters to lowercase for consistent comparisons
 sky_sub = strlowcase(sky_sub)
 ps = strlowcase(ps)
 comps = strlowcase(comps)
 sky_lev = strlowcase(sky_lev)
 plot_reff = strlowcase(plot_reff)
 plot_rmask = strlowcase(plot_rmask)
 
 ; Validate yes/no parameters
 if(sky_sub ne 'y' and sky_sub ne 'n') then begin
   print, 'COG_INPUT: Error - sky_sub must be "y" or "n"'
   ierr=1L
 endif
 
 if(ps ne 'y' and ps ne 'n') then begin
   print, 'COG_INPUT: Error - ps must be "y" or "n"'
   ierr=1L
 endif
 
 if(comps ne 'y' and comps ne 'n') then begin
   print, 'COG_INPUT: Error - comps must be "y" or "n"'
   ierr=1L
 endif
 
 if(plot_log ne 0L and plot_log ne 1L) then begin
   print,'COG_INPUT: Error - plot_log must be equal to 0 or 1'
   ierr=1L
 endif
 
 ; File existence check
 if ~file_test(fname) then begin
   print, 'COG_INPUT: Error - GALFIT output file '+fname+' not found'
   ierr=1L
 endif
 
 if ~file_test(mname) then begin
   print, 'COG_INPUT: Warning - Mask file '+mname+' not found'
 endif
 
 return
end
