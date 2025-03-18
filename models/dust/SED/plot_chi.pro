pro plot_chi,bfit_para,i,j,error,chimin,sig_meas,meas_noise,psol,bd,thet_gal,chi,thet_gal_app,$
 thet_gal_ran,lam_bd,thet_inc,d_gal_ran,q_uv,q_opt,data,s_model,sigma,bd_c,bd_app,sfr_ran,old_ran,ibound_sfr,$
 ibound_old,itest,dpara,plot_range_min,plot_range_max,n_par,gal_name,para_est,nofs,stze,paratitle,ierr
;
; Plots the chi-squared topology for the free parameters i and j. Defines an array of parameter 
; values and calculates chi-squared for each row or parameter combination, with the other 
; parameters (free or fixed) fixed at their best fitting value. The parameters are plotted 
; in the user specified range PLOT_RANGE_MIN to PLOT_RANGE_MAX, and with the user specified 
; sampling points N_PAR. Contour plots are generated for the parameter combination. The plots 
; show the solution, the estimate of the solution from the solve_tau_f output, and the search 
; area used in quadfit_chi. 
;
; Input parameters
;       BFIT_PARA: Vector of best fitting parameter values.
;       I: Integer identifty of the free parameter to be plotted on the x axis in the
;             BFIT_PARA vector.         
;       J: Integer identifty of the free parameter to be plotted on the y axis in the
;             BFIT_PARA vector.         
;       ERROR: Vector of free parameter uncertainties
;       DPARA: Vector of parameter precisions for the solution.
;       PLOT_RANGE_MIN: Vector of minimum plotted parameter values.
;       PLOT_RANGE_MAX: Vector of maximum plotted parameter values.
;       N_PAR: Vector of the number of plotted points for each parameter.
;       PARA_EST: Vector of parameter estimates returned by solve_tau_f.
;       NOFS: The number of offsets sampled (must be at least 2).
;       STZE: A multiple of the search precision defined in sed_input. This defines
;             the size of the steps away from the solve_tau_f output. 
;       PARATITLE: String vector of parameter names included in the file name.
;       GAL_NAME: String identifier of the galaxy included in the file name.
;       CHIMIN: Value of chi-squared for the solution.
; 
; Output
;       Contour plot file of the form GAL_NAME_PARI_PARJ_CHIPLOT.PS.
;
 set_plot,'PS'
;
 chi_plot=dblarr(n_par(i),n_par(j))
 plot_para_x=dblarr(n_par(i),n_par(j))
 plot_para_y=dblarr(n_par(i),n_par(j))
;
; define the sampling of parameter values, as used in quadfit_chi
;
 nn=2*nofs+1L
;
; define steps, as a multiple of the desired search range,
; as defined in sed_input and used in quadfit_chi.
;
 delta_para_x=stze*dpara(i)
 delta_para_y=stze*dpara(j)
;
; parameter values for chi-squared determination
;
; step lengths for parameters
; 
 delta_x=(plot_range_max(i)-plot_range_min(i))/n_par(i)
 delta_y=(plot_range_max(j)-plot_range_min(j))/n_par(j)
;
; trial values for parameter i, as used in quadfit_chi
;
 trial_para_x=para_est(i)+delta_para_x*(-nofs+findgen(nn))
 trial_para_y=para_est(j)+delta_para_y*(-nofs+findgen(nn))
;
 para_x=plot_range_min(i)+delta_x*(findgen(n_par(i)+1L))
 para_y=plot_range_min(j)+delta_y*(findgen(n_par(j)+1L))
;
; initialise para vector to solution
;
  para=bfit_para
;
; x axis loop
;
for nx=0,n_par(i)-1L do begin
;
; increment in x direction
;
  para(i)=para_x(nx)
;
; y axis loop
;
for ny=0,n_par(j)-1L do begin
;
; increment in y direction
;
 para(j)=para_y(ny)
;
; calculate chi-squared for para vector
; 
findchi_all,para(0),para(1),para(2),para(3),sig_meas,meas_noise,$
 psol,bd,thet_gal,chi,thet_gal_app,thet_gal_ran,lam_bd,thet_inc,$
 d_gal_ran,q_uv,q_opt,data,s_model,sigma,bd_c,sfr_ran,old_ran,$
 bd_app,ibound_sfr,ibound_old,itest,ierr
;
 plot_para_x(nx,ny)=para(i)
 plot_para_y(nx,ny)=para(j)
 chi_plot(nx,ny)=chi
;
endfor ; ny=0,n_par(j)-1L
endfor ; nx=0,n_par(i)-1L
;
; set up contours
; linear contours, straddling 0.0 with negative contour
;
 nc=41L ; number of contours
 ncneg=1L
 lev=max(chi_plot)*(findgen(nc)-float(ncneg-1L)-0.5)/float(nc-1) ; contour levels
;
 levstart=0.2d0
 nc=11L ; number of contours
 ncneg=1L ; must be 1 for log contours
 g=(1.0d0/(float(nc)-2.0d0))*alog10(max(chi_plot)/levstart)
 g=10.0d0^g
 levdum=levstart*g^findgen(nc) ; contour levels
 lev=[-levstart,levdum[0:nc-2]] ; contour levels from negative levstart 
;
; assign the file name to the contour plot
;
 fname=gal_name+'_'+string(paratitle(i))+'_'+string(paratitle(j))+'_chiplot'
 device,file=fname+'.ps'
;
; contour plot
;
 contour,chi_plot,plot_para_x,plot_para_y,xst=1,yst=1,xtit=paratitle(i),ytit=paratitle(j),$
  c_linest=1.,lev=lev[ncneg:nc-1]
;
; plot search area as used in quadfit_chi call
;
 plots,[trial_para_x(0),trial_para_x(nn-1),trial_para_x(nn-1),trial_para_x(0),trial_para_x(0)],$
  [trial_para_y(0),trial_para_y(0),trial_para_y(nn-1),trial_para_y(nn-1),trial_para_y(0)]
;
; plot the estimate of the solution from solve_tau_f output
;
 plots,para_est(i),para_est(j),psym=4
;
; plot best fitting parameters
;
 plots,bfit_para(i),bfit_para(j),psym=1
;
; draw error bars on best fit solution
;
 plots,[bfit_para(i)-error(i),bfit_para(i)+error(i)],[bfit_para(j),bfit_para(j)]
 plots,[bfit_para(i),bfit_para(i)],[bfit_para(j)-error(j),bfit_para(j)+error(j)]
;
 device,/close
;
c999:
;
return 
;
end
