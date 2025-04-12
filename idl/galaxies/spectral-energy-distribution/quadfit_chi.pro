pro quadfit_chi,chi_input,fpara,para,para_ran_min,para_ran_max,dpara,eps,thet_gal_ran,thet_gal_app,lam_bd,$
    bd_app,thet_inc,d_gal_ran,sfr_ran,old_ran,nofs,stze,tol,para_sol,pol_coeff_qa,ierr
;
; Defines a grid of values in chi-squared for a single free parameter by offsetting the value of the 
; parameter from the estimated value of the parameter at the minimum. A quadratic fit of the form 
; y=a+bx+cx^2 is performed over the parameter range, giving a more precise best fitting value at the 
; minimum of the quadratic. A test is performed over the range max_para to min_para to see if the 
; chi-squared topology is quadratic. To check that this is the case two further quadratics are fitted 
; to the parameter space: one to the three central data points and a second to the minimum, central 
; and maximum data points. The quadratic approximation is considered valid if the c coefficients agree 
; to within the tolerance specified by TOL.
;
; Input parameters
;       FPARA: Value of the free parameter at the estimated minimum chi-squared
;       PARA: Values of fixed parameters at the estimated minimum chi-squared
;       PARA_RAN_MIN: Minimum allowed value in the search range of FPARA
;       PARA_RAN_MAX: Maximum allowed value in the search range of FPARA
;       DPARA: desired precision of the solution of FPARA
;       CHI_INPUT: Value of chi-squared at input values of FPARA and PARA 
;       NOFS: The number of offsets sampled (must be at least 2)
;       STZE: A multiple of the search precision defined in sed_input. This defines
;             the size of the steps away from the solve_tau_f output. 
;       TOL: A value above which quadratic approximation in parameter space is invalid. 
;       EPS: Machine precision.
;
; Output parameters
;       PARA_SOL: best fitting FPARA value 
;       POL_COEFF_QA: Coefficents of the quadratic, which has been fit to all sampled 
;                     FPARA space
;
; check free parameter constraints on input
;
if(fpara lt para_ran_min or fpara gt para_ran_max) then begin
 print,'quadfit_chi: parameter out of range on input - aborting'
 ierr=1.
goto, c999
endif; fpara lt par_ran_min or fpara gt par_ran_max
;
; check that there are sufficient offsets
;
if(nofs lt 2) then begin
 print,'quadfit_chi: insufficient offsets - aborting'
 ierr=1.
goto, c999
endif; nofs lt 2
;
 nn=2*nofs+1L
 chitest=dblarr(nn)
 chitest(nofs)=chi_input
;
; define steps, as a multiple of the desired search range
; as defined in sed_input.
;
 delta_para=stze*dpara
;
; define parameters offsets for sampling chi-squared
;
 trial_para=fpara+delta_para*(-nofs+findgen(nn))
; 
; ensure constraints
;
 max_para=max(trial_para)
 min_para=min(trial_para)
;
; check if para offsets are inside allowed search range and, if not, 
; adjust offset
;
if(max_para gt para_ran_max) then begin
 cpara=(para_ran_max-min_para)/2.
 delta_para=(cpara-para_ran_min)/nofs
 trial_para=cpara+delta_para*(-nofs+findgen(nn))
 chitest(nofs)=0L
 print,'sedfit: max offset out of range, redefining max_para'
endif ; max_para gt para_ran_max 
;
if(min_para lt para_ran_min) then begin
 cpara=(para_ran_min-max_para)/2.
 delta_para=(cpara-para_ran_min)/nofs
 trial_para=cpara+delta_para*(-nofs+findgen(nn))
 chitest(nofs)=0L
 print,'sedfit: min offset out of range, redefining min_para'
endif ; min_para lt para_ran_min
;
; check step length is larger than machine precision
;
if(delta_para lt eps) then begin 
 print,'quadfit_chi: trial steps are smaller than machine precision - aborting'
 ierr=1.
 goto, c999
endif ; delta_para lt eps
;
 inchi=where(chitest eq 0, nchi) ; 
 fp=where(para eq 0) ; element for free parameter
;
for i=0,nchi-1 do begin
;
 n=inchi(i)
;
 para(fp)=trial_para(n)
;
findchi_all,para(0),para(1),para(2),para(3),bd,thet_gal,chi,$
 thet_gal_app,thet_gal_ran,lam_bd,thet_inc,d_gal_ran,$
 q_uv,q_opt,data,s_model,sigma,bd_c,sfr_ran,old_ran,bd_app,$
 ibound_sfr,ibound_old,ierr
;
 chitest(n)=chi
;
endfor ; n=0,nchi-1
;
; fit three quadratics (a,b,c) to the parameter space
;
; fit quadratic "a" to all the parameter space
;
 pol_coeff_qa=poly_fit(trial_para,chitest,2)
;
; create a new array in chi-squared and FPARA,
; considering 3 points, one either side and adjacent to the minimum
;
 chitest_qb=[chitest(nofs-1),chitest(nofs),chitest(nofs+1)]
 trial_para=[trial_para(nofs-1),trial_para(nofs),trial_para(nofs+1)]
;
; fit quadratic "b" to this array
;
 pol_coeff_qb=poly_fit(trial_para,chitest_qb,2)
;
; create another array in chi-squared and FPARA, 
; considering 3 points, two at the extremes and one 
; at the minimum
;
 chitest_qc=[chitest(0),chitest(nofs),chitest(nchi)]
 trial_para=[min_para,trial_para(nofs),max_para]
;
; fit quadratic "c" to this array
;
 pol_coeff_qc=poly_fit(trial_para,chitest_qc,2)
;
; compare the b,c quadratics to check approximation is valid
;
 qud=(pol_coeff_qb(2)-pol_coeff_qc(2))^2/pol_coeff_qc(2)^2
;
; find the minimum of quadratic a 
;
 min=(-pol_coeff_qa(1))/(2*pol_coeff_qa(2))
 para_sol=min
;
if(qud gt tol) then begin
 print,'quadfit_chi: quadratic approximation invalid - aborting'
 ierr=1.
goto, c999
endif; qa gt tol
;
c999:
;
return    
end
