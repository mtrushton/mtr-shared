pro sedfit,tau_b,f,sfr_mpc,old_mpc,bd,thet_gal,chi,$
   thet_gal_app,thet_gal_ran,lam_bd,thet_inc,q_uv,$
   q_opt,d_gal_ran,bd_c,tau_b_ran,dtau_b,f_ran,df,$
   sfr_ran,dsfr,old_ran,dold,bd_app,d_gal,nofs,stze,$
   tol,ibound_sfr,ibound_old,ierr
;
; Optimisation through a search in the free parameters defined in the set
; (TAU_B, F, SFR_MPC AND OLD_MPC) across the SED (UVO/NIR/MIR/SUB-MM). A 
; 4d parameter space in chi-squared is defined by offsetting the values 
; of these parameters from the approximate solution of the parameters 
; returned by solve_tau_f and a test is performed to see if the topology 
; in chi-squared is quadratic. In the case that the approximate solution is 
; sufficiently close to the true solution that the topology in chi-squared 
; can be regarded as quadratic, the value at the minimum of the topology 
; is returned, giving the best fitting value of the parameter. 
; 
; Quadratic fits of the form chi^2=a+bx+cx^2, where x is a free parameter, 
; are performed in the quadfit_chi routine and the minimum of the function 
; is returned to give the best-fitting value. 
;
; Once the minimum has been found, the uncertainties in the free parameters 
; are determined from the Covariance Matrix, which is the inverse square root 
; of the diagonal elements of the Fisher (Curvature) Matrix. The Curvature 
; Matrix is given by F_xy = 1/2 *d^2chi^2/dxdy, where x and y are the free 
; model parameters. The second order partial derivatives of chi-squared wrt
; each parameter combination are determined by offsetting the best-fitting
; values using the precision of the fits specified in the input parameters
; tau_b_ran, f_ran, sfr_ran and old_ran.
;
; For the case of 4 free model parameters, the Curvature Matrix is a 4x4 array:
;
;        F=1/2* 
; |  F11 F12 F13 F14 |
; |  F21 F22 F23 F24 |
; |  F31 F32 F33 F34 |
; |  F41 F42 F43 F44 |
;
; where, for example, F12 is the second order partial derivature of chi-squared
; wrt the free parameters 1 and 2. The diagonal terms (F11, F22, F33, F44) are
; simply given by 2 times the "c" coefficent of the quadratic fit, which is an 
; output of the quadfit_chi routine. The standard errors are given by the 
; Covariance Matrix: covar = ([F_xx]^-1)1/2).
;
; Input parameters
;       TAU_B: estimate for value of TAU_B
;       F: estimate for value of F
;       SFR_MPC: estimate for value of SFR
;       OLD_MPC: estimate for value of OLD
;       DTAU: desired precision for the solution in TAU_B
;       DF: desired precision for the solution in F
;       DSFR:  desired precision for the solution in SFR
;       DOLD: desired precision for the solution in OLD
;       TAU_B_RAN(2): search range in TAU_B. If elements equal, then
;            no search in TAU_B
;       F_RAN(2): search range in F. If elements equal, then no
;              search in F 
;       SFR_RAN(2): search range in SFR. If elements equal, then 
;              no search in SFR
;       OLD_RAN(2): search range in OLD. If elements equal, then no 
;              search in OLD 
;       TOL: The tolerance value above which the quadratic approximation
;            is invalid. 
;
; Output parameters
;       TAU_B_SOL: the scalar best-fitting tau value 
;       F_SOL: the scalar best-fitting f value 
;       SFR_MPC_SOL: the scalar best-fitting sfr value 
;       OLD_MPC_SOL: the scalar best-fitting old value 
;       COVAR: The covariance matrix or standard errors in the free parameters.
;       CHIMIN: the chi-squared minimum for the UVO/NIR/MIR/SUB-MM range.
;
; set machine precision
;
 res=machar(double=double)
 eps=sqrt(res.eps)
;
; Check constraints
; 
if(tau_b lt tau_b_ran(0) or tau_b gt tau_b_ran(1)) then begin
 print,'sedfit: TAU_B out of range on input - aborting'
 ierr=1.
goto,c999
endif; tau_b lt tau_b_ran(0) or tau_b gt tau_b_ran(1) 
;
if(f lt f_ran(0) or f gt f_ran(1)) then begin
 print,'sedfit: F out of range on input - aborting'
 ierr=1.
goto,c999
endif; f lt f_ran(0) or f gt f_ran(1) 
;
if(sfr_mpc lt sfr_ran(0) or sfr_mpc gt sfr_ran(1)) then begin
 print,'sedfit: SFR out of range on input - aborting'
 ierr=1.
goto,c999
endif; sfr_mpc lt sfr_ran(0) or sfr_mpc gt sfr_ran(1) 
;
if(old_mpc lt old_ran(0) or old_mpc gt old_ran(1)) then begin
 print,'sedfit: OLD out of range on input - aborting'
 ierr=1.
goto,c999
endif; old_mpc lt old_ran(0) or old_mpc gt old_ran(1) 
;
; find chi-squared with the input parameters 
; from solve_tau_f for UVO/NIR/FIR/SUB-MM 
;
findchi_all,tau_b,f,sfr_mpc,old_mpc,bd,thet_gal,chi,$
 thet_gal_app,thet_gal_ran,lam_bd,thet_inc,d_gal_ran,$
 q_uv,q_opt,data,s_model,sigma,bd_c,sfr_ran,old_ran,$
 bd_app,ibound_sfr,ibound_old,ierr
;
; print chi-squared with solve_tau output parameters for UVO/NIR/MIR/FIR/SUB-MM SED
; Note that the output below should show chi-squared less than (or equal to) this 
; value
;
print,'sedfit input:'
print,'chi-squared =',chi
;
 s_model_sol=s_model
;
 a=[tau_b,f,sfr_mpc,old_mpc] ;sed_input?
 npara=n_elements(a) ; number of parameters
 bfit_para=fltarr(npara) ; array containing solutions
 para_ran_min=[tau_b_ran(0),f_ran(0),sfr_ran(0),old_ran(0)]
 para_ran_max=[tau_b_ran(1),f_ran(1),sfr_ran(1),old_ran(1)]
 dpara=[dtau_b,df,dsfr*sfr_mpc,dold*old_mpc] ; sed_input?
 para=fltarr(npara) ; array containing fixed parameters on quad call 
;
; identify the free parameters
;
 parafind=indgen(npara) ; array to determine free parameters
;
 if(tau_b_ran(0) eq tau_b_ran(1)) then parafind(0)=0L else parafind(0)=1L
 if(f_ran(0) eq f_ran(1))  then parafind(1)=0L else parafind(1)=1L
 if(sfr_ran(0) eq sfr_ran(1)) then parafind(2)=0L else parafind(2)=1L
 if(old_ran(0) eq old_ran(1)) then parafind(3)=0L else parafind(3)=1L
;
 ifix=where(parafind eq 0, nfix) ; locate fixed parameters
;
; check for free parameters 
;
if(nfix eq npara) then begin
 print,'sedfit: no free parameter - aborting'
 goto, c999
endif
;
 if(ifix gt 0) then bfit_para(ifix)=a(ifix) ;!!!EDIT JULY 2019
; bfit_para(ifix)=a(ifix) ; add fixed parameters to the solutions array
;
 ifree=where(parafind eq 1, nfree) ; locate free parameters
 F=fltarr(npara,npara) ; Curvature Matrix arrays
 diag_F=fltarr(npara) ; array for diagonal elements of the Fisher Matrix
;
for i=0,npara-1 do begin
;
if(parafind(i) eq 0) then goto, c444
;
 para=a ; assign all the parameters to input values
 para(i)=0L ; set the free parameter within the array to zero
;
quadfit_chi,chi,a(i),para,para_ran_min(i),para_ran_max(i),dpara(i),eps,$
 thet_gal_ran,thet_gal_app,lam_bd,bd_app,thet_inc,d_gal_ran,sfr_ran,old_ran,$
 nofs,stze,tol,para_sol,pol_coeff,ierr
; 
if(ierr ne 0) then goto, c999
;
 diag_F(i)=2.0d0*pol_coeff(2)
;
 bfit_para(i)=para_sol
;
c444:
;
endfor ; i=0,npara-1
;
; chi-squared at the solution
;
findchi_all,bfit_para(0),bfit_para(1),bfit_para(2),bfit_para(3),bd,thet_gal,$
 chi,thet_gal_app,thet_gal_ran,lam_bd,thet_inc,d_gal_ran,q_uv,q_opt,data,$
 s_model,sigma,bd_c,sfr_ran,old_ran,bd_app,ibound_sfr,ibound_old,ierr
;
 chimin=chi
;
; calculate second order partial derivatives
; wrt ij free parameter combinations
;
for i=0,npara-1 do begin
;
if(parafind(i) eq 0) then goto, c555
;
for j=0,npara-1 do begin
;
if(parafind(j) eq 0) then goto, c666
;
 if(i eq j) then goto, c666
;
; positive offsets
;
 para(i)=bfit_para(i)+dpara(i)
;
 para(j)=bfit_para(j)+dpara(j)
;
; calculate chi-squared at offsets in parameter
; combinations (positive shift from minimum)
;
findchi_all,para(0),para(1),para(2),para(3),bd,thet_gal,chi,$
 thet_gal_app,thet_gal_ran,lam_bd,thet_inc,d_gal_ran,$
 q_uv,q_opt,data,s_model,sigma,bd_c,sfr_ran,old_ran,bd_app,$
 ibound_sfr,ibound_old,ierr
;
 chi1=chi
 dchi1=chi1-chimin
;
; negative offsets
;
 para(i)=bfit_para(i)-dpara(i)
;
 para(j)=bfit_para(j)-dpara(j)
;
; calculate chi-squared at offsets in parameter
; combinations (negative shift from minimum)
;
findchi_all,para(0),para(1),para(2),para(3),bd,thet_gal,chi,$
 thet_gal_app,thet_gal_ran,lam_bd,thet_inc,d_gal_ran,$
 q_uv,q_opt,data,s_model,sigma,bd_c,sfr_ran,old_ran,bd_app,$
 ibound_sfr,ibound_old,ierr
;
 chi2=chi
 dchi2=chi2-chimin
;
 F(i,j)=(dchi1+dchi2)/(4.0d0*dpara(i)*dpara(j))
;
c666:
;
endfor; j=0,npara-1
c555:
endfor ; i=0,npara-1
;
; define and assign the diagonal elements 
; of the curvature matrix
;
 diag=lindgen(npara)*(npara+1) ;
;
 F(diag)=diag_F
;
; invert and square root diagonal elements of the 
; curvature matrix for 1-sigma errors in free parameters
;
 cmat=invert(F)
 covar=sqrt(abs(cmat(diag)))
;
print,'sedfit output:'
print,'chi-squared =',chimin
print,'Best fit parameters:'
print,bfit_para
print,'Parameter errors:'
print,covar
;;
c999:
;
return    
end
