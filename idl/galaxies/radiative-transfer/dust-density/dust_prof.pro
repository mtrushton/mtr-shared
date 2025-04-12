pro dust_prof, nposr, nposz
  ; Define the dust distribution and calculate opacities
  ;
  ; Input below the dust geometry for the galaxy as defined  
  ; in the GEOMETRY.IN file
  ;
  ; Define parameters for dust components
  params = {$
    ; First component (thick disk)
    xid0: 0.7, hd: 12e3, hdin: 2.0e3, hdtin: 0.0, hdsolar: 12.e3, $
    zd: 0.16e3, zdin: 0.16e3, zdsolar: 0.16e3, rtrun: 17.e3, tau0: 1.44, $
    
    ; Second component (thin disk)
    xid1: 0.7, hd1: 3.2e3, hd1in: 2.0e3, hd1tin: 0.0, hd1solar: 3.2e3, $
    zd1: 0.09e3, zd1in: 0.09e3, zd1solar: 0.09e3, rtrun1: 14.e3, tau1: 0.56, $
    
    ; Third component (thin ring)
    xid2: 0.0001, hd2: 0.2e3, hd2in: 0.5e3, hd2solar: 0.2e3, $
    zd2: 0.09e3, zd2in: 0.09e3, zd2solar: 0.09e3, rtrun2: 20.e3, tau2: 0.0, $
    
    ; General parameters
    sha: 0.02, dr: 10., dz: 50., nposr: 2100, nposz: 12, $
    kappa: 8.07, z_t: 84e3, R_out: 5.8e3}
  
  nposr = params.nposr
  nposz = params.nposz
  
  ; Calculate opacity coefficients
  a0 = params.tau0 / (2.0d0 * params.hd)
  a1 = params.tau1 / (2.0d0 * params.hd1)
  a2 = params.tau2 / (2.0d0 * params.hd2)
  
  ; Create output files
  fileNames = ['opacity.dat', 'opacity0.dat', 'opacity1.dat', 'opacity2.dat']
  openw, unit, fileNames[0], /get_lun
  openw, unit0, fileNames[1], /get_lun
  openw, unit1, fileNames[2], /get_lun
  openw, unit2, fileNames[3], /get_lun
  
  ; Calculate sigtruncate values for each component
  sigtruncate = params.sha * params.rtrun
  sigtruncate1 = params.sha * params.rtrun1
  sigtruncate2 = params.sha * params.rtrun2
  
  for i = 0, params.nposr-1 do begin
    rho = i * params.dr
    
    for n = 0, params.nposz-1 do begin
      z = n * params.dz
      
      ; Component 1: Thick disk
      truncate = 0.5 * (1.0 - erf((rho - params.rtrun) / sigtruncate))
      
      case 1 of
        (rho lt params.hdtin): adum = 0
        (rho lt params.hdin): adum = a0 * (params.zd/params.zd) * $
                             ((rho/params.hdin) * (1. - params.xid0) + params.xid0) * $
                             exp(-params.hdin/params.hd - abs(z)/params.zd)
        else: adum = a0 * (params.zd/params.zd) * exp(-rho/params.hd - abs(z)/params.zd)
      endcase
      adum *= truncate
      
      ; Component 2: Thin disk
      truncate1 = 0.5 * (1.0 - erf((rho - params.rtrun1) / sigtruncate1))
      
      case 1 of
        (rho lt params.hd1tin): adum1 = 0
        (rho lt params.hd1in): adum1 = a1 * (params.zd1/params.zd1) * $
                              ((rho/params.hd1in) * (1. - params.xid1) + params.xid1) * $
                              exp(-params.hd1in/params.hd1 - abs(z)/params.zd1)
        else: adum1 = a1 * (params.zd1/params.zd1) * exp(-rho/params.hd1 - abs(z)/params.zd1)
      endcase
      adum1 *= truncate1
      
      ; Component 3: Thin ring
      truncate2 = 0.5 * (1.0 - erf((rho - params.rtrun2) / sigtruncate2))
      
      if (rho lt params.hd2in) then begin
        adum2 = a2 * (params.zd2/params.zd2) * $
               ((rho/params.hd2in) * (1. - params.xid2) + params.xid2) * $
               exp(-params.hd2in/params.hd2 - abs(z)/params.zd2)
      endif else begin
        adum2 = a2 * (params.zd2/params.zd2) * exp(-rho/params.hd2 - abs(z)/params.zd2)
      endelse
      adum2 *= truncate2
      
      ; Total opacity
      a = adum + adum1 + adum2
      
      printf, unit, rho, z, a
      printf, unit0, rho, z, adum
      printf, unit1, rho, z, adum1
      printf, unit2, rho, z, adum2
    endfor ; n = 0, params.nposz-1 
  endfor ; i = 0, params.nposr-1
  
  free_lun, unit
  free_lun, unit0
  free_lun, unit1
  free_lun, unit2
  
  ; Calculate global parameters
  outfile = 'global_param.dat'
  openw, ounit, outfile, /get_lun
  
  M = dblarr(3)
  tau_in = dblarr(3)
  tau_out = dblarr(3)
  
  for j = 0, 2 do begin
    ; Set component-specific parameters
    if (j eq 0) then begin
      tau_c = params.tau0
      z = params.zd
      h = params.hd
      R_in = params.hdin
      R_trun = params.rtrun
      R_tin = params.hdtin
      xi = params.xid0
    endif else if (j eq 1) then begin
      tau_c = params.tau1
      z = params.zd1
      h = params.hd1
      R_in = params.hd1in
      R_trun = params.rtrun1
      R_tin = params.hd1tin
      xi = params.xid1
    endif else begin
      tau_c = params.tau2
      z = params.zd2
      h = params.hd2
      R_in = params.hd2in
      R_trun = params.rtrun2
      R_tin = 0.
      xi = params.xid2
    endelse
    
    rho_c = tau_c / (2 * params.kappa * z)
    
    ; Calculate mass and opacity
    T_R = exp(-(R_in/h)) - exp(-(R_trun/h)) + $
          (R_in/h) * exp(-(R_in/h)) - (R_trun/h) * exp(-(R_trun/h))
    T_z = 1 - exp(-(params.z_t/z))
    
    M[j] = 4. * !pi * rho_c * $
          ((1./3.) * ((1 + (xi/2)) * R_in^2 - (1-xi) * (R_tin^3/R_in) - ((3./2.) * xi * R_tin^2)) * $
          exp(-(R_in/h)) + (h^2 * T_R)) * z * T_z
    
    tau_in[j] = tau_c * exp(-R_in/h)
    tau_out[j] = tau_c * exp(-params.R_out/h)
    
    ; Write results to file
    printf, ounit, 'component'
    printf, ounit, j+1
    printf, ounit, 'dust mass M_sol'
    printf, ounit, M[j]
    printf, ounit, 'dust opacity at R_in'
    printf, ounit, tau_in[j]
    printf, ounit, ''
    printf, ounit, 'dust opacity at R_out'
    printf, ounit, tau_out[j]
    printf, ounit, ''
  endfor ; j = 0, 2 
  
  printf, ounit, 'total dust mass'
  printf, ounit, total(M)
  
  free_lun, ounit
end
