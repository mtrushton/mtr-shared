pro write_global_params, comps
; Physical constants
 kappa = 8.07    ; opacity in pc^2/M_sol
 z_t = 84*10^3   ; total z-extent in pc
 R_out = 5.8e3   ; outer radius in pc
  
; Initialize arrays
 M = dblarr(6)           ; Component masses
 tau_in = dblarr(6)      ; Inner opacity
 tau_out = dblarr(6)     ; Outer opacity
  
; Open output file
 openw, ounit, 'global_param.dat', /get_lun
  
; Calculate parameters for each component
 for j=0, 5 do begin
   ; Current component params
    tau_c = comps[j].tau
    z = comps[j].zd
    h = comps[j].hd
    R_in = comps[j].hdin
    R_trun = comps[j].rtrun
    R_tin = comps[j].hdtin
    xi = comps[j].xid
    
    ; Central space density
    rho_c = tau_c/(2*kappa*z)
    
    ; Calculate radial and vertical integrals
    T_R = exp(-(R_in/h))-exp(-(R_trun/h))+(R_in/h)*exp(-(R_in/h))-(R_trun/h)*exp(-(R_trun/h))
    T_z = 1-exp(-(z_t/z))
    
    ; Calculate dust mass
    M[j] = 4.*!pi*rho_c*((1./3.)*((1+(xi/2))*R_in^2.-(1-xi)*(R_tin^3./R_in)-((3./2.)*xi*R_tin^2.))*exp(-(R_in/h))+(h^2.*T_R))*z*T_z
    
    ; Calculate opacities
    tau_in[j] = tau_c*exp(-R_in/h)
    tau_out[j] = tau_c*exp(-R_out/h)
    
    ; Write component information
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
  endfor ; j=0, 5 
  
  ; Write total dust mass
  printf, ounit, 'total dust mass'
  printf, ounit, total(M)
  
  ; Close file
  free_lun, ounit
end

pro dust_prof, nposr, nposz
;
; Define the dust distribution and calculate opacities
  
; Grid parameters
 dr = 15.                ; radial step size
 dz = 50.                ; vertical step size
 if n_elements(nposr) eq 0 then nposr = 2500  ; number of radial positions (with default)
 if n_elements(nposz) eq 0 then nposz = 12    ; number of vertical positions (with default)
 sha = 0.02              ; truncation sharpness factor
  
; Component parameters (using structures for cleaner organization)
; Format: [disk_params, thick disk, thin disk, thin ring1, thin ring2, outer disk1, outer disk2]
 comp = {tau: 0.0, $      ; central face-on opacity
          xid: 0.0, $      ; radial profile parameter
          hd: 0.0, $       ; radial scale length (pc)
          hdin: 0.0, $     ; inner radius (pc)
          hdtin: 0.0, $    ; truncation inner radius (pc)
          hdsolar: 0.0, $  ; solar radius parameter (pc)
          zd: 0.0, $       ; vertical scale height (pc)
          zdin: 0.0, $     ; inner vertical parameter (pc)
          zdsolar: 0.0, $  ; solar vertical parameter (pc)
          rtrun: 0.0 }     ; truncation radius (pc)
  
; Define 6 dust components (depends on galaxy)
 comps = replicate(comp, 6)
  
; Component 0: Thick disk
 comps[0].tau = 2.30
 comps[0].xid = 0.4
 comps[0].hd = 8e3
 comps[0].hdin = 2.5e3
 comps[0].hdtin = 0.0
 comps[0].hdsolar = 20.e3
 comps[0].zd = 0.27e3
 comps[0].zdin = 0.16e3
 comps[0].zdsolar = 0.27e3
 comps[0].rtrun = 30.e3
  
; Component 1: Thin disk
 comps[1].tau = 0.12
 comps[1].xid = -0.7
 comps[1].hd = 6.2e3
 comps[1].hdin = 2.5e3
 comps[1].hdtin = 1.029e3
 comps[1].hdsolar = 3.2e3
 comps[1].zd = 0.09e3
 comps[1].zdin = 0.09e3
 comps[1].zdsolar = 8.97e3
 comps[1].rtrun = 30.e3
 
; Component 2: Thin ring 1
 comps[2].tau = 1.91
 comps[2].xid = 0.000001
 comps[2].hd = 0.9e3
 comps[2].hdin = 0.000001
 comps[2].hdtin = 0.0
 comps[2].hdsolar = 8.972e3
 comps[2].zd = 0.27e3
 comps[2].zdin = 0.27e3
 comps[2].zdsolar = 0.27e3
 comps[2].rtrun = 5.e3
  
; Component 3: Thin ring 2
 comps[3].tau = 0.21
 comps[3].xid = 0.000001
 comps[3].hd = 0.4e3
 comps[3].hdin = 0.000001
 comps[3].hdtin = 0.000001
 comps[3].hdsolar = 10.5e3
 comps[3].zd = 0.09e3
 comps[3].zdin = 0.09e3
 comps[3].zdsolar = 0.09e3
 comps[3].rtrun = 5.e3
  
; Component 4: Outer disk 1
 comps[4].tau = 1e-10
 comps[4].xid = -20.5
 comps[4].hd = 7.9e3
 comps[4].hdin = 16.0e3
 comps[4].hdtin = 15.256e3
 comps[4].hdsolar = 10.2e3
 comps[4].zd = 0.27e3
 comps[4].zdin = 0.27e3
 comps[4].zdsolar = 0.27e3
 comps[4].rtrun = 30.0e3
  
; Component 5: Outer disk 2
 comps[5].tau = 1e-10
 comps[5].xid = -20.5
 comps[5].hd = 9.6e3
 comps[5].hdin = 16.0e3
 comps[5].hdtin = 15.256e3
 comps[5].hdsolar = 10.5e3
 comps[5].zd = 0.09e3
 comps[5].zdin = 0.09e3
 comps[5].zdsolar = 0.09e3
 comps[5].rtrun = 30.0e3
 
; Calculate central space density parameters
 a = dblarr(6)
 for i=0, 5 do a[i] = comps[i].tau/(2.0d0*comps[i].hd)
  
  ; Parameters for vertical flaring (only needed for comp 1)
  zd01 = comps[0].zd
  flared02 = alog10(comps[1].hdsolar/comps[1].hdin)
  flared01 = alog10((comps[0].zdsolar-comps[1].zd)/(comps[1].zdin-comps[1].zd))
  flared11 = flared01/flared02
  
  ; Calculate truncation parameters
  sigtruncate = fltarr(6)
  for i=0, 5 do sigtruncate[i] = sha*comps[i].rtrun
  
  ; Output filenames
  filenames = ['opacity.dat', 'opacity0.dat', 'opacity1.dat', 'opacity2.dat', $
               'opacity3.dat', 'opacity4.dat', 'opacity5.dat']
  
  ; Open output files
  units = lonarr(7) 
  for i=0, 6 do begin
    openw, lun, filenames[i], /get_lun
    units[i] = lun
  endfor; i=0, 6 
    
  ; Calculate dust density at each position
  rho = 0.
  for i=0, nposr-1 do begin
    z = 0.
    for n=0, nposz-1 do begin
      ; Initialise dust densities for each component
      dust_dens = fltarr(6)
      
      ; Calculate dust density for each component
      for c=0, 5 do begin
        ; Calculate truncation factor
        trunc = 0.5*(1.0-erf((rho-comps[c].rtrun)/sigtruncate[c]))
        
        ; Component-specific z-scale height
        zd1 = comps[c].zd
        
        ; Calculate dust density based on radial region
        if (rho lt comps[c].hdtin) then begin
          dust_dens[c] = 0
        endif else if (rho ge comps[c].hdtin and rho lt comps[c].hdin) then begin
          dust_dens[c] = a[c]*(comps[c].zd/zd1)* $
              ((rho/comps[c].hdin)*(1.-comps[c].xid)+comps[c].xid)* $
              exp(-comps[c].hdin/comps[c].hd-abs(z)/zd1)
        endif else begin
          dust_dens[c] = a[c]*(comps[c].zd/zd1)* $
              exp(-rho/comps[c].hd-abs(z)/zd1)
        endelse
        
        ; Apply truncation and enforce non-negative values
        dust_dens[c] = dust_dens[c]*trunc > 0
      endfor ; c=0, 5 
      
      ; Total dust density
      total_dens = total(dust_dens)
      
      ; Write to files
      printf, units[0], rho, z, total_dens
      for c=0, 5 do printf, units[c+1], rho, z, dust_dens[c]
      
      ; Move to next z position
      z += dz
    endfor ; n=0, nposz-1 
    
    ; Move to next radial position
    rho += dr
  endfor ; i=0, nposr-1
  
  ; Close output files
  for i=0, 6 do free_lun, units[i]
  
  ; Calculate global parameters and dust masses
  write_global_params, comps
end
