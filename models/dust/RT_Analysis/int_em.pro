pro int_em,npos,nposr,nposz,ncomp,file,total_energy
; Integrates component emissivities over galaxy.
; Values stored in TOTAL_ENERGY vector.
;---------------------------------------------------
; files with emissivities in arbitrary units
;file = ['emissivity_bulge.dat', 'emissivity_d.dat', $
;         'emissivity_td.dat', 'emissivity_tdi.dat']
file = ['emissivity_bulge.dat', 'emissivity_d.dat', $
         'emissivity_td.dat']
; define no. of positions in r and z in file
nposr = 2100
nposz = 12
npos = nposr * nposz

ncomp = n_elements(file)
total_energy = dblarr(ncomp) ; integrated emissivities

rho = dblarr(npos)
z = dblarr(npos)
eta = dblarr(npos)
;
; read emissvities for each component
;
for i = 0, ncomp - 1 do begin
 openr, unit, file(i), /get_lun

for n = 0, npos - 1L do begin
 readf,unit,dum,dum1,dum2
 rho(n) = dum
 z(n)= dum1
 eta(n) = dum2
endfor ; n = 0, npos -1L
;
; integrate component emisssivities over galaxy
;
etadouble = dblarr(nposr,nposz)
rrr = dblarr(nposr,nposz)
zzz = dblarr(nposr,nposz)
k = 0
for kr = 0L, nposr - 1 do begin
    for kz = 0L, nposz - 1 do begin
        etadouble(kr, kz) = eta(k)
        rrr(kr, kz) = rho(k)
        zzz(kr, kz) = z(k)
        k = k + 1
    endfor
endfor

pass = 0.
pass_energy = 0.
for kz = 0L,nposz-2 do begin
  for kr = 0L,nposr-2 do begin
    pass = !pi * (rrr(kr+1,kz)^2-rrr(kr,kz)^2) * (zzz(kr,kz+1)-zzz(kr,kz))$
           * (etadouble(kr,kz) + etadouble(kr+1,kz+1))/2.
    pass_energy = pass_energy + pass  ;in W/Hz
  endfor
endfor
total_energy(i) = pass_energy * 2

free_lun, unit

endfor ; i = 0, ncomp - 1 do begin

end
