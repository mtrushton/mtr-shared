pro fname,model,qyear,tau,sfr,old,bd,scaabs,slam,namep
; Creates file names from the best fitting parameters
rootdir='/nfs/d58/vlk/sedmodel/'
dir4='urad_ngc0628_sca/indata/'

;define some strings for file names
suv=strcompress(string(fix(sfr*100)),/remove_all)
sbd=strcompress(string(fix(bd*100)),/remove_all)
sold=strcompress(string(fix(old*100)),/remove_all)
stau=strcompress(string(fix(tau*10)),/remove_all)

;read file with geometry
ss = ''
ifinns1 = 0L
ifinns2 = 0L
ifinnd1 = 0L
ifinnd2 = 0L
idisk1 = 0L
idisk2 = 0L

filename = 'geometry.in'
name = rootdir+dir4+filename
openr,unit,name,/get_lun
print, 'read ', name

; Read all parameters 
readf, unit, ss & readf, unit, tau1
readf, unit, ss & readf, unit, tau2
readf, unit, ss & readf, unit, hd
readf, unit, ss & readf, unit, zd
readf, unit, ss & readf, unit, hdin
readf, unit, ss & readf, unit, hdtin
readf, unit, ss & readf, unit, zdin
readf, unit, ss & readf, unit, hdsolar
readf, unit, ss & readf, unit, zdsolar
readf, unit, ss & readf, unit, hd1
readf, unit, ss & readf, unit, zd1
readf, unit, ss & readf, unit, hd1in
readf, unit, ss & readf, unit, hd1tin
readf, unit, ss & readf, unit, zd1in
readf, unit, ss & readf, unit, h_d1solar
readf, unit, ss & readf, unit, z_d1solar
readf, unit, ss & readf, unit, h_uSDSSdisk
readf, unit, ss & readf, unit, h_udisk
readf, unit, ss & readf, unit, h_bdisk
readf, unit, ss & readf, unit, h_gdisk
readf, unit, ss & readf, unit, h_vdisk
readf, unit, ss & readf, unit, h_rSDSSdisk
readf, unit, ss & readf, unit, h_iSDSSdisk
readf, unit, ss & readf, unit, h_idisk
readf, unit, ss & readf, unit, h_zdisk
readf, unit, ss & readf, unit, h_jdisk
readf, unit, ss & readf, unit, h_hdisk
readf, unit, ss & readf, unit, h_kdisk
readf, unit, ss & readf, unit, zs
readf, unit, ss & readf, unit, hsin
readf, unit, ss & readf, unit, hstin
readf, unit, ss & readf, unit, zsin
readf, unit, ss & readf, unit, hssolar
readf, unit, ss & readf, unit, zssolar
readf, unit, ss & readf, unit, hs1
readf, unit, ss & readf, unit, zs1
readf, unit, ss & readf, unit, hs1in
readf, unit, ss & readf, unit, hs1tin
readf, unit, ss & readf, unit, zs1in
readf, unit, ss & readf, unit, hs1solar
readf, unit, ss & readf, unit, zs1solar
readf, unit, ss & readf, unit, rtrun
readf, unit, ss & readf, unit, sharp
readf, unit, ss & readf, unit, rtrun1
readf, unit, ss & readf, unit, sharp1
readf, unit, ss & readf, unit, reff_uSDSSbulge
readf, unit, ss & readf, unit, reff_ubulge
readf, unit, ss & readf, unit, reff_bbulge
readf, unit, ss & readf, unit, reff_gbulge
readf, unit, ss & readf, unit, reff_vbulge
readf, unit, ss & readf, unit, reff_rSDSSbulge
readf, unit, ss & readf, unit, reff_iSDSSbulge
readf, unit, ss & readf, unit, reff_ibulge
readf, unit, ss & readf, unit, reff_zbulge
readf, unit, ss & readf, unit, reff_jbulge
readf, unit, ss & readf, unit, reff_hbulge
readf, unit, ss & readf, unit, reff_kbulge
readf, unit, ss & readf, unit, ellipt
readf, unit, ss & readf, unit, nsersic
readf, unit, ss & readf, unit, ifinns1
readf, unit, ss & readf, unit, ifinns2
readf, unit, ss & readf, unit, ifinnd1
readf, unit, ss & readf, unit, ifinnd2
readf, unit, ss & readf, unit, idisk1
readf, unit, ss & readf, unit, idisk2

free_lun, unit

; Set derived parameters
hs = h_bdisk
reff = reff_bbulge

; Format geometry parameters for filename
shd=strcompress(string(fix(hd*1000.)),/remove_all)
szd=strcompress(string(fix(zd*1000.)),/remove_all)
shd1=strcompress(string(fix(hd1*1000.)),/remove_all)
szd1=strcompress(string(fix(zd1*1000.)),/remove_all)
shs=strcompress(string(fix(hs*1000.)),/remove_all)
szs=strcompress(string(fix(zs*1000.)),/remove_all)
shs1=strcompress(string(fix(hs1*1000.)),/remove_all)
szs1=strcompress(string(fix(zs1*1000.)),/remove_all)
sreff=strcompress(string(fix(reff*1000.)),/remove_all)
sellipt=strcompress(string(fix(ellipt*100.)),/remove_all)

; Create filename
namep=model+'_q'+qyear+$
      '_t'+stau+'_s'+suv+'_no'+sold+'_bd'+sbd+$
      '_hd'+shd+'_zd'+szd+'_hd1_'+shd1+'_zd1_'+szd1+$
      '_hs'+shs+'_zs'+szs+'_hs1_'+shs1+'_zs1_'+szs1+$
      '_reff'+sreff+'_ell'+sellipt+'_'+scaabs+'_l'+slam+'um'

return
end
