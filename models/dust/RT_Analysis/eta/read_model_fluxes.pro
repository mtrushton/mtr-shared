pro read_model_fluxes, ndat, freq, flux_model_bulge, flux_model_old_stellar_disk, flux_model_young_stellar_disk
;
; Read the files containing the intrinsic stellar fluxes at each frequency, used for scaling the models to the
; galaxy
;
; Set number of data points
 ndat = 16
  
; Initialise arrays
 freq = dblarr(ndat)
 flux_model_bulge = dblarr(ndat)
 flux_model_old_stellar_disk = dblarr(ndat)
 flux_model_young_stellar_disk = dblarr(ndat)
  
; Define input directory and files
 dir = '../intrinsic_luminosities/'
 files = [dir + 'old_stellar_bulge.dat', $
     dir + 'old_stellar_disk.dat', $
     dir + 'young_stellar_disk_clumpy.dat']
           
; Read bulge data (includes frequency values)
 data = dblarr(2, ndat)
 openr, unit, files[0], /get_lun
 readf, unit, data
 free_lun, unit
  
; Store frequency values (only need to read once)
 freq = data[0,*]
 flux_model_bulge = data[1,*]
  
; Read old stellar disk data
 openr, unit, files[1], /get_lun
 readf, unit, data
 free_lun, unit
 flux_model_old_stellar_disk = data[1,*]
  
; Read young stellar disk data
 openr, unit, files[2], /get_lun
 readf, unit, data
 free_lun, unit
 flux_model_young_stellar_disk = data[1,*]
end
