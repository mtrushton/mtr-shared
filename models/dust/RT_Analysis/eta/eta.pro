pro eta
;
; Creates radial profiles of the surface densities of SFR and stellar mass for an RT 
; modelled galaxy. 
; 
; Stellar emissivities of the bulge, disk, and young stellar disk(s) are calculated from
; the optimised morphological parameters (EMISSIVITY_*) , which are then scaled to the 
; derived intrinsic model fluxes of the galaxy (supplied as ascii files for each 
; morphological component).
;
; The radial plots are surface densities (i.e., integrated over z) in INT_Z.
;
; SFR Plot (MAKE_FIG_UV):
; The scaled young stellar disk emissivity is integrated from the UV to the B band for the 
; creation of an sfr radial plot from 'L_young_unit' 
; see Eqn 16 in Popescu et al., 2011, A&A, 527, 109 
; Data saved in 'sfr.xdr'
;
; Stellar mass plot (MAKE_FIG_STELLAR_MASS):
; The intrinsic L band model luminosity is used to create a radial profile of the stellar 
; mass distribution. Alternatively, the stellar mass can be used as an input for 
; normalising the plot. 
; Data saved in 'stellar_mass.xdr'
;
; The main inputs are handled in ETA_INPUT. Note that the model fluxes used in the scaling
; should be provided in separate files which are read by READ_MODEL_FLUXES 
;
; Add Coyote library to path
 !PATH = Expand_Path('+/nfs/d58/vlk/sedmodel/MTR/N0628/gn/PROFILE_RV/sca/coyote/') + ':' + !PATH 
  
; Read input parameters for all components
 eta_input, hstin, hs, hsin, hssolar, zs, zsin, zssolar, xis0, zs01, sha, rtrun, $
            hs1tin, hs1, hs1in, hs1solar, zs1, zs1solar, zs1in, xis1, zs11, sha1, rtrun1, $
            hs2tin, hs2, hs2in, hs2solar, zs2, zs2solar, zs2in, xis2, zs21, sha2, rtrun2, $
            reff, ellipt, nsersic, pixsize, refeta, shab, rtrunb, bsersic, $
            dust_mass, tau0, tau1, xid0, xid1, hd, hdin, hdsolar, zd, zdin, zdsolar, $
            hd1, hd1in, hd1solar, zd1, zd1in, zd1solar, rtrund, rtrund1, shad, hd1tin, $
            tau2, hd2, zd2, hd2in, hd2tin, rtrund2, xid2, $
            nposr, nposz, dr, dz, Mval, scale_mass_prof
  
; Calculate emissivities of galaxy components
 emissivity_bulge, nposr, nposz, dr, dz, reff, ellipt, nsersic, pixsize, refeta, sha, rtrun, bsersic
 emissivity_disk, nposr, nposz, dr, dz, hstin, hs, hsin, hssolar, zs, zsin, zssolar, xis0, zs01, sha, rtrun
 emissivity_tdisk, nposr, nposz, dr, dz, hs1tin, hs1, hs1in, hs1solar, zs1, zs1in, zs1solar, xis1, zs11, sha1, rtrun1
  
; Integrate models
 int_model, nposr, nposz, total_model_tdisk, total_model_tdisk2, total_model_disk, total_model_bulge
  
; Integrate emissivity along z axis
 int_z, nposr, nposz, 'emissivity_tdisk.dat', r, eta_tdisk_intz
  
; Read model fluxes for different components
 read_model_fluxes, ndat, freq, flux_model_bulge, flux_model_old_stellar_disk, flux_model_young_stellar_disk
  
; Integrate the young stellar components between the UV and the B band for SFR plot 
 scale_int_freq, nposr, ndat, freq, flux_model_young_stellar_disk, 0.45, r, total_model_tdisk, $
     eta_tdisk_intz, scaled_model_tdisk
  
; Create SFR plot
 make_fig_uv, r, scaled_model_tdisk, scaled_model_tdisk2, flux_tot_uv, 'uv.ps', 'uv_sfr_units.ps'
  
; Integrate other components along z axis
 int_z, nposr, nposz, 'emissivity_disk.dat', r, eta_disk_intz
 int_z, nposr, nposz, 'emissivity_bulge.dat', r, eta_bulge_intz
  
; Scale each component to the model fluxes at 3.5 microns for use of the L band-stellar mass relation
 scale_freq, nposr, ndat, freq, flux_model_old_stellar_disk, 3.5, r, total_model_disk, $
     eta_disk_intz, scaled_model_disk
 scale_freq, nposr, ndat, freq, flux_model_bulge, 3.5, r, total_model_bulge, $
     eta_bulge_intz, scaled_model_bulge
 scale_freq, nposr, ndat, freq, flux_model_young_stellar_disk, 3.5, r, total_model_tdisk, $
     eta_tdisk_intz, scaled_model_tdisk
  
; Create stellar mass plot
 make_fig_stellar_mass, r, scaled_model_bulge, scaled_model_disk, scaled_model_tdisk, $
     scaled_model_tdisk2, 'stellar_mass.ps', stellar_mass, scale_factor_mass, Mval, scale_mass_prof
  
end
