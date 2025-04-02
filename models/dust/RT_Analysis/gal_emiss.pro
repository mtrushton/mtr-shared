pro gal_emiss
;
; Creates radial profiles of the stellar emissivity at z = 0
; at a reference wavelength, usually at about 4400 Angstrom, for 
; a galaxy modelled with the RT analysis.
;
; INT_EM takes the radial emissivity profiles for each component
; and integrates them over wavelength. This is then used by 
; SCALE_FACTOR for scaling to the physical units for the galaxy.

; PLOT_EMISS_PROF plots the scaled radial profiles of the stellar
; emissivities for each component. A zoomed plot is also included,
; showing detail in the centre (this is defined as 0-2 kpc in the
; program).  
;
; INPUTS can be defined in the following routines as follows:
; INT_EM: file names for unscaled emissivity radial profiles, along
;     with the grid parameters
; SCALE_FACTOR: IDL save file containing fluxes for each component.
;    These fluxes should be given for each component at a specfic 
;     waveband (B is a good choice). 
; PLOT_EMISS_PROF: file name for the output 
;

int_em,npos,nposr,nposz,ncomp,file,total_energy
scale_factor,ncomp,total_energy,scale
plot_emiss_prof,npos,nposr,nposz,ncomp,file,scale

end
