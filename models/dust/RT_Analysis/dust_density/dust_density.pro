pro dust_density
;
; Plots a radial profile of the dust density along the galactic
; plane for a PT11 RT modelled galaxy.
;
; Three routines called as follows:
; DUST_PROF : defines the dust distribution and derives opacities
; INT : Scales the dust opacity distribution from the derived dust masses
; PLOT_DUST_DENSITY : Creates the plot
;
; OUTPUTS:
; DUST_DENSITY_R.PS : Plot of the dust density radial distribution at z=0
; DUST_DENSITY.SAV : IDL save file containing dust density at z=0 for each 
;    radial position in the galaxy

  dust_prof, nposr, nposz
  int, nposr, nposz, npos, scale, ncomp, file
  plot_dust_density, npos, ncomp, file, nposr, nposz, scale
end
