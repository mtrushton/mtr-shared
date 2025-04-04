pro dust_mass
; Main program to calculate and plot the radial profile of the 
; surface density of the dust mass distribution for a galaxy
; modelled with the RT analysis
  
; Define dust density distribution
 dust_prof, nposr, nposz
  
; Integrate over z-axis and calculate scaling factors
 int_z, npos, nposr, nposz, ncomp, files, r, dust_dens_z, scale
  
; Plot dust mass distribution
 plot_dust_mass, nposr, ncomp, files, nposr, nposz, dust_dens_z, scale, r
end
