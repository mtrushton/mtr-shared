pro scale_int_freq, nposr, nwav, freq, flux, wav_ref, r, total_model, pass_energy_z, energy_rho
;
; Scale model by flux/total_model ratio
;  
 scale = flux / total_model
  
; Create scaled model for all positions and wavelengths using array operations
 scaled_model = rebin(pass_energy_z, nposr, nwav) * rebin(transpose(scale), nposr, nwav)
  
; Total energy (integration over radial dimension)
 r_diffs = r[1:nposr-1]^2 - r[0:nposr-2]^2
 avg_energy_z = (pass_energy_z[0:nposr-2] + pass_energy_z[1:nposr-1]) / 2.0
 pass_energy = !pi * total(r_diffs * avg_energy_z)
  
 wav = 2.998d8 / freq
  
; Select frequencies below reference wavelength
 ind = where(wav lt wav_ref * 1d-6)
 freq_samp = freq[ind]
 model = scaled_model[*, ind]
  
; Energy density for each radius
 energy_rho = dblarr(nposr)
  
; Integrate over frequency for each radial position
 for n = 0, nposr - 1 do begin
    energy_rho[n] = int_tabulated(freq_samp, reform(model[n,*]))
 endfor ; n = 0, nposr - 1
end
