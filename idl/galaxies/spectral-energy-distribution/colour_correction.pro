
pro colour_correction
;
;  Calculates colour corrections for IRAC, MIPS, and WISE filters for user supplied
;  model in NAME_MODEL_FILE, or for testing against three function types:
;  blackbody (see input TEFF), modified blackbody (see input TEFF, BETA_INDEX), 
;  and power law (see input BETA_INDEX).
;----------------------------------------------------------------------------
;  INPUTS
;
;  NAME_MODEL_FILE: String name of file for model fluxes. The file should
;       be in ascii format, containing two columns of data: wavelength (um) 
;       and flux (frequency units).
;  INSTRUMENT: String with value 'IRAC', 'MIPS', or 'WISE' to select instrument.
;  TEFF: Float value of effective temperature in the blackbody function (BB_nu) 
;        used for model spectra, if desired.
;  BETA_INDEX: Index for power law or modified blackbody model spectra.
;  BBTEST: Integer of value 1 or 0. If set to '1', a blackbody of temperature 
;        TEFF is used as an input model. 
;  MODBBTEST: Integer. Same as BBTEST, but for a modified blackbody, with a 
;        beta index of BETA_INDEX.
;  POWERLAWTEST: Integer of value 1 or 0. If set to '1', a power law of index
;        BETA_INDEX is used as an input model. 
;----------------------------------------------------------------------------
;  OUTPUTS
;  
;  LAMBDA_EFF: Float effective wavelength of the selected filter.
;  WIDTH: Float band width of the selected filter.
;  CORR: Float colour correction. To convert the observed fluxes into 
;        monochromatic values, divide by this factor.
; ---------------------------------------------------------------------------
;  Requires ascii files of instrument response curves in the appropriate directory.
;
;  See: 
;  IRAC: https://irsa.ipac.caltech.edu/data/SPITZER/docs/irac/iracinstrumenthandbook/18/
;  MIPS: https://irsa.ipac.caltech.edu/data/SPITZER/docs/mips/mipsinstrumenthandbook/51/
;  WISE: http://svo2.cab.inta-csic.es/svo/theory/fps/index.php?mode=browse&gname=WISE
;       and Wright et al. 2010, AJ, 140, 1868.
;
print,'Unified Colour Correction Program'
print,''

; USER INPUTS ----------------------------------------------------------------
name_model_file = 'IR_model_fluxes.dat'

; Select instrument: 'IRAC', 'MIPS', or 'WISE'
instrument = 'IRAC'  ; Change as needed

Teff = 1000 ; K

beta_index = -2.0

; Set the integer if a function is required for the input model. The function is 
;    used if the corresponding value is set to 1 and unused if the value is 0.
;    CORR values for one function type are calculated only. If none are equal to 
;    1, the values are calculated from the input model file.

BBTest = 0L      ; Input blackbody 1 or 0. If '1', set Teff.
ModBBTest = 0L   ; Input modified blackbody 1 or 0. If '1', set Teff & beta_index.
PowerLawTest = 0L ; Input power law 1 or 0. If '1', set beta_index.
;------------------------------------------------------------------------------

; Constants for functions
c = 2.998d8  ; m/s
h = 6.626d-34 
k = 1.381d-23
hckt = h * c / (k * Teff)
cnst = 2 * h * c^2

; Setup instrument-specific parameters
CASE STRUPCASE(instrument) OF
  'IRAC': BEGIN
    print,'Colour Corrections for IRAC Channels 1-4.'
    lambda_ref = [3.550, 4.493, 5.731, 7.872]
    name_response_file = ['irac3.6_response.dat', 'irac4.5_response.dat', $
                          'irac5.7_response.dat', 'irac8_response.dat']
    unit_conv = 1.0  ; Wavelength in microns
    reference_temp = 0  ; No specific reference temp for IRAC
  END
  
  'MIPS': BEGIN
    print,'Colour Corrections for MIPS 24/70/160 filters.'
    lambda_ref = [23.68, 71.42, 155.9]
    name_response_file = ['mips24_response.dat', 'mips70_response.dat', 'mips160_response.dat']
    unit_conv = 1.0  ; Wavelength in microns
    reference_temp = 10000.0  ; 10,000 K blackbody reference
    hcktreff = h * c / (k * reference_temp)
  END
  
  'WISE': BEGIN
    print,'Colour Corrections for WISE filters 1-4.'
    lambda_ref = [3.3526, 4.6028, 11.5608, 22.0883]
    name_response_file = ['wise34_response.dat', 'wise46_response.dat', $
                         'wise11_response.dat', 'wise23_response.dat']
    unit_conv = 1.0e-4  ; Convert Angstroms to microns
    reference_temp = 14454.0  ; Vega-like reference temperature
    hcktref = h * c / (k * reference_temp)
  END
  
  ELSE: BEGIN
    print,'Unknown instrument. Please specify IRAC, MIPS, or WISE'
    RETURN
  END
ENDCASE

; Validate input parameters
if(BBTest + ModBBTest + PowerLawTest eq 0L) then print,'Input model spectrum' + $
  ' will be read from file '+name_model_file+'.'

if(BBTest + ModBBTest + PowerLawTest gt 1L) then begin 
  print,'Integers set incorrectly for input model functions. Values should be' + $  
  ' 0 or 1, and only one value can be equal to 1.'
  stop
endif

if(BBTest eq 1) then print,'Blackbody function used as input model, with an' + $
  ' effective temperature of'+ strcompress(string(fix(Teff))) +' K.'
if(ModBBTest eq 1) then print,'Modified blackbody function used as input model' + $
  ' spectrum, with an effective temperature of'+strcompress(string(fix(Teff)))+' K' + $
  ' and a beta index of'+strcompress(string(fix(beta_index)))+'.'
if(PowerLawTest eq 1) then print,'Power law used as input model spectrum, with' + $
  ' beta index of'+strcompress(string(fix(beta_index)))+'.'

; Read the model file if needed
ncorr = n_elements(name_response_file)
nflux = numlines(name_model_file)

if(BBTest + ModBBTest + PowerLawTest eq 0L) then begin
  openr, unit1, name_model_file, /get_lun
  
  lambda = dblarr(nflux)
  flux = dblarr(nflux)
  for i = 0, nflux - 1L do begin
    readf, unit1, dum, dum1
    lambda(i) = dum  ; microns
    flux(i) = dum1   ; Jy or frequency units
    flux(i) = flux(i) * c/lambda(i)^2  ; wavelength units
  endfor 
  free_lun, unit1
endif

; Process each filter in the selected instrument
for n = 0, ncorr - 1L do begin
  print,''
  CASE STRUPCASE(instrument) OF
    'IRAC': print,'IRAC Channel' + strcompress(string(fix(n+1)))
    'MIPS': print,'MIPS ' + strcompress(string(fix(lambda_ref[n])))
    'WISE': print,'WISE filter' + strcompress(string(fix(n+1)))
  ENDCASE

  ; Read response function file
  nr = numlines(name_response_file[n])
  openr, unit2, name_response_file[n], /get_lun
  
  wav_res = dblarr(nr)
  response = dblarr(nr)
  for i = 0, nr - 1L do begin
    readf, unit2, dum, dum1
    wav_res(i) = dum * unit_conv  ; Convert to microns if needed
    response(i) = dum1 
  endfor
  
  ; Normalise response function for bandwidth calculation
  response = response/max(response)
  
  ; Calculate response integral and band properties
  R = 0.0d0
  delta_wav = dblarr(nr)
  for i = 0, nr - 1L do begin
    if(i eq 0) then begin
      delta_wav(i) = wav_res(1) - wav_res(0)
    endif else begin
      delta_wav(i) = wav_res(i) - wav_res(i-1)
    endelse
    R = R + response(i) * delta_wav(i)
  endfor
  
  ; Calculate effective wavelength and bandwidth
  lambda_eff_calc = 0.0d0
  delta = 0.0d0
  
  for i = 0, nr - 1L do begin
    lambda_eff_calc = lambda_eff_calc + wav_res(i) * response(i) * delta_wav(i)/R
    delta = delta + response(i) * delta_wav(i)/wav_res(i)
  endfor
  
  print,'effective wavelength (calculated) = ', lambda_eff_calc, ' um'
  
  ; Use reference wavelength from instrument documentation
  lambda_eff = lambda_ref[n]
  print,'effective wavelength (reference) = ', lambda_eff, ' um'
  
  width = lambda_eff * delta
  print,'band width = ', width, ' um'
  
  free_lun, unit2
  
  ; Interpolate model flux if using input file
  if(BBTest + ModBBTest + PowerLawTest eq 0L) then begin
    lambda_max = max(wav_res)
    lambda_min = min(wav_res)
    
    nlower = where(lambda lt lambda_min)
    nhigher = where(lambda gt lambda_max)
    down = max(nlower)
    up = min(nhigher)
    
    lambda_cut = lambda(down:up)
    flux_cut = flux(down:up)
    
    flux_int = interpol(flux_cut, lambda_cut, wav_res)
  endif
  
  ; Calculate test functions
  bb = dblarr(nr)
  mod_bb = dblarr(nr)
  pl = dblarr(nr)
  
  ; Initialise reference function arrays based on instrument
  CASE STRUPCASE(instrument) OF
    'MIPS': BEGIN
      bb_num = dblarr(nr)
      bb_dem = dblarr(nr)
    END
    'WISE': BEGIN
      bb_ref = dblarr(nr)
      vega = dblarr(nr)
    END
  ENDCASE
  
  ; Calculate the various test functions
  for i = 0, nr - 1L do begin
    ; Blackbody
    bb(i) = cnst/(wav_res(i) * 1d-6)^5 * (1 / (exp(hckt/(wav_res(i) * 1d-6) - 1)))
    
    ; Modified blackbody
    mod_bb(i) = bb(i)/(wav_res(i) * 1d-6)^beta_index
    
    ; Power law - format differs between instruments
    CASE STRUPCASE(instrument) OF
      'IRAC': pl(i) = wav_res(i)^(beta_index + 2L)
      'MIPS': pl(i) = wav_res(i)^(beta_index + 2L)
      'WISE': pl(i) = 1L /(wav_res(i)^(beta_index + 2L))
    ENDCASE
    
    ; Calculate reference functions for specific instruments
    CASE STRUPCASE(instrument) OF
      'MIPS': BEGIN
        bb_num(i) = cnst/(wav_res(i) * 1d-6)^5 * (1 / (exp(hcktreff/(wav_res(i) * 1d-6) - 1)))
        bb_dem(i) = cnst/(lambda_eff * 1d-6)^5 * (1 / (exp(hcktreff/(lambda_eff * 1d-6) - 1)))
      END
      'WISE': BEGIN
        bb_ref(i) = cnst/(wav_res(i) * 1d-6)^5 * (1 / (exp(hcktref/(wav_res(i) * 1d-6) - 1)))
        vega(i) = 1.0158d-16 * (1 - 0.0083 * alog(wav_res(i)/8.891)^2) * bb_ref(i)
      END
    ENDCASE
  endfor
  
  ; Select model function based on test flags
  if(BBTest eq 1L) then flux_int = bb
  if(ModBBTest eq 1L) then flux_int = mod_bb
  if(PowerLawTest eq 1L) then flux_int = pl
  
  ; Find the reference wavelength index
  freff = abs(wav_res - lambda_eff)
  index = where(freff eq min(freff))
  F0 = flux_int(index)
  
  ; Calculate colour correction based on instrument method
  col_num = 0.0
  col_dem = 0.0
  
  for i = 0, nr - 1L do begin
    CASE STRUPCASE(instrument) OF
      'IRAC': BEGIN
        col_num = col_num + flux_int(i)/F0 * (lambda_eff/wav_res(i))^(-1) * response(i) * delta_wav(i)
        col_dem = col_dem + (lambda_eff/wav_res(i))^(-2) * response(i) * delta_wav(i)
      END
      'MIPS': BEGIN
        col_num = col_num + flux_int(i)/F0 * response(i) * wav_res(i) * delta_wav(i)
        col_dem = col_dem + bb_num(i)/bb_dem(i) * wav_res(i) * response(i) * delta_wav(i)
      END
      'WISE': BEGIN
        col_num = col_num + flux_int(i)/F0 * response(i) * delta_wav(i)
        col_dem = col_dem + (lambda_eff/wav_res(i))^4 * response(i) * delta_wav(i)
      END
    ENDCASE
  endfor
  
  corr = col_num/col_dem
  print,'colour correction = ', corr
  
endfor

end
