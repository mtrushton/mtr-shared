      program specfit_dusty_grid
      implicit none
c-----------------------------------------------------------------------
c     Program to fit IR spectra from a grid of Dusty models
c-----------------------------------------------------------------------
      integer ndata, nmodel, ntemp, ntau, nlines
      parameter (ndata = 357)
      parameter (nmodel = 420)
      parameter (ntemp = 21)
      parameter (ntau = 5)
      parameter (nlines = 15)
      
c     Arrays for observed data
      real lambda(ndata), flux(ndata), error1(ndata), error3(ndata)
      real sed_flux(ndata)
      
c     Arrays for model data
      real lambda_model(nmodel), flux_model(nmodel)
      real lambda_model_trim(ndata), flux_model_trim(ndata)
      real flux_scaled_model(ndata), t(ndata)
      
c     Arrays for line data and calculations
      real line_ls(nlines), line_ll(nlines)
      logical exclude_line(ndata)
      
c     Model grid parameters
      character*3 temp(ntemp), best_temp
      character*4 tau(ntau), best_tau
      
c     File names and paths
      character*17 model_dir
      character*6 model_prefix
      character*2 model_infix
      character*5 model_suffix
      character*80 filename
      
c     Variables for calculations
      integer i, j, k, jk, il, istat, istart, xi, xii, is
      real max_lambda, min_lambda
      real S, Sx, Sy, Stt, sum1, a, b, chisq, best_chisq
      real scaled_flux
      
c     File units
      integer unit_input, unit_lines, unit_sed, unit_model, unit_output
      parameter (unit_input = 10)
      parameter (unit_lines = 11)
      parameter (unit_sed = 12)
      parameter (unit_model = 20)
      parameter (unit_output = 30)
      
c-----------------------------------------------------------------------
c     Initialize model grid parameters
c-----------------------------------------------------------------------
      temp(1) = '400'
      temp(2) = '410'
      temp(3) = '420'
      temp(4) = '430'
      temp(5) = '440'
      temp(6) = '450'
      temp(7) = '460'
      temp(8) = '470'
      temp(9) = '480'
      temp(10) = '490'
      temp(11) = '500'
      temp(12) = '510'
      temp(13) = '520'
      temp(14) = '530'
      temp(15) = '540'
      temp(16) = '550'
      temp(17) = '560'
      temp(18) = '570'
      temp(19) = '580'
      temp(20) = '590'
      temp(21) = '600'
      
      tau(1) = '0018'
      tau(2) = '0020'
      tau(3) = '0022'
      tau(4) = '0024'
      tau(5) = '0025'
      
      model_dir = 'models/rsoph0709/'
      model_prefix = 'rsoph_'
      model_infix = 'k_'
      model_suffix = '.s001'
      
c-----------------------------------------------------------------------
c     Read observation data files
c-----------------------------------------------------------------------
      open(unit_input, file='cassis_tbl_spcf_17681152_interp.dat',
     &     status='old', err=901)
      
      open(unit_lines, file='lines.dat', status='old', err=902)
      
      open(unit_sed, file='sed_sub_interp.dat', status='old', err=903)
      
c     Read spectral data, error values, and SED data
      do 10 i = 1, ndata
          read(unit_input, 100, err=904) lambda(i), flux(i), 
     &         error1(i), error3(i)
          read(unit_sed, 150, err=905) lambda(i), sed_flux(i)
          
c         Subtract SED from flux
          flux(i) = flux(i) - 1.0 * sed_flux(i)
          
c         Initialize exclusion array
          exclude_line(i) = .false.
10    continue

c     Read line data
      do 20 il = 1, nlines
          read(unit_lines, 120, err=906) line_ls(il), line_ll(il)
20    continue

c     Close input files
      close(unit_input)
      close(unit_lines)
      close(unit_sed)
      
c     Mark points that fall within spectral lines
      do 30 i = 1, ndata
          do 35 is = 1, nlines
              if (lambda(i) .gt. line_ls(is) .and. 
     &            lambda(i) .lt. line_ll(is)) then
                  exclude_line(i) = .true.
                  goto 39
              endif
35        continue
39        continue
30    continue
      
c-----------------------------------------------------------------------
c     Open output file for results
c-----------------------------------------------------------------------
      open(unit_output, file='specfit_results.dat', 
     &     status='unknown', err=907)
      
      write(unit_output, 500) 'Temperature  Tau       Chi-square'
      
c-----------------------------------------------------------------------
c     Initialize best fit tracking
c-----------------------------------------------------------------------
      best_chisq = 1.0e32
      best_temp = '   '
      best_tau = '    '
      
c-----------------------------------------------------------------------
c     Main loop through temperature and optical depth grid
c-----------------------------------------------------------------------
      do 80 k = 1, ntemp
          do 90 jk = 1, ntau
c             Create model filename
              filename = model_dir // model_prefix // temp(k) //
     &                  model_infix // tau(jk) // model_suffix
              
c             Read model file
              open(unit_model, file=filename, status='old', err=908)
              
c             Skip header lines
              read(unit_model, *)
              read(unit_model, *)
              read(unit_model, *)
              read(unit_model, *)
              
c             Read model data
              do 15 xi = 1, nmodel
                  read(unit_model, 200, err=909) lambda_model(xi), 
     &                 flux_model(xi)
15            continue
              
c             Close model file
              close(unit_model)
              
c             Find wavelength range of data
              max_lambda = lambda(1)
              min_lambda = lambda(1)
              do 25 i = 2, ndata
                  if (lambda(i) .gt. max_lambda) then
                      max_lambda = lambda(i)
                  endif
                  if (lambda(i) .lt. min_lambda) then
                      min_lambda = lambda(i)
                  endif
25            continue
              
c             Add margin to wavelength range
              max_lambda = max_lambda + 0.1
              min_lambda = min_lambda - 0.1
              
c             Find starting index in model that covers data range
              istat = 0
              istart = 0
              do 45 xi = 1, nmodel
                  if (lambda_model(xi) .ge. min_lambda .and. 
     &                lambda_model(xi) .le. max_lambda) then
                      if (istat .eq. 0) then
                          istart = xi
                      endif
                      istat = istat + 1
                  endif
45            continue
              
c             Check if model covers enough of the data range
              if (istat .lt. ndata) then
                  write(*, *) 'Warning: Model ', filename
                  write(*, *) 'does not cover entire data range'
              endif
              
c             Create trimmed model arrays
              do 55 xii = 1, ndata
                  if (xii + istart - 1 .le. nmodel) then
                      lambda_model_trim(xii) = lambda_model(xii+istart-1)
                      flux_model_trim(xii) = flux_model(xii+istart-1) *
     &                    lambda_model_trim(xii) / 2.988e8
                  else
                      lambda_model_trim(xii) = 0.0
                      flux_model_trim(xii) = 0.0
                  endif
55            continue
              
c-----------------------------------------------------------------------
c             Calculate scaling parameters excluding spectral lines
c-----------------------------------------------------------------------
              S = 0.0
              Sx = 0.0
              Sy = 0.0
              
              do 60 i = 1, ndata
                  if (.not. exclude_line(i) .and. error1(i) .gt. 0.0) then
                      S = S + (1.0 / error1(i)**2)
                      Sx = Sx + (flux_model_trim(i) / error1(i)**2)
                      Sy = Sy + (flux(i) / error1(i)**2)
                  endif
60            continue
              
              Stt = 0.0
              sum1 = 0.0
              
              do 65 i = 1, ndata
                  if (.not. exclude_line(i) .and. error1(i) .gt. 0.0) then
                      t(i) = (1.0/error1(i)) * (flux_model_trim(i)-(Sx/S))
                      Stt = Stt + t(i)**2
                      sum1 = sum1 + ((t(i) * flux(i)) / error1(i))
                  else
                      t(i) = 0.0
                  endif
65            continue
              
c             Calculate linear regression coefficients
              if (Stt .gt. 0.0) then
                  b = (1.0 / Stt) * sum1
              else
                  b = 0.0
              endif
              
              if (S .gt. 0.0) then
                  a = (Sy - (Sx * b)) / S
              else
                  a = 0.0
              endif
              
c-----------------------------------------------------------------------
c             Calculate chi-square value excluding spectral lines
c-----------------------------------------------------------------------
              chisq = 0.0
              
              do 70 j = 1, ndata
                  if (.not. exclude_line(j) .and. error1(j) .gt. 0.0) then
                      scaled_flux = a + b * flux_model_trim(j)
                      chisq = chisq + ((flux(j) - scaled_flux)**2 / 
     &                        error1(j)**2)
                  endif
70            continue
              
c             Write results to output file
              write(unit_output, 550) temp(k), tau(jk), chisq
              
c             Update best fit if better
              if (chisq .lt. best_chisq) then
                  best_chisq = chisq
                  best_temp = temp(k)
                  best_tau = tau(jk)
              endif
              
90        continue
80    continue
      
c     Write best fit results
      write(unit_output, *)
      write(unit_output, 500) 'Best fit:'
      write(unit_output, 600) 'Temperature = ', best_temp, 
     &      ' K, Tau = ', best_tau, ', Chi-square = ', best_chisq
      
c     Close output file
      close(unit_output)
      
      write(*, *) 'Analysis complete. Results written to specfit_results.dat'
      stop
      
c-----------------------------------------------------------------------
c     Format statements
c-----------------------------------------------------------------------
100   format(f16.7,f16.8,2f16.10)
120   format(f8.5, f11.5)
150   format(2f16.7)
200   format(2f11.7)
500   format(A)
550   format(A3, 2X, A4, 2X, E12.5)
600   format(A14, A3, A8, A4, A15, E12.5)
      
c-----------------------------------------------------------------------
c     Error handlers
c-----------------------------------------------------------------------
901   write(*, *) 'Error opening input spectrum file'
      stop
902   write(*, *) 'Error opening lines file'
      stop
903   write(*, *) 'Error opening SED file'
      stop
904   write(*, *) 'Error reading from spectrum file'
      stop
905   write(*, *) 'Error reading from SED file'
      stop
906   write(*, *) 'Error reading from lines file'
      stop
907   write(*, *) 'Error opening output file'
      stop
908   write(*, *) 'Error opening model file: ', filename
      stop
909   write(*, *) 'Error reading from model file: ', filename
      stop
      
      end
