PROGRAM MAKE_DUSTY_GRID
C-------------------------------------------------------------
C  Calculates a grid of Dusty models for range of dust temperatures
C  and dust opacities.
C  Dusty Manual:
C     https://arxiv.org/pdf/astro-ph/9910475.pdf
C-------------------------------------------------------------
C  INPUTS
C  MAX_TEMP: Float maximum model dust temperature in Kelvin.
C  MIN_TEMP: Float minium model dust temperature in Kelvin.     
C  DTEMP: Integer step size in dust temperature in Kelvin.
C  MAX_TAU: Float maximum model dust opacity.
C  MIN_TAU: Float minimum model dust opacity.
C  DTAU: Float step size in dust opacity.
C  DIR: String directory for the output models. 
C  FILE_SED_INPUT: String input SED (two column file with
C   lambda in microns and lambda*F_lambda).
C 
C  OUTPUTS
C  FILE_SPEC: Dusty output file with the emerging flux in 
C   lambda*f_lambda (see section 4.2.2 in the Dusty manual).
C   File name in the form model_TEMPk_TAU.s001.
C  FILE_OUT: Dusty output file with parameters for each run
C   (see section 4.1 in the Dusty manual). File name in the 
C   form model_TEMPk_TAU.out.
C ----------------------------------------------------------
      IMPLICIT NONE
      
C     Constants - defined at the beginning
      INTEGER UNIT_INP
      PARAMETER (UNIT_INP = 10)
      REAL MAX_TEMP, MIN_TEMP
      PARAMETER (MAX_TEMP = 700.0, MIN_TEMP = 300.0)
      INTEGER DTEMP
      PARAMETER (DTEMP = 20)
      REAL MAX_TAU, MIN_TAU, DTAU
      PARAMETER (MAX_TAU = 0.050, MIN_TAU = 0.001, DTAU = 0.002)
      
C     Directory and file names (fixed length in F77)
      CHARACTER*7 DIR
      PARAMETER (DIR = 'models/')
      CHARACTER*13 FILE_SED_INPUT
      PARAMETER (FILE_SED_INPUT = 'SED_INPUT.dat')
      CHARACTER*7 DUSTY_EXE
      PARAMETER (DUSTY_EXE = './dusty')
      CHARACTER*9 DUSTY_INP
      PARAMETER (DUSTY_INP = 'multi.inp')
      CHARACTER*10 DUSTY_SPEC
      PARAMETER (DUSTY_SPEC = 'multi.s001')
      CHARACTER*9 DUSTY_OUT
      PARAMETER (DUSTY_OUT = 'multi.out')
      
C     Variables
      INTEGER NTEMP, NTAU, N, I, IOS, STATUS
      REAL TEMP, TAU
      CHARACTER*4 STEMP
      CHARACTER*6 STAU
      CHARACTER*100 FILE_SPEC, FILE_OUT, COMMAND

C     Calculate the number of temperature and tau steps
      NTEMP = INT((MAX_TEMP - MIN_TEMP) / DTEMP)
      NTAU = INT((MAX_TAU - MIN_TAU) / DTAU)

C     Check if output directory exists, create if not
      COMMAND = 'if [ ! -d ' // DIR // ' ]; then mkdir -p ' 
     &          // DIR // '; fi'
      CALL SYSTEM(COMMAND)

C     Process each temperature in the range
      DO 100 N = 0, NTEMP
         TEMP = MIN_TEMP + N * DTEMP
         
C        Process each opacity in the range
         DO 200 I = 0, NTAU
            TAU = MIN_TAU + I * DTAU
            
C           Format temperature and tau strings for filenames
            WRITE(STEMP, 10) NINT(TEMP)
10          FORMAT(I3)
            WRITE(STAU, 20) TAU
20          FORMAT(F5.3)
            
            WRITE(*,*) 'Processing: temp = ', STEMP, 'K  ',
     &                'tau = ', STAU

C           Construct output filenames
            FILE_SPEC = DIR // 'model_' // STEMP // 
     &                  'k_' // STAU // '.s001'
            FILE_OUT = DIR // 'model_' // STEMP // 
     &                 'k_' // STAU // '.out'
            
            WRITE(*,*) 'Output spec file: ', FILE_SPEC
            WRITE(*,*) 'Output parameters file: ', FILE_OUT

C           Create DUSTY input file for this run
            CALL CREATE_DUSTY_INPUT(TEMP, TAU, STEMP, STAU, 
     &                              FILE_SED_INPUT)
            
C           Run DUSTY and handle outputs
            CALL RUN_DUSTY_MODEL(DUSTY_SPEC, DUSTY_OUT,
     &                           FILE_SPEC, FILE_OUT)
            
            WRITE(*,*) 'Completed run for temp=', STEMP, 
     &                 ', tau=', STAU
200      CONTINUE
100   CONTINUE

      WRITE(*,*) 'All model calculations completed successfully.'
      END

C-------------------------------------------------------------
C  Subroutine to create the DUSTY input file for a specific run
C-------------------------------------------------------------
      SUBROUTINE CREATE_DUSTY_INPUT(TEMP, TAU, STEMP, STAU, 
     &                              FILE_SED_INPUT)
      IMPLICIT NONE
      
C     Input parameters
      REAL TEMP, TAU
      CHARACTER*4 STEMP
      CHARACTER*6 STAU
      CHARACTER*13 FILE_SED_INPUT
      
C     Open the input file for writing
      OPEN(UNIT=10, FILE='multi.inp', STATUS='UNKNOWN')

C     Write DUSTY input parameters
      WRITE(10, 1000) 
1000  FORMAT('   I PHYSICAL PARAMETERS')
      WRITE(10, 1010) 
1010  FORMAT('      1) External radiation:')
      WRITE(10, 1020) 
1020  FORMAT('         Spectrum =  4')
      WRITE(10, 1030) FILE_SED_INPUT
1030  FORMAT('        ',A)
      WRITE(10, 1040) 
1040  FORMAT(' ')
      WRITE(10, 1050) 
1050  FORMAT('      2) Dust Properties')
      WRITE(10, 1060) 
1060  FORMAT('  ')
      WRITE(10, 1070) 
1070  FORMAT('        2.1 Chemical composition')
      WRITE(10, 1080) 
1080  FORMAT('            optical properties index = 1')
      WRITE(10, 1090) 
1090  FORMAT('            Abundances for supported grain types:')
      WRITE(10, 1100) 
1100  FORMAT('      Sil-Ow  Sil-Oc  Sil-DL  grf-DL  amC-Hn  SiC-Pg')
      WRITE(10, 1110) 
1110  FORMAT('    x =  1.00    0.00    0.00    0.00    0.00    0.00')
      WRITE(10, 1120) 
1120  FORMAT(' ')
      WRITE(10, 1130) 
1130  FORMAT('        2.2 Grain size distribution')
      WRITE(10, 1140) 
1140  FORMAT('')
      WRITE(10, 1150) 
1150  FORMAT('         - size distribution = 1  % standard MRN')
      WRITE(10, 1160) 
1160  FORMAT(' ')
      WRITE(10, 1170) 
1170  FORMAT('        2.3 Dust temperature on inner boundary:')
      WRITE(10, 1180) 
1180  FORMAT(' ')
      WRITE(10, 1190) STEMP
1190  FORMAT('         - temperature = ',A,'K ')
      WRITE(10, 1200) 
1200  FORMAT(' ')
      WRITE(10, 1210) 
1210  FORMAT(' ')
      WRITE(10, 1220) 
1220  FORMAT('      3) Density Distribution')
      WRITE(10, 1230) 
1230  FORMAT('                density type = 3;     Y = 1.e4')
      WRITE(10, 1240) 
1240  FORMAT(' ')
      WRITE(10, 1250) 
1250  FORMAT('  ')
      WRITE(10, 1260) 
1260  FORMAT(' ')
      WRITE(10, 1270) 
1270  FORMAT('      4) Optical Depth')
      WRITE(10, 1280) 
1280  FORMAT('         - grid type = 1            %  linear grid')
      WRITE(10, 1290) 
1290  FORMAT('         - lambda0 = 0.55 micron    % optical depth')
      WRITE(10, 1300) STAU, STAU
1300  FORMAT('         - tau(min) = ',A,'; tau(max) = ',A,' %visual')
      WRITE(10, 1310) 
1310  FORMAT('         - number of models = 1            ')
      WRITE(10, 1320) 
1320  FORMAT('---------------------------------------------------')
      WRITE(10, 1330) 
1330  FORMAT(' ')
      WRITE(10, 1340) 
1340  FORMAT('   II NUMERICS ')
      WRITE(10, 1350) 
1350  FORMAT('       ')
      WRITE(10, 1360) 
1360  FORMAT('      - accuracy for flux conservation = 0.05  ')
      WRITE(10, 1370) 
1370  FORMAT('   -----------------------------------------------')
      WRITE(10, 1380) 
1380  FORMAT(' ')
      WRITE(10, 1390) 
1390  FORMAT('   III OUTPUT PARAMETERS       ')
      WRITE(10, 1400) 
1400  FORMAT(' ')
      WRITE(10, 1410) 
1410  FORMAT(' ')
      WRITE(10, 1420) 
1420  FORMAT('         FILE DESCRIPTION                FLAG      ')
      WRITE(10, 1430) 
1430  FORMAT(' -----------------------------------------------')
      WRITE(10, 1440) 
1440  FORMAT('  - verbosity flag;                       verbose = 1')
      WRITE(10, 1450) 
1450  FORMAT('        - properties of emerging spectra; fname.spp = 0')
      WRITE(10, 1460) 
1460  FORMAT('        - detailed spectra for each model; fname.s### = 2')
      WRITE(10, 1470) 
1470  FORMAT('        - images at specified wavelengths; fname.i### = 0')
      WRITE(10, 1480) 
1480  FORMAT('        - radial profiles for each model;  fname.r### = 0')
      WRITE(10, 1490) 
1490  FORMAT('        - detailed run-time messages;      fname.m### = 0')
      WRITE(10, 1500) 
1500  FORMAT('  --------------------------------------------------')
      WRITE(10, 1510) 
1510  FORMAT('  ')
      WRITE(10, 1520) 
1520  FORMAT(' ')
      WRITE(10, 1530) 
1530  FORMAT('   The end of the input parameters listing.     ')
      WRITE(10, 1540) 
1540  FORMAT('  ************************************************')
      WRITE(10, 1550) 
1550  FORMAT('  --------------------------------------------------')
      WRITE(10, 1560) 
1560  FORMAT('      The end of the input file.                    ')
      WRITE(10, 1570) 
1570  FORMAT('  ***************************************')
      
      CLOSE(10)
      END

C-------------------------------------------------------------
C  Subroutine to run DUSTY and copy the output files
C-------------------------------------------------------------
      SUBROUTINE RUN_DUSTY_MODEL(DUSTY_SPEC, DUSTY_OUT,
     &                           FILE_SPEC, FILE_OUT)
      IMPLICIT NONE
      
C     Input parameters
      CHARACTER*10 DUSTY_SPEC
      CHARACTER*9 DUSTY_OUT
      CHARACTER*100 FILE_SPEC, FILE_OUT
      
C     Local variables
      CHARACTER*200 COMMAND
      
C     Run DUSTY
      CALL SYSTEM('./dusty')
      
C     Copy the output files to their destination
      COMMAND = 'cp ' // DUSTY_SPEC // ' ' // FILE_SPEC
      CALL SYSTEM(COMMAND)
      
      COMMAND = 'cp ' // DUSTY_OUT // ' ' // FILE_OUT
      CALL SYSTEM(COMMAND)
      
      END
