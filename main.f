      PROGRAM DARK_MATTER_PM_SIMULATION
      IMPLICIT NONE

      ! Variable Declarations
      REAL T0, R0, U0, N0, M0
      INTEGER NUM_PARTICLES, NX, NY, NZ, IERR, NSTEP, FILE_CNT,NW
      REAL G
      REAL MASS_PARTICLE
      REAL DT, R_C, BETA, DL
      REAL, ALLOCATABLE :: POS(:, :), VEL(:, :), ACC(:, :)
      REAL, ALLOCATABLE ::  DENSITY(:, :, :), POTENTIAL(:, :, :)
      REAL CLUSTER1_FRACTION
      INTEGER STEP, I
      CHARACTER*20 FNAME_DENS, FNAME_POT

      ! Parameter Initialization (CGS Unit System)
      PARAMETER (G = 6.67430D-8)                  ! cm^3 g^-1 s^-2
      PARAMETER (T0 = 3.1557D16 * 1.0D0 * 5.2D0)  ! 5.2gy s (시간 단위)
      PARAMETER (R0 = 3.086D24 * 4.0D0 * 5.0D0)          ! 20Mpc cm (거리 단위)
      PARAMETER (U0 = R0/T0)      ! cm/s (속도 단위)
      PARAMETER (M0= 1.0D15*1.989D33)  ! solar mass
c      PARAMETER (N0 = 1.0D-3)                     ! cm^-3 (밀도 단위)
c      PARAMETER (PMASS = 1.6726D-24)              ! g (양성자 질량)
      PARAMETER (NX=128,NY=NX,NZ=NX,NW=NX)
      PARAMETER (NUM_PARTICLES=10000000)
      PARAMETER (MASS_PARTICLE=1.0D0/DBLE(NUM_PARTICLES))
      PARAMETER (DL=1.0D0/DBLE(NW))
      ! 무차원화된 단위 계산
      PARAMETER (NSTEP=100000)
      PARAMETER (DT=1.0D0/DBLE(NSTEP))
      PARAMETER (R_C=0.5D0)
      PARAMETER (BETA=0.5D0)
      PARAMETER (CLUSTER1_FRACTION = 0.5D0)
c      PARAMETER (DENS_SCAL=M0*1.989D48/R0**3)
      FILE_CNT = 0

      ALLOCATE(POS(3, NUM_PARTICLES))
      ALLOCATE(VEL(3, NUM_PARTICLES))
      ALLOCATE(ACC(3, NUM_PARTICLES))
      ALLOCATE(DENSITY(NX, NY, NZ))
      ALLOCATE(POTENTIAL(NX, NY, NZ))

      ! Initialize Clusters
      CALL INIT_CLUSTERS(POS, VEL, NUM_PARTICLES, 
     &     CLUSTER1_FRACTION, R_C, BETA)
      PRINT *, 'Clusters initialized'

      ! Compute Initial Density
      CALL COMPUTE_DENSITY(POS, NUM_PARTICLES, DENSITY, DL, 
     &     NX, NY, NZ, MASS_PARTICLE)
      CALL SAVE_DENSITY(DENSITY, NX, NY, NZ, 'initial_density.dat')
      PRINT *, 'Initial density computation done'

      ! Compute Initial Potential
      CALL COMPUTE_POTENTIAL(DENSITY, POTENTIAL, NX, NY, NZ,
     &     G)
      PRINT *, "Potential computation end"
      CALL SAVE_POTENTIAL(POTENTIAL, NX, NY, NZ,
     &  'initial_potential.dat')
      PRINT *, 'Initial potential computation done'
      ! Main Simulation Loop
      DO STEP = 1, NSTEP
c         PRINT *, 'Step: ', STEP
         DO I = 1, NUM_PARTICLES
            CALL COMPUTE_ACCEL(ACC(1:3, I), POS(1:3, I), DENSITY, 
     &           POTENTIAL, NX,NY,NZ)
         ENDDO
         DO I = 1, NUM_PARTICLES
            CALL UPDATE_POS_VEL(POS(1:3, I), VEL(1:3, I), 
     &           ACC(1:3, I), DT)
         ENDDO

         CALL COMPUTE_DENSITY(POS, NUM_PARTICLES, DENSITY,DL,
     &     NX, NY, NZ, MASS_PARTICLE)
c         PRINT *, 'Density recomputed for step', STEP
         
         CALL COMPUTE_POTENTIAL(DENSITY, POTENTIAL, NX, NY, NZ,
     &     G)

         IF (MOD(STEP, 10000) .EQ. 0) THEN
            FILE_CNT = FILE_CNT + 1
            WRITE(FNAME_DENS, '(A,I4.4)') 'dens', FILE_CNT
            FNAME_DENS = FNAME_DENS // '.dat'
            CALL SAVE_DENSITY(DENSITY, NX, NY, NZ, FNAME_DENS)

            WRITE(FNAME_POT, '(A,I4.4)') 'poten', FILE_CNT
            FNAME_POT = FNAME_POT // '.dat'
            CALL SAVE_DENSITY(POTENTIAL, NX, NY, NZ, FNAME_POT)
            PRINT *, "STEP : ", STEP
         ENDIF
      ENDDO

      PRINT *, 'Simulation completed'

      STOP
      END

