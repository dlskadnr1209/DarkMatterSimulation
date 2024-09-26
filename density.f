      SUBROUTINE COMPUTE_DENSITY(POS, NUM_PARTICLES, DENSITY, 
     &     DL,NX, NY, NZ, MASS_PARTICLE)
      IMPLICIT NONE
      REAL POS(3, *), DENSITY(NX, NY, NZ), DL
      INTEGER NUM_PARTICLES, NX, NY, NZ
      INTEGER I, IX, IY, IZ
      REAL DX, DY, DZ, MASS_PARTICLE,VOLUME


      DX = DL
      DY = DL
      DZ = DL
      VOLUME = DX * DY * DZ
c      PRINT *, VOLUME_NORM,MASS_PARTICLE_NORM
      ! Initialize density array to zero
      DENSITY=0.0D0
      ! Compute density
      DO I = 1, NUM_PARTICLES
         ! Convert position to grid index with periodic boundary conditions
         IX = MOD(INT(POS(1, I) / DX), NX) + 1
         IY = MOD(INT(POS(2, I) / DY), NY) + 1
         IZ = MOD(INT(POS(3, I) / DZ), NZ) + 1

         ! Ensure indices are within bounds and update density
         IF (IX >= 1 .AND. IX <= NX .AND. IY >= 1 .AND. IY <= NY 
     &        .AND. IZ >= 1 .AND. IZ <= NZ) THEN
            DENSITY(IX, IY, IZ) = DENSITY(IX, IY, IZ) 
     &                        +MASS_PARTICLE / VOLUME
c            PRINT *, IX,IY,IZ,MASS_PARTICLE/VOLUME
         ENDIF
      ENDDO

      RETURN
      END

      SUBROUTINE SAVE_DENSITY(DENSITY, NX, NY, NZ, FILENAME)
      IMPLICIT NONE
      REAL DENSITY(NX, NY, NZ)
      INTEGER NX, NY, NZ,K,J,I
      CHARACTER(*) FILENAME

      INTEGER FILE_UNIT, IERR

      OPEN(UNIT=15, FILE=FILENAME, STATUS='REPLACE', IOSTAT=IERR)
      IF (IERR .NE. 0) THEN
         PRINT *, 'Error opening file: ', FILENAME
         STOP
      ENDIF

      DO K = 1, NZ
         DO J = 1, NY
            DO I = 1, NX
               WRITE(15, '(E15.6E3)') DENSITY(I, J, K)
            ENDDO
         ENDDO
      ENDDO
      CLOSE(15)

      RETURN
      END

