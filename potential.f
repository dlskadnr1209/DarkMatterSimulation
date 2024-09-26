      SUBROUTINE COMPUTE_POTENTIAL(DENSITY, POTENTIAL, NX, NY, NZ,G)
      IMPLICIT NONE
      REAL DENSITY(NX, NY, NZ), POTENTIAL(NX, NY, NZ)
      REAL G,POTENTIAL_THRESHOLD
      INTEGER NX, NY, NZ
      INTEGER I, J, K, KX_INDEX, KY_INDEX, KZ_INDEX
      REAL KX, KY, KZ, K_SQUARED
      REAL PI
      COMPLEX DENS_FFT(NX/2+1, NY, NZ), POT_FFT(NX/2+1, NY, NZ)
      COMPLEX, ALLOCATABLE :: FULL_POT_FFT(:,:,:)
      PARAMETER (PI = 3.141592653589793D0)
      REAL :: min_threshold, max_threshold

      ALLOCATE(FULL_POT_FFT(NX, NY, NZ))
      DO K = 1, NZ
        DO J = 1, NY
          DO I = 1, NX/2+1
              DENS_FFT(I, J, K) = CMPLX(DENSITY(I, J, K), 0.0D0)
c              PRINT *, DENS_FFT(I,J,K)
          ENDDO
        ENDDO
      ENDDO
c      PRINT *, DENSITY(32,64,64)
      ! Forward FFT to convert density to Fourier space
      CALL FFT3RKs(DENS_FFT, NX, NY, NZ)
c      PRINT *, DENS_FFT(32,64,64)
c      PRINT *, DENS_FFT
c      PRINT *, DENS_FFT(16,64,64)
      ! Compute potential in Fourier space
      DO K = 1, NZ
         KZ_INDEX = K - 1
         IF (KZ_INDEX > NZ/2) KZ_INDEX = KZ_INDEX - NZ
         DO J = 1, NY
            KY_INDEX = J - 1
            IF (KY_INDEX > NY/2) KY_INDEX = KY_INDEX - NY

            DO I = 1, NX/2+1
               KX_INDEX = I - 1
               IF (KX_INDEX > NX/2) KX_INDEX = KX_INDEX - NX
                
               KX = 2.0D0 * PI * KX_INDEX
               KY = 2.0D0 * PI * KY_INDEX
               KZ = 2.0D0 * PI * KZ_INDEX

               K_SQUARED = KX*KX + KY*KY + KZ*KZ
               IF (K_SQUARED .NE. 0.0D0) THEN
                  POT_FFT(I, J, K) = -4.0D0 * PI * G *
     &          DENS_FFT(I,J,K) / K_SQUARED
c                PRINT *, POT_FFT(I,J,K)
               ELSE
                  POT_FFT(I, J, K) = 0.0D0  ! k = 0 mode is set to zero
               END IF
c               IF (K .eq. 64) THEN
c                 IF (J .eq. 64) THEN
c                   IF (I .eq. 32) THEN
c                     PRINT *, KX,KY,KZ,K_SQUARED
c                     PRINT *, POT_FFT(I,J,K),DENS_FFT(I,J,K)
c                     
c                   ENDIF
c                 ENDIF
c               ENDIF
               
            ENDDO
         ENDDO
      ENDDO
      ! Inverse FFT to convert potential back to real space
      CALL FFT3KRs(POT_FFT, NX, NY, NZ)
        DO K = 1, NZ
          DO J = 1, NY
            DO I = 1, NX/2
              FULL_POT_FFT(I, J, K) = POT_FFT(I, J, K)
            ENDDO
            DO I = NX/2+1, NX
              FULL_POT_FFT(I, J, K) = CONJG(POT_FFT(NX-I+1, J, K))
            ENDDO
          ENDDO
        ENDDO
c        min_threshold = 1.0D-10
c        max_threshold = 1.0D+10
        DO K = 1, NZ
          DO J = 1, NY
            DO I = 1, NX
              POTENTIAL(I, J, K) = REAL(FULL_POT_FFT(I, J, K))
c              PRINT *, POTENTIAL(I,J,K)
c              IF (ABS(POTENTIAL(I, J, K)) < min_threshold .OR.
c     &             ABS(POTENTIAL(I, J, K)) > max_threshold) THEN
c                POTENTIAL(I, J, K) = 0.0D0
c              END IF
            ENDDO
          ENDDO
        ENDDO
        DEALLOCATE(FULL_POT_FFT)
c        PRINT *, SHAPE(POTENTIAL)
        RETURN
      
      RETURN
      END SUBROUTINE COMPUTE_POTENTIAL
      SUBROUTINE SAVE_POTENTIAL(POTENTIAL, NX, NY, NZ, FILENAME)
      IMPLICIT NONE
      REAL POTENTIAL(NX, NY, NZ)
      INTEGER NX, NY, NZ,K,J,I
      CHARACTER(*) FILENAME

      INTEGER FILE_UNIT, IERR

      OPEN(UNIT=16, FILE=FILENAME, STATUS='REPLACE', IOSTAT=IERR)
      IF (IERR .NE. 0) THEN
         PRINT *, 'Error opening file: ', FILENAME
         STOP
      ENDIF

      DO K = 1, NZ
         DO J = 1, NY
            DO I = 1, NX
               WRITE(16, '(E15.6E3)') POTENTIAL(I, J, K)
            ENDDO
         ENDDO
      ENDDO
      CLOSE(16)

      RETURN
      END

