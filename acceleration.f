      SUBROUTINE COMPUTE_ACCEL(ACC, POS, DENSITY, POTENTIAL,
     & NX,NY,NZ)
      IMPLICIT NONE
      REAL ACC(3), POS(3), DENSITY(NX, NY, NZ), POTENTIAL(NX, NY, NZ)
      INTEGER NX, NY, NZ
      INTEGER IX, IY, IZ
      REAL DX, DY, DZ

      DX = 1.0D0/DBLE(NX)
      DY = 1.0D0/DBLE(NY)
      DZ = 1.0D0/DBLE(NZ)

      ! Convert position to grid indices
      IX = MOD(INT(POS(1) / DX), NX) + 1
      IY = MOD(INT(POS(2) / DY), NY) + 1
      IZ = MOD(INT(POS(3) / DZ), NZ) + 1

      ! Compute acceleration from the gradient of the potential
      IF (IX >= 1 .AND. IX <= NX .AND. IY >= 1 .AND. 
     &       IY <= NY .AND. IZ >= 1 .AND. IZ <= NZ) THEN
        IF (IX == 1) THEN
            ACC(1) = -(POTENTIAL(IX+1, IY, IZ) - POTENTIAL(IX, IY, IZ))
     &    / DX
        ELSE IF (IX == NX) THEN
            ACC(1) = -(POTENTIAL(IX, IY, IZ) - POTENTIAL(IX-1, IY, IZ))
     &   / DX
        ELSE
            ACC(1) = -(POTENTIAL(IX+1, IY, IZ) 
     &  - POTENTIAL(IX-1, IY, IZ)) / (2.0D0 * DX)
        END IF

        IF (IY == 1) THEN
            ACC(2) = -(POTENTIAL(IX, IY+1, IZ) - POTENTIAL(IX, IY, IZ))
     &  / DY
        ELSE IF (IY == NY) THEN
            ACC(2) = -(POTENTIAL(IX, IY, IZ) - POTENTIAL(IX, IY-1, IZ))
     &  / DY
        ELSE
            ACC(2) = -(POTENTIAL(IX, IY+1, IZ) 
     & - POTENTIAL(IX, IY-1, IZ)) / (2.0D0 * DY)
        END IF

        IF (IZ == 1) THEN
            ACC(3) = -(POTENTIAL(IX, IY, IZ+1) - POTENTIAL(IX, IY, IZ))
     & / DZ
        ELSE IF (IZ == NZ) THEN
            ACC(3) = -(POTENTIAL(IX, IY, IZ) - POTENTIAL(IX, IY, IZ-1))
     & / DZ
        ELSE
            ACC(3) = -(POTENTIAL(IX, IY, IZ+1)
     & - POTENTIAL(IX, IY, IZ-1)) / (2.0D0 * DZ)
        END IF
      ELSE
        ACC(1) = 0.0D0
        ACC(2) = 0.0D0
        ACC(3) = 0.0D0
      END IF
c      PRINT *, ACC(1),ACC(2),ACC(3)
      RETURN
      END SUBROUTINE COMPUTE_ACCEL


