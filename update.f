      SUBROUTINE UPDATE_POS_VEL(POS, VEL, ACC, DT)
      IMPLICIT NONE
      REAL POS(3), VEL(3), ACC(3), DT

      ! Update velocity
      VEL(1) = VEL(1) + ACC(1) * DT
      VEL(2) = VEL(2) + ACC(2) * DT
      VEL(3) = VEL(3) + ACC(3) * DT

      ! Update position
      POS(1) = POS(1) + VEL(1) * DT
      POS(2) = POS(2) + VEL(2) * DT
      POS(3) = POS(3) + VEL(3) * DT

      RETURN
      END SUBROUTINE UPDATE_POS_VEL

