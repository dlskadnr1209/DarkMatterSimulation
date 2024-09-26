      SUBROUTINE INIT_CLUSTERS(POS, VEL, NUM_PARTICLES,
     &    CLUSTER1_FRACTION, R_C, BETA)
        IMPLICIT NONE
        REAL POS(3, *), VEL(3, *), R_C, BETA,V_MAG,V_THETA,V_PHI
        INTEGER NUM_PARTICLES
        REAL CLUSTER1_FRACTION
        INTEGER I
        REAL R, THETA, PHI, X, Y, Z,V_MIN,V_MAX
        REAL CLUSTER1_CENTER(3), CLUSTER2_CENTER(3)
        REAL RAND_ARRAY(3)
        INTEGER SEED(1)

        ! Define cluster centers in terms of normalized box size
        CLUSTER1_CENTER(1) = 0.25D0
        CLUSTER1_CENTER(2) = 0.5D0
        CLUSTER1_CENTER(3) = 0.5D0

        CLUSTER2_CENTER(1) = 0.75D0
        CLUSTER2_CENTER(2) = 0.5D0
        CLUSTER2_CENTER(3) = 0.5D0
        V_MIN=-0.1D0
        V_MAX=1.0D0
        ! Seed for random number generator
        CALL RANDOM_NUMBER(RAND_ARRAY)
        SEED(1) = INT(SUM(RAND_ARRAY * 100000.0D0)) + 1
        CALL RANDOM_SEED(PUT=SEED)

        ! Initialize positions and velocities for Cluster 1
        DO I = 1, INT(NUM_PARTICLES * CLUSTER1_FRACTION)
          CALL RANDOM_NUMBER(RAND_ARRAY)
          R = R_C * (RAND_ARRAY(1) ** (1.0D0 / BETA))
          THETA = ACOS(2.0D0 * RAND_ARRAY(2) - 1.0D0)
          PHI = 2.0D0 * 3.141592653589793D0 * RAND_ARRAY(3)
c          PRINT *,R_C,BETA
          X = R * SIN(THETA) * COS(PHI)
          Y = R * SIN(THETA) * SIN(PHI)
          Z = R * COS(THETA)
c          print *, R,X,Y,Z
          POS(1, I) = CLUSTER1_CENTER(1) + X
          POS(2, I) = CLUSTER1_CENTER(2) + Y
          POS(3, I) = CLUSTER1_CENTER(3) + Z
          CALL RANDOM_NUMBER(RAND_ARRAY)
 
          V_MAG = V_MIN + (V_MAX - V_MIN) * RAND_ARRAY(1)  ! Random velocity between V_MIN and V_MAX
          V_THETA = 2.0D0 * 3.141592653589793D0 * RAND_ARRAY(1)
          V_PHI = ACOS(2.0D0 * RAND_ARRAY(2) - 1.0D0)

          VEL(1, I) = V_MAG * SIN(V_PHI) * COS(V_THETA) + 1.0D0
          VEL(2, I) = V_MAG * SIN(V_PHI) * SIN(V_THETA)
          VEL(3, I) = V_MAG * COS(V_PHI)
        ENDDO

        ! Initialize positions and velocities for Cluster 2
        DO I = INT(NUM_PARTICLES * CLUSTER1_FRACTION) + 1, NUM_PARTICLES
          CALL RANDOM_NUMBER(RAND_ARRAY)
          R = R_C * (RAND_ARRAY(1) ** (1.0D0 / BETA))
          THETA = ACOS(2.0D0 * RAND_ARRAY(2) - 1.0D0)
          PHI = 2.0D0 * 3.141592653589793D0 * RAND_ARRAY(3)

          X = R * SIN(THETA) * COS(PHI)
          Y = R * SIN(THETA) * SIN(PHI)
          Z = R * COS(THETA)

          POS(1, I) = CLUSTER2_CENTER(1) + X
          POS(2, I) = CLUSTER2_CENTER(2) + Y
          POS(3, I) = CLUSTER2_CENTER(3) + Z
          CALL RANDOM_NUMBER(RAND_ARRAY)

          V_MAG = V_MIN + (V_MAX - V_MIN) * RAND_ARRAY(1)
          V_THETA = 2.0D0 * 3.141592653589793D0 * RAND_ARRAY(1)
          V_PHI = ACOS(2.0D0 * RAND_ARRAY(2) - 1.0D0)

          VEL(1, I) = V_MAG * SIN(V_PHI) * COS(V_THETA) - 1.0D0
          VEL(2, I) = V_MAG * SIN(V_PHI) * SIN(V_THETA)
          VEL(3, I) = V_MAG * COS(V_PHI)

        ENDDO

      RETURN
      END SUBROUTINE INIT_CLUSTERS

