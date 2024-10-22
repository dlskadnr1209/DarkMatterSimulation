      subroutine init

      implicit none
      include "com"
      real*8 vmag, vtheta, vphi, r, theta, phii
      real*8 x, y, z, vmin, vmax
      integer ii

      cluster1_center(1) = 0.25D0
      cluster1_center(2) = 0.5D0
      cluster1_center(3) = 0.5D0

      cluster2_center(1) = 0.75D0
      cluster2_center(2) = 0.5D0
      cluster2_center(3) = 0.5D0
      vmin = -0.1D0
      vmax = 0.1D0
      call random_number(rand_array)
      seed(1) = int(sum(rand_array * 10000.0D0)) +1
      call random_seed(put=seed)
      do ii=1, int(numprt*cluster_fraction)
        call random_number(rand_array)
        r = rc * (rand_array(1) ** (1.0D0/BETA))
        theta = acos(2.0D0 * rand_array(2) - 1.0D0)
        phii = 2.0D0 * pi * rand_array(3)
        x = r * sin(theta) * cos(phii)
        y = r * sin(theta) * sin(phii)
        z = r * cos(theta)
        pos(1,ii) = cluster1_center(1)+x
        pos(2,ii) = cluster1_center(2)+y
        pos(3,ii) = cluster1_center(3)+z
        call random_number(rand_array)
        vmag = vmin + (vmax - vmin) * rand_array(1)
        vtheta = 2.0D0 * pi * rand_array(1)
        vphi = acos(2.0D0 * rand_array(2) - 1.0D0)
        vel(1,ii) = vmag * sin(vphi) * cos(vtheta) + 1.0D0
        vel(2,ii) = vmag * sin(vphi) * sin(vtheta)
        vel(3,ii) = vmag * cos(vphi)
      enddo

      do ii=int(numprt*cluster_fraction) + 1, numprt
        call random_number(rand_array)
        r = rc * (rand_array(1) ** (1.0D0/BETA))
        theta = acos(2.0D0 * rand_array(2) - 1.0D0)
        phii = 2.0D0 * pi * rand_array(3)
        x = r * sin(theta) * cos(phii)
        y = r * sin(theta) * sin(phii)
        z = r * cos(theta)
        pos(1,ii) = cluster2_center(1)+x
        pos(2,ii) = cluster2_center(2)+y
        pos(3,ii) = cluster2_center(3)+z
        call random_number(rand_array)
        vmag = vmin + (vmax - vmin) * rand_array(1)
        vtheta = 2.0D0 * pi * rand_array(1)
        vphi = acos(2.0D0 * rand_array(2) - 1.0D0)
        vel(1,ii) = vmag * sin(vphi) * cos(vtheta) - 1.0D0
        vel(2,ii) = vmag * sin(vphi) * sin(vtheta)
        vel(3,ii) = vmag * cos(vphi)
      enddo
      return
      end      
