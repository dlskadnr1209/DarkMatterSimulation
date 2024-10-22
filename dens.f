      subroutine dens

      implicit none
      include "com"      
      integer ii,ix,iy,iz
      rho=0.0D0
c      print *, massprt,v 
      do ii=1, numprt
        ix = mod(int(pos(1,ii)/dl),nx) +1
        iy = mod(int(pos(2,ii)/dl),ny) +1
        iz = mod(int(pos(3,ii)/dl),nz) +1
        if (ix >=1 .and. ix <= nx .and. iy >= 1 .and. iy <= ny
     & .and. iz>= 1 .and. iz <= nz) then
          rho(ix,iy,iz)= rho(ix,iy,iz)+massprt/v
        endif
      enddo
      return
      end

