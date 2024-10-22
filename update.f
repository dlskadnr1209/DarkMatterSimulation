      subroutine update
      implicit none
      include "com"
      integer i

      do i = 1, numprt
      vel(1,i) = vel(1,i) + acc(1,i) * dt
      vel(2,i) = vel(2,i) + acc(2,i) * dt
      vel(3,i) = vel(3,i) + acc(3,i) * dt

      pos(1,i) = pos(1,i) + vel(1,i) * dt
      pos(2,i) = pos(2,i) + vel(2,i) * dt
      pos(3,i) = pos(3,i) + vel(3,i) * dt
      enddo
      return
      end
