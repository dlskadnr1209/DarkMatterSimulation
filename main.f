      program main
c
      implicit none
c
      include "com"

c
      fcnt = 0
      nstep = 0
      call init
c
  100 continue
      nstep = nstep +1
c      print *, nstep
c
c
      call dens
c
      call poten
c
      call update
c
c
      if(mod(nstep,dout).eq. 1) then
        call out
      endif

      if(nstep.le.mstep) then
        goto 100
      endif
c
c
c
      stop
      end
