      subroutine out
      implicit none
c
      include 'com'
      integer ix, iy, iz      

      write(dens_dat, '("dens", I4.4, ".dat")') fcnt
      write(poten_dat, '("poten", I4.4, ".dat")') fcnt
      open(unit=10, file=dens_dat)
      open(unit=11, file=poten_dat)
      write(10,'(E15.6E3)') (((rho(ix,iy,iz),ix=1,nx),iy=1,ny),iz=1,nz)
      write(11,'(E15.6E3)') (((phi(ix,iy,iz),ix=1,nx),iy=1,ny),iz=1,nz)
      close(10)
      close(11)
      fcnt = fcnt + 1
      print *, fcnt
      return
      end
