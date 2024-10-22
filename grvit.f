      subroutine grvit
c
c
      implicit none
c
c      include "par"
      include "com"
c
      integer*8 ix, iy, iz, ix1, iy1, iz1
      real twopi, fps, xfac, xk
c
c
      twopi = 8.*atan(1.)
      fps   = twopi*twopi
      fknyq = real(nx)/2.
c
c     read density field 
c
c      do iz=1,nz
c         do iy=1,ny
c            do ix=1,nx
c               tempa5(ix+1,iy+1,iz+1) = (-1.0 + tempa4(ix+1,iy+1,iz+1)
c     +         	                      + tempa5(ix+1,iy+1,iz+1))/aaa
c            enddo
c         enddo
c      enddo
c
c
c      do iz=1,nz
c         do iy=1,ny
c            do ix=1,nx
c               tempa4(ix,iy,iz) = tempa5(ix+1,iy+1,iz+1)
c            enddo
c         enddo
c      enddo
c
c
      do iz=1,nz
         do iy=1,ny
            do ix=1,nx
               tempa4(ix,iy,iz) = rho(ix,iy,iz)
               tempa5(ix,iy,iz) = 0.
            enddo
         enddo
      enddo
c
c     fourier transform
c
      do iz=1,nz
         do iy=1,ny
            do ix=1,nx
               cw4(ix) = tempa4(ix,iy,iz)-tempa5(ix,iy,iz)*(0.,1.)
            enddo
            call fftcc (cw4,nw,iwk2,wk2)
            do ix=1,nx
               tempa4(ix,iy,iz) = real(cw4(ix))
               tempa5(ix,iy,iz) = -aimag(cw4(ix))
            enddo
         enddo
      enddo
c
      do ix=1,nx
         do iz=1,nz
            do iy=1,ny
               cw4(iy) = tempa4(ix,iy,iz)-tempa5(ix,iy,iz)*(0.,1.)
            enddo
            call fftcc (cw4,nw,iwk2,wk2)
            do iy=1,ny
               tempa4(ix,iy,iz) = real(cw4(iy))
               tempa5(ix,iy,iz) = -aimag(cw4(iy))
            enddo
         enddo
      enddo
c
      do iy=1,ny
         do ix=1,nx
            do iz=1,nz
               cw4(iz) = tempa4(ix,iy,iz)-tempa5(ix,iy,iz)*(0.,1.)
            enddo
            call fftcc (cw4,nw,iwk2,wk2)
            do iz=1,nz
               tempa4(ix,iy,iz) = real(cw4(iz))
               tempa5(ix,iy,iz) = -aimag(cw4(iz))
            enddo
         enddo
      enddo
c
c     gravitational potential in fourier space
c      
      do iz=1,nz
         iz1=iz-1
         if(iz.gt.nz/2+1) iz1=iz1-nz
         do iy=1,ny
            iy1=iy-1
            if(iy.gt.ny/2+1) iy1=iy1-ny
            do ix=1,nx
               ix1=ix-1
               if(ix.gt.nx/2+1) ix1=ix1-nx
c
               xfac  = 0.
               xk = sqrt(ix1*ix1*fps +iy1*iy1*fps +iz1*iz1*fps)
               if(xk.gt.0.0 .and. xk .le. (twopi*fknyq)) 
     +           xfac = -1./xk/xk
               tempa4(ix,iy,iz) = tempa4(ix,iy,iz)*xfac*v
               tempa5(ix,iy,iz) = tempa5(ix,iy,iz)*xfac*v
            enddo
         enddo
      enddo
c
c     inverse fourier transform
c
      do iz=1,nz
         do iy=1,ny
            do ix=1,nx
               cw4(ix) = tempa4(ix,iy,iz)+tempa5(ix,iy,iz)*(0.,1.)
            enddo
            call fftcc (cw4,nw,iwk2,wk2)
            do ix=1,nx
               tempa4(ix,iy,iz) = real(cw4(ix))
               tempa5(ix,iy,iz) = aimag(cw4(ix))
            enddo
         enddo
      enddo
c
      do ix=1,nx
         do iz=1,nz
            do iy=1,ny
               cw4(iy) = tempa4(ix,iy,iz)+tempa5(ix,iy,iz)*(0.,1.)
            enddo
            call fftcc (cw4,nw,iwk2,wk2)
            do iy=1,ny
               tempa4(ix,iy,iz) = real(cw4(iy))
               tempa5(ix,iy,iz) = aimag(cw4(iy))
            enddo
         enddo
      enddo
c
      do iy=1,ny
         do ix=1,nx
            do iz=1,nz
               cw4(iz) = tempa4(ix,iy,iz)+tempa5(ix,iy,iz)*(0.,1.)
            enddo
            call fftcc (cw4,nw,iwk2,wk2)
            do iz=1,nz
               phi(ix,iy,iz) = real(cw4(iz))
            enddo
         enddo
      enddo
c
c     put one point in each side along three directions
c
c      do iz=nz+1,2,-1
c         do iy=ny+1,2,-1
c            do ix=nx+1,2,-1
c               phi(ix,iy,iz) = phi(ix-1,iy-1,iz-1)
c            enddo
c         enddo
c      enddo
c
c      do iz=2,nz+1
c         do iy=2,ny+1
c            phi(1,iy,iz) = phi(nx+1,iy,iz)
c            phi(nx+2,iy,iz) = phi(2,iy,iz)
c         enddo
c      enddo
c
c      do iz=2,nz+1
c         do ix=1,nx+2
c            phi(ix,1,iz) = phi(ix,ny+1,iz)
c            phi(ix,ny+2,iz) = phi(ix,2,iz)
c         enddo
c      enddo
c
c      do iy=1,ny+2
c         do ix=1,nx+2
c            phi(ix,iy,1) = phi(ix,iy,nz+1)
c            phi(ix,iy,nz+2) = phi(ix,iy,2)
c         enddo
c      enddo
c
c
      return
      end
