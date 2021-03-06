c
c
c     ############################################################
c     ##  COPYRIGHT (C) 1995 by Yong Kong & Jay William Ponder  ##
c     ##                  All Rights Reserved                   ##
c     ############################################################
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine rotpole  --  rotate multipoles to global frame  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "rotpole" constructs the set of atomic multipoles in the global
c     frame by applying the correct rotation matrix for each site
c
c
      subroutine rotpole
      use sizes
      use mpole
      use kpolr_aniso
      use polar_aniso
      implicit none
      integer i
      real*8 a(3,3)
      integer j
c
c
c     rotate the atomic multipoles at each site in turn
c
      do i = 1, npole
c         if(i.ne.4) then
         call rotmat (i,a)
c         print*,"i=",i,"After rotmat"
         call rotsite (i,a)
c         print*,"i=",i,"After rotsite"
c         end if
      end do
c      print*,"After rotpole"

c   FOR TESTING PURPOSES USING ISOTROPIC POLARIZ
c      do i=1,npole
c         do j=1,9
c            ranisopolarity(j,i)=anisopolarity(j,i)
c         end do
c      end do
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine rotmat  --  find global frame rotation matrix  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "rotmat" finds the rotation matrix that converts from the local
c     coordinate system to the global frame at a multipole site
c
c
      subroutine rotmat (i,a)
      use sizes
      use atoms
      use mpole
      use atomid
      implicit none
      integer i,ii
      integer ix,iy,iz
      real*8 r,dot
      real*8 random
      real*8 xi,yi,zi
      real*8 dx,dy,dz
      real*8 dx1,dy1,dz1
      real*8 dx2,dy2,dz2
      real*8 dx3,dy3,dz3
      real*8 a(3,3)
c
c
c     get coordinates and frame definition for the multipole site
c
      ii = ipole(i)
      xi = x(ii)
      yi = y(ii)
      zi = z(ii)
      ix = xaxis(i)
      iy = yaxis(i)
      iz = zaxis(i)
c
c     use the identity matrix as the default rotation matrix
c
      a(1,1) = 1.0d0
      a(2,1) = 0.0d0
      a(3,1) = 0.0d0
      a(1,3) = 0.0d0
      a(2,3) = 0.0d0
      a(3,3) = 1.0d0
c
c     Z-Only method rotation matrix elements for z-axis only
c
      if (polaxe(i) .eq. 'Z-Only') then
c         print*,i,"Atom ",name(i)," Polaxe ",polaxe(i)
c         print*,i,"Atom ",name(i),ix," Atom X-axis ",name(ix)
c         print*,i,"Atom ",name(i),iy," Atom Y-axis ",name(iy)
c         print*,i,"Atom ",name(i),iz," Atom Z-axis ",name(iz)
         dx = x(iz) - xi
         dy = y(iz) - yi
         dz = z(iz) - zi
         r = sqrt(dx*dx + dy*dy + dz*dz)
         a(1,3) = dx / r
         a(2,3) = dy / r
         a(3,3) = dz / r
         dx = random ()
         dy = random ()
         dz = random ()
         dot = dx*a(1,3) + dy*a(2,3) + dz*a(3,3)
         dx = dx - dot*a(1,3)
         dy = dy - dot*a(2,3)
         dz = dz - dot*a(3,3)
         r = sqrt(dx*dx + dy*dy + dz*dz)
         a(1,1) = dx / r
         a(2,1) = dy / r
         a(3,1) = dz / r
c
c     Z-then-X method rotation matrix elements for z- and x-axes
c
      else if (polaxe(i) .eq. 'Z-then-X') then
c         print*,i,"Atom ",name(i)," Polaxe ",polaxe(i)
c         print*,i,"Atom ",name(i),ix," Atom X-axis ",name(ix)
c         print*,i,"Atom ",name(i),iy," Atom Y-axis ",name(iy)
c         print*,i,"Atom ",name(i),iz," Atom Z-axis ",name(iz)
         dx = x(iz) - xi
         dy = y(iz) - yi
         dz = z(iz) - zi
         r = sqrt(dx*dx + dy*dy + dz*dz)
         a(1,3) = dx / r
         a(2,3) = dy / r
         a(3,3) = dz / r
         dx = x(ix) - xi
         dy = y(ix) - yi
         dz = z(ix) - zi
         dot = dx*a(1,3) + dy*a(2,3) + dz*a(3,3)
         dx = dx - dot*a(1,3)
         dy = dy - dot*a(2,3)
         dz = dz - dot*a(3,3)
         r = sqrt(dx*dx + dy*dy + dz*dz)
         a(1,1) = dx / r
         a(2,1) = dy / r
         a(3,1) = dz / r
c
c     Bisector method rotation matrix elements for z- and x-axes
c
      else if (polaxe(i) .eq. 'Bisector') then
c         print*,i,"Atom ",name(i)," Polaxe ",polaxe(i)
c         print*,i,"Atom ",name(i),ix," Atom X-axis ",name(ix)
c         print*,i,"Atom ",name(i),iy," Atom Y-axis ",name(iy)
c         print*,i,"Atom ",name(i),iz," Atom Z-axis ",name(iz)
         dx = x(iz) - xi
         dy = y(iz) - yi
         dz = z(iz) - zi
         r = sqrt(dx*dx + dy*dy + dz*dz)
         dx1 = dx / r
         dy1 = dy / r
         dz1 = dz / r
         dx = x(ix) - xi
         dy = y(ix) - yi
         dz = z(ix) - zi
         r = sqrt(dx*dx + dy*dy + dz*dz)
         dx2 = dx / r
         dy2 = dy / r
         dz2 = dz / r
         dx = dx1 + dx2
         dy = dy1 + dy2
         dz = dz1 + dz2
         r = sqrt(dx*dx + dy*dy + dz*dz)
         a(1,3) = dx / r
         a(2,3) = dy / r
         a(3,3) = dz / r
         dot = dx2*a(1,3) + dy2*a(2,3) + dz2*a(3,3)
         dx = dx2 - dot*a(1,3)
         dy = dy2 - dot*a(2,3)
         dz = dz2 - dot*a(3,3)
         r = sqrt(dx*dx + dy*dy + dz*dz)
         a(1,1) = dx / r
         a(2,1) = dy / r
         a(3,1) = dz / r
c
c     Z-Bisect method rotation matrix elements for z- and x-axes
c
      else if (polaxe(i) .eq. 'Z-Bisect') then
c         print*,i,"Atom ",name(i)," Polaxe ",polaxe(i)
c         print*,i,"Atom ",name(i),ix," Atom X-axis ",name(ix)
c         print*,i,"Atom ",name(i),iy," Atom Y-axis ",name(iy)
c         print*,i,"Atom ",name(i),iz," Atom Z-axis ",name(iz)
         dx = x(iz) - xi
         dy = y(iz) - yi
         dz = z(iz) - zi
         r = sqrt(dx*dx + dy*dy + dz*dz)
         a(1,3) = dx / r
         a(2,3) = dy / r
         a(3,3) = dz / r
         dx = x(ix) - xi
         dy = y(ix) - yi
         dz = z(ix) - zi
         r = sqrt(dx*dx + dy*dy + dz*dz)
         dx1 = dx / r
         dy1 = dy / r
         dz1 = dz / r
         dx = x(iy) - xi
         dy = y(iy) - yi
         dz = z(iy) - zi
         r = sqrt(dx*dx + dy*dy + dz*dz)
         dx2 = dx / r
         dy2 = dy / r
         dz2 = dz / r
         dx = dx1 + dx2
         dy = dy1 + dy2
         dz = dz1 + dz2
         r = sqrt(dx*dx + dy*dy + dz*dz)
         dx = dx / r
         dy = dy / r
         dz = dz / r
         dot = dx*a(1,3) + dy*a(2,3) + dz*a(3,3)
         dx = dx - dot*a(1,3)
         dy = dy - dot*a(2,3)
         dz = dz - dot*a(3,3)
         r = sqrt(dx*dx + dy*dy + dz*dz)
         a(1,1) = dx / r
         a(2,1) = dy / r
         a(3,1) = dz / r
c
c     3-Fold method rotation matrix elements for z- and x-axes
c
      else if (polaxe(i) .eq. '3-Fold') then
c         print*,i,"Atom ",name(i)," Polaxe ",polaxe(i)
c         print*,i,"Atom ",name(i),ix," Atom X-axis ",name(ix)
c         print*,i,"Atom ",name(i),iy," Atom Y-axis ",name(iy)
c         print*,i,"Atom ",name(i),iz," Atom Z-axis ",name(iz)
         dx = x(iz) - xi
         dy = y(iz) - yi
         dz = z(iz) - zi
         r = sqrt(dx*dx + dy*dy + dz*dz)
         dx1 = dx / r
         dy1 = dy / r
         dz1 = dz / r
         dx = x(ix) - xi
         dy = y(ix) - yi
         dz = z(ix) - zi
         r = sqrt(dx*dx + dy*dy + dz*dz)
         dx2 = dx / r
         dy2 = dy / r
         dz2 = dz / r
         dx = x(iy) - xi
         dy = y(iy) - yi
         dz = z(iy) - zi
         r = sqrt(dx*dx + dy*dy + dz*dz)
         dx3 = dx / r
         dy3 = dy / r
         dz3 = dz / r
         dx = dx1 + dx2 + dx3
         dy = dy1 + dy2 + dy3
         dz = dz1 + dz2 + dz3
         r = sqrt(dx*dx + dy*dy + dz*dz)
         a(1,3) = dx / r
         a(2,3) = dy / r
         a(3,3) = dz / r
         dot = dx2*a(1,3) + dy2*a(2,3) + dz2*a(3,3)
         dx = dx2 - dot*a(1,3)
         dy = dy2 - dot*a(2,3)
         dz = dz2 - dot*a(3,3)
         r = sqrt(dx*dx + dy*dy + dz*dz)
         a(1,1) = dx / r
         a(2,1) = dy / r
         a(3,1) = dz / r
      end if
c
c     finally, find rotation matrix elements for the y-axis
c
      a(1,2) = a(3,1)*a(2,3) - a(2,1)*a(3,3)
      a(2,2) = a(1,1)*a(3,3) - a(3,1)*a(1,3)
      a(3,2) = a(2,1)*a(1,3) - a(1,1)*a(2,3)
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine rotsite  --  rotate multipoles at single site  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "rotsite" computes the atomic multipoles at a specified site
c     in the global coordinate frame by applying a rotation matrix
c
c
      subroutine rotsite (isite,a)
      use sizes
      use atoms
      use mpole
      use kpolr_aniso
      use polar_aniso
      implicit none
      integer i,j,k,m
      integer isite
      real*8 a(3,3)
      real*8 m2(3,3)
      real*8 r2(3,3)
c
c
c     monopoles have the same value in any coordinate frame
c
      rpole(1,isite) = pole(1,isite)
c
c     rotate the dipoles to the global coordinate frame
c
      do i = 2, 4
         rpole(i,isite) = 0.0d0
         do j = 2, 4
            rpole(i,isite) = rpole(i,isite) + pole(j,isite)*a(i-1,j-1)
         end do
      end do
c
c     rotate the quadrupoles to the global coordinate frame
c
      k = 5
      do i = 1, 3
         do j = 1, 3
            m2(i,j) = pole(k,isite)
            r2(i,j) = 0.0d0
            k = k + 1
         end do
      end do
      do i = 1, 3
         do j = 1, 3
            if (j .lt. i) then
               r2(i,j) = r2(j,i)
            else
               do k = 1, 3
                  do m = 1, 3
                     r2(i,j) = r2(i,j) + a(i,k)*a(j,m)*m2(k,m)
                  end do
               end do
            end if
         end do
      end do
      k = 5
      do i = 1, 3
         do j = 1, 3
            rpole(k,isite) = r2(i,j)
            k = k + 1
         end do
      end do

c
c     rotate the anisotropic polarizabilities to the global coordinate frame
c

      if(use_anisopolz) then
c      print*,"In rotsite using aniso"
      k = 1
      do i = 1, 3
         do j = 1, 3
            m2(i,j) = anisopolarity(k,isite)
            r2(i,j) = 0.0d0
            k = k + 1
         end do
      end do
      do i = 1, 3
         do j = 1, 3
            if (j .lt. i) then
               r2(i,j) = r2(j,i)
            else
               do k = 1, 3
                  do m = 1, 3
                     r2(i,j) = r2(i,j) + a(i,k)*a(j,m)*m2(k,m)
                  end do
               end do
            end if
         end do
      end do
      k = 1
      do i = 1, 3
         do j = 1, 3
            ranisopolarity(k,isite) = r2(i,j)
            k = k + 1
         end do
      end do

      end if

      return
      end
