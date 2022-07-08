c
c
c     #############################################################
c     ##  COPYRIGHT (C) 1999 by Pengyu Ren & Jay William Ponder  ##
c     ##                   All Rights Reserved                   ##
c     #############################################################
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine induce  --  evaluate induced dipole moments  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "induce" computes the induced dipole moments at polarizable
c     sites due to direct or mutual polarization
c
c     assumes multipole components have already been rotated into
c     the global coordinate frame; computes induced dipoles based
c     on full system, use of active or inactive atoms is ignored
c
c
      subroutine induce_aniso
      use sizes
      use inform
      use iounit
      use limits
      use mpole
      use polar
      use potent
      use solute
      use units
      use uprior
      implicit none
      integer i,j,k
      real*8 norm
      logical header
c
c
c     choose the method for computation of induced dipoles
c
c      if (solvtyp(1:2) .eq. 'PB') then
c         call induce0e
c      else if (solvtyp(1:2) .eq. 'GK') then
c         call induce0d
c      else
         call induce0a_aniso
c      end if
c
c     update the lists of previous induced dipole values
c
      if (use_pred) then
         nualt = min(nualt+1,maxualt)
         do i = 1, npole
            do j = 1, 3
               do k = nualt, 2, -1
                  udalt(k,j,i) = udalt(k-1,j,i)
                  upalt(k,j,i) = upalt(k-1,j,i)
               end do
               udalt(1,j,i) = uind(j,i)
               upalt(1,j,i) = uinp(j,i)
               if (use_solv) then
                  do k = nualt, 2, -1
                     usalt(k,j,i) = usalt(k-1,j,i)
                     upsalt(k,j,i) = upsalt(k-1,j,i)
                  end do
                  usalt(1,j,i) = uinds(j,i)
                  upsalt(1,j,i) = uinps(j,i)
               end if
            end do
         end do
      end if
c
c     print out a list of the final induced dipole moments
c
      if (debug) then
         header = .true.
         do i = 1, npole
            if (polarity(i) .ne. 0.0d0) then
               if (header) then
                  header = .false.
                  if (solvtyp.eq.'GK' .or. solvtyp.eq.'PB') then
                     write (iout,10)
   10                format (/,' Vacuum Induced Dipole Moments',
     &                          ' (Debyes) :')
                  else
                     write (iout,20)
   20                format (/,' Induced Dipole Moments (Debyes) :')
                  end if
                  if (digits .ge. 8) then
                     write (iout,30)
   30                format (/,4x,'Atom',14x,'X',15x,'Y',15x,'Z',
     &                          15x,'Total',/)
                  else if (digits .ge. 6) then
                     write (iout,40)
   40                format (/,4x,'Atom',14x,'X',13x,'Y',13x,'Z',
     &                          12x,'Total',/)
                  else
                     write (iout,50)
   50                format (/,4x,'Atom',14x,'X',11x,'Y',11x,'Z',
     &                          9x,'Total',/)
                  end if
               end if
               k = ipole(i)
               norm = sqrt(uind(1,i)**2+uind(2,i)**2+uind(3,i)**2)
               if (digits .ge. 8) then
                  write (iout,60)  k,(debye*uind(j,i),j=1,3),debye*norm
   60             format (i8,3x,4f16.8)
               else if (digits .ge. 6) then
                  write (iout,70)  k,(debye*uind(j,i),j=1,3),debye*norm
   70             format (i8,4x,4f14.6)
               else
                  write (iout,80)  k,(debye*uind(j,i),j=1,3),debye*norm
   80             format (i8,5x,4f12.4)
               end if
            end if
         end do
c         header = .true.
c         if (solvtyp.eq.'GK' .or. solvtyp.eq.'PB') then
c            do i = 1, npole
c               if (polarity(i) .ne. 0.0d0) then
c                  if (header) then
c                     header = .false.
c                     write (iout,90)
c   90                format (/,' SCRF Induced Dipole Moments',
c     &                          ' (Debyes) :')
c                     if (digits .ge. 8) then
c                        write (iout,100)
c  100                   format (/,4x,'Atom',14x,'X',15x,'Y',15x,'Z',
c     &                             15x,'Total',/)
c                     else if (digits .ge. 6) then
c                        write (iout,110)
c  110                   format (/,4x,'Atom',14x,'X',13x,'Y',13x,'Z',
c     &                             12x,'Total',/)
c                     else
c                        write (iout,120)
c  120                   format (/,4x,'Atom',14x,'X',11x,'Y',11x,'Z',
c     &                             9x,'Total',/)
c                     end if
c                  end if
c                  k = ipole(i)
c                  norm = sqrt(uinds(1,i)**2+uinds(2,i)**2+uinds(3,i)**2)
c                  if (digits .ge. 8) then
c                     write (iout,130)  k,(debye*uinds(j,i),j=1,3),
c     &                                 debye*norm
c  130                format (i8,3x,4f16.8)
c                  else if (digits .ge. 6) then
c                     write (iout,140)  k,(debye*uinds(j,i),j=1,3),
c     &                                 debye*norm
c  140                format (i8,4x,4f14.6)
c                  else
c                     write (iout,150)  k,(debye*uinds(j,i),j=1,3),
c     &                                 debye*norm
c  150                format (i8,5x,4f12.4)
c                  end if
c               end if
c            end do
c         end if
c      end if
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine induce0a_aniso  --  conjugate gradient dipole solver  ##
c     ##  for use with anisotropic polarizabilities     
c     ##                                                             ##
c     #################################################################
c
c
c     "induce0a_aniso" computes the induced dipole moments at polarizable
c     sites using a preconditioned conjugate gradient solver
c     Here, anisotropic polarizabilities calculated with the
c     Williams-Stone-Misquitta Method are used 
c     At present, neither Ewald nor neighbor lists are enabled.
c
      subroutine induce0a_aniso
      use sizes
      use atoms
      use inform
      use iounit
      use limits
      use mpole
      use polar
      use polpot
      use potent
      use units
      use uprior
      use polar_aniso
      implicit none
      integer i,j,k,iter
      integer maxiter
      real*8 polmin
      real*8 eps,epsold
      real*8 epsd,epsp
      real*8 udsum,upsum
      real*8 a,ap,b,bp
      real*8 sum,sump
      real*8, allocatable :: anisopoli(:,:)
      real*8, allocatable :: field(:,:)
      real*8, allocatable :: fieldp(:,:)
      real*8, allocatable :: udir(:,:)
      real*8, allocatable :: udirp(:,:)
      real*8, allocatable :: rsd(:,:)
      real*8, allocatable :: rsdp(:,:)
      real*8, allocatable :: zrsd(:,:)
      real*8, allocatable :: zrsdp(:,:)
      real*8, allocatable :: conj(:,:)
      real*8, allocatable :: conjp(:,:)
      real*8, allocatable :: vec(:,:)
      real*8, allocatable :: vecp(:,:)
      logical done
      character*6 mode
c
c
c     zero out the induced dipoles at each site
c
      do i = 1, npole
         do j = 1, 3
            uind(j,i) = 0.0d0
            uinp(j,i) = 0.0d0
         end do
      end do
      if (.not. use_polar)  return
c
c     perform dynamic allocation of some local arrays
c
      allocate (field(3,npole))
      allocate (fieldp(3,npole))
      allocate (udir(3,npole))
      allocate (udirp(3,npole))
c
c     get the electrostatic field due to permanent multipoles
c
      if (use_ewald) then
             write (iout,11)
   11      format (/,' Anisotropic polarizability is not yet available',
     &              ' for Particle Mesh Ewald')
             call fatal
      else if (use_mlist) then
             write (iout,21)
   21      format (/,' Anisotropic polarizability is not yet available',
     &              ' for Neighbor Lists without Ewald')
             call fatal
      else
         call dfield0a (field,fieldp)
      end if
c
c     set induced dipoles to polarizability times direct field
c
      do i = 1, npole
            udir(1,i) = ranisopolarity(1,i)*field(1,i) +
     &      ranisopolarity(2,i)*field(2,i) 
     &      + ranisopolarity(3,i)*field(3,i)

            udir(2,i) = ranisopolarity(4,i)*field(1,i) + 
     &      ranisopolarity(5,i)*field(2,i) 
     &      + ranisopolarity(6,i)*field(3,i)

            udir(3,i) = ranisopolarity(7,i)*field(1,i) +
     &      ranisopolarity(8,i)*field(2,i)
     &      + ranisopolarity(9,i)*field(3,i)

            udirp(1,i) = ranisopolarity(1,i)*fieldp(1,i) + 
     &      ranisopolarity(2,i)*fieldp(2,i) 
     &      + ranisopolarity(3,i)*fieldp(3,i)
            
            udirp(2,i) = ranisopolarity(4,i)*fieldp(1,i) +
     &      ranisopolarity(5,i)*fieldp(2,i)
     &      + ranisopolarity(6,i)*fieldp(3,i)

            udirp(3,i) = ranisopolarity(7,i)*fieldp(1,i) +
     &      ranisopolarity(8,i)*fieldp(2,i)
     &      + ranisopolarity(9,i)*fieldp(3,i)

            do j=1,3
            uind(j,i) = udir(j,i)
            uinp(j,i) = udirp(j,i)
            end do
      end do
c
c     set tolerances for computation of mutual induced dipoles
c
      if (poltyp .eq. 'MUTUAL') then
         done = .false.
         maxiter = 500
         iter = 0
         polmin = 0.00000001d0
         eps = 100.0d0
c
c     estimated induced dipoles from polynomial predictor
c
         if (use_pred .and. nualt.eq.maxualt) then
            call ulspred
            do i = 1, npole
               do j = 1, 3
                  udsum = 0.0d0
                  upsum = 0.0d0
                  do k = 1, nualt-1
                     udsum = udsum + bpred(k)*udalt(k,j,i)
                     upsum = upsum + bpredp(k)*upalt(k,j,i)
                  end do
                  uind(j,i) = udsum
                  uinp(j,i) = upsum
               end do
            end do
         end if
c
c     perform dynamic allocation of some local arrays
c
         allocate (anisopoli(9,npole))
         allocate (rsd(3,npole))
         allocate (rsdp(3,npole))
         allocate (zrsd(3,npole))
         allocate (zrsdp(3,npole))
         allocate (conj(3,npole))
         allocate (conjp(3,npole))
         allocate (vec(3,npole))
         allocate (vecp(3,npole))
c
c     get the electrostatic field due to induced dipoles
c
         if (use_ewald) then
             write (iout,31)
   31      format (/,' Anisotropic polarizability is not yet available',
     &              ' for Particle Mesh Ewald')
             call fatal

         else if (use_mlist) then
             write (iout,41)
   41      format (/,' Anisotropic polarizability is not yet available',
     &              ' for Neighbor Lists without Ewald')
             call fatal
         else
            call ufield0a (field,fieldp)
         end if

c LEFT OFF HERE A;SLDKFJ;ADKSJ

c
c     set initial conjugate gradient residual and conjugate vector
c
         do i = 1, npole
            do j=1,9
            anisopoli(j,i) = max(polmin,ranisopolarity(j,i))
            end do
            do j = 1, 3
               rsd(j,i) = (udir(j,i)-uind(j,i))/poli(i)
     &                       + field(j,i)
               rsdp(j,i) = (udirp(j,i)-uinp(j,i))/poli(i)
     &                       + fieldp(j,i)
            end do
         end do
         mode = 'BUILD'
         if (use_mlist) then
             write (iout,51)
   51      format (/,' Anisotropic polarizability is not yet available',
     &              ' for Neighbor Lists without Ewald')
             call fatal
         else
            call uscale0a (mode,rsd,rsdp,zrsd,zrsdp)
            mode = 'APPLY'
            call uscale0a (mode,rsd,rsdp,zrsd,zrsdp)
         end if
         do i = 1, npole
            do j = 1, 3
               conj(j,i) = zrsd(j,i)
               conjp(j,i) = zrsdp(j,i)
            end do
         end do
c
c     conjugate gradient iteration of the mutual induced dipoles
c
         do while (.not. done)
            iter = iter + 1
            do i = 1, npole
               do j = 1, 3
                  vec(j,i) = uind(j,i)
                  vecp(j,i) = uinp(j,i)
                  uind(j,i) = conj(j,i)
                  uinp(j,i) = conjp(j,i)
               end do
            end do
            if (use_ewald) then
             write (iout,61)
   61      format (/,' Anisotropic polarizability is not yet available',
     &              ' for Particle Mesh Ewald')
             call fatal
            else if (use_mlist) then
             write (iout,71)
   71      format (/,' Anisotropic polarizability is not yet available',
     &              ' for Neighbor Lists without Ewald')
             call fatal
            else
               call ufield0a (field,fieldp)
            end if
            do i = 1, npole
               do j = 1, 3
                  uind(j,i) = vec(j,i)
                  uinp(j,i) = vecp(j,i)
                  vec(j,i) = conj(j,i)/poli(i) - field(j,i)
                  vecp(j,i) = conjp(j,i)/poli(i) - fieldp(j,i)
               end do
            end do
            a = 0.0d0
            ap = 0.0d0
            sum = 0.0d0
            sump = 0.0d0
            do i = 1, npole
               do j = 1, 3
                  a = a + conj(j,i)*vec(j,i)
                  ap = ap + conjp(j,i)*vecp(j,i)
                  sum = sum + rsd(j,i)*zrsd(j,i)
                  sump = sump + rsdp(j,i)*zrsdp(j,i)
               end do
            end do
            if (a .ne. 0.0d0)  a = sum / a
            if (ap .ne. 0.0d0)  ap = sump / ap
            do i = 1, npole
               do j = 1, 3
                  uind(j,i) = uind(j,i) + a*conj(j,i)
                  uinp(j,i) = uinp(j,i) + ap*conjp(j,i)
                  rsd(j,i) = rsd(j,i) - a*vec(j,i)
                  rsdp(j,i) = rsdp(j,i) - ap*vecp(j,i)
               end do
            end do
            if (use_mlist) then
             write (iout,81)
   81      format (/,' Anisotropic polarizability is not yet available',
     &              ' for Neighbor Lists without Ewald')
             call fatal
            else
               call uscale0a (mode,rsd,rsdp,zrsd,zrsdp)
            end if
            b = 0.0d0
            bp = 0.0d0
            do i = 1, npole
               do j = 1, 3
                  b = b + rsd(j,i)*zrsd(j,i)
                  bp = bp + rsdp(j,i)*zrsdp(j,i)
               end do
            end do
            if (sum .ne. 0.0d0)  b = b / sum
            if (sump .ne. 0.0d0)  bp = bp / sump
            epsd = 0.0d0
            epsp = 0.0d0
            do i = 1, npole
               do j = 1, 3
                  conj(j,i) = zrsd(j,i) + b*conj(j,i)
                  conjp(j,i) = zrsdp(j,i) + bp*conjp(j,i)
                  epsd = epsd + rsd(j,i)*rsd(j,i)
                  epsp = epsp + rsdp(j,i)*rsdp(j,i)
               end do
            end do
c
c     check the convergence of the mutual induced dipoles
c
            epsold = eps
            eps = max(epsd,epsp)
            eps = debye * sqrt(eps/dble(npolar))
            if (debug) then
               if (iter .eq. 1) then
                  write (iout,10)
   10             format (/,' Determination of Induced Dipole',
     &                       ' Moments :',
     &                    //,4x,'Iter',8x,'RMS Change (Debyes)',/)
               end if
               write (iout,20)  iter,eps
   20          format (i8,7x,f16.10)
            end if
            if (eps .lt. poleps)  done = .true.
            if (eps .gt. epsold)  done = .true.
            if (iter .ge. politer)  done = .true.
         end do
c
c     perform deallocation of some local arrays
c
         deallocate (poli)
         deallocate (rsd)
         deallocate (rsdp)
         deallocate (zrsd)
         deallocate (zrsdp)
         deallocate (conj)
         deallocate (conjp)
         deallocate (vec)
         deallocate (vecp)
c
c     print the results from the conjugate gradient iteration
c
         if (debug) then
            write (iout,30)  iter,eps
   30       format (/,' Induced Dipoles :',6x,'Iterations',i5,
     &                 6x,'RMS Change',f15.10)
         end if
c
c     terminate the calculation if dipoles failed to converge
c
         if (iter.ge.maxiter .or. eps.gt.epsold) then
            write (iout,40)
   40       format (/,' INDUCE  --  Warning, Induced Dipoles',
     &                 ' are not Converged')
            call prterr
            call fatal
         end if
      end if
c
c     perform deallocation of some local arrays
c
      deallocate (field)
      deallocate (fieldp)
      deallocate (udir)
      deallocate (udirp)
      return
      end
c
c
c
c
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine uscale0a  --  dipole preconditioner via loop  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "uscale0a" builds and applies a preconditioner for the conjugate
c     gradient induced dipole solver using a double loop
c
c
      subroutine uscale0a_aniso (mode,rsd,rsdp,zrsd,zrsdp)
      use sizes
      use atoms
      use limits
      use mpole
      use polar
      use polgrp
      use polpot
      use usolve
      implicit none
      integer i,j,k,m
      integer ii,kk
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 r,r2,rr3,rr5
      real*8 pdi,pti
      real*8 polmin
      real*8 poli,polik
      real*8 damp,expdamp
      real*8 pgamma,off2
      real*8 scale3,scale5
      real*8 m1,m2,m3
      real*8 m4,m5,m6
      real*8, allocatable :: dscale(:)
      real*8 rsd(3,*)
      real*8 rsdp(3,*)
      real*8 zrsd(3,*)
      real*8 zrsdp(3,*)
      character*6 mode
c
c
c     apply the preconditioning matrix to the current residual
c
      if (mode .eq. 'APPLY') then
c
c     use diagonal preconditioner elements as first approximation
c
         polmin = 0.00000001d0
         do i = 1, npole
            poli = udiag * max(polmin,polarity(i))
            do j = 1, 3
               zrsd(j,i) = poli * rsd(j,i)
               zrsdp(j,i) = poli * rsdp(j,i)
            end do
         end do
c
c     use the off-diagonal preconditioner elements in second phase
c
         off2 = usolvcut * usolvcut
         j = 0
         do i = 1, npole-1
            ii = ipole(i)
            do k = i+1, npole
               kk = ipole(k)
               xr = x(kk) - x(ii)
               yr = y(kk) - y(ii)
               zr = z(kk) - z(ii)
               call image (xr,yr,zr)
               r2 = xr*xr + yr* yr + zr*zr
               if (r2 .le. off2) then
                  m1 = minv(j+1)
                  m2 = minv(j+2)
                  m3 = minv(j+3)
                  m4 = minv(j+4)
                  m5 = minv(j+5)
                  m6 = minv(j+6)
                  j = j + 6
                  zrsd(1,i) = zrsd(1,i) + m1*rsd(1,k) + m2*rsd(2,k)
     &                           + m3*rsd(3,k)
                  zrsd(2,i) = zrsd(2,i) + m2*rsd(1,k) + m4*rsd(2,k)
     &                           + m5*rsd(3,k)
                  zrsd(3,i) = zrsd(3,i) + m3*rsd(1,k) + m5*rsd(2,k)
     &                           + m6*rsd(3,k)
                  zrsd(1,k) = zrsd(1,k) + m1*rsd(1,i) + m2*rsd(2,i)
     &                           + m3*rsd(3,i)
                  zrsd(2,k) = zrsd(2,k) + m2*rsd(1,i) + m4*rsd(2,i)
     &                           + m5*rsd(3,i)
                  zrsd(3,k) = zrsd(3,k) + m3*rsd(1,i) + m5*rsd(2,i)
     &                           + m6*rsd(3,i)
                  zrsdp(1,i) = zrsdp(1,i) + m1*rsdp(1,k) + m2*rsdp(2,k)
     &                            + m3*rsdp(3,k)
                  zrsdp(2,i) = zrsdp(2,i) + m2*rsdp(1,k) + m4*rsdp(2,k)
     &                            + m5*rsdp(3,k)
                  zrsdp(3,i) = zrsdp(3,i) + m3*rsdp(1,k) + m5*rsdp(2,k)
     &                            + m6*rsdp(3,k)
                  zrsdp(1,k) = zrsdp(1,k) + m1*rsdp(1,i) + m2*rsdp(2,i)
     &                            + m3*rsdp(3,i)
                  zrsdp(2,k) = zrsdp(2,k) + m2*rsdp(1,i) + m4*rsdp(2,i)
     &                            + m5*rsdp(3,i)
                  zrsdp(3,k) = zrsdp(3,k) + m3*rsdp(1,i) + m5*rsdp(2,i)
     &                            + m6*rsdp(3,i)
               end if
            end do
         end do
c
c     construct off-diagonal elements of preconditioning matrix
c
      else if (mode .eq. 'BUILD') then
c
c     perform dynamic allocation of some local arrays
c
         allocate (dscale(n))
c
c     set array needed to scale connected atom interactions
c
         do i = 1, n
            dscale(i) = 1.0d0
         end do
c
c     determine the off-diagonal elements of the preconditioner
c
         off2 = usolvcut * usolvcut
         m = 0
         do i = 1, npole-1
            ii = ipole(i)
            xi = x(ii)
            yi = y(ii)
            zi = z(ii)
            pdi = pdamp(i)
            pti = thole(i)
            poli = polarity(i)
            do j = i+1, npole
               dscale(ipole(j)) = 1.0d0
            end do
            do j = 1, np11(ii)
               dscale(ip11(j,ii)) = u1scale
            end do
            do j = 1, np12(ii)
               dscale(ip12(j,ii)) = u2scale
            end do
            do j = 1, np13(ii)
               dscale(ip13(j,ii)) = u3scale
            end do
            do j = 1, np14(ii)
               dscale(ip14(j,ii)) = u4scale
            end do
            do k = i+1, npole
               kk = ipole(k)
               xr = x(kk) - xi
               yr = y(kk) - yi
               zr = z(kk) - zi
               call image (xr,yr,zr)
               r2 = xr*xr + yr* yr + zr*zr
               if (r2 .le. off2) then
                  r = sqrt(r2)
                  scale3 = dscale(kk)
                  scale5 = dscale(kk)
                  damp = pdi * pdamp(k)
                  if (damp .ne. 0.0d0) then
                     pgamma = min(pti,thole(k))
                     damp = -pgamma * (r/damp)**3
                     if (damp .gt. -50.0d0) then
                        expdamp = exp(damp)
                        scale3 = scale3 * (1.0d0-expdamp)
                        scale5 = scale5 * (1.0d0-expdamp*(1.0d0-damp))
                     end if
                  end if
                  polik = poli * polarity(k)
                  rr3 = scale3 * polik / (r*r2)
                  rr5 = 3.0d0 * scale5 * polik / (r*r2*r2)
                  minv(m+1) = rr5*xr*xr - rr3
                  minv(m+2) = rr5*xr*yr
                  minv(m+3) = rr5*xr*zr
                  minv(m+4) = rr5*yr*yr - rr3
                  minv(m+5) = rr5*yr*zr
                  minv(m+6) = rr5*zr*zr - rr3
                  m = m + 6
               end if
            end do
c
c     reset interaction scaling coefficients for connected atoms
c
            do j = 1, np11(ii)
               dscale(ip11(j,ii)) = 1.0d0
            end do
            do j = 1, np12(ii)
               dscale(ip12(j,ii)) = 1.0d0
            end do
            do j = 1, np13(ii)
               dscale(ip13(j,ii)) = 1.0d0
            end do
            do j = 1, np14(ii)
               dscale(ip14(j,ii)) = 1.0d0
            end do
         end do
c
c     perform deallocation of some local arrays
c
         deallocate (dscale)
      end if
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
