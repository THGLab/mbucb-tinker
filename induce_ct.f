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
      subroutine induce_ct
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
         call induce0a_ct
c
c     update the lists of previous induced dipole values
c
c      if (use_pred) then
c         nualt = min(nualt+1,maxualt)
c         do i = 1, npole
c            do j = 1, 3
c               do k = nualt, 2, -1
c                  udalt(k,j,i) = udalt(k-1,j,i)
c                  upalt(k,j,i) = upalt(k-1,j,i)
c               end do
c               udalt(1,j,i) = uind(j,i)
c               upalt(1,j,i) = uinp(j,i)
c               if (use_solv) then
c                  do k = nualt, 2, -1
c                     usalt(k,j,i) = usalt(k-1,j,i)
c                     upsalt(k,j,i) = upsalt(k-1,j,i)
c                  end do
c                  usalt(1,j,i) = uinds(j,i)
c                  upsalt(1,j,i) = uinps(j,i)
c               end if
c            end do
c         end do
c      end if
c
c     print out a list of the final induced dipole moments
c
c      if (debug) then
c         header = .true.
c         do i = 1, npole
c            if (polarity(i) .ne. 0.0d0) then
c               if (header) then
c                  header = .false.
c                  if (solvtyp.eq.'GK' .or. solvtyp.eq.'PB') then
c                     write (iout,10)
c   10                format (/,' Vacuum Induced Dipole Moments',
c     &                          ' (Debyes) :')
c                  else
c                     write (iout,20)
c   20                format (/,' Induced Dipole Moments (Debyes) :')
c                  end if
c                  if (digits .ge. 8) then
c                     write (iout,30)
c   30                format (/,4x,'Atom',14x,'X',15x,'Y',15x,'Z',
c     &                          15x,'Total',/)
c                  else if (digits .ge. 6) then
c                     write (iout,40)
c   40                format (/,4x,'Atom',14x,'X',13x,'Y',13x,'Z',
c     &                          12x,'Total',/)
c                  else
c                     write (iout,50)
c   50                format (/,4x,'Atom',14x,'X',11x,'Y',11x,'Z',
c     &                          9x,'Total',/)
c                  end if
c               end if
c               k = ipole(i)
c               norm = sqrt(uind(1,i)**2+uind(2,i)**2+uind(3,i)**2)
c               if (digits .ge. 8) then
c                  write (iout,60)  k,(debye*uind(j,i),j=1,3),debye*norm
c   60             format (i8,3x,4f16.8)
c               else if (digits .ge. 6) then
c                  write (iout,70)  k,(debye*uind(j,i),j=1,3),debye*norm
c   70             format (i8,4x,4f14.6)
c               else
c                  write (iout,80)  k,(debye*uind(j,i),j=1,3),debye*norm
c   80             format (i8,5x,4f12.4)
c               end if
c            end if
c         end do
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
c      return
      end
cc
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine induce0a  --  conjugate gradient dipole solver  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "induce0a" computes the induced dipole moments at polarizable
c     sites using a preconditioned conjugate gradient solver
c
c
      subroutine induce0a_ct
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
      use chargetransfer
      implicit none
      integer i,j,k,iter
      integer maxiter
      real*8 polmin
      real*8 eps,epsold
      real*8 epsd,epsp
      real*8 udsum,upsum
      real*8 a,ap,b,bp
      real*8 sum,sump
      real*8, allocatable :: poli(:)
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
c      if (use_ewald) then
c         call dfield0c (field,fieldp)
c      else if (use_mlist) then
c         call dfield0b (field,fieldp)
c      else
         call dfield0a_ct (field,fieldp)
c      end if
c
c     set induced dipoles to polarizability times direct field
c
      do i = 1, npole
         do j = 1, 3
            udir(j,i) = polarity(i) * field(j,i)
            udirp(j,i) = polarity(i) * fieldp(j,i)
            uind(j,i) = udir(j,i)
            uinp(j,i) = udirp(j,i) 
         end do
      end do
c
c     set tolerances for computation of mutual induced dipoles
c
c      if (poltyp .eq. 'MUTUAL') then
c         done = .false.
c         maxiter = 500
c         iter = 0
c         polmin = 0.00000001d0
c         eps = 100.0d0
c
c     estimated induced dipoles from polynomial predictor
c
c         if (use_pred .and. nualt.eq.maxualt) then
c            call ulspred
c            do i = 1, npole
c               do j = 1, 3
c                  udsum = 0.0d0
c                  upsum = 0.0d0
c                  do k = 1, nualt-1
c                     udsum = udsum + bpred(k)*udalt(k,j,i)
c                     upsum = upsum + bpredp(k)*upalt(k,j,i)
c                  end do
c                  uind(j,i) = udsum
c                  uinp(j,i) = upsum
c               end do
c            end do
c         end if
c
c     perform dynamic allocation of some local arrays
c
         allocate (poli(npole))
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
c         if (use_ewald) then
c            call ufield0c (field,fieldp)
c         else if (use_mlist) then
c            call ufield0b (field,fieldp)
c         else
c            call ufield0a (field,fieldp)
c         end if
c
c     set initial conjugate gradient residual and conjugate vector
c
c         do i = 1, npole
c            poli(i) = max(polmin,polarity(i))
c            do j = 1, 3
c               rsd(j,i) = (udir(j,i)-uind(j,i))/poli(i)
c     &                       + field(j,i)
c               rsdp(j,i) = (udirp(j,i)-uinp(j,i))/poli(i)
c     &                       + fieldp(j,i)
c            end do
c         end do
c         mode = 'BUILD'
c         if (use_mlist) then
c            call uscale0b (mode,rsd,rsdp,zrsd,zrsdp)
c            mode = 'APPLY'
c            call uscale0b (mode,rsd,rsdp,zrsd,zrsdp)
c         else
c            call uscale0a (mode,rsd,rsdp,zrsd,zrsdp)
c            mode = 'APPLY'
c            call uscale0a (mode,rsd,rsdp,zrsd,zrsdp)
c         end if
c         do i = 1, npole
c            do j = 1, 3
c               conj(j,i) = zrsd(j,i)
c               conjp(j,i) = zrsdp(j,i)
c            end do
c         end do
cc
cc     conjugate gradient iteration of the mutual induced dipoles
cc
c         do while (.not. done)
c            iter = iter + 1
c            do i = 1, npole
c               do j = 1, 3
c                  vec(j,i) = uind(j,i)
c                  vecp(j,i) = uinp(j,i)
c                  uind(j,i) = conj(j,i)
c                  uinp(j,i) = conjp(j,i)
c               end do
c            end do
c            if (use_ewald) then
c               call ufield0c (field,fieldp)
c            else if (use_mlist) then
c               call ufield0b (field,fieldp)
c            else
c               call ufield0a (field,fieldp)
c            end if
c            do i = 1, npole
c               do j = 1, 3
c                  uind(j,i) = vec(j,i)
c                  uinp(j,i) = vecp(j,i)
c                  vec(j,i) = conj(j,i)/poli(i) - field(j,i)
c                  vecp(j,i) = conjp(j,i)/poli(i) - fieldp(j,i)
c               end do
c            end do
c            a = 0.0d0
c            ap = 0.0d0
c            sum = 0.0d0
c            sump = 0.0d0
c            do i = 1, npole
c               do j = 1, 3
c                  a = a + conj(j,i)*vec(j,i)
c                  ap = ap + conjp(j,i)*vecp(j,i)
c                  sum = sum + rsd(j,i)*zrsd(j,i)
c                  sump = sump + rsdp(j,i)*zrsdp(j,i)
c               end do
c            end do
c            if (a .ne. 0.0d0)  a = sum / a
c            if (ap .ne. 0.0d0)  ap = sump / ap
c            do i = 1, npole
c               do j = 1, 3
c                  uind(j,i) = uind(j,i) + a*conj(j,i)
c                  uinp(j,i) = uinp(j,i) + ap*conjp(j,i)
c                  rsd(j,i) = rsd(j,i) - a*vec(j,i)
c                  rsdp(j,i) = rsdp(j,i) - ap*vecp(j,i)
c               end do
c            end do
c            if (use_mlist) then
c               call uscale0b (mode,rsd,rsdp,zrsd,zrsdp)
c            else
c               call uscale0a (mode,rsd,rsdp,zrsd,zrsdp)
c            end if
c            b = 0.0d0
c            bp = 0.0d0
c            do i = 1, npole
c               do j = 1, 3
c                  b = b + rsd(j,i)*zrsd(j,i)
c                  bp = bp + rsdp(j,i)*zrsdp(j,i)
c               end do
c            end do
c            if (sum .ne. 0.0d0)  b = b / sum
c            if (sump .ne. 0.0d0)  bp = bp / sump
c            epsd = 0.0d0
c            epsp = 0.0d0
c            do i = 1, npole
c               do j = 1, 3
c                  conj(j,i) = zrsd(j,i) + b*conj(j,i)
c                  conjp(j,i) = zrsdp(j,i) + bp*conjp(j,i)
c                  epsd = epsd + rsd(j,i)*rsd(j,i)
c                  epsp = epsp + rsdp(j,i)*rsdp(j,i)
c               end do
c            end do
cc
cc     check the convergence of the mutual induced dipoles
cc
c            epsold = eps
c            eps = max(epsd,epsp)
c            eps = debye * sqrt(eps/dble(npolar))
c            if (debug) then
c               if (iter .eq. 1) then
c                  write (iout,10)
c   10             format (/,' Determination of Induced Dipole',
c     &                       ' Moments :',
c     &                    //,4x,'Iter',8x,'RMS Change (Debyes)',/)
c               end if
c               write (iout,20)  iter,eps
c   20          format (i8,7x,f16.10)
c            end if
c            if (eps .lt. poleps)  done = .true.
c            if (eps .gt. epsold)  done = .true.
c            if (iter .ge. politer)  done = .true.
c         end do
cc
cc     perform deallocation of some local arrays
cc
c         deallocate (poli)
c         deallocate (rsd)
c         deallocate (rsdp)
c         deallocate (zrsd)
c         deallocate (zrsdp)
c         deallocate (conj)
c         deallocate (conjp)
c         deallocate (vec)
c         deallocate (vecp)
cc
cc     print the results from the conjugate gradient iteration
cc
c         if (debug) then
c            write (iout,30)  iter,eps
c   30       format (/,' Induced Dipoles :',6x,'Iterations',i5,
c     &                 6x,'RMS Change',f15.10)
c         end if
cc
cc     terminate the calculation if dipoles failed to converge
cc
c         if (iter.ge.maxiter .or. eps.gt.epsold) then
c            write (iout,40)
c   40       format (/,' INDUCE  --  Warning, Induced Dipoles',
c     &                 ' are not Converged')
c            call prterr
c            call fatal
c         end if
c      end if
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
c     #################################################################
c     ##                                                             ##
c     ##  subroutine dfield0a  --  direct induction via double loop  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "dfield0a" computes the direct electrostatic field due to
c     permanent multipole moments via a double loop
c
c
      subroutine dfield0a_ct (field,fieldp)
      use sizes
      use atoms
      use bound
      use cell
      use couple
      use group
      use mpole
      use polar
      use polgrp
      use polpot
      use shunt
      use chargetransfer
      implicit none
      integer i,j,k,m
      integer ii,kk
      real*8 xr,yr,zr
      real*8 fgrp,r,r2
      real*8 rr3,rr5,rr7
      real*8 ci,dix,diy,diz
      real*8 qixx,qixy,qixz
      real*8 qiyy,qiyz,qizz
      real*8 ck,dkx,dky,dkz
      real*8 qkxx,qkxy,qkxz
      real*8 qkyy,qkyz,qkzz
      real*8 dir,dkr
      real*8 qix,qiy,qiz,qir
      real*8 qkx,qky,qkz,qkr
      real*8 damp,expdamp
      real*8 scale3,scale5
      real*8 scale7
      real*8 pdi,pti,pgamma
      real*8 fid(3),fkd(3)
      real*8 fip(3),fkp(3)
      real*8, allocatable :: dscale(:)
      real*8, allocatable :: pscale(:)
      real*8 field(3,*)
      real*8 fieldp(3,*)
      logical proceed
      character*6 mode
c
c
c     zero out the value of the field at each site
c
      do i = 1, npole
         do j = 1, 3
            field(j,i) = 0.0d0
            fieldp(j,i) = 0.0d0
         end do
      end do
c
c     set the switching function coefficients
c
      mode = 'MPOLE'
      call switch (mode)
c
c     perform dynamic allocation of some local arrays
c
      allocate (dscale(n))
      allocate (pscale(n))
c
c     find the electrostatic field due to permanent multipoles
c
      do i = 1, npole-1
         ii = ipole(i)
         pdi = pdamp(i)
c         pti = thole(i)
         pti = charge_trans(i)
         ci = rpole(1,i)
         dix = rpole(2,i)
         diy = rpole(3,i)
         diz = rpole(4,i)
         qixx = rpole(5,i)
         qixy = rpole(6,i)
         qixz = rpole(7,i)
         qiyy = rpole(9,i)
         qiyz = rpole(10,i)
         qizz = rpole(13,i)
         do j = i+1, npole
            dscale(ipole(j)) = 1.0d0
            pscale(ipole(j)) = 1.0d0
         end do
         do j = 1, n12(ii)
            pscale(i12(j,ii)) = p2scale
         end do
         do j = 1, n13(ii)
            pscale(i13(j,ii)) = p3scale
         end do
         do j = 1, n14(ii)
            pscale(i14(j,ii)) = p4scale
            do k = 1, np11(ii)
               if (i14(j,ii) .eq. ip11(k,ii))
     &            pscale(i14(j,ii)) = p4scale * p41scale
            end do
         end do
         do j = 1, n15(ii)
            pscale(i15(j,ii)) = p5scale
         end do
         do j = 1, np11(ii)
            dscale(ip11(j,ii)) = d1scale
         end do
         do j = 1, np12(ii)
            dscale(ip12(j,ii)) = d2scale
         end do
         do j = 1, np13(ii)
            dscale(ip13(j,ii)) = d3scale
         end do
         do j = 1, np14(ii)
            dscale(ip14(j,ii)) = d4scale
         end do
         do k = i+1, npole
            kk = ipole(k)
            proceed = .true.
            if (use_intra)  call groups (proceed,fgrp,ii,kk,0,0,0,0)
            if (proceed) then
               xr = x(kk) - x(ii)
               yr = y(kk) - y(ii)
               zr = z(kk) - z(ii)
               if (use_bounds)  call image (xr,yr,zr)
               r2 = xr*xr + yr* yr + zr*zr
               if (r2 .le. off2) then
                  r = sqrt(r2)
                  ck = rpole(1,k)
                  dkx = rpole(2,k)
                  dky = rpole(3,k)
                  dkz = rpole(4,k)
                  qkxx = rpole(5,k)
                  qkxy = rpole(6,k)
                  qkxz = rpole(7,k)
                  qkyy = rpole(9,k)
                  qkyz = rpole(10,k)
                  qkzz = rpole(13,k)
                  scale3 = 1.0d0
                  scale5 = 1.0d0
                  scale7 = 1.0d0
                  damp = pdi * pdamp(k)
                  if (damp .ne. 0.0d0) then
                     pgamma = min(pti,charge_trans(k))
                     damp = -pgamma * (r/damp)**3
                     if (damp .gt. -50.0d0) then
                        expdamp = exp(damp)
                        scale3 = expdamp
c                        scale5 = 1.0d0 - expdamp
                        scale5 = (1.0d0 -damp)*expdamp
c                        scale7 = 1.0d0 - expdamp*(1.0d0-damp)
                           scale7 = (1.0d0-damp+0.6d0*damp**2)
     &                                             *expdamp
                     end if
                  end if
                  rr3 = scale3 / (r*r2)
                  rr5 = 3.0d0 * scale5 / (r*r2*r2)
                  rr7 = 15.0d0 * scale7 / (r*r2*r2*r2)
                  dir = dix*xr + diy*yr + diz*zr
                  qix = qixx*xr + qixy*yr + qixz*zr
                  qiy = qixy*xr + qiyy*yr + qiyz*zr
                  qiz = qixz*xr + qiyz*yr + qizz*zr
                  qir = qix*xr + qiy*yr + qiz*zr
                  dkr = dkx*xr + dky*yr + dkz*zr
                  qkx = qkxx*xr + qkxy*yr + qkxz*zr
                  qky = qkxy*xr + qkyy*yr + qkyz*zr
                  qkz = qkxz*xr + qkyz*yr + qkzz*zr
                  qkr = qkx*xr + qky*yr + qkz*zr
                  fid(1) = -xr*(rr3*ck-rr5*dkr+rr7*qkr)
     &                        - rr3*dkx + 2.0d0*rr5*qkx
                  fid(2) = -yr*(rr3*ck-rr5*dkr+rr7*qkr)
     &                        - rr3*dky + 2.0d0*rr5*qky
                  fid(3) = -zr*(rr3*ck-rr5*dkr+rr7*qkr)
     &                        - rr3*dkz + 2.0d0*rr5*qkz
                  fkd(1) = xr*(rr3*ci+rr5*dir+rr7*qir)
     &                        - rr3*dix - 2.0d0*rr5*qix
                  fkd(2) = yr*(rr3*ci+rr5*dir+rr7*qir)
     &                        - rr3*diy - 2.0d0*rr5*qiy
                  fkd(3) = zr*(rr3*ci+rr5*dir+rr7*qir)
     &                        - rr3*diz - 2.0d0*rr5*qiz
                  do j = 1, 3
                     field(j,i) = field(j,i) + fid(j)*dscale(kk)
                     field(j,k) = field(j,k) + fkd(j)*dscale(kk)
                     fieldp(j,i) = fieldp(j,i) + fid(j)*pscale(kk)
                     fieldp(j,k) = fieldp(j,k) + fkd(j)*pscale(kk)
                  end do
               end if
            end if
         end do
      end do
c
c     periodic boundary for large cutoffs via replicates method
c
      if (use_replica) then
         do i = 1, npole
            ii = ipole(i)
            pdi = pdamp(i)
            pti = charge_trans(i)
            ci = rpole(1,i)
            dix = rpole(2,i)
            diy = rpole(3,i)
            diz = rpole(4,i)
            qixx = rpole(5,i)
            qixy = rpole(6,i)
            qixz = rpole(7,i)
            qiyy = rpole(9,i)
            qiyz = rpole(10,i)
            qizz = rpole(13,i)
            do j = i, npole
               dscale(ipole(j)) = 1.0d0
               pscale(ipole(j)) = 1.0d0
            end do
            do j = 1, n12(ii)
               pscale(i12(j,ii)) = p2scale
            end do
            do j = 1, n13(ii)
               pscale(i13(j,ii)) = p3scale
            end do
            do j = 1, n14(ii)
               pscale(i14(j,ii)) = p4scale
            end do
            do j = 1, n15(ii)
               pscale(i15(j,ii)) = p5scale
            end do
            do j = 1, np11(ii)
               dscale(ip11(j,ii)) = d1scale
            end do
            do j = 1, np12(ii)
               dscale(ip12(j,ii)) = d2scale
            end do
            do j = 1, np13(ii)
               dscale(ip13(j,ii)) = d3scale
            end do
            do j = 1, np14(ii)
               dscale(ip14(j,ii)) = d4scale
            end do
            do k = i, npole
               kk = ipole(k)
               ck = rpole(1,k)
               dkx = rpole(2,k)
               dky = rpole(3,k)
               dkz = rpole(4,k)
               qkxx = rpole(5,k)
               qkxy = rpole(6,k)
               qkxz = rpole(7,k)
               qkyy = rpole(9,k)
               qkyz = rpole(10,k)
               qkzz = rpole(13,k)
               do m = 1, ncell
                  xr = x(kk) - x(ii)
                  yr = y(kk) - y(ii)
                  zr = z(kk) - z(ii)
                  call imager (xr,yr,zr,m)
                  r2 = xr*xr + yr* yr + zr*zr
                  if (r2 .le. off2) then
                     r = sqrt(r2)
                     scale3 = 1.0d0
                     scale5 = 1.0d0
                     scale7 = 1.0d0
                     damp = pdi * pdamp(k)
                     if (damp .ne. 0.0d0) then
                        pgamma = min(pti,charge_trans(k))
                        damp = -pgamma * (r/damp)**3
                        if (damp .gt. -50.0d0) then
                           expdamp = exp(damp)
                           scale3 = expdamp
                           scale5 = 1.0d0 - expdamp
                           scale7 = 1.0d0 - expdamp*(1.0d0-damp)
                        end if
                     end if
                     rr3 = scale3 / (r*r2)
                     rr5 = 3.0d0 * scale5 / (r*r2*r2)
                     rr7 = 15.0d0 * scale7 / (r*r2*r2*r2)
                     dir = dix*xr + diy*yr + diz*zr
                     qix = qixx*xr + qixy*yr + qixz*zr
                     qiy = qixy*xr + qiyy*yr + qiyz*zr
                     qiz = qixz*xr + qiyz*yr + qizz*zr
                     qir = qix*xr + qiy*yr + qiz*zr
                     dkr = dkx*xr + dky*yr + dkz*zr
                     qkx = qkxx*xr + qkxy*yr + qkxz*zr
                     qky = qkxy*xr + qkyy*yr + qkyz*zr
                     qkz = qkxz*xr + qkyz*yr + qkzz*zr
                     qkr = qkx*xr + qky*yr + qkz*zr
                     fid(1) = -xr*(rr3*ck-rr5*dkr+rr7*qkr)
     &                        - rr3*dkx + 2.0d0*rr5*qkx
                     fid(2) = -yr*(rr3*ck-rr5*dkr+rr7*qkr)
     &                        - rr3*dky + 2.0d0*rr5*qky
                     fid(3) = -zr*(rr3*ck-rr5*dkr+rr7*qkr)
     &                        - rr3*dkz + 2.0d0*rr5*qkz
                     fkd(1) = xr*(rr3*ci+rr5*dir+rr7*qir)
     &                        - rr3*dix - 2.0d0*rr5*qix
                     fkd(2) = yr*(rr3*ci+rr5*dir+rr7*qir)
     &                        - rr3*diy - 2.0d0*rr5*qiy
                     fkd(3) = zr*(rr3*ci+rr5*dir+rr7*qir)
     &                        - rr3*diz - 2.0d0*rr5*qiz
                     do j = 1, 3
                        fip(j) = fid(j)
                        fkp(j) = fkd(j)
                     end do
                     if (use_polymer .and. r2 .le. polycut2) then
                        do j = 1, 3
                           fid(j) = fid(j) * dscale(kk)
                           fip(j) = fip(j) * pscale(kk)
                           fkd(j) = fkd(j) * dscale(kk)
                           fkp(j) = fkp(j) * pscale(kk)
                        end do
                     end if
                     do j = 1, 3
                        field(j,i) = field(j,i) + fid(j)
                        fieldp(j,i) = fieldp(j,i) + fip(j)
                        if (ii .ne. kk) then
                           field(j,k) = field(j,k) + fkd(j)
                           fieldp(j,k) = fieldp(j,k) + fkp(j)
                        end if
                     end do
                  end if
               end do
            end do
         end do
      end if
c
c     perform deallocation of some local arrays
c
      deallocate (dscale)
      deallocate (pscale)
      return
      end
