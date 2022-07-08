c
c
c     #############################################################
c     ##  COPYRIGHT (C) 1999 by Pengyu Ren & Jay William Ponder  ##
c     ##                   All Rights Reserved                   ##
c     #############################################################
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine empole3  --  mpole/polar energy & analysis  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "empole3" calculates the electrostatic energy due to
c     atomic multipole and dipole polarizability interactions,
c     and partitions the energy among the atoms
c
c
      subroutine empole3_ct
      use sizes
      use analyz
      use energi
      use limits
      use mpole
      use potent
      use iounit
      use chargetransfer
      implicit none
      integer i,ii
c
c
c     choose the method for summing over multipole interactions
c
      if (use_ewald) then
         if (use_mlist) then
            call empole3d_ct
         else
            call empole3c_ct
         end if
      else
         if (use_mlist) then
            call empole3b_ct
         else
            call empole3a_ct
         end if
      end if
c
c     zero out energy terms and analysis which are not in use
c
      if (.not. use_chargetranfer) then
         ep_ct = 0.0d0
         do i = 1, npole
            ii = ipole(i)
            aep_ct(ii) = 0.0d0
         end do
      end if
      return
      end


c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine empole3a  --  double loop multipole analysis  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "empole3a" calculates the atomic multipole and dipole
c     polarizability interaction energy using a double loop,
c     and partitions the energy among the atoms
c
c
      subroutine empole3a_ct
      use sizes
      use action
      use analyz
      use atomid
      use atoms
      use bound
      use boxes
      use cell
      use chgpot
      use couple
      use energi
      use group
      use inform
      use inter
      use iounit
      use math
      use molcul
      use mplpot
      use mpole
      use polar
      use polgrp
      use polpot
      use potent
      use shunt
      use usage
      use chargetransfer
      implicit none
      integer i,j,k
      integer ii,kk
      integer ix,iy,iz
      integer kx,ky,kz
      real*8 e,ei_ct,fgrp
      real*8 f,fm,fp
      real*8 r,r2,xr,yr,zr
      real*8 damp_ct,expdamp_ct,damp_beta_ct
      real*8 pdi_ct,pti_ct,pgamma_ct,pti_beta_ct
      real*8 rr1,rr3,rr5
      real*8 rr7,rr9
      real*8 ci,dix,diy,diz
      real*8 uix,uiy,uiz
      real*8 qixx,qixy,qixz
      real*8 qiyy,qiyz,qizz
      real*8 ck,dkx,dky,dkz
      real*8 ukx,uky,ukz
      real*8 qkxx,qkxy,qkxz
      real*8 qkyy,qkyz,qkzz
      real*8 qix,qiy,qiz
      real*8 qkx,qky,qkz
      real*8 scale3_ct,scale5_ct
      real*8 scale7_ct
      real*8 sc(10),sci(8)
      real*8 gl(0:4),gli(3)
      real*8, allocatable :: mscale(:)
      real*8, allocatable :: pscale(:)
      logical proceed
      logical header,huge
      logical usei,usek
      logical muse,puse
      character*6 mode
c
c
c     zero out multipole and polarization energy and partitioning
c
c      nem = 0
      nep_ct = 0
c      em = 0.0d0
      ep_ct = 0.0d0
      do i = 1, n
c         aem(i) = 0.0d0
         aep_ct(i) = 0.0d0
      end do
      header = .true.
      if (npole .eq. 0)  return
c
c     check the sign of multipole components at chiral sites
c
      call chkpole
c
c     rotate the multipole components into the global frame
c
      call rotpole
c
c     compute the induced dipoles at each polarizable atom
c
       call induce_ct_akd
c
c     perform dynamic allocation of some local arrays
c
c      allocate (mscale(n))
      allocate (pscale(n))
c
c     set arrays needed to scale connected atom interactions
c
      do i = 1, n
c         mscale(i) = 1.0d0
         pscale(i) = 1.0d0
c         print*,"i=",i,"charge_trans=",charge_trans(i)
      end do
c
c     set conversion factor, cutoff and switching coefficients
c
      f = electric / dielec
      mode = 'MPOLE'
      call switch (mode)
c
c     calculate the multipole interaction energy term
c
      do i = 1, npole-1
         ii = ipole(i)
         iz = zaxis(i)
         ix = xaxis(i)
         iy = yaxis(i)
         pdi_ct = pdamp_ct(i)
         pti_ct = thole_ct(i)
         pti_beta_ct = beta_damp_ct(i) 
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
         uix = uind_ct(1,i)
         uiy = uind_ct(2,i)
         uiz = uind_ct(3,i)
         usei = (use(ii) .or. use(iz) .or. use(ix) .or. use(iy))
c
c     set interaction scaling coefficients for connected atoms
c
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
     &            pscale(i14(j,ii)) = p4scale * p41scale_ct
            end do
         end do
         do j = 1, n15(ii)
            pscale(i15(j,ii)) = p5scale
         end do
c
c     decide whether to compute the current interaction
c
         do k = i+1, npole
            kk = ipole(k)
            kz = zaxis(k)
            kx = xaxis(k)
            ky = yaxis(k)
            usek = (use(kk) .or. use(kz) .or. use(kx) .or. use(ky))
            proceed = .true.
            if (use_group)  call groups (proceed,fgrp,ii,kk,0,0,0,0)
            if (.not. use_intra)  proceed = .true.
            if (proceed)  proceed = (usei .or. usek)
c
c     compute the energy contribution for this interaction
c
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
                  ukx = uind_ct(1,k)
                  uky = uind_ct(2,k)
                  ukz = uind_ct(3,k)
c
c     construct some intermediate quadrupole values
c
                  qix = qixx*xr + qixy*yr + qixz*zr
                  qiy = qixy*xr + qiyy*yr + qiyz*zr
                  qiz = qixz*xr + qiyz*yr + qizz*zr
                  qkx = qkxx*xr + qkxy*yr + qkxz*zr
                  qky = qkxy*xr + qkyy*yr + qkyz*zr
                  qkz = qkxz*xr + qkyz*yr + qkzz*zr
c
c     calculate the scalar products for permanent multipoles
c
                  sc(2) = dix*dkx + diy*dky + diz*dkz
                  sc(3) = dix*xr + diy*yr + diz*zr
                  sc(4) = dkx*xr + dky*yr + dkz*zr
                  sc(5) = qix*xr + qiy*yr + qiz*zr
                  sc(6) = qkx*xr + qky*yr + qkz*zr
                  sc(7) = qix*dkx + qiy*dky + qiz*dkz
                  sc(8) = qkx*dix + qky*diy + qkz*diz
                  sc(9) = qix*qkx + qiy*qky + qiz*qkz
                  sc(10) = 2.0d0*(qixy*qkxy+qixz*qkxz+qiyz*qkyz)
     &                        + qixx*qkxx + qiyy*qkyy + qizz*qkzz
c
c     calculate the scalar products for polarization components
c
                  sci(2) = uix*dkx + dix*ukx + uiy*dky
     &                        + diy*uky + uiz*dkz + diz*ukz
                  sci(3) = uix*xr + uiy*yr + uiz*zr
                  sci(4) = ukx*xr + uky*yr + ukz*zr
                  sci(7) = qix*ukx + qiy*uky + qiz*ukz
                  sci(8) = qkx*uix + qky*uiy + qkz*uiz
c
c     calculate the gl functions for permanent multipoles
c
c                  gl(0) = ci*ck
c                  gl(1) = ck*sc(3) - ci*sc(4) + sc(2)
c                  gl(2) = ci*sc(6) + ck*sc(5) - sc(3)*sc(4)
c     &                       + 2.0d0*(sc(7)-sc(8)+sc(10))
c                  gl(3) = sc(3)*sc(6) - sc(4)*sc(5) - 4.0d0*sc(9)
c                  gl(4) = sc(5)*sc(6)
c
c     calculate the gl functions for polarization components
c
                  gli(1) = ck*sci(3) - ci*sci(4) + sci(2)
                  gli(2) = 2.0d0*(sci(7)-sci(8)) - sci(3)*sc(4)
     &                        - sc(3)*sci(4)
                  gli(3) = sci(3)*sc(6) - sci(4)*sc(5)
c
c     compute the energy contributions for this interaction
c
                  rr1 = 1.0d0 / r
                  rr3 = rr1 / r2
                  rr5 = 3.0d0 * rr3 / r2
                  rr7 = 5.0d0 * rr5 / r2
                  rr9 = 7.0d0 * rr7 / r2
                  scale3_ct = 1.0d0
                  scale5_ct = 1.0d0
                  scale7_ct = 1.0d0
                  damp_ct = pdi_ct * pdamp_ct(k)
                  damp_beta_ct = pti_beta_ct * beta_damp_ct(k)
                  if (damp_ct .ne. 0.0d0) then
                     pgamma_ct = min(pti_ct,thole_ct(k))
                     damp_ct = -pgamma_ct * (r/damp_ct)**3
                     if (damp_ct .gt. -50.0d0) then
                        expdamp_ct = exp(damp_ct)
             scale3_ct = 1.0d0 - expdamp_ct*damp_beta_ct
             scale5_ct = 1.0d0 - (1.0d0-damp_ct)*expdamp_ct*damp_beta_ct
             scale7_ct = 1.0d0 - (1.0d0-damp_ct+0.6d0*damp_ct**2)
     &                                          *expdamp_ct*damp_beta_ct
                     end if
                  end if
c                  e = gl(0)*rr1 + gl(1)*rr3 + gl(2)*rr5
c     &                   + gl(3)*rr7 + gl(4)*rr9
                  ei_ct = gli(1)*rr3*scale3_ct + gli(2)*rr5*scale5_ct
     &                    + gli(3)*rr7*scale7_ct
c
c     apply the energy adjustments for scaled interactions
c
                  fp = f * pscale(kk)
c                  e = fm * e
                  ei_ct = 0.5d0 * fp * ei_ct
c
c     scale the interaction based on its group membership;
c     polarization cannot be group scaled as it is not pairwise
c
                  if (use_group) then
c                     e = e * fgrp
                    ei_ct = ei_ct * fgrp
                  end if
c
c     increment the overall multipole and polarization energies
c
                  puse = (use_polar .and. pscale(kk).ne.0.0d0)
c                  if (muse)  nem = nem + 1
                  if (puse)  nep_ct = nep_ct + 1
c                  em = em + e
                  ep_ct = ep_ct + ei_ct
                  aep_ct(ii) = aep_ct(ii) + 0.5d0*ei_ct
                  aep_ct(kk) = aep_ct(kk) + 0.5d0*ei_ct 
c
c     increment the total intermolecular energy
c
                  if (molcule(ii) .ne. molcule(kk)) then
                     einter = einter + ei_ct
                  end if
c
c     print a message if the energy of this interaction is large
c
                  huge = (abs(ei_ct) .gt. 100.0d0)
                  if (debug .or. (verbose.and.huge)) then
                     if (puse) then
                        if (header) then
                           header = .false.
                           write (iout,10)
   10                      format (/,' Individual Multipole and',
     &                                ' Charge transfer Interactions :',
     &                             //,' Type',14x,'Atom Names',
     &                                15x,'Distance',6x,'Energies',
     &                                ' (MPol,Polar)',/)
                        end if
                        write (iout,20)  ii,name(ii),kk,name(kk),r,ei_ct
   20                   format (' M-Pole',4x,2(i7,'-',a3),9x,
     &                             f10.4,2x, f12.4)
                     end if
                  end if
               end if
            end if
c                  write(*,*)i, k,"charge_transfer energy", ep_ct, einter
         end do
c
c     reset interaction scaling coefficients for connected atoms
c
         do j = 1, n12(ii)
            pscale(i12(j,ii)) = 1.0d0
         end do
         do j = 1, n13(ii)
            pscale(i13(j,ii)) = 1.0d0
         end do
         do j = 1, n14(ii)
            pscale(i14(j,ii)) = 1.0d0
         end do
         do j = 1, n15(ii)
            pscale(i15(j,ii)) = 1.0d0
         end do
      end do
c
c     for periodic boundary conditions with large cutoffs
c     neighbors must be found by the replicates method
c
      if (.not. use_replica)  return
c
c     calculate interaction energy with other unit cells
c
      do i = 1, npole
         ii = ipole(i)
         iz = zaxis(i)
         ix = xaxis(i)
         iy = yaxis(k)
         pdi_ct = pdamp_ct(i)
         pti_ct = thole_ct(i)
         pti_beta_ct = beta_damp_ct(i) 
         usei = (use(ii) .or. use(iz) .or. use(ix) .or. use(iy))
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
         uix = uind_ct(1,i)
         uiy = uind_ct(2,i)
         uiz = uind_ct(3,i)
c
c     set interaction scaling coefficients for connected atoms
c
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
c
c     decide whether to compute the current interaction
c
         do k = i, npole
            kk = ipole_ct(k)
            kz = zaxis(k)
            kx = xaxis(k)
            ky = yaxis(k)
            usek = (use(kk) .or. use(kz) .or. use(kx) .or. use(ky))
            if (use_group)  call groups (proceed,fgrp,ii,kk,0,0,0,0)
            proceed = .true.
            if (proceed)  proceed = (usei .or. usek)
c
c     compute the energy contribution for this interaction
c
            if (proceed) then
               do j = 1, ncell
                  xr = x(kk) - x(ii)
                  yr = y(kk) - y(ii)
                  zr = z(kk) - z(ii)
                  call imager (xr,yr,zr,j)
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
                     ukx = uind_ct(1,k)
                     uky = uind_ct(2,k)
                     ukz = uind_ct(3,k)
c
c     construct some intermediate quadrupole values
c
                     qix = qixx*xr + qixy*yr + qixz*zr
                     qiy = qixy*xr + qiyy*yr + qiyz*zr
                     qiz = qixz*xr + qiyz*yr + qizz*zr
                     qkx = qkxx*xr + qkxy*yr + qkxz*zr
                     qky = qkxy*xr + qkyy*yr + qkyz*zr
                     qkz = qkxz*xr + qkyz*yr + qkzz*zr
c
c     calculate the scalar products for permanent multipoles
c
                     sc(2) = dix*dkx + diy*dky + diz*dkz
                     sc(3) = dix*xr + diy*yr + diz*zr
                     sc(4) = dkx*xr + dky*yr + dkz*zr
                     sc(5) = qix*xr + qiy*yr + qiz*zr
                     sc(6) = qkx*xr + qky*yr + qkz*zr
                     sc(7) = qix*dkx + qiy*dky + qiz*dkz
                     sc(8) = qkx*dix + qky*diy + qkz*diz
                     sc(9) = qix*qkx + qiy*qky + qiz*qkz
                     sc(10) = 2.0d0*(qixy*qkxy+qixz*qkxz+qiyz*qkyz)
     &                           + qixx*qkxx + qiyy*qkyy + qizz*qkzz
c
c     calculate the scalar products for polarization components
c
                     sci(2) = uix*dkx + dix*ukx + uiy*dky
     &                           + diy*uky + uiz*dkz + diz*ukz
                     sci(3) = uix*xr + uiy*yr + uiz*zr
                     sci(4) = ukx*xr + uky*yr + ukz*zr
                     sci(7) = qix*ukx + qiy*uky + qiz*ukz
                     sci(8) = qkx*uix + qky*uiy + qkz*uiz
c
c     calculate the gl functions for permanent multipoles
c
c                     gl(0) = ci*ck
c                     gl(1) = ck*sc(3) - ci*sc(4) + sc(2)
c                     gl(2) = ci*sc(6) + ck*sc(5) - sc(3)*sc(4)
c     &                          + 2.0d0*(sc(7)-sc(8)+sc(10))
c                     gl(3) = sc(3)*sc(6) - sc(4)*sc(5) - 4.0d0*sc(9)
c                     gl(4) = sc(5)*sc(6)
c
c     calculate the gl functions for polarization components
c
                     gli(1) = ck*sci(3) - ci*sci(4) + sci(2)
                     gli(2) = 2.0d0*(sci(7)-sci(8)) - sci(3)*sc(4)
     &                           - sc(3)*sci(4)
                     gli(3) = sci(3)*sc(6) - sci(4)*sc(5)
c
c     compute the energy contributions for this interaction
c
                     rr1 = 1.0d0 / r
                     rr3 = rr1 / r2
                     rr5 = 3.0d0 * rr3 / r2
                     rr7 = 5.0d0 * rr5 / r2
                     rr9 = 7.0d0 * rr7 / r2
                     scale3_ct = 1.0d0
                     scale5_ct = 1.0d0
                     scale7_ct = 1.0d0
                     damp_ct = pdi_ct * pdamp_ct(k)
                     damp_beta_ct = pti_beta_ct * beta_damp_ct(k)
                     if (damp_ct .ne. 0.0d0) then
                        pgamma_ct = min(pti_ct,thole_ct(k))
                        damp_ct = -pgamma_ct * (r/damp_ct)**3
                        if (damp_ct .gt. -50.0d0) then
                           expdamp_ct = exp(damp_ct)
                   scale3_ct = 1.0d0 - expdamp_ct*damp_beta_ct
                   scale5_ct = 1.0d0 - (1.0d0-damp_ct)
     &                            *expdamp_ct*damp_beta_ct
                   scale7_ct = 1.0d0 - (1.0d0-damp_ct+0.6d0*damp_ct**2)
     &                                          *expdamp_ct*damp_beta_ct
                        end if
                     end if
c                     e = gl(0)*rr1 + gl(1)*rr3 + gl(2)*rr5
c     &                      + gl(3)*rr7 + gl(4)*rr9
                     ei_ct = gli(1)*rr3*scale3_ct + gli(2)*rr5*scale5_ct
     &                       + gli(3)*rr7*scale7_ct
c
c     apply the energy adjustments for scaled interactions
c
c                     fm = f
                     fp = f
                     if (use_polymer) then
                        if (r2 .le. polycut2) then
                           fp = fp * pscale(kk)
                        end if
                     end if
c                     e = fm * e
                     ei_ct = 0.5d0 * fp * ei_ct
c
c     scale the interaction based on its group membership;
c     polarization cannot be group scaled as it is not pairwise
c
                     if (use_group) then
c                        e = e * fgrp
                       ei_ct = ei_ct * fgrp
                     end if
c
c     increment the overall multipole and polarization energies
c
                     if (ii .eq. kk) then
c                        e = 0.5d0 * e
                        ei_ct = 0.5d0 * ei_ct
                     end if
c                     nem = nem + 1
                     nep_ct = nep_ct + 1
c                     em = em + e
                     ep_ct = ep_ct + ei_ct
c                     aem(ii) = aem(ii) + 0.5d0*e
c                     aem(kk) = aem(kk) + 0.5d0*e
                     aep_ct(ii) = aep_ct(ii) + 0.5d0*ei_ct
                     aep_ct(kk) = aep_ct(kk) + 0.5d0*ei_ct
c
c     increment the total intermolecular energy
c
                     einter = einter + ei_ct
c
c     print a message if the energy of this interaction is large
c
                     huge = ((abs(ei_ct)) .gt. 100.0d0)
                     if (debug .or. (verbose.and.huge)) then
                        if (header) then
                           header = .false.
                           write (iout,30)
   30                      format (/,' Individual Multipole and',
     &                                ' Polarization Interactions :',
     &                             //,' Type',14x,'Atom Names',
     &                                15x,'Distance',6x,'Energies',
     &                                ' (MPol,Polar)',/)
                        end if
                     write (iout,40)  ii,name(ii),kk,name(kk),r,ei_ct
   40                   format (' M-Pole',4x,2(i7,'-',a3),1x,
     &                             '(X)',5x,f10.4,2x,f12.4)
                     end if
                  end if
               end do
            end if
         end do
c
c     reset interaction scaling coefficients for connected atoms
c
         do j = 1, n12(ii)
            pscale(i12(j,ii)) = 1.0d0
         end do
         do j = 1, n13(ii)
            pscale(i13(j,ii)) = 1.0d0
         end do
         do j = 1, n14(ii)
            pscale(i14(j,ii)) = 1.0d0
         end do
         do j = 1, n15(ii)
            pscale(i15(j,ii)) = 1.0d0
         end do
      end do

c     perform deallocation of some local arrays
c
      deallocate (pscale)
      return
      end

c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine empole3b  --  neighbor list multipole analysis  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "empole3b" calculates the atomic multipole and dipole
c     polarizability interaction energy using a neighbor list,
c     and partitions the energy among the atoms
c
c
      subroutine empole3b_ct
      use sizes
      use action
      use analyz
      use atomid
      use atoms
      use bound
      use boxes
      use chgpot
      use couple
      use energi
      use group
      use inform
      use inter
      use iounit
      use math
      use molcul
      use mplpot
      use mpole
      use neigh
      use polar
      use polgrp
      use polpot
      use potent
      use shunt
      use usage
      use chargetransfer
      implicit none
      integer i,j,k
      integer ii,kk,kkk
      integer ix,iy,iz
      integer kx,ky,kz
      integer nemo,nepo
      real*8 e,ei_ct,fgrp
      real*8 f,fm,fp
      real*8 r,r2,xr,yr,zr
      real*8 damp_ct,expdamp_ct,damp_beta_ct
      real*8 pdi_ct,pti_ct,pgamma_ct,pti_beta_ct
      real*8 rr1,rr3,rr5
      real*8 rr7,rr9
      real*8 ci,dix,diy,diz
      real*8 uix,uiy,uiz
      real*8 qixx,qixy,qixz
      real*8 qiyy,qiyz,qizz
      real*8 ck,dkx,dky,dkz
      real*8 ukx,uky,ukz
      real*8 qkxx,qkxy,qkxz
      real*8 qkyy,qkyz,qkzz
      real*8 qix,qiy,qiz
      real*8 qkx,qky,qkz
      real*8 scale3_ct,scale5_ct
      real*8 scale7_ct
      real*8 emo,epo
      real*8 eintero
      real*8 sc(10),sci(8)
      real*8 gl(0:4),gli(3)
      real*8, allocatable :: mscale(:)
      real*8, allocatable :: pscale(:)
      real*8, allocatable :: aemo(:)
      real*8, allocatable :: aepo(:)
      logical proceed
      logical header,huge
      logical usei,usek
      logical muse,puse
      character*6 mode
c
c
c     zero out multipole and polarization energy and partitioning
c
c      nem = 0
      nep_ct = 0
c      em = 0.0d0
      ep_ct = 0.0d0
      do i = 1, n
         aep_ct(i) = 0.0d0
      end do
      header = .true.
c
c     check the sign of multipole components at chiral sites
c
      call chkpole
c
c     rotate the multipole components into the global frame
c
      call rotpole
c
c     compute the induced dipoles at each polarizable atom
c
      call induce_ct_akd
c
c     perform dynamic allocation of some local arrays
c
      allocate (pscale(n))
c      allocate (aemo(n))
      allocate (aepo(n))
c
c     set arrays needed to scale connected atom interactions
c
      do i = 1, n
         pscale(i) = 1.0d0
      end do
c
c     set conversion factor, cutoff and switching coefficients
c
      f = electric / dielec
      mode = 'MPOLE'
      call switch (mode)
c
c     initialize local variables for OpenMP calculation
c
c      emo = 0.0d0
      epo = 0.0d0
      eintero = einter
c      nemo = nem
      nepo = nep_ct
      do i = 1, n
c         aemo(i) = aem(i)
         aepo(i) = aep_ct(i)
      end do
c
c
c     set OpenMP directives for the major loop structure
c
c!$OMP PARALLEL default(shared) firstprivate(f)
c!$OMP&
c!private(i,j,k,ii,ix,iy,iz,usei,kk,kx,ky,kz,usek,kkk,proceed,ei_ct,
c!$OMP& damp,expdamp_ct,pdi_ct,pti_ct,pgamma,scale3_ct,scale5_ct,scale7_ct,xr,yr,zr,r,r2,
c!$OMP& rr1,rr3,rr5,rr7,rr9,ci,dix,diy,diz,qix,qiy,qiz,qixx,qixy,qixz,
c!$OMP& qiyy,qiyz,qizz,uix,uiy,uiz,ck,dkx,dky,dkz,qkx,qky,qkz,qkxx,qkxy,
c!$OMP& qkxz,qkyy,qkyz,qkzz,ukx,uky,ukz,fgrp,fp,sc,sci,gli)
c!$OMP& firstprivate(pscale)
c!$OMP DO reduction(+:epo,eintero,nepo,aepo)
c!$OMP& schedule(guided)
c
c     calculate the multipole interaction energy term
c
      do i = 1, npole-1
         ii = ipole(i)
         iz = zaxis(i)
         ix = xaxis(i)
         iy = yaxis(i)
         pdi_ct = pdamp_ct(i)
         pti_ct = thole_ct(i)
         pti_beta_ct = beta_damp_ct(i) 
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
         uix = uind_ct(1,i)
         uiy = uind_ct(2,i)
         uiz = uind_ct(3,i)
         usei = (use(ii) .or. use(iz) .or. use(ix) .or. use(iy))
c
c     set interaction scaling coefficients for connected atoms
c
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
c
c     decide whether to compute the current interaction
c
         do kkk = 1, nelst(i)
            k = elst(kkk,i)
            kk = ipole(k)
            kz = zaxis(k)
            kx = xaxis(k)
            ky = yaxis(k)
            usek = (use(kk) .or. use(kz) .or. use(kx) .or. use(ky))
            proceed = .true.
            if (use_group)  call groups (proceed,fgrp,ii,kk,0,0,0,0)
            if (.not. use_intra)  proceed = .true.
            if (proceed)  proceed = (usei .or. usek)
c
c     compute the energy contribution for this interaction
c
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
                  ukx = uind_ct(1,k)
                  uky = uind_ct(2,k)
                  ukz = uind_ct(3,k)
c
c     construct some intermediate quadrupole values
c
                  qix = qixx*xr + qixy*yr + qixz*zr
                  qiy = qixy*xr + qiyy*yr + qiyz*zr
                  qiz = qixz*xr + qiyz*yr + qizz*zr
                  qkx = qkxx*xr + qkxy*yr + qkxz*zr
                  qky = qkxy*xr + qkyy*yr + qkyz*zr
                  qkz = qkxz*xr + qkyz*yr + qkzz*zr
c
c     calculate the scalar products for permanent multipoles
c
                  sc(2) = dix*dkx + diy*dky + diz*dkz
                  sc(3) = dix*xr + diy*yr + diz*zr
                  sc(4) = dkx*xr + dky*yr + dkz*zr
                  sc(5) = qix*xr + qiy*yr + qiz*zr
                  sc(6) = qkx*xr + qky*yr + qkz*zr
                  sc(7) = qix*dkx + qiy*dky + qiz*dkz
                  sc(8) = qkx*dix + qky*diy + qkz*diz
                  sc(9) = qix*qkx + qiy*qky + qiz*qkz
                  sc(10) = 2.0d0*(qixy*qkxy+qixz*qkxz+qiyz*qkyz)
     &                        + qixx*qkxx + qiyy*qkyy + qizz*qkzz
c
c     calculate the scalar products for polarization components
c
                  sci(2) = uix*dkx + dix*ukx + uiy*dky
     &                        + diy*uky + uiz*dkz + diz*ukz
                  sci(3) = uix*xr + uiy*yr + uiz*zr
                  sci(4) = ukx*xr + uky*yr + ukz*zr
                  sci(7) = qix*ukx + qiy*uky + qiz*ukz
                  sci(8) = qkx*uix + qky*uiy + qkz*uiz
c
c     calculate the gl functions for permanent multipoles
c
c                  gl(0) = ci*ck
c                  gl(1) = ck*sc(3) - ci*sc(4) + sc(2)
c                  gl(2) = ci*sc(6) + ck*sc(5) - sc(3)*sc(4)
c     &                       + 2.0d0*(sc(7)-sc(8)+sc(10))
c                  gl(3) = sc(3)*sc(6) - sc(4)*sc(5) - 4.0d0*sc(9)
c                  gl(4) = sc(5)*sc(6)
c
c     calculate the gl functions for polarization components
c
                  gli(1) = ck*sci(3) - ci*sci(4) + sci(2)
                  gli(2) = 2.0d0*(sci(7)-sci(8)) - sci(3)*sc(4)
     &                        - sc(3)*sci(4)
                  gli(3) = sci(3)*sc(6) - sci(4)*sc(5)
c
c     compute the energy contributions for this interaction
c
                  rr1 = 1.0d0 / r
                  rr3 = rr1 / r2
                  rr5 = 3.0d0 * rr3 / r2
                  rr7 = 5.0d0 * rr5 / r2
                  rr9 = 7.0d0 * rr7 / r2
                  scale3_ct = 1.0d0
                  scale5_ct = 1.0d0
                  scale7_ct = 1.0d0
                  damp_ct = pdi_ct * pdamp_ct(k)
                  damp_beta_ct = pti_beta_ct * beta_damp_ct(k)
                  if (damp_ct .ne. 0.0d0) then
                     pgamma_ct = min(pti_ct,thole_ct(k))
                     damp_ct = -pgamma_ct * (r/damp_ct)**3
                     if (damp_ct .gt. -50.0d0) then
                        expdamp_ct = exp(damp_ct)
             scale3_ct = 1.0d0 - expdamp_ct*damp_beta_ct
             scale5_ct = 1.0d0 - (1.0d0-damp_ct)*expdamp_ct*damp_beta_ct
             scale7_ct = 1.0d0 - (1.0d0-damp_ct+0.6d0*damp_ct**2)
     &                                          *expdamp_ct*damp_beta_ct
                     end if
                  end if
c                  e = gl(0)*rr1 + gl(1)*rr3 + gl(2)*rr5
c     &                   + gl(3)*rr7 + gl(4)*rr9
                  ei_ct = gli(1)*rr3*scale3_ct + gli(2)*rr5*scale5_ct
     &                    + gli(3)*rr7*scale7_ct
c
c     apply the energy adjustments for scaled interactions
c
                  fp = f * pscale(kk)
c                  e = fm * e
                  ei_ct = 0.5d0 * fp * ei_ct
c
c     scale the interaction based on its group membership;
c     polarization cannot be group scaled as it is not pairwise
c
                  if (use_group) then
c                     e = e * fgrp
                    ei_ct = ei_ct * fgrp
                  end if
c
c     increment the overall multipole and polarization energies
c
                  puse = (use_polar .and. pscale(kk).ne.0.0d0)
c                  if (muse)  nemo = nemo + 1
                  if (puse)  nepo = nepo + 1
c                  emo = emo + e
                  epo = epo + ei_ct
c                  aemo(ii) = aemo(ii) + 0.5d0*e
c                  aemo(kk) = aemo(kk) + 0.5d0*e
                  aepo(ii) = aepo(ii) + 0.5d0*ei_ct
                  aepo(kk) = aepo(kk) + 0.5d0*ei_ct
c
c     increment the total intermolecular energy
c
                  if (molcule(ii) .ne. molcule(kk)) then
                     eintero = eintero + ei_ct
                  end if
c
c     print a message if the energy of this interaction is large
c
                  huge = ((abs(ei_ct)) .gt. 100.0d0)
                  if (debug .or. (verbose.and.huge)) then
                     if (muse .or. puse) then
                        if (header) then
                           header = .false.
                           write (iout,10)
   10                      format (/,' Individual Multipole and',
     &                                ' Polarization Interactions :',
     &                             //,' Type',14x,'Atom Names',
     &                                15x,'Distance',6x,'Energies',
     &                                ' (MPol,Polar)',/)
                        end if
                       write (iout,20)  ii,name(ii),kk,name(kk),r,ei_ct
   20                   format (' M-Pole',4x,2(i7,'-',a3),9x,
     &                             f10.4,2x,f12.4)
                     end if
                  end if
               end if
            end if
         end do
c
c     reset interaction scaling coefficients for connected atoms
c
         do j = 1, n12(ii)
            pscale(i12(j,ii)) = 1.0d0
         end do
         do j = 1, n13(ii)
            pscale(i13(j,ii)) = 1.0d0
         end do
         do j = 1, n14(ii)
            pscale(i14(j,ii)) = 1.0d0
         end do
         do j = 1, n15(ii)
            pscale(i15(j,ii)) = 1.0d0
         end do
      end do
c
c     end OpenMP directives for the major loop structure
c
c!$OMP END DO
c!$OMP END PARALLEL
c
c     add local copies to global variables for OpenMP calculation
c
c      em = em + emo
      ep_ct = ep_ct + epo
      einter = eintero
c      nem = nemo
      nep_ct = nepo
      do i = 1, n
c         aem(i) = aemo(i)
         aep_ct(i) = aepo(i)
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (pscale)
c      deallocate (aemo)
      deallocate (aepo)
      return
      end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine empole3c  --  Ewald multipole analysis via loop  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "empole3c" calculates the atomic multipole and dipole
c     polarizability interaction energy using a particle mesh
c     Ewald summation and double loop, and partitions the energy
c     among the atoms
c
c
      subroutine empole3c_ct
      use sizes
      use action
      use analyz
      use atoms
      use boxes
      use chgpot
      use energi
      use ewald
      use inter
      use math
      use mpole
      use polar
      use chargetransfer
      implicit none
      integer i,ii
      real*8 e,ei_ct,eintra
      real*8 f,term,fterm
      real*8 cii,dii,qii,uii
      real*8 xd,yd,zd
      real*8 xu,yu,zu
      real*8 ci,dix,diy,diz
      real*8 uix,uiy,uiz
      real*8 qixx,qixy,qixz
      real*8 qiyy,qiyz,qizz
c
c
c     zero out the multipole and polarization energies
c
c      nem = 0
      nep_ct = 0
c      em = 0.0d0
      ep_ct = 0.0d0
      do i = 1, n
c         aem(i) = 0.0d0
         aep_ct(i) = 0.0d0
      end do
      if (npole .eq. 0)  return
c
c     set the energy unit conversion factor
c
      f = electric / dielec
c
c     check the sign of multipole components at chiral sites
c
      call chkpole
c
c     rotate the multipole components into the global frame
c
      call rotpole
c
c     compute the induced dipoles at each polarizable atom
c
      call induce_ct_akd
c
c     compute the reciprocal space part of the Ewald summation
c
      call emrecip_ct
c
c     compute the real space part of the Ewald summation
c
      call ereal3c_ct (eintra)
c
c     compute the self-energy part of the Ewald summation
c
      term = 2.0d0 * aewald * aewald
      fterm = -f * aewald / sqrtpi
      do i = 1, npole
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
         uix = uind_ct(1,i)
         uiy = uind_ct(2,i)
         uiz = uind_ct(3,i)
         cii = ci*ci
         dii = dix*dix + diy*diy + diz*diz
         qii = qixx*qixx + qiyy*qiyy + qizz*qizz
     &            + 2.0d0*(qixy*qixy+qixz*qixz+qiyz*qiyz)
         uii = dix*uix + diy*uiy + diz*uiz
c         e = fterm * (cii + term*(dii/3.0d0+2.0d0*term*qii/5.0d0))
         ei_ct = fterm * term * uii / 3.0d0
c         nem = nem + 1
         nep_ct = nep_ct + 1
c         em = em + e
         ep_ct = ep_ct + ei_ct
c         aem(i) = aem(i) + e
         aep_ct(i) = aep_ct(i) + ei_ct
      end do
c
c     compute the cell dipole boundary correction term
c
      if (boundary .eq. 'VACUUM') then
         xd = 0.0d0
         yd = 0.0d0
         zd = 0.0d0
         xu = 0.0d0
         yu = 0.0d0
         zu = 0.0d0
         do i = 1, npole
            ii = ipole(i)
            dix = rpole(2,i)
            diy = rpole(3,i)
            diz = rpole(4,i)
            uix = uind_ct(1,i)
            uiy = uind_ct(2,i)
            uiz = uind_ct(3,i)
            xd = xd + dix + rpole(1,i)*x(ii)
            yd = yd + diy + rpole(1,i)*y(ii)
            zd = zd + diz + rpole(1,i)*z(ii)
            xu = xu + uix
            yu = yu + uiy
            zu = zu + uiz
         end do
         term = (2.0d0/3.0d0) * f * (pi/volbox)
c         nem = nem + 1
         nep_ct = nep_ct + 1
c         em = em + term*(xd*xd+yd*yd+zd*zd)
         ep_ct = ep_ct + term*(xd*xu+yd*yu+zd*zu)
      end if
c
c     intermolecular energy is total minus intramolecular part
c
      einter = einter + ep_ct - eintra
      return
      end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine ereal3c  --  real space mpole analysis via loop  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "ereal3c" evaluates the real space portion of the Ewald sum
c     energy due to atomic multipole interactions and dipole
c     polarizability and partitions the energy among the atoms
c
c     literature reference:
c
c     W. Smith, "Point Multipoles in the Ewald Summation (Revisited)",
c     CCP5 Newsletter, 46, 18-30, 1998  (see http://www.ccp5.org/)
c
c
      subroutine ereal3c_ct (eintra)
      use sizes
      use action
      use analyz
      use atomid
      use atoms
      use bound
      use boxes
      use cell
      use chgpot
      use couple
      use energi
      use ewald
      use inform
      use iounit
      use math
      use molcul
      use mplpot
      use mpole
      use polar
      use polgrp
      use polpot
      use potent
      use shunt
      use chargetransfer
      implicit none
      integer i,j,k,m
      integer ii,kk
      real*8 e,ei_ct,eintra
      real*8 f,erfc
      real*8 r,r2,xr,yr,zr
      real*8 bfac,exp2a
      real*8 efix,eifix
      real*8 ralpha
      real*8 damp_ct,expdamp_ct,damp_beta_ct
      real*8 pdi_ct,pti_ct,pgamma_ct,pti_beta_ct
      real*8 alsq2,alsq2n
      real*8 rr1,rr3,rr5
      real*8 rr7,rr9
      real*8 ci,dix,diy,diz
      real*8 uix,uiy,uiz
      real*8 qixx,qixy,qixz
      real*8 qiyy,qiyz,qizz
      real*8 ck,dkx,dky,dkz
      real*8 ukx,uky,ukz
      real*8 qkxx,qkxy,qkxz
      real*8 qkyy,qkyz,qkzz
      real*8 qix,qiy,qiz
      real*8 qkx,qky,qkz
      real*8 scale3_ct,scale5_ct
      real*8 scale7_ct
      real*8 sc(10),sci(8)
      real*8 gl(0:4),gli(3)
      real*8 bn(0:4)
      real*8, allocatable :: mscale_ct(:)
      real*8, allocatable :: pscale(:)
      logical header,huge
      logical muse,puse
      character*6 mode
      external erfc
c
c
c     zero out the intramolecular portion of the Ewald energy
c
      eintra = 0.0d0
      header = .true.
c
c     perform dynamic allocation of some local arrays
c
      allocate (pscale(n))
c
c     set arrays needed to scale connected atom interactions
c
      do i = 1, n
         pscale(i) = 1.0d0
      end do
c
c     set conversion factor, cutoff and switching coefficients
c
      f = electric / dielec
      mode = 'EWALD'
      call switch (mode)
c
c     compute the real space portion of the Ewald summation
c
      do i = 1, npole-1
         ii = ipole(i)
         pdi_ct = pdamp_ct(i)
         pti_ct = thole_ct(i)
         pti_beta_ct = beta_damp_ct(i) 
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
         uix = uind_ct(1,i)
         uiy = uind_ct(2,i)
         uiz = uind_ct(3,i)
c
c     set interaction scaling coefficients for connected atoms
c
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
     &            pscale(i14(j,ii)) = p4scale * p41scale_ct
            end do
         end do
         do j = 1, n15(ii)
            pscale(i15(j,ii)) = p5scale
         end do
         do k = i+1, npole
            kk = ipole(k)
            xr = x(kk) - x(ii)
            yr = y(kk) - y(ii)
            zr = z(kk) - z(ii)
            call image (xr,yr,zr)
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
               ukx = uind_ct(1,k)
               uky = uind_ct(2,k)
               ukz = uind_ct(3,k)
c
c     calculate the error function damping terms
c
               ralpha = aewald * r
               bn(0) = erfc(ralpha) / r
               alsq2 = 2.0d0 * aewald**2
               alsq2n = 0.0d0
               if (aewald .gt. 0.0d0)
     &            alsq2n = 1.0d0 / (sqrtpi*aewald)
               exp2a = exp(-ralpha**2)
               do m = 1, 4
                  bfac = dble(m+m-1)
                  alsq2n = alsq2 * alsq2n
                  bn(m) = (bfac*bn(m-1)+alsq2n*exp2a) / r2
               end do
c
c     construct some intermediate quadrupole values
c
               qix = qixx*xr + qixy*yr + qixz*zr
               qiy = qixy*xr + qiyy*yr + qiyz*zr
               qiz = qixz*xr + qiyz*yr + qizz*zr
               qkx = qkxx*xr + qkxy*yr + qkxz*zr
               qky = qkxy*xr + qkyy*yr + qkyz*zr
               qkz = qkxz*xr + qkyz*yr + qkzz*zr
c
c     calculate the scalar products for permanent multipoles
c
               sc(2) = dix*dkx + diy*dky + diz*dkz
               sc(3) = dix*xr + diy*yr + diz*zr
               sc(4) = dkx*xr + dky*yr + dkz*zr
               sc(5) = qix*xr + qiy*yr + qiz*zr
               sc(6) = qkx*xr + qky*yr + qkz*zr
               sc(7) = qix*dkx + qiy*dky + qiz*dkz
               sc(8) = qkx*dix + qky*diy + qkz*diz
               sc(9) = qix*qkx + qiy*qky + qiz*qkz
               sc(10) = 2.0d0*(qixy*qkxy+qixz*qkxz+qiyz*qkyz)
     &                     + qixx*qkxx + qiyy*qkyy + qizz*qkzz
c
c     calculate the scalar products for polarization components
c
               sci(2) = uix*dkx + dix*ukx + uiy*dky
     &                     + diy*uky + uiz*dkz + diz*ukz
               sci(3) = uix*xr + uiy*yr + uiz*zr
               sci(4) = ukx*xr + uky*yr + ukz*zr
               sci(7) = qix*ukx + qiy*uky + qiz*ukz
               sci(8) = qkx*uix + qky*uiy + qkz*uiz
c
c     calculate the gl functions for permanent multipoles
c
c               gl(0) = ci*ck
c               gl(1) = ck*sc(3) - ci*sc(4) + sc(2)
c               gl(2) = ci*sc(6) + ck*sc(5) - sc(3)*sc(4)
c     &                    + 2.0d0*(sc(7)-sc(8)+sc(10))
c               gl(3) = sc(3)*sc(6) - sc(4)*sc(5) - 4.0d0*sc(9)
c               gl(4) = sc(5)*sc(6)
c
c     calculate the gl functions for polarization components
c
               gli(1) = ck*sci(3) - ci*sci(4) + sci(2)
               gli(2) = 2.0d0*(sci(7)-sci(8)) - sci(3)*sc(4)
     &                     - sc(3)*sci(4)
               gli(3) = sci(3)*sc(6) - sci(4)*sc(5)
c
c     compute the energy contributions for this interaction
c
c               e = gl(0)*bn(0) + gl(1)*bn(1) + gl(2)*bn(2)
c     &                + gl(3)*bn(3) + gl(4)*bn(4)
               ei_ct = gli(1)*bn(1) + gli(2)*bn(2) + gli(3)*bn(3)
c
c     full real space energies needed for scaled interactions
c
               rr1 = 1.0d0 / r
               rr3 = rr1 / r2
               rr5 = 3.0d0 * rr3 / r2
               rr7 = 5.0d0 * rr5 / r2
               rr9 = 7.0d0 * rr7 / r2
               scale3_ct = pscale(kk)
               scale5_ct = pscale(kk)
               scale7_ct = pscale(kk)
               damp_ct = pdi_ct * pdamp_ct(k)
               damp_beta_ct = pti_beta_ct * beta_damp_ct(k)
               if (damp_ct .ne. 0.0d0) then
                  pgamma_ct = min(pti_ct,thole_ct(k))
                  damp_ct = -pgamma_ct * (r/damp_ct)**3
                  if (damp_ct .gt. -50.0d0) then
                     expdamp_ct = exp(damp_ct)
                scale3_ct = scale3_ct * (1.0d0-expdamp_ct*damp_beta_ct)
                scale5_ct = scale5_ct *
     &                  (1.0d0-(1.0d0-damp_ct)*expdamp_ct*damp_beta_ct)
                scale7_ct = scale7_ct * (1.0d0-(1.0d0-damp_ct
     &
     &                      +0.6d0*damp_ct**2)*expdamp_ct*damp_beta_ct)
                  end if
               end if
c               efix = gl(0)*rr1 + gl(1)*rr3 + gl(2)*rr5
c     &                   + gl(3)*rr7 + gl(4)*rr9
               eifix = gli(1)*rr3*(1.0d0-scale3_ct)
     &                    + gli(2)*rr5*(1.0d0-scale5_ct)
     &                    + gli(3)*rr7*(1.0d0-scale7_ct)
c
c     apply the energy adjustments for scaled interactions
c
               ei_ct = ei_ct - eifix
c               e = f * e
               ei_ct = 0.5d0 * f * ei_ct
c
c     increment the overall multipole and polarization energies
c
c               muse = use_mpole
               puse = use_polar
c              if (muse)  nem = nem + 1
               if (puse)  nep_ct = nep_ct + 1
c               em = em + e
               ep_ct = ep_ct + ei_ct
c               aem(ii) = aem(ii) + 0.5d0*e
c               aem(kk) = aem(kk) + 0.5d0*e
               aep_ct(ii) = aep_ct(ii) + 0.5d0*ei_ct
               aep_ct(kk) = aep_ct(kk) + 0.5d0*ei_ct
c
c     increment the total intramolecular energy
c
               eifix = gli(1)*rr3*scale3_ct + gli(2)*rr5*scale5_ct
     &                    + gli(3)*rr7*scale7_ct
               eifix = 0.5d0 * f * eifix
               if (molcule(ii) .eq. molcule(kk)) then
                  eintra = eintra + eifix
               end if
c
c     print a message if the energy of this interaction is large
c
               huge = ((abs(eifix)) .gt. 100.0d0)
               if (debug .or. (verbose.and.huge)) then
                  if (muse .or. puse) then
                     if (header) then
                        header = .false.
                        write (iout,10)
   10                   format (/,' Real Space Multipole and',
     &                             ' Polarization Interactions :',
     &                          //,' Type',14x,'Atom Names',
     &                             15x,'Distance',6x,'Energies',
     &                             ' (MPol,Polar)',/)
                     end if
                     write (iout,20)  ii,name(ii),kk,name(kk),r,
     &                                eifix
   20                format (' M-Pole',4x,2(i7,'-',a3),9x,
     &                          f10.4,2x,f12.4)
                  end if
               end if
            end if
         end do
c
c     reset interaction scaling coefficients for connected atoms
c
         do j = 1, n12(ii)
            pscale(i12(j,ii)) = 1.0d0
         end do
         do j = 1, n13(ii)
            pscale(i13(j,ii)) = 1.0d0
         end do
         do j = 1, n14(ii)
            pscale(i14(j,ii)) = 1.0d0
         end do
         do j = 1, n15(ii)
            pscale(i15(j,ii)) = 1.0d0
         end do
      end do
c
c     for periodic boundary conditions with large cutoffs
c     neighbors must be found by the replicates method
c
      if (.not. use_replica)  return
c
c     calculate interaction energy with other unit cells
c
      do i = 1, npole
         ii = ipole(i)
         pdi_ct = pdamp_ct(i)
         pti_ct = thole_ct(i)
         pti_beta_ct = beta_damp_ct(i) 
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
         uix = uind_ct(1,i)
         uiy = uind_ct(2,i)
         uiz = uind_ct(3,i)
c
c     set interaction scaling coefficients for connected atoms
c
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
         do k = i, npole
            kk = ipole(k)
            do j = 1, ncell
               xr = x(kk) - x(ii)
               yr = y(kk) - y(ii)
               zr = z(kk) - z(ii)
               call imager (xr,yr,zr,j)
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
                  ukx = uind_ct(1,k)
                  uky = uind_ct(2,k)
                  ukz = uind_ct(3,k)
c
c     calculate the error function damping terms
c
                  ralpha = aewald * r
                  bn(0) = erfc(ralpha) / r
                  alsq2 = 2.0d0 * aewald**2
                  alsq2n = 0.0d0
                  if (aewald .gt. 0.0d0)
     &               alsq2n = 1.0d0 / (sqrtpi*aewald)
                  exp2a = exp(-ralpha**2)
                  do m = 1, 4
                     bfac = dble(m+m-1)
                     alsq2n = alsq2 * alsq2n
                     bn(m) = (bfac*bn(m-1)+alsq2n*exp2a) / r2
                  end do
c
c     construct some intermediate quadrupole values
c
                  qix = qixx*xr + qixy*yr + qixz*zr
                  qiy = qixy*xr + qiyy*yr + qiyz*zr
                  qiz = qixz*xr + qiyz*yr + qizz*zr
                  qkx = qkxx*xr + qkxy*yr + qkxz*zr
                  qky = qkxy*xr + qkyy*yr + qkyz*zr
                  qkz = qkxz*xr + qkyz*yr + qkzz*zr
c
c     calculate the scalar products for permanent multipoles
c
                  sc(2) = dix*dkx + diy*dky + diz*dkz
                  sc(3) = dix*xr + diy*yr + diz*zr
                  sc(4) = dkx*xr + dky*yr + dkz*zr
                  sc(5) = qix*xr + qiy*yr + qiz*zr
                  sc(6) = qkx*xr + qky*yr + qkz*zr
                  sc(7) = qix*dkx + qiy*dky + qiz*dkz
                  sc(8) = qkx*dix + qky*diy + qkz*diz
                  sc(9) = qix*qkx + qiy*qky + qiz*qkz
                  sc(10) = 2.0d0*(qixy*qkxy+qixz*qkxz+qiyz*qkyz)
     &                        + qixx*qkxx + qiyy*qkyy + qizz*qkzz
c
c     calculate the scalar products for polarization components
c
                  sci(2) = uix*dkx + dix*ukx + uiy*dky
     &                        + diy*uky + uiz*dkz + diz*ukz
                  sci(3) = uix*xr + uiy*yr + uiz*zr
                  sci(4) = ukx*xr + uky*yr + ukz*zr
                  sci(7) = qix*ukx + qiy*uky + qiz*ukz
                  sci(8) = qkx*uix + qky*uiy + qkz*uiz
c
c     calculate the gl functions for permanent multipoles
c
c                  gl(0) = ci*ck
c                  gl(1) = ck*sc(3) - ci*sc(4) + sc(2)
c                  gl(2) = ci*sc(6) + ck*sc(5) - sc(3)*sc(4)
c     &                       + 2.0d0*(sc(7)-sc(8)+sc(10))
c                  gl(3) = sc(3)*sc(6) - sc(4)*sc(5) - 4.0d0*sc(9)
c                  gl(4) = sc(5)*sc(6)
c
c     calculate the gl functions for polarization components
c
                  gli(1) = ck*sci(3) - ci*sci(4) + sci(2)
                  gli(2) = 2.0d0*(sci(7)-sci(8)) - sci(3)*sc(4)
     &                        - sc(3)*sci(4)
                  gli(3) = sci(3)*sc(6) - sci(4)*sc(5)
c
c     compute the energy contributions for this interaction
c
c                  e = gl(0)*bn(0) + gl(1)*bn(1) + gl(2)*bn(2)
c     &                   + gl(3)*bn(3) + gl(4)*bn(4)
                  ei_ct = gli(1)*bn(1) + gli(2)*bn(2) + gli(3)*bn(3)
c
c     full real space energies needed for scaled interactions
c
                  rr1 = 1.0d0 / r
                  rr3 = rr1 / r2
                  rr5 = 3.0d0 * rr3 / r2
                  rr7 = 5.0d0 * rr5 / r2
                  rr9 = 7.0d0 * rr7 / r2
                  scale3_ct = 1.0d0
                  scale5_ct = 1.0d0
                  scale7_ct = 1.0d0
                  damp_ct = pdi_ct * pdamp_ct(k)
                  damp_beta_ct = pti_beta_ct * beta_damp_ct(k)
                  if (damp_ct .ne. 0.0d0) then
                     pgamma_ct = min(pti_ct,thole_ct(k))
                     damp_ct = -pgamma_ct * (r/damp_ct)**3
                     if (damp_ct .gt. -50.0d0) then
                        expdamp_ct = exp(damp_ct)
             scale3_ct = 1.0d0 - expdamp_ct*damp_beta_ct
             scale5_ct = 1.0d0 - (1.0d0-damp_ct)*expdamp_ct*damp_beta_ct
             scale7_ct = 1.0d0 - (1.0d0-damp_ct+0.6d0*damp_ct**2)
     &                                          *expdamp_ct*damp_beta_ct
                        if (use_polymer .and. r2.le.polycut2) then
                           scale3_ct = scale3_ct * pscale(kk)
                           scale5_ct = scale5_ct * pscale(kk)
                           scale7_ct = scale7_ct * pscale(kk)
                        end if
                     end if
                  end if
c                  efix = gl(0)*rr1 + gl(1)*rr3 + gl(2)*rr5
c     &                      + gl(3)*rr7 + gl(4)*rr9
                  eifix = gli(1)*rr3*(1.0d0-scale3_ct)
     &                       + gli(2)*rr5*(1.0d0-scale5_ct)
     &                       + gli(3)*rr7*(1.0d0-scale7_ct)
c
c     apply the energy adjustments for scaled interactions
c
                  if (use_polymer .and. r2.le.polycut2)
     &             ei_ct = ei_ct - eifix
c
c     increment the overall multipole and polarization energies
c
c                  e = f * e
                  ei_ct = 0.5d0 * f * ei_ct
                  if (ii .eq. kk) then
c                     e = 0.5d0 * e
                     ei_ct = 0.5d0 * ei_ct
                  end if
c                  nem = nem + 1
                  nep_ct = nep_ct + 1
c                  em = em + e
                  ep_ct = ep_ct + ei_ct
c                  aem(ii) = aem(ii) + 0.5d0*e
c                  aem(kk) = aem(kk) + 0.5d0*e
                  aep_ct(ii) = aep_ct(ii) + 0.5d0*ei_ct
                  aep_ct(kk) = aep_ct(kk) + 0.5d0*ei_ct
c
c     print a message if the energy of this interaction is large
c
                  eifix = gli(1)*rr3*scale3_ct + gli(2)*rr5*scale5_ct
     &                       + gli(3)*rr7*scale7_ct
                  eifix = 0.5d0 * f * eifix
                  huge = ((abs(eifix)) .gt. 100.0d0)
                  if (debug .or. (verbose.and.huge)) then
                     if (header) then
                        header = .false.
                        write (iout,30)
   30                   format (/,' Real Space Multipole and',
     &                             ' Polarization Interactions :',
     &                          //,' Type',14x,'Atom Names',
     &                             15x,'Distance',6x,'Energies',
     &                             ' (MPol,Polar)',/)
                     end if
                     write (iout,40)  ii,name(ii),kk,name(kk),r,
     &                                eifix
   40                format (' M-Pole',4x,2(i7,'-',a3),1x,
     &                          '(X)',5x,f10.4,2x,f12.4)
                  end if
               end if
            end do
         end do
c
c     reset interaction scaling coefficients for connected atoms
c
         do j = 1, n12(ii)
            pscale(i12(j,ii)) = 1.0d0
         end do
         do j = 1, n13(ii)
            pscale(i13(j,ii)) = 1.0d0
         end do
         do j = 1, n14(ii)
            pscale(i14(j,ii)) = 1.0d0
         end do
         do j = 1, n15(ii)
            pscale(i15(j,ii)) = 1.0d0
         end do
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (pscale)
      return
      end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine empole3d  --  Ewald multipole analysis via list  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "empole3d" calculates the atomic multipole and dipole
c     polarizability interaction energy using a particle mesh
c     Ewald summation and a neighbor list, and partitions the
c     energy among the atoms
c
c
      subroutine empole3d_ct
      use sizes
      use action
      use analyz
      use atoms
      use boxes
      use chgpot
      use energi
      use ewald
      use inter
      use math
      use mpole
      use polar
      use chargetransfer
c      use kpolr_aniso
      implicit none
      integer i,ii
      real*8 e,ei_ct,eintra
      real*8 f,term,fterm
      real*8 cii,dii,qii,uii
      real*8 xd,yd,zd
      real*8 xu,yu,zu
      real*8 ci,dix,diy,diz
      real*8 uix,uiy,uiz
      real*8 qixx,qixy,qixz
      real*8 qiyy,qiyz,qizz
c
c
c     zero out the multipole and polarization energies
c
c      nem = 0
      nep_ct = 0
c      em = 0.0d0
      ep_ct = 0.0d0
      do i = 1, n
c         aem(i) = 0.0d0
         aep_ct(i) = 0.0d0
      end do
      if (npole .eq. 0)  return
c
c     set the energy unit conversion factor
c
      f = electric / dielec
c
c     check the sign of multipole components at chiral sites
c
      call chkpole
c
c     rotate the multipole components into the global frame
c
      call rotpole
c
c     compute the induced dipoles at each polarizable atom
c
      call induce_ct_akd
c
c     compute the reciprocal space part of the Ewald summation
c
      call emrecip_ct
c
c     compute the real space part of the Ewald summation
c
      call ereal3d_ct (eintra)
c
c     compute the self-energy part of the Ewald summation
c
      term = 2.0d0 * aewald * aewald
      fterm = -f * aewald / sqrtpi
      do i = 1, npole
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
         uix = uind_ct(1,i)
         uiy = uind_ct(2,i)
         uiz = uind_ct(3,i)
         cii = ci*ci
         dii = dix*dix + diy*diy + diz*diz
         qii = qixx*qixx + qiyy*qiyy + qizz*qizz
     &            + 2.0d0*(qixy*qixy+qixz*qixz+qiyz*qiyz)
         uii = dix*uix + diy*uiy + diz*uiz
c         e = fterm * (cii + term*(dii/3.0d0+2.0d0*term*qii/5.0d0))
         ei_ct = fterm * term * uii / 3.0d0
c         nem = nem + 1
         nep_ct = nep_ct + 1
c         em = em + e
         ep_ct = ep_ct + ei_ct
c         aem(i) = aem(i) + e
         aep_ct(i) = aep_ct(i) + ei_ct
      end do
c
c     compute the cell dipole boundary correction term
c
      if (boundary .eq. 'VACUUM') then
         xd = 0.0d0
         yd = 0.0d0
         zd = 0.0d0
         xu = 0.0d0
         yu = 0.0d0
         zu = 0.0d0
         do i = 1, npole
            ii = ipole(i)
            dix = rpole(2,i)
            diy = rpole(3,i)
            diz = rpole(4,i)
            uix = uind_ct(1,i)
            uiy = uind_ct(2,i)
            uiz = uind_ct(3,i)
            xd = xd + dix + rpole(1,i)*x(ii)
            yd = yd + diy + rpole(1,i)*y(ii)
            zd = zd + diz + rpole(1,i)*z(ii)
            xu = xu + uix
            yu = yu + uiy
            zu = zu + uiz
         end do
         term = (2.0d0/3.0d0) * f * (pi/volbox)
c         nem = nem + 1
         nep_ct = nep_ct + 1
c         em = em + term*(xd*xd+yd*yd+zd*zd)
         ep_ct = ep_ct + term*(xd*xu+yd*yu+zd*zu)
      end if
c
c     intermolecular energy is total minus intramolecular part
c
      einter = einter + ep_ct - eintra
      return
      end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine ereal3d  --  real space mpole analysis via list  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "ereal3d" evaluates the real space portion of the Ewald sum
c     energy due to atomic multipole interactions and dipole
c     polarizability and partitions the energy among the atoms
c     using a pairwise neighbor list
c
c     literature reference:
c
c     W. Smith, "Point Multipoles in the Ewald Summation (Revisited)",
c     CCP5 Newsletter, 46, 18-30, 1998  (see http://www.ccp5.org/)
c
c
      subroutine ereal3d_ct (eintra)
      use sizes
      use action
      use analyz
      use atomid
      use atoms
      use bound
      use boxes
      use chgpot
      use couple
      use energi
      use ewald
      use inform
      use iounit
      use math
      use molcul
      use mplpot
      use mpole
      use neigh
      use polar
      use polgrp
      use polpot
      use potent
      use shunt
      use chargetransfer
      implicit none
      integer i,j,k,m
      integer ii,kk,kkk
      integer nemo,nepo
      real*8 e,ei_ct,eintra
      real*8 f,erfc,r,r2
      real*8 xr,yr,zr
      real*8 bfac,exp2a
      real*8 efix,eifix
      real*8 ralpha
      real*8 damp_ct,expdamp_ct,damp_beta_ct
      real*8 pdi_ct,pti_ct,pgamma_ct,pti_beta_ct
      real*8 alsq2,alsq2n
      real*8 rr1,rr3,rr5
      real*8 rr7,rr9
      real*8 ci,dix,diy,diz
      real*8 uix,uiy,uiz
      real*8 qixx,qixy,qixz
      real*8 qiyy,qiyz,qizz
      real*8 ck,dkx,dky,dkz
      real*8 ukx,uky,ukz
      real*8 qkxx,qkxy,qkxz
      real*8 qkyy,qkyz,qkzz
      real*8 qix,qiy,qiz
      real*8 qkx,qky,qkz
      real*8 scale3_ct,scale5_ct
      real*8 scale7_ct
      real*8 emo,epo
      real*8 eintrao
      real*8 sc(10),sci(8)
      real*8 gl(0:4),gli(3)
      real*8 bn(0:4)
      real*8, allocatable :: mscale(:)
      real*8, allocatable :: pscale(:)
      real*8, allocatable :: aemo(:)
      real*8, allocatable :: aepo(:)
      logical header,huge
      logical muse,puse
      character*6 mode
      external erfc
c
c
c     zero out the intramolecular portion of the Ewald energy
c
      eintra = 0.0d0
      header = .true.
c
c     perform dynamic allocation of some local arrays
c
      allocate (pscale(n))
c      allocate (aemo(n))
      allocate (aepo(n))
c
c     set arrays needed to scale connected atom interactions
c
      do i = 1, n
         pscale(i) = 1.0d0
      end do
c
c     set conversion factor, cutoff and switching coefficients
c
      f = electric / dielec
      mode = 'EWALD'
      call switch (mode)
c
c     initialize local variables for OpenMP calculation
c
c      emo = 0.0d0
      epo = 0.0d0
      eintrao = eintra
c      nemo = nem
      nepo = nep_ct
      do i = 1, n
c         aemo(i) = aem(i)
         aepo(i) = aep_ct(i)
      end do

c
c     set OpenMP directives for the major loop structure
c
c!$OMP PARALLEL default(shared) firstprivate(f)
c!$OMP& private(i,j,k,ii,kk,kkk,e,ei_ct,efix,eifix,bfac,damp_ct,expdamp_ct,
c!$OMP& pdi_ct,pti_ct,pgamma_ct,scale3_ct,scale5_ct,scale7_ct,alsq2,alsq2n,
c!$OMP& exp2a,ralpha,xr,yr,zr,r,r2,rr1,rr3,rr5,rr7,rr9,
c!$OMP& ci,dix,diy,diz,qixx,qixy,qixz,qiyy,qiyz,qizz,uix,uiy,uiz,
c!$OMP& ck,dkx,dky,dkz,qkxx,qkxy,qkxz,qkyy,qkyz,qkzz,ukx,uky,ukz,
c!$OMP& bn,sc,gl,sci,gli)
c!$OMP& firstprivate(mscale,pscale)
c!$OMP DO reduction(+:emo,epo,eintrao,nemo,nepo,aemo,aepo)
c!$OMP& schedule(guided)
c
c     compute the real space portion of the Ewald summation
c
      do i = 1, npole
         ii = ipole(i)
         pdi_ct = pdamp_ct(i)
         pti_ct = thole_ct(i)
         pti_beta_ct = beta_damp_ct(i) 
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
         uix = uind_ct(1,i)
         uiy = uind_ct(2,i)
         uiz = uind_ct(3,i)
c
c     set interaction scaling coefficients for connected atoms
c
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
         do kkk = 1, nelst(i)
            k = elst(kkk,i)
            kk = ipole(k)
            xr = x(kk) - x(ii)
            yr = y(kk) - y(ii)
            zr = z(kk) - z(ii)
            call image (xr,yr,zr)
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
               ukx = uind_ct(1,k)
               uky = uind_ct(2,k)
               ukz = uind_ct(3,k)
c
c     calculate the error function damping terms
c
               ralpha = aewald * r
               bn(0) = erfc(ralpha) / r
               alsq2 = 2.0d0 * aewald**2
               alsq2n = 0.0d0
               if (aewald .gt. 0.0d0)
     &            alsq2n = 1.0d0 / (sqrtpi*aewald)
               exp2a = exp(-ralpha**2)
               do m = 1, 4
                  bfac = dble(m+m-1)
                  alsq2n = alsq2 * alsq2n
                  bn(m) = (bfac*bn(m-1)+alsq2n*exp2a) / r2
               end do
c
c     construct some intermediate quadrupole values
c
               qix = qixx*xr + qixy*yr + qixz*zr
               qiy = qixy*xr + qiyy*yr + qiyz*zr
               qiz = qixz*xr + qiyz*yr + qizz*zr
               qkx = qkxx*xr + qkxy*yr + qkxz*zr
               qky = qkxy*xr + qkyy*yr + qkyz*zr
               qkz = qkxz*xr + qkyz*yr + qkzz*zr
c
c     calculate the scalar products for permanent multipoles
c
               sc(2) = dix*dkx + diy*dky + diz*dkz
               sc(3) = dix*xr + diy*yr + diz*zr
               sc(4) = dkx*xr + dky*yr + dkz*zr
               sc(5) = qix*xr + qiy*yr + qiz*zr
               sc(6) = qkx*xr + qky*yr + qkz*zr
               sc(7) = qix*dkx + qiy*dky + qiz*dkz
               sc(8) = qkx*dix + qky*diy + qkz*diz
               sc(9) = qix*qkx + qiy*qky + qiz*qkz
               sc(10) = 2.0d0*(qixy*qkxy+qixz*qkxz+qiyz*qkyz)
     &                     + qixx*qkxx + qiyy*qkyy + qizz*qkzz
c
c     calculate the scalar products for polarization components
c
               sci(2) = uix*dkx + dix*ukx + uiy*dky
     &                     + diy*uky + uiz*dkz + diz*ukz
               sci(3) = uix*xr + uiy*yr + uiz*zr
               sci(4) = ukx*xr + uky*yr + ukz*zr
               sci(7) = qix*ukx + qiy*uky + qiz*ukz
               sci(8) = qkx*uix + qky*uiy + qkz*uiz
c
c     calculate the gl functions for permanent multipoles
c
c               gl(0) = ci*ck
c               gl(1) = ck*sc(3) - ci*sc(4) + sc(2)
c               gl(2) = ci*sc(6) + ck*sc(5) - sc(3)*sc(4)
c     &                    + 2.0d0*(sc(7)-sc(8)+sc(10))
c               gl(3) = sc(3)*sc(6) - sc(4)*sc(5) - 4.0d0*sc(9)
c               gl(4) = sc(5)*sc(6)
c
c     calculate the gl functions for polarization components
c
               gli(1) = ck*sci(3) - ci*sci(4) + sci(2)
               gli(2) = 2.0d0*(sci(7)-sci(8)) - sci(3)*sc(4)
     &                     - sc(3)*sci(4)
               gli(3) = sci(3)*sc(6) - sci(4)*sc(5)
c
c     compute the energy contributions for this interaction
c
c               e = gl(0)*bn(0) + gl(1)*bn(1) + gl(2)*bn(2)
c     &                + gl(3)*bn(3) + gl(4)*bn(4)
               ei_ct = gli(1)*bn(1) + gli(2)*bn(2) + gli(3)*bn(3)
c
c     full real space energies needed for scaled interactions
c
               rr1 = 1.0d0 / r
               rr3 = rr1 / r2
               rr5 = 3.0d0 * rr3 / r2
               rr7 = 5.0d0 * rr5 / r2
               rr9 = 7.0d0 * rr7 / r2
               scale3_ct = pscale(kk)
               scale5_ct = pscale(kk)
               scale7_ct = pscale(kk)
               damp_ct = pdi_ct * pdamp_ct(k)
               damp_beta_ct = pti_beta_ct * beta_damp_ct(k)
               if (damp_ct .ne. 0.0d0) then
                  pgamma_ct = min(pti_ct,thole_ct(k))
                  damp_ct = -pgamma_ct * (r/damp_ct)**3
                  if (damp_ct .gt. -50.0d0) then
                     expdamp_ct = exp(damp_ct)
                 scale3_ct = scale3_ct * (1.0d0-expdamp_ct*damp_beta_ct)
                 scale5_ct = scale5_ct *
     &                  (1.0d0-(1.0d0-damp_ct)*expdamp_ct*damp_beta_ct)
                 scale7_ct = scale7_ct * (1.0d0-(1.0d0-damp_ct
     &
     &                      +0.6d0*damp_ct**2)*expdamp_ct*damp_beta_ct)
                  end if
               end if
c               efix = gl(0)*rr1 + gl(1)*rr3 + gl(2)*rr5
c     &                   + gl(3)*rr7 + gl(4)*rr9
               eifix = gli(1)*rr3*(1.0d0-scale3_ct)
     &                    + gli(2)*rr5*(1.0d0-scale5_ct)
     &                    + gli(3)*rr7*(1.0d0-scale7_ct)
c
c     apply the energy adjustments for scaled interactions
c
               ei_ct = ei_ct - eifix
c
c     increment the overall multipole and polarization energies
c
c               e = f * e
               ei_ct = 0.5d0 * f * ei_ct
c               muse = use_mpole
               puse = use_polar
c               if (muse)  nemo = nemo + 1
               if (puse)  nepo = nepo + 1
c               emo = emo + e
               epo = epo + ei_ct
c               aemo(ii) = aemo(ii) + 0.5d0*e
c               aemo(kk) = aemo(kk) + 0.5d0*e
               aepo(ii) = aepo(ii) + 0.5d0*ei_ct
               aepo(kk) = aepo(kk) + 0.5d0*ei_ct
c
c     increment the total intramolecular energy
c
               eifix = gli(1)*rr3*scale3_ct + gli(2)*rr5*scale5_ct
     &                    + gli(3)*rr7*scale7_ct
               eifix = 0.5d0 * f * eifix
               if (molcule(ii) .eq. molcule(kk)) then
                  eintrao = eintrao + eifix
               end if
c
c     print a message if the energy of this interaction is large
c
               huge = ((abs(eifix)) .gt. 100.0d0)
               if (debug .or. (verbose.and.huge)) then
                  if (muse .or. puse) then
                     if (header) then
                        header = .false.
                        write (iout,10)
   10                   format (/,' Real Space Multipole and',
     &                             ' Polarization Interactions :',
     &                          //,' Type',14x,'Atom Names',
     &                             15x,'Distance',6x,'Energies',
     &                             ' (MPol,Polar)',/)
                     end if
                     write (iout,20)  ii,name(ii),kk,name(kk),r,
     &                                eifix
   20                format (' M-Pole',4x,2(i7,'-',a3),9x,
     &                          f10.4,2x,f12.4)
                  end if
               end if
            end if
         end do
c
c     reset interaction scaling coefficients for connected atoms
c
         do j = 1, n12(ii)
            pscale(i12(j,ii)) = 1.0d0
         end do
         do j = 1, n13(ii)
            pscale(i13(j,ii)) = 1.0d0
         end do
         do j = 1, n14(ii)
            pscale(i14(j,ii)) = 1.0d0
         end do
         do j = 1, n15(ii)
            pscale(i15(j,ii)) = 1.0d0
         end do
      end do
c
c     end OpenMP directives for the major loop structure
c
c!$OMP END DO
c!$OMP END PARALLEL
c
c     add local copies to global variables for OpenMP calculation
c
c      em = em + emo
      ep_ct = ep_ct + epo
      eintra = eintrao
c      nem = nemo
      nep_ct = nepo
      do i = 1, n
c         aem(i) = aemo(i)
         aep_ct(i) = aepo(i)
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (pscale)
c      deallocate (aemo)
      deallocate (aepo)
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine emrecip  --  PME recip space multipole energy  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "emrecip" evaluates the reciprocal space portion of the particle
c     mesh Ewald energy due to atomic multipole interactions and
c     dipole polarizability
c
c     literature reference:
c
c     C. Sagui, L. G. Pedersen and T. A. Darden, "Towards an Accurate
c     Representation of Electrostatics in Classical Force Fields:
c     Efficient Implementation of Multipolar Interactions in
c     Biomolecular Simulations", Journal of Chemical Physics, 120,
c     73-87 (2004)
c
c     modifications for nonperiodic systems suggested by Tom Darden
c     during May 2007
c
c
      subroutine emrecip_ct
      use sizes
      use bound
      use boxes
      use chgpot
      use energi
      use ewald
      use math
      use mpole
      use pme
      use polar
      use potent
      use chargetransfer
      implicit none
      integer i,j,k,ntot
      integer k1,k2,k3
      integer m1,m2,m3
      integer nff,nf1,nf2,nf3
      real*8 e,r1,r2,r3
      real*8 h1,h2,h3
      real*8 volterm,denom
      real*8 hsq,expterm
      real*8 term,pterm
      real*8 a(3,3)
      real*8, allocatable :: fuind(:,:)
      real*8, allocatable :: cmp(:,:)
      real*8, allocatable :: fmp(:,:)
      real*8, allocatable :: fphi(:,:)
c
c
c     return if the Ewald coefficient is zero
c
      if (aewald .lt. 1.0d-6)  return
c
c     perform dynamic allocation of some local arrays
c
      allocate (fuind(3,npole))
      allocate (cmp(10,npole))
      allocate (fmp(10,npole))
      allocate (fphi(20,npole))
c
c     copy the multipole moments into local storage areas
c
      do i = 1, npole
         cmp(1,i) = rpole(1,i)
         cmp(2,i) = rpole(2,i)
         cmp(3,i) = rpole(3,i)
         cmp(4,i) = rpole(4,i)
         cmp(5,i) = rpole(5,i)
         cmp(6,i) = rpole(9,i)
         cmp(7,i) = rpole(13,i)
         cmp(8,i) = 2.0d0 * rpole(6,i)
         cmp(9,i) = 2.0d0 * rpole(7,i)
         cmp(10,i) = 2.0d0 * rpole(10,i)
      end do
c
c     compute B-spline coefficients and spatial decomposition
c
      if (.not. use_chargetranfer) then
c      if (.not. use_polar) then
         call bspline_fill
         call table_fill
      end if
c
c     convert Cartesian multipoles to fractional coordinates
c
      call cmp_to_fmp (cmp,fmp)
c
c     assign PME grid and perform 3-D FFT forward transform
c
      call grid_mpole (fmp)
      call fftfront
c
c     make the scalar summation over reciprocal lattice
c
      ntot = nfft1 * nfft2 * nfft3
      pterm = (pi/aewald)**2
      volterm = pi * volbox
      nff = nfft1 * nfft2
      nf1 = (nfft1+1) / 2
      nf2 = (nfft2+1) / 2
      nf3 = (nfft3+1) / 2
      do i = 1, ntot-1
         k3 = i/nff + 1
         j = i - (k3-1)*nff
         k2 = j/nfft1 + 1
         k1 = j - (k2-1)*nfft1 + 1
         m1 = k1 - 1
         m2 = k2 - 1
         m3 = k3 - 1
         if (k1 .gt. nf1)  m1 = m1 - nfft1
         if (k2 .gt. nf2)  m2 = m2 - nfft2
         if (k3 .gt. nf3)  m3 = m3 - nfft3
         r1 = dble(m1)
         r2 = dble(m2)
         r3 = dble(m3)
         h1 = recip(1,1)*r1 + recip(1,2)*r2 + recip(1,3)*r3
         h2 = recip(2,1)*r1 + recip(2,2)*r2 + recip(2,3)*r3
         h3 = recip(3,1)*r1 + recip(3,2)*r2 + recip(3,3)*r3
         hsq = h1*h1 + h2*h2 + h3*h3
         term = -pterm * hsq
         expterm = 0.0d0
         if (term .gt. -50.0d0) then
            denom = volterm*hsq*bsmod1(k1)*bsmod2(k2)*bsmod3(k3)
            expterm = exp(term) / denom
            if (.not. use_bounds) then
               expterm = expterm * (1.0d0-cos(pi*xbox*sqrt(hsq)))
            else if (octahedron) then
               if (mod(m1+m2+m3,2) .ne. 0)  expterm = 0.0d0
            end if
         end if
         qfac(k1,k2,k3) = expterm
      end do
c
c     account for the zeroth grid point for a finite system
c
      qfac(1,1,1) = 0.0d0
      if (.not. use_bounds) then
         expterm = 0.5d0 * pi / xbox
         qfac(1,1,1) = expterm
      end if
c
c     complete the transformation of the charge grid
c
      do k = 1, nfft3
         do j = 1, nfft2
            do i = 1, nfft1
               term = qfac(i,j,k)
               qgrid(1,i,j,k) = term * qgrid(1,i,j,k)
               qgrid(2,i,j,k) = term * qgrid(2,i,j,k)
            end do
         end do
      end do
c
c     perform 3-D FFT backward transform and get potential
c
      call fftback
      call fphi_mpole (fphi)
c
c     sum over multipoles and increment total multipole energy
c
c      e = 0.0d0
c      do i = 1, npole
c         do k = 1, 10
c            e = e + fmp(k,i)*fphi(k,i)
c         end do
c      end do
c      e = 0.5d0 * electric * e
c      em = em + e
c
c     convert Cartesian induced dipoles to fractional coordinates
c
      if (use_polar) then
         do i = 1, 3
            a(1,i) = dble(nfft1) * recip(i,1)
            a(2,i) = dble(nfft2) * recip(i,2)
            a(3,i) = dble(nfft3) * recip(i,3)
         end do
         do i = 1, npole
            do k = 1, 3
               fuind(k,i) = a(k,1)*uind_ct(1,i) + a(k,2)*uind_ct(2,i)
     &                         + a(k,3)*uind_ct(3,i)
            end do
         end do
c
c     sum over induced dipoles and increment total induced energy
c
         e = 0.0d0
         do i = 1, npole
            do k = 1, 3
               e = e + fuind(k,i)*fphi(k+1,i)
            end do
         end do
         e = 0.5d0 * electric * e
         ep_ct = ep_ct + e
      end if
c
c     perform deallocation of some local arrays
c
      deallocate (fuind)
      deallocate (cmp)
      deallocate (fmp)
      deallocate (fphi)
      return
      end
