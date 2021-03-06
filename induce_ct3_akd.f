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
      subroutine induce_ct_akd
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
      use chargetransfer
      implicit none
      integer i,j,k
      real*8 norm
      logical header
c
c
c     choose the method for computation of induced dipoles
c
c      if (solvtyp(1:2) .eq. 'PB') then
c         call induce0e_ct
c      else if (solvtyp(1:2) .eq. 'GK') then
c         call induce0d_ct
c      else
         call induce0a_ct
c      end if

c
c     update the lists of previous induced dipole values
c
      if (use_pred_ct) then
         nualt_ct = min(nualt_ct+1,maxualt_ct)
         do i = 1, npole
            do j = 1, 3
               do k = nualt_ct, 2, -1
                  udalt_ct(k,j,i) = udalt_ct(k-1,j,i)
                  upalt_ct(k,j,i) = upalt_ct(k-1,j,i)
               end do
               udalt_ct(1,j,i) = uind_ct(j,i)
               upalt_ct(1,j,i) = uinp_ct(j,i)
               if (use_solv) then
                  do k = nualt, 2, -1
                     usalt_ct(k,j,i) = usalt_ct(k-1,j,i)
                     upsalt_ct(k,j,i) = upsalt_ct(k-1,j,i)
                  end do
                  usalt_ct(1,j,i) = uinds_ct(j,i)
                  upsalt_ct(1,j,i) = uinps_ct(j,i)
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
            if (charge_trans(i) .ne. 0.0d0) then
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
            norm = sqrt(uind_ct(1,i)**2+uind_ct(2,i)**2+uind_ct(3,i)**2)
               if (digits .ge. 8) then
                write (iout,60)  k,(debye*uind_ct(j,i),j=1,3),debye*norm
   60             format (i8,3x,4f16.8)
               else if (digits .ge. 6) then
                write (iout,70)  k,(debye*uind_ct(j,i),j=1,3),debye*norm
   70             format (i8,4x,4f14.6)
               else
                write (iout,80)  k,(debye*uind_ct(j,i),j=1,3),debye*norm
   80             format (i8,5x,4f12.4)
               end if
            end if
         end do
         header = .true.
         if (solvtyp.eq.'GK' .or. solvtyp.eq.'PB') then
            do i = 1, npole
               if (charge_trans(i) .ne. 0.0d0) then
                  if (header) then
                     header = .false.
                     write (iout,90)
   90                format (/,' SCRF Induced Dipole Moments',
     &                          ' (Debyes) :')
                     if (digits .ge. 8) then
                        write (iout,100)
  100                   format (/,4x,'Atom',14x,'X',15x,'Y',15x,'Z',
     &                             15x,'Total',/)
                     else if (digits .ge. 6) then
                        write (iout,110)
  110                   format (/,4x,'Atom',14x,'X',13x,'Y',13x,'Z',
     &                             12x,'Total',/)
                     else
                        write (iout,120)
  120                   format (/,4x,'Atom',14x,'X',11x,'Y',11x,'Z',
     &                             9x,'Total',/)
                     end if
                  end if
                  k = ipole(i)
         norm = sqrt(uinds_ct(1,i)**2+uinds_ct(2,i)**2+uinds_ct(3,i)**2)
                  if (digits .ge. 8) then
                     write (iout,130)  k,(debye*uinds_ct(j,i),j=1,3),
     &                                 debye*norm
  130                format (i8,3x,4f16.8)
                  else if (digits .ge. 6) then
                     write (iout,140)  k,(debye*uinds_ct(j,i),j=1,3),
     &                                 debye*norm
  140                format (i8,4x,4f14.6)
                  else
                     write (iout,150)  k,(debye*uinds_ct(j,i),j=1,3),
     &                                 debye*norm
  150                format (i8,5x,4f12.4)
                  end if
               end if
            end do
         end if
      end if
      return
      end
c
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
      integer i,j,k,iter_ct
      integer maxiter_ct
      real*8 polmin_ct
      real*8 eps_ct,epsold_ct
      real*8 epsd_ct,epsp_ct
      real*8 udsum_ct,upsum_ct
      real*8 a_ct,ap_ct,b_ct,bp_ct
      real*8 sum_ct,sump_ct
      real*8, allocatable :: poli_ct(:)
      real*8, allocatable :: field_ct(:,:)
      real*8, allocatable :: fieldp_ct(:,:)
      real*8, allocatable :: udir_ct(:,:)
      real*8, allocatable :: udirp_ct(:,:)
      real*8, allocatable :: rsd_ct(:,:)
      real*8, allocatable :: rsdp_ct(:,:)
      real*8, allocatable :: zrsd_ct(:,:)
      real*8, allocatable :: zrsdp_ct(:,:)
      real*8, allocatable :: conj_ct(:,:)
      real*8, allocatable :: conjp_ct(:,:)
      real*8, allocatable :: vec_ct(:,:)
      real*8, allocatable :: vecp_ct(:,:)
      logical done
      character*6 mode
c
c
c     zero out the induced dipoles at each site
c
      do i = 1, npole
         do j = 1, 3
            uind_ct(j,i) = 0.0d0
            uinp_ct(j,i) = 0.0d0
         end do
      end do
      if (.not. use_polar)  return
c
c     perform dynamic allocation of some local arrays
c
      allocate (field_ct(3,npole))
      allocate (fieldp_ct(3,npole))
      allocate (udir_ct(3,npole))
      allocate (udirp_ct(3,npole))
c
c     get the electrostatic field due to permanent multipoles
c
      if (use_ewald) then
         call dfield0c_ct (field_ct,fieldp_ct)
c      else if (use_mlist) then
c        call dfield0b_ct (field,fieldp)
      else
         call dfield0a_ct (field_ct,fieldp_ct)
      end if
          do i =1,n
        write(*,*)'akdctperm',field_ct(1,i),field_ct(2,i),field_ct(3,i)
          end do
c
c     set induced dipoles to polarizability times direct field
c
      do i = 1, npole
         do j = 1, 3
            udir_ct(j,i) = charge_trans(i) * field_ct(j,i)
            udirp_ct(j,i) = charge_trans(i) * fieldp_ct(j,i)
            uind_ct(j,i) = udir_ct(j,i)
            uinp_ct(j,i) = udirp_ct(j,i)
         end do
      end do
c
c     set tolerances for computation of mutual induced dipoles
c
      if (poltyp .eq. 'MUTUAL') then
         done = .false.
         maxiter_ct = 500
         iter_ct = 0
         polmin_ct = 0.00000001d0
         eps_ct = 100.0d0
c
c     estimated induced dipoles from polynomial predictor
c
         if (use_pred_ct .and. nualt_ct.eq.maxualt_ct) then
            call ulspred_ct
            do i = 1, npole
               do j = 1, 3
                  udsum_ct = 0.0d0
                  upsum_ct = 0.0d0
                  do k = 1, nualt_ct-1
                     udsum_ct = udsum_ct + bpred_ct(k)*udalt_ct(k,j,i)
                     upsum_ct = upsum_ct + bpredp_ct(k)*upalt_ct(k,j,i)
                  end do
                  uind_ct(j,i) = udsum_ct
                  uinp_ct(j,i) = upsum_ct
               end do
            end do
         end if
c
c     perform dynamic allocation of some local arrays
c
         allocate (poli_ct(npole))
         allocate (rsd_ct(3,npole))
         allocate (rsdp_ct(3,npole))
         allocate (zrsd_ct(3,npole))
         allocate (zrsdp_ct(3,npole))
         allocate (conj_ct(3,npole))
         allocate (conjp_ct(3,npole))
         allocate (vec_ct(3,npole))
         allocate (vecp_ct(3,npole))
        do i=1,n
       write(*,*)'before',uind_ct(1,i),uind_ct(2,i),uind_ct(3,i)
        end do
c
c     get the electrostatic field due to induced dipoles
c
         if (use_ewald) then
            call ufield0c_ct (field_ct,fieldp_ct)
c         else if (use_mlist) then
c            call ufield0b_ct (field,fieldp)
         else
            call ufield0a_ct (field_ct,fieldp_ct)
         end if
c
c     set initial conjugate gradient residual and conjugate vector
c
         do i = 1, npole
            poli_ct(i) = max(polmin_ct,charge_trans(i))
            do j = 1, 3
               rsd_ct(j,i) = (udir_ct(j,i)-uind_ct(j,i))/poli_ct(i)
     &                       + field_ct(j,i)
               rsdp_ct(j,i) = (udirp_ct(j,i)-uinp_ct(j,i))/poli_ct(i)
     &                       + fieldp_ct(j,i)
            end do
         end do
         mode = 'BUILD'
c         if (use_mlist) then
c            call uscale0b_ct (mode,rsd,rsdp,zrsd,zrsdp)
c            mode = 'APPLY'
c            call uscale0b_ct (mode,rsd,rsdp,zrsd,zrsdp)
c         else
            call uscale0a_ct (mode,rsd_ct,rsdp_ct,zrsd_ct,zrsdp_ct)
            mode = 'APPLY'
            call uscale0a_ct (mode,rsd_ct,rsdp_ct,zrsd_ct,zrsdp_ct)
c         end if
         do i = 1, npole
            do j = 1, 3
               conj_ct(j,i) = zrsd_ct(j,i)
               conjp_ct(j,i) = zrsdp_ct(j,i)
            end do
         end do
c
c     conjugate gradient iteration of the mutual induced dipoles
c
         do while (.not. done)
            iter_ct = iter_ct + 1
            do i = 1, npole
               do j = 1, 3
                  vec_ct(j,i) = uind_ct(j,i)
                  vecp_ct(j,i) = uinp_ct(j,i)
                  uind_ct(j,i) = conj_ct(j,i)
                  uinp_ct(j,i) = conjp_ct(j,i)
               end do
            end do
            if (use_ewald) then
               call ufield0c_ct (field_ct,fieldp_ct)
c            else if (use_mlist) then
c               call ufield0b_ct (field,fieldp)
            else
               call ufield0a_ct (field_ct,fieldp_ct)
            end if
            do i = 1, npole
               do j = 1, 3
                  uind_ct(j,i) = vec_ct(j,i)
                  uinp_ct(j,i) = vecp_ct(j,i)
                vec_ct(j,i) = conj_ct(j,i)/poli_ct(i) - field_ct(j,i)
                vecp_ct(j,i) = conjp_ct(j,i)/poli_ct(i) - fieldp_ct(j,i)
               end do
            end do
            a_ct = 0.0d0
            ap_ct = 0.0d0
            sum_ct = 0.0d0
            sump_ct = 0.0d0
            do i = 1, npole
               do j = 1, 3
                  a_ct = a_ct + conj_ct(j,i)*vec_ct(j,i)
                  ap_ct = ap_ct + conjp_ct(j,i)*vecp_ct(j,i)
                  sum_ct = sum_ct + rsd_ct(j,i)*zrsd_ct(j,i)
                  sump_ct = sump_ct + rsdp_ct(j,i)*zrsdp_ct(j,i)
               end do
            end do
            if (a_ct .ne. 0.0d0)  a_ct = sum_ct / a_ct
            if (ap_ct .ne. 0.0d0)  ap_ct = sump_ct / ap_ct
            do i = 1, npole
               do j = 1, 3
                  uind_ct(j,i) = uind_ct(j,i) + a_ct*conj_ct(j,i)
                  uinp_ct(j,i) = uinp_ct(j,i) + ap_ct*conjp_ct(j,i)
                  rsd_ct(j,i) = rsd_ct(j,i) - a_ct*vec_ct(j,i)
                  rsdp_ct(j,i) = rsdp_ct(j,i) - ap_ct*vecp_ct(j,i)
               end do
            end do
c            if (use_mlist) then
c               call uscale0b_ct (mode,rsd,rsdp,zrsd,zrsdp)
c            else
              call uscale0a_ct (mode,rsd_ct,rsdp_ct,zrsd_ct,zrsdp_ct)
c            end if
            b_ct = 0.0d0
            bp_ct = 0.0d0
            do i = 1, npole
               do j = 1, 3
                  b_ct = b_ct + rsd_ct(j,i)*zrsd_ct(j,i)
                  bp_ct = bp_ct + rsdp_ct(j,i)*zrsdp_ct(j,i)
               end do
            end do
            if (sum_ct .ne. 0.0d0)  b_ct = b_ct / sum_ct
            if (sump_ct .ne. 0.0d0)  bp_ct = bp_ct / sump_ct
            epsd_ct = 0.0d0
            epsp_ct = 0.0d0
            do i = 1, npole
               do j = 1, 3
                  conj_ct(j,i) = zrsd_ct(j,i) + b_ct*conj_ct(j,i)
                  conjp_ct(j,i) = zrsdp_ct(j,i) + bp_ct*conjp_ct(j,i)
                  epsd_ct = epsd_ct + rsd_ct(j,i)*rsd_ct(j,i)
                  epsp_ct = epsp_ct + rsdp_ct(j,i)*rsdp_ct(j,i)
               end do
            end do
c
c     check the convergence of the mutual induced dipoles
c
            epsold_ct = eps_ct
            eps_ct = max(epsd_ct,epsp_ct)
            eps_ct = debye * sqrt(eps_ct/dble(npolar))
            if (debug) then
               if (iter_ct .eq. 1) then
                  write (iout,10)
   10             format (/,' Determination of Induced Dipole',
     &                       ' Moments :',
     &                    //,4x,'Iter',8x,'RMS Change (Debyes)',/)
               end if
               write (iout,20)  iter_ct,eps_ct
   20          format (i8,7x,f16.10)
            end if
            if (eps_ct .lt. poleps) done = .true.
            if (eps_ct .gt. epsold_ct) done = .true.
            if (iter_ct .ge. politer_ct) done = .true.
         end do
c
c     perform deallocation of some local arrays
c
         deallocate (poli_ct)
         deallocate (rsd_ct)
         deallocate (rsdp_ct)
         deallocate (zrsd_ct)
         deallocate (zrsdp_ct)
         deallocate (conj_ct)
         deallocate (conjp_ct)
         deallocate (vec_ct)
         deallocate (vecp_ct)
c
c     print the results from the conjugate gradient iteration
c
         if (debug) then
            write (iout,30)  iter_ct,eps_ct
   30       format (/,' Induced Dipoles :',6x,'Iterations',i5,
     &                 6x,'RMS Change',f15.10)
         end if
c
c     terminate the calculation if dipoles failed to converge
c
         if (iter_ct.ge.maxiter_ct .or. eps_ct.gt.epsold_ct) then
            write (iout,40)
   40       format (/,' INDUCE  --  Warning, Induced Charge Transfer',
     &                 ' are not Converged')
            call prterr
            call fatal
         end if
      end if
c
c     perform deallocation of some local arrays
c
      deallocate (field_ct)
      deallocate (fieldp_ct)
      deallocate (udir_ct)
      deallocate (udirp_ct)
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
      subroutine dfield0a_ct (field_ct,fieldp_ct)
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
      real*8 damp_ct,expdamp_ct,damp_beta_ct
      real*8 scale3_ct,scale5_ct
      real*8 scale7_ct
      real*8 pdi_ct,pti_ct,pgamma_ct,pti_beta_ct
      real*8 fid_ct(3),fkd_ct(3)
      real*8 fip_ct(3),fkp_ct(3)
      real*8, allocatable :: dscale_ct(:)
      real*8, allocatable :: pscale_ct(:)
      real*8 field_ct(3,*)
      real*8 fieldp_ct(3,*)
      logical proceed
      character*6 mode
c
c
c     zero out the value of the field at each site
c
      do i = 1, npole
         do j = 1, 3
            field_ct(j,i) = 0.0d0
            fieldp_ct(j,i) = 0.0d0
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
      allocate (dscale_ct(n))
      allocate (pscale_ct(n))
c
c     find the electrostatic field due to permanent multipoles
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
         do j = i+1, npole
            dscale_ct(ipole(j)) = 1.0d0
            pscale_ct(ipole(j)) = 1.0d0
         end do
         do j = 1, n12(ii)
            pscale_ct(i12(j,ii)) = p2scale_ct
         end do
         do j = 1, n13(ii)
            pscale_ct(i13(j,ii)) = p3scale_ct
         end do
         do j = 1, n14(ii)
            pscale_ct(i14(j,ii)) = p4scale_ct
            do k = 1, np11(ii)
               if (i14(j,ii) .eq. ip11(k,ii))
     &            pscale_ct(i14(j,ii)) = p4scale_ct * p41scale_ct
            end do
         end do
         do j = 1, n15(ii)
            pscale_ct(i15(j,ii)) = p5scale_ct
         end do
         do j = 1, np11(ii)
            dscale_ct(ip11(j,ii)) = d1scale_ct
         end do
         do j = 1, np12(ii)
            dscale_ct(ip12(j,ii)) = d2scale_ct
         end do
         do j = 1, np13(ii)
            dscale_ct(ip13(j,ii)) = d3scale_ct
         end do
         do j = 1, np14(ii)
            dscale_ct(ip14(j,ii)) = d4scale_ct
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
             scale5_ct = 1.0d0 - expdamp_ct*damp_beta_ct*(1.0d0-damp_ct)
             scale7_ct = 1.0d0 - expdamp_ct*damp_beta_ct
     &                              *(1.0d0-damp_ct+0.6d0*damp_ct**2)
                     end if
                  end if
                  rr3 = scale3_ct / (r*r2)
                  rr5 = 3.0d0 * scale5_ct / (r*r2*r2)
                  rr7 = 15.0d0 * scale7_ct / (r*r2*r2*r2)
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
                  fid_ct(1) = -xr*(rr3*ck-rr5*dkr+rr7*qkr)
     &                        - rr3*dkx + 2.0d0*rr5*qkx
                  fid_ct(2) = -yr*(rr3*ck-rr5*dkr+rr7*qkr)
     &                        - rr3*dky + 2.0d0*rr5*qky
                  fid_ct(3) = -zr*(rr3*ck-rr5*dkr+rr7*qkr)
     &                        - rr3*dkz + 2.0d0*rr5*qkz
                  fkd_ct(1) = xr*(rr3*ci+rr5*dir+rr7*qir)
     &                        - rr3*dix - 2.0d0*rr5*qix
                  fkd_ct(2) = yr*(rr3*ci+rr5*dir+rr7*qir)
     &                        - rr3*diy - 2.0d0*rr5*qiy
                  fkd_ct(3) = zr*(rr3*ci+rr5*dir+rr7*qir)
     &                        - rr3*diz - 2.0d0*rr5*qiz
                  do j = 1, 3
                 field_ct(j,i) = field_ct(j,i) + fid_ct(j)*dscale_ct(kk)
                 field_ct(j,k) = field_ct(j,k) + fkd_ct(j)*dscale_ct(kk)
               fieldp_ct(j,i) = fieldp_ct(j,i) + fid_ct(j)*pscale_ct(kk)
               fieldp_ct(j,k) = fieldp_ct(j,k) + fkd_ct(j)*pscale_ct(kk)
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
            do j = i, npole
               dscale_ct(ipole(j)) = 1.0d0
               pscale_ct(ipole(j)) = 1.0d0
            end do
            do j = 1, n12(ii)
               pscale_ct(i12(j,ii)) = p2scale_ct
            end do
            do j = 1, n13(ii)
               pscale_ct(i13(j,ii)) = p3scale_ct
            end do
            do j = 1, n14(ii)
               pscale_ct(i14(j,ii)) = p4scale_ct
            end do
            do j = 1, n15(ii)
               pscale_ct(i15(j,ii)) = p5scale_ct
            end do
            do j = 1, np11(ii)
               dscale_ct(ip11(j,ii)) = d1scale_ct
            end do
            do j = 1, np12(ii)
               dscale_ct(ip12(j,ii)) = d2scale_ct
            end do
            do j = 1, np13(ii)
               dscale_ct(ip13(j,ii)) = d3scale_ct
            end do
            do j = 1, np14(ii)
               dscale_ct(ip14(j,ii)) = d4scale_ct
            end do
            do k = i, npole
               kk = ipole_ct(k)
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
             scale5_ct = 1.0d0 - expdamp_ct*damp_beta_ct*(1.0d0-damp_ct)
             scale7_ct = 1.0d0 - expdamp_ct*damp_beta_ct
     &                              *(1.0d0-damp_ct+0.6d0*damp_ct**2)
                        end if
                     end if
                     rr3 = scale3_ct / (r*r2)
                     rr5 = 3.0d0 * scale5_ct / (r*r2*r2)
                     rr7 = 15.0d0 * scale7_ct / (r*r2*r2*r2)
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
                     fid_ct(1) = -xr*(rr3*ck-rr5*dkr+rr7*qkr)
     &                        - rr3*dkx + 2.0d0*rr5*qkx
                     fid_ct(2) = -yr*(rr3*ck-rr5*dkr+rr7*qkr)
     &                        - rr3*dky + 2.0d0*rr5*qky
                     fid_ct(3) = -zr*(rr3*ck-rr5*dkr+rr7*qkr)
     &                        - rr3*dkz + 2.0d0*rr5*qkz
                     fkd_ct(1) = xr*(rr3*ci+rr5*dir+rr7*qir)
     &                        - rr3*dix - 2.0d0*rr5*qix
                     fkd_ct(2) = yr*(rr3*ci+rr5*dir+rr7*qir)
     &                        - rr3*diy - 2.0d0*rr5*qiy
                     fkd_ct(3) = zr*(rr3*ci+rr5*dir+rr7*qir)
     &                        - rr3*diz - 2.0d0*rr5*qiz
                     do j = 1, 3
                        fip_ct(j) = fid_ct(j)
                        fkp_ct(j) = fkd_ct(j)
                     end do
                     if (use_polymer .and. r2 .le. polycut2) then
                        do j = 1, 3
                           fid_ct(j) = fid_ct(j) * dscale_ct(kk)
                           fip_ct(j) = fip_ct(j) * pscale_ct(kk)
                           fkd_ct(j) = fkd_ct(j) * dscale_ct(kk)
                           fkp_ct(j) = fkp_ct(j) * pscale_ct(kk)
                        end do
                     end if
                     do j = 1, 3
                        field_ct(j,i) = field_ct(j,i) + fid_ct(j)
                        fieldp_ct(j,i) = fieldp_ct(j,i) + fip_ct(j)
                        if (ii .ne. kk) then
                           field_ct(j,k) = field_ct(j,k) + fkd_ct(j)
                           fieldp_ct(j,k) = fieldp_ct(j,k) + fkp_ct(j)
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
      deallocate (dscale_ct)
      deallocate (pscale_ct)
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine ufield0a  --  mutual induction via double loop  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "ufield0a" computes the mutual electrostatic field due to
c     induced dipole moments via a double loop
c
c
      subroutine ufield0a_ct (field_ct,fieldp_ct)
      use sizes
      use atoms
      use bound
      use cell
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
      real*8 rr3,rr5
      real*8 duix,duiy,duiz
      real*8 puix,puiy,puiz
      real*8 dukx,duky,dukz
      real*8 pukx,puky,pukz
      real*8 duir,puir
      real*8 dukr,pukr
      real*8 damp_ct,expdamp_ct,damp_beta_ct
      real*8 scale3_ct,scale5_ct
      real*8 pdi_ct,pti_ct,pgamma_ct,pti_beta_ct
      real*8 fid_ct(3),fkd_ct(3)
      real*8 fip_ct(3),fkp_ct(3)
      real*8, allocatable :: dscale_ct(:)
      real*8 field_ct(3,*)
      real*8 fieldp_ct(3,*)
      logical proceed
      character*6 mode
c
c
c     zero out the value of the field at each site
c
      do i = 1, npole
         do j = 1, 3
            field_ct(j,i) = 0.0d0
            fieldp_ct(j,i) = 0.0d0
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
      allocate (dscale_ct(n))
c
c     find the electrostatic field due to mutual induced dipoles
c
      do i = 1, npole-1
         ii = ipole(i)
         pdi_ct = pdamp_ct(i)
         pti_ct = thole_ct(i)
         pti_beta_ct = beta_damp_ct(i) 
         duix = uind_ct(1,i)
         duiy = uind_ct(2,i)
         duiz = uind_ct(3,i)
         puix = uinp_ct(1,i)
         puiy = uinp_ct(2,i)
         puiz = uinp_ct(3,i)
         do j = i+1, npole
            dscale_ct(ipole(j)) = 1.0d0
         end do
         do j = 1, np11(ii)
            dscale_ct(ip11(j,ii)) = u1scale_ct
         end do
         do j = 1, np12(ii)
            dscale_ct(ip12(j,ii)) = u2scale_ct
         end do
         do j = 1, np13(ii)
            dscale_ct(ip13(j,ii)) = u3scale_ct
         end do
         do j = 1, np14(ii)
            dscale_ct(ip14(j,ii)) = u4scale_ct
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
                  dukx = uind_ct(1,k)
                  duky = uind_ct(2,k)
                  dukz = uind_ct(3,k)
                  pukx = uinp_ct(1,k)
                  puky = uinp_ct(2,k)
                  pukz = uinp_ct(3,k)
                  scale3_ct = dscale_ct(kk)
                  scale5_ct = dscale_ct(kk)
                  damp_ct = pdi_ct * pdamp_ct(k)
                  damp_beta_ct = pti_beta_ct * beta_damp_ct(k)
                  if (damp_ct .ne. 0.0d0) then
                     pgamma_ct = min(pti_ct,thole_ct(k))
                     damp_ct = -pgamma_ct * (r/damp_ct)**3
                     if (damp_ct .gt. -50.0d0) then
                        expdamp_ct = exp(damp_ct)
                 scale3_ct = scale3_ct * (1.0d0-expdamp_ct*damp_beta_ct)
                 scale5_ct = scale5_ct * 
     &                  (1.0d0-expdamp_ct*damp_beta_ct*(1.0d0-damp_ct))
                     end if
                  end if
                  rr3 = scale3_ct / (r*r2)
                  rr5 = 3.0d0 * scale5_ct / (r*r2*r2)
                  duir = xr*duix + yr*duiy + zr*duiz
                  dukr = xr*dukx + yr*duky + zr*dukz
                  puir = xr*puix + yr*puiy + zr*puiz
                  pukr = xr*pukx + yr*puky + zr*pukz
                  fid_ct(1) = -rr3*dukx + rr5*dukr*xr
                  fid_ct(2) = -rr3*duky + rr5*dukr*yr
                  fid_ct(3) = -rr3*dukz + rr5*dukr*zr
                  fkd_ct(1) = -rr3*duix + rr5*duir*xr
                  fkd_ct(2) = -rr3*duiy + rr5*duir*yr
                  fkd_ct(3) = -rr3*duiz + rr5*duir*zr
                  fip_ct(1) = -rr3*pukx + rr5*pukr*xr
                  fip_ct(2) = -rr3*puky + rr5*pukr*yr
                  fip_ct(3) = -rr3*pukz + rr5*pukr*zr
                  fkp_ct(1) = -rr3*puix + rr5*puir*xr
                  fkp_ct(2) = -rr3*puiy + rr5*puir*yr
                  fkp_ct(3) = -rr3*puiz + rr5*puir*zr
                  do j = 1, 3
                     field_ct(j,i) = field_ct(j,i) + fid_ct(j)
                     field_ct(j,k) = field_ct(j,k) + fkd_ct(j)
                     fieldp_ct(j,i) = fieldp_ct(j,i) + fip_ct(j)
                     fieldp_ct(j,k) = fieldp_ct(j,k) + fkp_ct(j)
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
            pdi_ct = pdamp_ct(i)
            pti_ct = thole_ct(i)
            pti_beta_ct = beta_damp_ct(i) 
            duix = uind_ct(1,i)
            duiy = uind_ct(2,i)
            duiz = uind_ct(3,i)
            puix = uinp_ct(1,i)
            puiy = uinp_ct(2,i)
            puiz = uinp_ct(3,i)
            do j = i, npole
               dscale_ct(ipole(j)) = 1.0d0
            end do
            do j = 1, np11(ii)
               dscale_ct(ip11(j,ii)) = u1scale_ct
            end do
            do j = 1, np12(ii)
               dscale_ct(ip12(j,ii)) = u2scale_ct
            end do
            do j = 1, np13(ii)
               dscale_ct(ip13(j,ii)) = u3scale_ct
            end do
            do j = 1, np14(ii)
               dscale_ct(ip14(j,ii)) = u4scale_ct
            end do
            do k = i, npole
               kk = ipole_ct(k)
               dukx = uind_ct(1,k)
               duky = uind_ct(2,k)
               dukz = uind_ct(3,k)
               pukx = uinp_ct(1,k)
               puky = uinp_ct(2,k)
               pukz = uinp_ct(3,k)
               proceed = .true.
               do m = 1, ncell
                  xr = x(kk) - x(ii)
                  yr = y(kk) - y(ii)
                  zr = z(kk) - z(ii)
                  call imager (xr,yr,zr,m)
                  r2 = xr*xr + yr* yr + zr*zr
                  if (r2 .le. off2) then
                     r = sqrt(r2)
                     scale3_ct = 1.0d0
                     scale5_ct = 1.0d0
                     damp_ct = pdi_ct * pdamp_ct(k)
                     damp_beta_ct = pti_beta_ct * beta_damp_ct(k)
                     if (damp_ct .ne. 0.0d0) then
                        pgamma_ct = min(pti_ct,thole_ct(k))
                        damp_ct = -pgamma_ct * (r/damp_ct)**3
                        if (damp_ct .gt. -50.0d0) then
                           expdamp_ct = exp(damp_ct)
                       scale3_ct = 1.0d0 - expdamp_ct*damp_beta_ct
                       scale5_ct = 1.0d0 - expdamp_ct*damp_beta_ct
     &                             *(1.0d0-damp_ct)
                        end if
                     end if
                     rr3 = scale3_ct / (r*r2)
                     rr5 = 3.0d0 * scale5_ct / (r*r2*r2)
                     duir = xr*duix + yr*duiy + zr*duiz
                     dukr = xr*dukx + yr*duky + zr*dukz
                     puir = xr*puix + yr*puiy + zr*puiz
                     pukr = xr*pukx + yr*puky + zr*pukz
                     fid_ct(1) = -rr3*dukx + rr5*dukr*xr
                     fid_ct(2) = -rr3*duky + rr5*dukr*yr
                     fid_ct(3) = -rr3*dukz + rr5*dukr*zr
                     fkd_ct(1) = -rr3*duix + rr5*duir*xr
                     fkd_ct(2) = -rr3*duiy + rr5*duir*yr
                     fkd_ct(3) = -rr3*duiz + rr5*duir*zr
                     fip_ct(1) = -rr3*pukx + rr5*pukr*xr
                     fip_ct(2) = -rr3*puky + rr5*pukr*yr
                     fip_ct(3) = -rr3*pukz + rr5*pukr*zr
                     fkp_ct(1) = -rr3*puix + rr5*puir*xr
                     fkp_ct(2) = -rr3*puiy + rr5*puir*yr
                     fkp_ct(3) = -rr3*puiz + rr5*puir*zr
                     if (use_polymer) then
                        if (r2 .le. polycut2) then
                           do j = 1, 3
                              fid_ct(j) = fid_ct(j) * dscale_ct(kk)
                              fkd_ct(j) = fkd_ct(j) * dscale_ct(kk)
                              fip_ct(j) = fip_ct(j) * dscale_ct(kk)
                              fkp_ct(j) = fkp_ct(j) * dscale_ct(kk)
                           end do
                        end if
                     end if
                     do j = 1, 3
                        field_ct(j,i) = field_ct(j,i) + fid_ct(j)
                        fieldp_ct(j,i) = fieldp_ct(j,i) + fip_ct(j)
                        if (ii .ne. kk) then
                           field_ct(j,k) = field_ct(j,k) + fkd_ct(j)
                           fieldp_ct(j,k) = fieldp_ct(j,k) + fkp_ct(j)
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
      deallocate (dscale_ct)
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine dfield0b  --  direct induction via pair list  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "dfield0b" computes the mutual electrostatic field due to
c     permanent multipole moments via a pair list
c
c
      subroutine dfield0b_ct (field,fieldp)
      use sizes
      use atoms
      use bound
      use couple
      use group
      use mpole
      use neigh
      use polar
      use polgrp
      use polpot
      use shunt
      use chargetransfer
      implicit none
      integer i,j,k
      integer ii,kk,kkk
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
      real*8 damp,expdamp,damp_beta
      real*8 scale3,scale5
      real*8 scale7
      real*8 pdi,pti,pgamma,pti_beta_ct
      real*8 fid(3),fkd(3)
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
      do i = 1, npole
         ii = ipole(i)
         pdi = pdamp_ct(i)
         pti = thole_ct(i)
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
         do j = 1, nelst(i)
            dscale(ipole(elst(j,i))) = 1.0d0
            pscale(ipole(elst(j,i))) = 1.0d0
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
         do kkk = 1, nelst(i)
            k = elst(kkk,i)
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
                  damp = pdi * pdamp_ct(k)
                  damp_beta = pti_beta_ct * beta_damp_ct(k)
                  if (damp .ne. 0.0d0) then
                     pgamma = min(pti,thole_ct(k))
                     damp = -pgamma * (r/damp)**3
                     if (damp .gt. -50.0d0) then
                        expdamp = exp(damp)
                        scale3 = 1.0d0 - expdamp*damp_beta
                        scale5 = 1.0d0 - expdamp*damp_beta*(1.0d0-damp)
                        scale7 = 1.0d0 - expdamp*damp_beta
     &                              *(1.0d0-damp+0.6d0*damp**2)
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
c     perform deallocation of some local arrays
c
      deallocate (dscale)
      deallocate (pscale)
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine ufield0b  --  mutual induction via pair list  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "ufield0b" computes the mutual electrostatic field due to
c     induced dipole moments via a pair list
c
c
      subroutine ufield0b_ct (field,fieldp)
      use sizes
      use atoms
      use bound
      use group
      use mpole
      use neigh
      use polar
      use polgrp
      use polpot
      use shunt
      use chargetransfer
      implicit none
      integer i,j,k
      integer ii,kk,kkk
      real*8 xr,yr,zr
      real*8 fgrp,r,r2
      real*8 rr3,rr5
      real*8 duix,duiy,duiz
      real*8 puix,puiy,puiz
      real*8 dukx,duky,dukz
      real*8 pukx,puky,pukz
      real*8 duir,puir
      real*8 dukr,pukr
      real*8 damp,expdamp,damp_beta
      real*8 scale3,scale5
      real*8 pdi,pti,pgamma,pti_beta_ct
      real*8 fid(3),fkd(3)
      real*8 fip(3),fkp(3)
      real*8, allocatable :: dscale(:)
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
c
c     find the electrostatic field due to mutual induced dipoles
c
      do i = 1, npole
         ii = ipole(i)
         pdi = pdamp_ct(i)
         pti = thole_ct(i)
         pti_beta_ct = beta_damp_ct(i) 
         duix = uind(1,i)
         duiy = uind(2,i)
         duiz = uind(3,i)
         puix = uinp(1,i)
         puiy = uinp(2,i)
         puiz = uinp(3,i)
         do j = 1, nelst(i)
            dscale(ipole(elst(j,i))) = 1.0d0
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
         do kkk = 1, nelst(i)
            k = elst(kkk,i)
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
                  dukx = uind(1,k)
                  duky = uind(2,k)
                  dukz = uind(3,k)
                  pukx = uinp(1,k)
                  puky = uinp(2,k)
                  pukz = uinp(3,k)
                  scale3 = dscale(kk)
                  scale5 = dscale(kk)
                  damp = pdi * pdamp_ct(k)
                  damp_beta = pti_beta_ct * beta_damp_ct(k)
                  if (damp .ne. 0.0d0) then
                     pgamma = min(pti,thole_ct(k))
                     damp = -pgamma * (r/damp)**3
                     if (damp .gt. -50.0d0) then
                        expdamp = exp(damp)
                        scale3 = scale3 * (1.0d0-expdamp*damp_beta)
                        scale5 = scale5 * 
     &                    (1.0d0-expdamp*damp_beta*(1.0d0-damp))
                     end if
                  end if
                  rr3 = scale3 / (r*r2)
                  rr5 = 3.0d0 * scale5 / (r*r2*r2)
                  duir = xr*duix + yr*duiy + zr*duiz
                  dukr = xr*dukx + yr*duky + zr*dukz
                  puir = xr*puix + yr*puiy + zr*puiz
                  pukr = xr*pukx + yr*puky + zr*pukz
                  fid(1) = -rr3*dukx + rr5*dukr*xr
                  fid(2) = -rr3*duky + rr5*dukr*yr
                  fid(3) = -rr3*dukz + rr5*dukr*zr
                  fkd(1) = -rr3*duix + rr5*duir*xr
                  fkd(2) = -rr3*duiy + rr5*duir*yr
                  fkd(3) = -rr3*duiz + rr5*duir*zr
                  fip(1) = -rr3*pukx + rr5*pukr*xr
                  fip(2) = -rr3*puky + rr5*pukr*yr
                  fip(3) = -rr3*pukz + rr5*pukr*zr
                  fkp(1) = -rr3*puix + rr5*puir*xr
                  fkp(2) = -rr3*puiy + rr5*puir*yr
                  fkp(3) = -rr3*puiz + rr5*puir*zr
                  do j = 1, 3
                     field(j,i) = field(j,i) + fid(j)
                     field(j,k) = field(j,k) + fkd(j)
                     fieldp(j,i) = fieldp(j,i) + fip(j)
                     fieldp(j,k) = fieldp(j,k) + fkp(j)
                  end do
               end if
            end if
         end do
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (dscale)
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine dfield0c  --  direct induction via Ewald sum  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "dfield0c" computes the mutual electrostatic field due to
c     permanent multipole moments via Ewald summation
c
c
      subroutine dfield0c_ct (field_ct,fieldp_ct)
      use sizes
      use atoms
      use boxes
      use ewald
      use limits
      use math
      use mpole
      use polar
      implicit none
      integer i,j,ii
      real*8 term
      real*8 ucell(3)
      real*8 field_ct(3,*)
      real*8 fieldp_ct(3,*)
c
c
c     zero out the value of the field at each site
c
      do i = 1, npole
         do j = 1, 3
            field_ct(j,i) = 0.0d0
            fieldp_ct(j,i) = 0.0d0
         end do
      end do
c
c     get the reciprocal space part of the electrostatic field
c
      call udirect1_ct (field_ct)
      do i = 1, npole
         do j = 1, 3
            fieldp_ct(j,i) = field_ct(j,i)
         end do
      end do
c
c     get the real space portion of the electrostatic field
c
      if (use_mlist) then
         call udirect2b_ct (field_ct,fieldp_ct)
      else
         call udirect2a_ct (field_ct,fieldp_ct)
      end if
c
c     get the self-energy portion of the electrostatic field
c
      term = (4.0d0/3.0d0) * aewald**3 / sqrtpi
      do i = 1, npole
         do j = 1, 3
            field_ct(j,i) = field_ct(j,i) + term*rpole(j+1,i)
            fieldp_ct(j,i) = fieldp_ct(j,i) + term*rpole(j+1,i)
         end do
      end do
c
c     compute the cell dipole boundary correction to field
c
      if (boundary .eq. 'VACUUM') then
         do i = 1, 3
            ucell(i) = 0.0d0
         end do
         do i = 1, npole
            ii = ipole(i)
            ucell(1) = ucell(1) + rpole(2,i) + rpole(1,i)*x(ii)
            ucell(2) = ucell(2) + rpole(3,i) + rpole(1,i)*y(ii)
            ucell(3) = ucell(3) + rpole(4,i) + rpole(1,i)*z(ii)
         end do
         term = (4.0d0/3.0d0) * pi/volbox
         do i = 1, npole
            do j = 1, 3
               field_ct(j,i) = field_ct(j,i) - term*ucell(j)
               fieldp_ct(j,i) = fieldp_ct(j,i) - term*ucell(j)
            end do
         end do
      end if
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine ufield0c  --  mutual induction via Ewald sum  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "ufield0c" computes the mutual electrostatic field due to
c     induced dipole moments via Ewald summation
c
c
      subroutine ufield0c_ct (field_ct,fieldp_ct)
      use sizes
      use atoms
      use boxes
      use ewald
      use limits
      use math
      use mpole
      use polar
      use chargetransfer
      implicit none
      integer i,j
      real*8 term
      real*8 ucell(3)
      real*8 ucellp(3)
      real*8 field_ct(3,*)
      real*8 fieldp_ct(3,*)
c
c
c     zero out the electrostatic field at each site
c
      do i = 1, npole
         do j = 1, 3
            field_ct(j,i) = 0.0d0
            fieldp_ct(j,i) = 0.0d0
         end do
      end do
c
c     get the reciprocal space part of the electrostatic field
c
      call umutual1_ct (field_ct,fieldp_ct)
c
c     get the real space portion of the electrostatic field
c
      if (use_mlist) then
         call umutual2b_ct (field_ct,fieldp_ct)
      else
         call umutual2a_ct (field_ct,fieldp_ct)
      end if
c
c     get the self-energy portion of the electrostatic field
c
      term = (4.0d0/3.0d0) * aewald**3 / sqrtpi
      do i = 1, npole
         do j = 1, 3
            field_ct(j,i) = field_ct(j,i) + term*uind_ct(j,i)
            fieldp_ct(j,i) = fieldp_ct(j,i) + term*uinp_ct(j,i)
         end do
      end do
c
c     compute the cell dipole boundary correction to the field
c
      if (boundary .eq. 'VACUUM') then
         do i = 1, 3
            ucell(i) = 0.0d0
            ucellp(i) = 0.0d0
         end do
         do i = 1, npole
            do j = 1, 3
               ucell(j) = ucell(j) + uind_ct(j,i)
               ucellp(j) = ucellp(j) + uinp_ct(j,i)
            end do
         end do
         term = (4.0d0/3.0d0) * pi/volbox
         do i = 1, npole
            do j = 1, 3
               field_ct(j,i) = field_ct(j,i) - term*ucell(j)
               fieldp_ct(j,i) = fieldp_ct(j,i) - term*ucellp(j)
            end do
         end do
      end if
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine udirect1  --  Ewald recip direct induced field  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "udirect1" computes the reciprocal space contribution of the
c     permanent atomic multipole moments to the field
c
c
      subroutine udirect1_ct (field_ct)
      use sizes
      use bound
      use boxes
      use ewald
      use math
      use mpole
      use pme
      implicit none
      integer i,j,k,ntot
      integer k1,k2,k3
      integer m1,m2,m3
      integer nff,nf1,nf2,nf3
      real*8 r1,r2,r3
      real*8 h1,h2,h3
      real*8 volterm,denom
      real*8 hsq,expterm
      real*8 term,pterm
      real*8 field_ct(3,*)
      real*8, allocatable :: cmp(:,:)
      real*8, allocatable :: fmp(:,:)
      real*8, allocatable :: cphi(:,:)
      real*8, allocatable :: fphi(:,:)
c
c
c     return if the Ewald coefficient is zero
c
      if (aewald .lt. 1.0d-6)  return
c
c     perform dynamic allocation of some local arrays
c
      allocate (cmp(10,npole))
      allocate (fmp(10,npole))
      allocate (cphi(10,npole))
      allocate (fphi(20,npole))
c
c     copy multipole moments and coordinates to local storage
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
      call bspline_fill
      call table_fill
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
      qfac(1,1,1) = 0.0d0
      pterm = (pi/aewald)**2
      volterm = pi * volbox
      nff = nfft1 * nfft2
      nf1 = (nfft1+1) / 2
      nf2 = (nfft2+1) / 2
      nf3 = (nfft3+1) / 2
      ntot = nfft1 * nfft2 * nfft3
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
c     complete the transformation of the PME grid
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
c     perform 3-D FFT backward transform and get field
c
      call fftback
      call fphi_mpole (fphi)
c
c     convert the field from fractional to Cartesian
c
      call fphi_to_cphi (fphi,cphi)
c
c     increment the field at each multipole site
c
      do i = 1, npole
         field_ct(1,i) = field_ct(1,i) - cphi(2,i)
         field_ct(2,i) = field_ct(2,i) - cphi(3,i)
         field_ct(3,i) = field_ct(3,i) - cphi(4,i)
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (cmp)
      deallocate (fmp)
      deallocate (cphi)
      deallocate (fphi)
      return
      end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine udirect2a  --  Ewald real direct field via loop  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "udirect2a" computes the real space contribution of the permanent
c     atomic multipole moments to the field via a double loop
c
c
      subroutine udirect2a_ct (field_ct,fieldp_ct)
      use sizes
      use atoms
      use boxes
      use bound
      use cell
      use couple
      use ewald
      use math
      use mpole
      use polar
      use polgrp
      use polpot
      use shunt
      use units
      use chargetransfer
      implicit none
      integer i,j,k,m
      integer ii,kk
      real*8 xr,yr,zr,r,r2
      real*8 rr1,rr2,rr3
      real*8 rr5,rr7
      real*8 erfc,bfac,exp2a
      real*8 ci,dix,diy,diz
      real*8 qixx,qiyy,qizz
      real*8 qixy,qixz,qiyz
      real*8 ck,dkx,dky,dkz
      real*8 qkxx,qkyy,qkzz
      real*8 qkxy,qkxz,qkyz
      real*8 dir,dkr
      real*8 qix,qiy,qiz,qir
      real*8 qkx,qky,qkz,qkr
      real*8 ralpha,aefac
      real*8 aesq2,aesq2n
      real*8 pdi_ct,pti_ct,pgamma_ct,pti_beta_ct
      real*8 damp_ct,expdamp_ct,damp_beta_ct
      real*8 scale3_ct,scale5_ct
      real*8 scale7_ct
      real*8 dsc3,dsc5,dsc7
      real*8 psc3,psc5,psc7
      real*8 bn(0:3),bcn(3)
      real*8 fimd_ct(3),fkmd_ct(3)
      real*8 fimp_ct(3),fkmp_ct(3)
      real*8, allocatable :: pscale_ct(:)
      real*8, allocatable :: dscale_ct(:)
      real*8 field_ct(3,*)
      real*8 fieldp_ct(3,*)
      character*6 mode
      external erfc
c
c
c     check for multipoles and set cutoff coefficients
c
      if (npole .eq. 0)  return
      mode = 'EWALD'
      call switch (mode)
      aesq2 = 2.0 * aewald * aewald
      aesq2n = 0.0d0
      if (aewald .gt. 0.0d0)  aesq2n = 1.0d0 / (sqrtpi*aewald)
c
c     perform dynamic allocation of some local arrays
c
      allocate (pscale_ct(n))
      allocate (dscale_ct(n))
c
c     set arrays needed to scale connected atom interactions
c
      do i = 1, n
         pscale_ct(i) = 1.0d0
         dscale_ct(i) = 1.0d0
      end do
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
         do j = 1, n12(ii)
            pscale_ct(i12(j,ii)) = p2scale_ct
         end do
         do j = 1, n13(ii)
            pscale_ct(i13(j,ii)) = p3scale_ct
         end do
         do j = 1, n14(ii)
            pscale_ct(i14(j,ii)) = p4scale_ct
            do k = 1, np11(ii)
               if (i14(j,ii) .eq. ip11(k,ii))
     &            pscale_ct(i14(j,ii)) = p4scale_ct * p41scale_ct
            end do
         end do
         do j = 1, n15(ii)
            pscale_ct(i15(j,ii)) = p5scale_ct
         end do
         do j = 1, np11(ii)
            dscale_ct(ip11(j,ii)) = d1scale_ct
         end do
         do j = 1, np12(ii)
            dscale_ct(ip12(j,ii)) = d2scale_ct
         end do
         do j = 1, np13(ii)
            dscale_ct(ip13(j,ii)) = d3scale_ct
         end do
         do j = 1, np14(ii)
            dscale_ct(ip14(j,ii)) = d4scale_ct
         end do
         do k = i+1, npole
            kk = ipole(k)
            xr = x(kk) - x(ii)
            yr = y(kk) - y(ii)
            zr = z(kk) - z(ii)
            call image (xr,yr,zr)
            r2 = xr*xr + yr* yr + zr*zr
            if (r2 .le. cut2) then
               r = sqrt(r2)
               rr1 = 1.0d0 / r
               rr2 = rr1 * rr1
               rr3 = rr2 * rr1
               rr5 = rr2 * rr3
               rr7 = rr2 * rr5
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
c
c     calculate the error function damping factors
c
               ralpha = aewald * r
               bn(0) = erfc(ralpha) * rr1
               exp2a = exp(-ralpha**2)
               aefac = aesq2n
               do j = 1, 3
                  bfac = dble(j+j-1)
                  aefac = aesq2 * aefac
                  bn(j) = (bfac*bn(j-1)+aefac*exp2a) * rr2
               end do
c
c     compute the polarization damping scale factors
c
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
             scale5_ct = 1.0d0 - expdamp_ct*damp_beta_ct*(1.0d0-damp_ct)
             scale7_ct = 1.0d0 - expdamp_ct*damp_beta_ct
     &                           *(1.0d0-damp_ct+0.6d0*damp_ct**2)
                  end if
               end if
c
c     find the field terms for the current interaction
c
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
             bcn(1) = bn(1) - (1.0d0-scale3_ct*dscale_ct(kk))*rr3
             bcn(2) = bn(2) - 3.0d0*(1.0d0-scale5_ct*dscale_ct(kk))*rr5
             bcn(3) = bn(3) - 15.0d0*(1.0d0-scale7_ct*dscale_ct(kk))*rr7
               fimd_ct(1) = -xr*(bcn(1)*ck-bcn(2)*dkr+bcn(3)*qkr)
     &                     - bcn(1)*dkx + 2.0d0*bcn(2)*qkx
               fimd_ct(2) = -yr*(bcn(1)*ck-bcn(2)*dkr+bcn(3)*qkr)
     &                     - bcn(1)*dky + 2.0d0*bcn(2)*qky
               fimd_ct(3) = -zr*(bcn(1)*ck-bcn(2)*dkr+bcn(3)*qkr)
     &                     - bcn(1)*dkz + 2.0d0*bcn(2)*qkz
               fkmd_ct(1) = xr*(bcn(1)*ci+bcn(2)*dir+bcn(3)*qir)
     &                     - bcn(1)*dix - 2.0d0*bcn(2)*qix
               fkmd_ct(2) = yr*(bcn(1)*ci+bcn(2)*dir+bcn(3)*qir)
     &                     - bcn(1)*diy - 2.0d0*bcn(2)*qiy
               fkmd_ct(3) = zr*(bcn(1)*ci+bcn(2)*dir+bcn(3)*qir)
     &                     - bcn(1)*diz - 2.0d0*bcn(2)*qiz
             bcn(1) = bn(1) - (1.0d0-scale3_ct*pscale_ct(kk))*rr3
             bcn(2) = bn(2) - 3.0d0*(1.0d0-scale5_ct*pscale_ct(kk))*rr5
             bcn(3) = bn(3) - 15.0d0*(1.0d0-scale7_ct*pscale_ct(kk))*rr7
               fimp_ct(1) = -xr*(bcn(1)*ck-bcn(2)*dkr+bcn(3)*qkr)
     &                     - bcn(1)*dkx + 2.0d0*bcn(2)*qkx
               fimp_ct(2) = -yr*(bcn(1)*ck-bcn(2)*dkr+bcn(3)*qkr)
     &                     - bcn(1)*dky + 2.0d0*bcn(2)*qky
               fimp_ct(3) = -zr*(bcn(1)*ck-bcn(2)*dkr+bcn(3)*qkr)
     &                     - bcn(1)*dkz + 2.0d0*bcn(2)*qkz
               fkmp_ct(1) = xr*(bcn(1)*ci+bcn(2)*dir+bcn(3)*qir)
     &                     - bcn(1)*dix - 2.0d0*bcn(2)*qix
               fkmp_ct(2) = yr*(bcn(1)*ci+bcn(2)*dir+bcn(3)*qir)
     &                     - bcn(1)*diy - 2.0d0*bcn(2)*qiy
               fkmp_ct(3) = zr*(bcn(1)*ci+bcn(2)*dir+bcn(3)*qir)
     &                     - bcn(1)*diz - 2.0d0*bcn(2)*qiz
c
c     increment the field at each site due to this interaction
c
               do j = 1, 3
                  field_ct(j,i) = field_ct(j,i) + fimd_ct(j)
                  field_ct(j,k) = field_ct(j,k) + fkmd_ct(j)
                  fieldp_ct(j,i) = fieldp_ct(j,i) + fimp_ct(j)
                  fieldp_ct(j,k) = fieldp_ct(j,k) + fkmp_ct(j)
               end do
            end if
         end do
c
c     reset interaction scaling coefficients for connected atoms
c
         do j = 1, n12(ii)
            pscale_ct(i12(j,ii)) = 1.0d0
         end do
         do j = 1, n13(ii)
            pscale_ct(i13(j,ii)) = 1.0d0
         end do
         do j = 1, n14(ii)
            pscale_ct(i14(j,ii)) = 1.0d0
         end do
         do j = 1, n15(ii)
            pscale_ct(i15(j,ii)) = 1.0d0
         end do
         do j = 1, np11(ii)
            dscale_ct(ip11(j,ii)) = 1.0d0
         end do
         do j = 1, np12(ii)
            dscale_ct(ip12(j,ii)) = 1.0d0
         end do
         do j = 1, np13(ii)
            dscale_ct(ip13(j,ii)) = 1.0d0
         end do
         do j = 1, np14(ii)
            dscale_ct(ip14(j,ii)) = 1.0d0
         end do
      end do
c
c     periodic boundary for large cutoffs via replicates method
c
      if (use_replica) then
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
            do j = 1, n12(ii)
               pscale_ct(i12(j,ii)) = p2scale_ct
            end do
            do j = 1, n13(ii)
               pscale_ct(i13(j,ii)) = p3scale_ct
            end do
            do j = 1, n14(ii)
               pscale_ct(i14(j,ii)) = p4scale_ct
            end do
            do j = 1, n15(ii)
               pscale_ct(i15(j,ii)) = p5scale_ct
            end do
            do j = 1, np11(ii)
               dscale_ct(ip11(j,ii)) = d1scale_ct
            end do
            do j = 1, np12(ii)
               dscale_ct(ip12(j,ii)) = d2scale_ct
            end do
            do j = 1, np13(ii)
               dscale_ct(ip13(j,ii)) = d3scale_ct
            end do
            do j = 1, np14(ii)
               dscale_ct(ip14(j,ii)) = d4scale_ct
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
c
c     calculate the error function damping factors
c
                  if (r2 .le. cut2) then
                     r = sqrt(r2)
                     rr1 = 1.0d0 / r
                     rr2 = rr1 * rr1
                     rr3 = rr2 * rr1
                     rr5 = rr2 * rr3
                     rr7 = rr2 * rr5
                     ralpha = aewald * r
                     bn(0) = erfc(ralpha) * rr1
                     exp2a = exp(-ralpha**2)
                     aefac = aesq2n
                     do j = 1, 3
                        bfac = dble(j+j-1)
                        aefac = aesq2 * aefac
                        bn(j) = (bfac*bn(j-1)+aefac*exp2a) * rr2
                     end do
c
c     compute the polarization damping scale factors
c
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
             scale5_ct = 1.0d0 - expdamp_ct*damp_beta_ct*(1.0d0-damp_ct)
             scale7_ct = 1.0d0 - expdamp_ct*damp_beta_ct
     &                           *(1.0d0-damp_ct+0.6d0*damp_ct**2)
                        end if
                     end if
                     dsc3 = scale3_ct
                     dsc5 = scale5_ct
                     dsc7 = scale7_ct
                     psc3 = scale3_ct
                     psc5 = scale5_ct
                     psc7 = scale7_ct
                     if (use_polymer) then
                        if (r2 .le. polycut2) then
                           dsc3 = scale3_ct * dscale_ct(kk)
                           dsc5 = scale5_ct * dscale_ct(kk)
                           dsc7 = scale7_ct * dscale_ct(kk)
                           psc3 = scale3_ct * pscale_ct(kk)
                           psc5 = scale5_ct * pscale_ct(kk)
                           psc7 = scale7_ct * pscale_ct(kk)
                        end if
                     end if
c
c     find the field terms for the current interaction
c
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
                     bcn(1) = bn(1) - (1.0d0-dsc3)*rr3
                     bcn(2) = bn(2) - 3.0d0*(1.0d0-dsc5)*rr5
                     bcn(3) = bn(3) - 15.0d0*(1.0d0-dsc7)*rr7
                     fimd_ct(1) = -xr*(bcn(1)*ck-bcn(2)*dkr+bcn(3)*qkr)
     &                           - bcn(1)*dkx + 2.0d0*bcn(2)*qkx
                     fimd_ct(2) = -yr*(bcn(1)*ck-bcn(2)*dkr+bcn(3)*qkr)
     &                           - bcn(1)*dky + 2.0d0*bcn(2)*qky
                     fimd_ct(3) = -zr*(bcn(1)*ck-bcn(2)*dkr+bcn(3)*qkr)
     &                           - bcn(1)*dkz + 2.0d0*bcn(2)*qkz
                     fkmd_ct(1) = xr*(bcn(1)*ci+bcn(2)*dir+bcn(3)*qir)
     &                           - bcn(1)*dix - 2.0d0*bcn(2)*qix
                     fkmd_ct(2) = yr*(bcn(1)*ci+bcn(2)*dir+bcn(3)*qir)
     &                           - bcn(1)*diy - 2.0d0*bcn(2)*qiy
                     fkmd_ct(3) = zr*(bcn(1)*ci+bcn(2)*dir+bcn(3)*qir)
     &                           - bcn(1)*diz - 2.0d0*bcn(2)*qiz
                     bcn(1) = bn(1) - (1.0d0-psc3)*rr3
                     bcn(2) = bn(2) - 3.0d0*(1.0d0-psc5)*rr5
                     bcn(3) = bn(3) - 15.0d0*(1.0d0-psc7)*rr7
                     fimp_ct(1) = -xr*(bcn(1)*ck-bcn(2)*dkr+bcn(3)*qkr)
     &                           - bcn(1)*dkx + 2.0d0*bcn(2)*qkx
                     fimp_ct(2) = -yr*(bcn(1)*ck-bcn(2)*dkr+bcn(3)*qkr)
     &                           - bcn(1)*dky + 2.0d0*bcn(2)*qky
                     fimp_ct(3) = -zr*(bcn(1)*ck-bcn(2)*dkr+bcn(3)*qkr)
     &                           - bcn(1)*dkz + 2.0d0*bcn(2)*qkz
                     fkmp_ct(1) = xr*(bcn(1)*ci+bcn(2)*dir+bcn(3)*qir)
     &                           - bcn(1)*dix - 2.0d0*bcn(2)*qix
                     fkmp_ct(2) = yr*(bcn(1)*ci+bcn(2)*dir+bcn(3)*qir)
     &                           - bcn(1)*diy - 2.0d0*bcn(2)*qiy
                     fkmp_ct(3) = zr*(bcn(1)*ci+bcn(2)*dir+bcn(3)*qir)
     &                           - bcn(1)*diz - 2.0d0*bcn(2)*qiz
c
c     increment the field at each site due to this interaction
c
                     do j = 1, 3
                        field_ct(j,i) = field_ct(j,i) + fimd_ct(j)
                        fieldp_ct(j,i) = fieldp_ct(j,i) + fimd_ct(j)
                        if (ii .ne. kk) then
                           field_ct(j,k) = field_ct(j,k) + fkmp_ct(j)
                           fieldp_ct(j,k) = fieldp_ct(j,k) + fkmp_ct(j)
                        end if
                     end do
                  end if
               end do
            end do
c
c     reset interaction scaling coefficients for connected atoms
c
            do j = 1, n12(ii)
               pscale_ct(i12(j,ii)) = 1.0d0
            end do
            do j = 1, n13(ii)
               pscale_ct(i13(j,ii)) = 1.0d0
            end do
            do j = 1, n14(ii)
               pscale_ct(i14(j,ii)) = 1.0d0
            end do
            do j = 1, n15(ii)
               pscale_ct(i15(j,ii)) = 1.0d0
            end do
            do j = 1, np11(ii)
               dscale_ct(ip11(j,ii)) = 1.0d0
            end do
            do j = 1, np12(ii)
               dscale_ct(ip12(j,ii)) = 1.0d0
            end do
            do j = 1, np13(ii)
               dscale_ct(ip13(j,ii)) = 1.0d0
            end do
            do j = 1, np14(ii)
               dscale_ct(ip14(j,ii)) = 1.0d0
            end do
         end do
      end if
c
c     perform deallocation of some local arrays
c
      deallocate (dscale_ct)
      deallocate (pscale_ct)
      return
      end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine udirect2b  --  Ewald real direct field via list  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "udirect2b" computes the real space contribution of the permanent
c     atomic multipole moments to the field via a neighbor list
c
c
      subroutine udirect2b_ct (field,fieldp)
      use sizes
      use atoms
      use boxes
      use bound
      use couple
      use ewald
      use math
      use mpole
      use neigh
      use openmp
      use polar
      use polgrp
      use polpot
      use shunt
      use tarray
      use units
      use chargetransfer
      implicit none
      integer i,j,k,m
      integer ii,kk,kkk
      integer nlocal,maxlocal
      integer tid,toffset0
!$    integer omp_get_thread_num
      integer, allocatable :: toffset(:)
      integer, allocatable :: ilocal(:,:)
      real*8 xr,yr,zr,r,r2
      real*8 rr1,rr2,rr3
      real*8 rr5,rr7
      real*8 erfc,bfac,exp2a
      real*8 ci,dix,diy,diz
      real*8 qixx,qiyy,qizz
      real*8 qixy,qixz,qiyz
      real*8 ck,dkx,dky,dkz
      real*8 qkxx,qkyy,qkzz
      real*8 qkxy,qkxz,qkyz
      real*8 dir,dkr
      real*8 qix,qiy,qiz,qir
      real*8 qkx,qky,qkz,qkr
      real*8 ralpha,aefac
      real*8 aesq2,aesq2n
      real*8 pdi,pti,pgamma,pti_beta_ct
      real*8 damp,expdamp,damp_beta
      real*8 scale3,scale5
      real*8 scale7
      real*8 bn(0:3),bcn(3)
      real*8 fimd(3),fkmd(3)
      real*8 fimp(3),fkmp(3)
      real*8, allocatable :: pscale(:)
      real*8, allocatable :: dscale(:)
      real*8, allocatable :: uscale(:)
      real*8 field(3,*)
      real*8 fieldp(3,*)
      real*8, allocatable :: fieldt(:,:)
      real*8, allocatable :: fieldtp(:,:)
      real*8, allocatable :: dlocal(:,:)
      character*6 mode
      external erfc
c
c
c     check for multipoles and set cutoff coefficients
c
      if (npole .eq. 0)  return
      mode = 'EWALD'
      call switch (mode)
      aesq2 = 2.0 * aewald * aewald
      aesq2n = 0.0d0
      if (aewald .gt. 0.0d0)  aesq2n = 1.0d0 / (sqrtpi*aewald)
      nlocal = 0
      toffset0 = 0
      maxlocal = int(dble(npole)*dble(maxelst)/dble(nthread))
c
c     perform dynamic allocation of some local arrays
c
      allocate (toffset(0:nthread-1))
      allocate (pscale(n))
      allocate (dscale(n))
      allocate (uscale(n))
      allocate (fieldt(3,npole))
      allocate (fieldtp(3,npole))
c
c     set arrays needed to scale connected atom interactions
c
      do i = 1, n
         pscale(i) = 1.0d0
         dscale(i) = 1.0d0
         uscale(i) = 1.0d0
      end do
c
c     set OpenMP directives for the major loop structure
c
!$OMP PARALLEL default(private) shared(n,npole,ipole,x,y,z,pdamp,thole,
!$OMP& rpole,p2scale,p3scale,p4scale,p41scale,p5scale,d1scale,d2scale,
!$OMP& d3scale,d4scale,u1scale,u2scale,u3scale,u4scale,n12,i12,n13,i13,
!$OMP& n14,i14,n15,i15,np11,ip11,np12,ip12,np13,ip13,np14,ip14,nelst,
!$OMP& elst,cut2,aewald,aesq2,aesq2n,poltyp,ntpair,tindex,tdipdip,
!$OMP& toffset,toffset0,field,fieldp,fieldt,fieldtp,maxlocal)
!$OMP& firstprivate(pscale,dscale,uscale,nlocal)
c
c     perform dynamic allocation of some local arrays
c
      if (poltyp .eq. 'MUTUAL') then
         allocate (ilocal(2,maxlocal))
         allocate (dlocal(6,maxlocal))
      end if
c
c     initialize local variables for OpenMP calculation
c
!$OMP DO collapse(2)
      do i = 1, npole
         do j = 1, 3
            fieldt(j,i) = 0.0d0
            fieldtp(j,i) = 0.0d0
         end do
      end do
!$OMP END DO
c
c     compute the real space portion of the Ewald summation
c
!$OMP DO reduction(+:fieldt,fieldtp) schedule(guided)
      do i = 1, npole
         ii = ipole(i)
         pdi = pdamp_ct(i)
         pti = thole_ct(i)
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
            uscale(ip11(j,ii)) = u1scale
         end do
         do j = 1, np12(ii)
            dscale(ip12(j,ii)) = d2scale
            uscale(ip12(j,ii)) = u2scale
         end do
         do j = 1, np13(ii)
            dscale(ip13(j,ii)) = d3scale
            uscale(ip13(j,ii)) = u3scale
         end do
         do j = 1, np14(ii)
            dscale(ip14(j,ii)) = d4scale
            uscale(ip14(j,ii)) = u4scale
         end do
         do kkk = 1, nelst(i)
            k = elst(kkk,i)
            kk = ipole(k)
            xr = x(kk) - x(ii)
            yr = y(kk) - y(ii)
            zr = z(kk) - z(ii)
            call image (xr,yr,zr)
            r2 = xr*xr + yr* yr + zr*zr
            if (r2 .le. cut2) then
               r = sqrt(r2)
               rr1 = 1.0d0 / r
               rr2 = rr1 * rr1
               rr3 = rr2 * rr1
               rr5 = rr2 * rr3
               rr7 = rr2 * rr5
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
c
c     calculate the error function damping factors
c
               ralpha = aewald * r
               bn(0) = erfc(ralpha) * rr1
               exp2a = exp(-ralpha**2)
               aefac = aesq2n
               do j = 1, 3
                  bfac = dble(j+j-1)
                  aefac = aesq2 * aefac
                  bn(j) = (bfac*bn(j-1)+aefac*exp2a) * rr2
               end do
c
c     compute the polarization damping scale factors
c
               scale3 = 1.0d0
               scale5 = 1.0d0
               scale7 = 1.0d0
               damp = pdi * pdamp_ct(k)
               damp_beta = pti_beta_ct * beta_damp_ct(k)
               if (damp .ne. 0.0d0) then
                  pgamma = min(pti,thole_ct(k))
                  damp = -pgamma * (r/damp)**3
                  if (damp .gt. -50.0d0) then
                     expdamp = exp(damp)
                     scale3 = 1.0d0 - expdamp*damp_beta
                     scale5 = 1.0d0 - expdamp*damp_beta*(1.0d0-damp)
                     scale7 = 1.0d0 - expdamp*damp_beta
     &                           *(1.0d0-damp+0.6d0*damp**2)
                  end if
               end if
c
c     find the field terms for the current interaction
c
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
               bcn(1) = bn(1) - (1.0d0-scale3*dscale(kk))*rr3
               bcn(2) = bn(2) - 3.0d0*(1.0d0-scale5*dscale(kk))*rr5
               bcn(3) = bn(3) - 15.0d0*(1.0d0-scale7*dscale(kk))*rr7
               fimd(1) = -xr*(bcn(1)*ck-bcn(2)*dkr+bcn(3)*qkr)
     &                     - bcn(1)*dkx + 2.0d0*bcn(2)*qkx
               fimd(2) = -yr*(bcn(1)*ck-bcn(2)*dkr+bcn(3)*qkr)
     &                     - bcn(1)*dky + 2.0d0*bcn(2)*qky
               fimd(3) = -zr*(bcn(1)*ck-bcn(2)*dkr+bcn(3)*qkr)
     &                     - bcn(1)*dkz + 2.0d0*bcn(2)*qkz
               fkmd(1) = xr*(bcn(1)*ci+bcn(2)*dir+bcn(3)*qir)
     &                     - bcn(1)*dix - 2.0d0*bcn(2)*qix
               fkmd(2) = yr*(bcn(1)*ci+bcn(2)*dir+bcn(3)*qir)
     &                     - bcn(1)*diy - 2.0d0*bcn(2)*qiy
               fkmd(3) = zr*(bcn(1)*ci+bcn(2)*dir+bcn(3)*qir)
     &                     - bcn(1)*diz - 2.0d0*bcn(2)*qiz
               bcn(1) = bn(1) - (1.0d0-scale3*pscale(kk))*rr3
               bcn(2) = bn(2) - 3.0d0*(1.0d0-scale5*pscale(kk))*rr5
               bcn(3) = bn(3) - 15.0d0*(1.0d0-scale7*pscale(kk))*rr7
               fimp(1) = -xr*(bcn(1)*ck-bcn(2)*dkr+bcn(3)*qkr)
     &                     - bcn(1)*dkx + 2.0d0*bcn(2)*qkx
               fimp(2) = -yr*(bcn(1)*ck-bcn(2)*dkr+bcn(3)*qkr)
     &                     - bcn(1)*dky + 2.0d0*bcn(2)*qky
               fimp(3) = -zr*(bcn(1)*ck-bcn(2)*dkr+bcn(3)*qkr)
     &                     - bcn(1)*dkz + 2.0d0*bcn(2)*qkz
               fkmp(1) = xr*(bcn(1)*ci+bcn(2)*dir+bcn(3)*qir)
     &                     - bcn(1)*dix - 2.0d0*bcn(2)*qix
               fkmp(2) = yr*(bcn(1)*ci+bcn(2)*dir+bcn(3)*qir)
     &                     - bcn(1)*diy - 2.0d0*bcn(2)*qiy
               fkmp(3) = zr*(bcn(1)*ci+bcn(2)*dir+bcn(3)*qir)
     &                     - bcn(1)*diz - 2.0d0*bcn(2)*qiz
c
c     find terms needed later to compute mutual polarization
c
               if (poltyp .eq. 'MUTUAL') then
                  bcn(1) = bn(1) - (1.0d0-scale3*uscale(kk))*rr3
                  bcn(2) = bn(2) - 3.0d0*(1.0d0-scale5*uscale(kk))*rr5
                  nlocal = nlocal + 1
                  ilocal(1,nlocal) = i
                  ilocal(2,nlocal) = k
                  dlocal(1,nlocal) = -bcn(1) + bcn(2)*xr*xr
                  dlocal(2,nlocal) = bcn(2)*xr*yr
                  dlocal(3,nlocal) = bcn(2)*xr*zr
                  dlocal(4,nlocal) = -bcn(1) + bcn(2)*yr*yr
                  dlocal(5,nlocal) = bcn(2)*yr*zr
                  dlocal(6,nlocal) = -bcn(1) + bcn(2)*zr*zr
               end if
c
c     increment the field at each site due to this interaction
c
               do j = 1, 3
                  fieldt(j,i) = fieldt(j,i) + fimd(j)
                  fieldt(j,k) = fieldt(j,k) + fkmd(j)
                  fieldtp(j,i) = fieldtp(j,i) + fimp(j)
                  fieldtp(j,k) = fieldtp(j,k) + fkmp(j)
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
         do j = 1, np11(ii)
            uscale(ip11(j,ii)) = 1.0d0
            dscale(ip11(j,ii)) = 1.0d0
         end do
         do j = 1, np12(ii)
            uscale(ip12(j,ii)) = 1.0d0
            dscale(ip12(j,ii)) = 1.0d0
         end do
         do j = 1, np13(ii)
            uscale(ip13(j,ii)) = 1.0d0
            dscale(ip13(j,ii)) = 1.0d0
         end do
         do j = 1, np14(ii)
            uscale(ip14(j,ii)) = 1.0d0
            dscale(ip14(j,ii)) = 1.0d0
         end do
      end do
!$OMP END DO
c
c     transfer the results from local to global arrays
c
!$OMP DO
      do i = 1, npole
         do j = 1, 3
            field(j,i) = fieldt(j,i) + field(j,i)
            fieldp(j,i) = fieldtp(j,i) + fieldp(j,i)
         end do
      end do
!$OMP END DO
c
c     store terms needed later to compute mutual polarization
c
!$OMP CRITICAL
      tid = 0
!$    tid = omp_get_thread_num ()
      toffset(tid) = toffset0
      toffset0 = toffset0 + nlocal
      ntpair = toffset0
!$OMP END CRITICAL
      if (poltyp .eq. 'MUTUAL') then
         k = toffset(tid)
         do i = 1, nlocal
            m = k + i
            tindex(1,m) = ilocal(1,i)
            tindex(2,m) = ilocal(2,i)
            do j = 1, 6
               tdipdip(j,m) = dlocal(j,i)
            end do
         end do
         deallocate (ilocal)
         deallocate (dlocal)
      end if
!$OMP END PARALLEL
c
c     perform deallocation of some local arrays
c
      deallocate (toffset)
      deallocate (pscale)
      deallocate (dscale)
      deallocate (uscale)
      deallocate (fieldt)
      deallocate (fieldtp)
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine umutual1  --  Ewald recip mutual induced field  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "umutual1" computes the reciprocal space contribution of the
c     induced atomic dipole moments to the field
c
c
      subroutine umutual1_ct (field_ct,fieldp_ct)
      use sizes
      use boxes
      use ewald
      use math
      use mpole
      use pme
      use polar
      use chargetransfer
      implicit none
      integer i,j,k
      real*8 term
      real*8 a(3,3)
      real*8 field_ct(3,*)
      real*8 fieldp_ct(3,*)
      real*8, allocatable :: fuind_ct(:,:)
      real*8, allocatable :: fuinp_ct(:,:)
      real*8, allocatable :: fdip_phi1(:,:)
      real*8, allocatable :: fdip_phi2(:,:)
      real*8, allocatable :: fdip_sum_phi(:,:)
      real*8, allocatable :: dipfield1(:,:)
      real*8, allocatable :: dipfield2(:,:)
c
c
c     return if the Ewald coefficient is zero
c
      if (aewald .lt. 1.0d-6)  return
c
c     perform dynamic allocation of some local arrays
c
      allocate (fuind_ct(3,npole))
      allocate (fuinp_ct(3,npole))
      allocate (fdip_phi1(10,npole))
      allocate (fdip_phi2(10,npole))
      allocate (fdip_sum_phi(20,npole))
      allocate (dipfield1(3,npole))
      allocate (dipfield2(3,npole))
c
c     convert Cartesian dipoles to fractional coordinates
c
      do i = 1, 3
         a(1,i) = dble(nfft1) * recip(i,1)
         a(2,i) = dble(nfft2) * recip(i,2)
         a(3,i) = dble(nfft3) * recip(i,3)
      end do
      do i = 1, npole
         do k = 1, 3
            fuind_ct(k,i) = a(k,1)*uind_ct(1,i) + a(k,2)*uind_ct(2,i)
     &                      + a(k,3)*uind_ct(3,i)
            fuinp_ct(k,i) = a(k,1)*uinp_ct(1,i) + a(k,2)*uinp_ct(2,i)
     &                      + a(k,3)*uinp_ct(3,i)
         end do
      end do
c
c     assign PME grid and perform 3-D FFT forward transform
c
      call grid_uind (fuind_ct,fuinp_ct)
      call fftfront
c
c     complete the transformation of the PME grid
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
c     perform 3-D FFT backward transform and get field
c
      call fftback
      call fphi_uind (fdip_phi1,fdip_phi2,fdip_sum_phi)
c
c     convert the dipole fields from fractional to Cartesian
c
      do i = 1, 3
         a(i,1) = dble(nfft1) * recip(i,1)
         a(i,2) = dble(nfft2) * recip(i,2)
         a(i,3) = dble(nfft3) * recip(i,3)
      end do
      do i = 1, npole
         do k = 1, 3
            dipfield1(k,i) = a(k,1)*fdip_phi1(2,i)
     &                          + a(k,2)*fdip_phi1(3,i)
     &                          + a(k,3)*fdip_phi1(4,i)
            dipfield2(k,i) = a(k,1)*fdip_phi2(2,i)
     &                          + a(k,2)*fdip_phi2(3,i)
     &                          + a(k,3)*fdip_phi2(4,i)
         end do
      end do
c
c     increment the field at each multipole site
c
      do i = 1, npole
         do k = 1, 3
            field_ct(k,i) = field_ct(k,i) - dipfield1(k,i)
            fieldp_ct(k,i) = fieldp_ct(k,i) - dipfield2(k,i)
         end do
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (fuind_ct)
      deallocate (fuinp_ct)
      deallocate (fdip_phi1)
      deallocate (fdip_phi2)
      deallocate (fdip_sum_phi)
      deallocate (dipfield1)
      deallocate (dipfield2)
      return
      end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine umutual2a  --  Ewald real mutual field via loop  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "umutual2a" computes the real space contribution of the induced
c     atomic dipole moments to the field via a double loop
c
c
      subroutine umutual2a_ct (field_ct,fieldp_ct)
      use sizes
      use atoms
      use boxes
      use bound
      use cell
      use couple
      use ewald
      use math
      use mpole
      use polar
      use polgrp
      use polpot
      use shunt
      use units
      use chargetransfer
      implicit none
      integer i,j,k,m
      integer ii,kk
      real*8 xr,yr,zr,r,r2
      real*8 rr1,rr2,rr3,rr5
      real*8 erfc,bfac,exp2a
      real*8 duir,dukr
      real*8 puir,pukr
      real*8 duix,duiy,duiz
      real*8 puix,puiy,puiz
      real*8 dukx,duky,dukz
      real*8 pukx,puky,pukz
      real*8 ralpha,aefac
      real*8 aesq2,aesq2n
      real*8 pdi_ct,pti_ct,pgamma_ct,pti_beta_ct 
      real*8 damp_ct,expdamp_ct,damp_beta_ct
      real*8 scale3_ct,scale5_ct
      real*8 bn(0:2)
      real*8 fimd_ct(3),fkmd_ct(3)
      real*8 fimp_ct(3),fkmp_ct(3)
      real*8, allocatable :: uscale_ct(:)
      real*8 field_ct(3,*)
      real*8 fieldp_ct(3,*)
      character*6 mode
      external erfc
c
c
c     check for multipoles and set cutoff coefficients
c
      if (npole .eq. 0)  return
      mode = 'EWALD'
      call switch (mode)
      aesq2 = 2.0 * aewald * aewald
      aesq2n = 0.0d0
      if (aewald .gt. 0.0d0)  aesq2n = 1.0d0 / (sqrtpi*aewald)
c
c     perform dynamic allocation of some local arrays
c
      allocate (uscale_ct(n))
c
c     set array needed to scale connected atom interactions
c
      do i = 1, n
         uscale_ct(i) = 1.0d0
      end do
c
c     compute the real space portion of the Ewald summation
c
      do i = 1, npole-1
         ii = ipole(i)
         pdi_ct = pdamp_ct(i)
         pti_ct = thole_ct(i)
         pti_beta_ct = beta_damp_ct(i)
         duix = uind_ct(1,i)
         duiy = uind_ct(2,i)
         duiz = uind_ct(3,i)
         puix = uinp_ct(1,i)
         puiy = uinp_ct(2,i)
         puiz = uinp_ct(3,i)
         do j = 1, np11(ii)
            uscale_ct(ip11(j,ii)) = u1scale_ct
         end do
         do j = 1, np12(ii)
            uscale_ct(ip12(j,ii)) = u2scale_ct
         end do
         do j = 1, np13(ii)
            uscale_ct(ip13(j,ii)) = u3scale_ct
         end do
         do j = 1, np14(ii)
            uscale_ct(ip14(j,ii)) = u4scale_ct
         end do
         do k = i+1, npole
            kk = ipole(k)
            xr = x(kk) - x(ii)
            yr = y(kk) - y(ii)
            zr = z(kk) - z(ii)
            call image (xr,yr,zr)
            r2 = xr*xr + yr* yr + zr*zr
            if (r2 .le. cut2) then
               r = sqrt(r2)
               rr1 = 1.0d0 / r
               rr2 = rr1 * rr1
               rr3 = rr2 * rr1
               rr5 = rr2 * rr3
               dukx = uind_ct(1,k)
               duky = uind_ct(2,k)
               dukz = uind_ct(3,k)
               pukx = uinp_ct(1,k)
               puky = uinp_ct(2,k)
               pukz = uinp_ct(3,k)
c
c     calculate the error function damping factors
c
               ralpha = aewald * r
               bn(0) = erfc(ralpha) * rr1
               exp2a = exp(-ralpha**2)
               aefac = aesq2n
               do j = 1, 2
                  bfac = dble(j+j-1)
                  aefac = aesq2 * aefac
                  bn(j) = (bfac*bn(j-1)+aefac*exp2a) * rr2
               end do
c
c     compute the polarization damping scale factors
c
               scale3_ct = uscale_ct(kk)
               scale5_ct = uscale_ct(kk)
               damp_ct = pdi_ct * pdamp_ct(k)
               damp_beta_ct = pti_beta_ct * beta_damp_ct(k)
               if (damp_ct .ne. 0.0d0) then
                  pgamma_ct = min(pti_ct,thole_ct(k))
                  damp_ct = -pgamma_ct * (r/damp_ct)**3
                  if (damp_ct .gt. -50.0d0) then
                     expdamp_ct = exp(damp_ct)
                 scale3_ct = scale3_ct * (1.0d0-expdamp_ct*damp_beta_ct)
                 scale5_ct = scale5_ct * 
     &                   (1.0d0-expdamp_ct*damp_beta_ct*(1.0d0-damp_ct))
                  end if
               end if
               bn(1) = -bn(1) + (1.0d0-scale3_ct)*rr3
               bn(2) = bn(2) - 3.0d0*(1.0d0-scale5_ct)*rr5
c
c     find the field terms for the current interaction
c
               duir = xr*duix + yr*duiy + zr*duiz
               dukr = xr*dukx + yr*duky + zr*dukz
               puir = xr*puix + yr*puiy + zr*puiz
               pukr = xr*pukx + yr*puky + zr*pukz
               fimd_ct(1) = bn(1)*dukx + bn(2)*dukr*xr
               fimd_ct(2) = bn(1)*duky + bn(2)*dukr*yr
               fimd_ct(3) = bn(1)*dukz + bn(2)*dukr*zr
               fkmd_ct(1) = bn(1)*duix + bn(2)*duir*xr
               fkmd_ct(2) = bn(1)*duiy + bn(2)*duir*yr
               fkmd_ct(3) = bn(1)*duiz + bn(2)*duir*zr
               fimp_ct(1) = bn(1)*pukx + bn(2)*pukr*xr
               fimp_ct(2) = bn(1)*puky + bn(2)*pukr*yr
               fimp_ct(3) = bn(1)*pukz + bn(2)*pukr*zr
               fkmp_ct(1) = bn(1)*puix + bn(2)*puir*xr
               fkmp_ct(2) = bn(1)*puiy + bn(2)*puir*yr
               fkmp_ct(3) = bn(1)*puiz + bn(2)*puir*zr
c
c     increment the field at each site due to this interaction
c
               do j = 1, 3
                  field_ct(j,i) = field_ct(j,i) + fimd_ct(j)
                  field_ct(j,k) = field_ct(j,k) + fkmd_ct(j)
                  fieldp_ct(j,i) = fieldp_ct(j,i) + fimp_ct(j)
                  fieldp_ct(j,k) = fieldp_ct(j,k) + fkmp_ct(j)
               end do
            end if
         end do
c
c     reset interaction scaling coefficients for connected atoms
c
         do j = 1, np11(ii)
            uscale_ct(ip11(j,ii)) = 1.0d0
         end do
         do j = 1, np12(ii)
            uscale_ct(ip12(j,ii)) = 1.0d0
         end do
         do j = 1, np13(ii)
            uscale_ct(ip13(j,ii)) = 1.0d0
         end do
         do j = 1, np14(ii)
            uscale_ct(ip14(j,ii)) = 1.0d0
         end do
      end do
c
c     periodic boundary for large cutoffs via replicates method
c
      if (use_replica) then
         do i = 1, npole
            ii = ipole(i)
            pdi_ct = pdamp_ct(i)
            pti_ct = thole_ct(i)
            pti_beta_ct = beta_damp_ct(i)
            duix = uind(1,i)
            duiy = uind(2,i)
            duiz = uind(3,i)
            puix = uinp(1,i)
            puiy = uinp(2,i)
            puiz = uinp(3,i)
            do j = 1, np11(ii)
               uscale_ct(ip11(j,ii)) = u1scale_ct
            end do
            do j = 1, np12(ii)
               uscale_ct(ip12(j,ii)) = u2scale_ct
            end do
            do j = 1, np13(ii)
               uscale_ct(ip13(j,ii)) = u3scale_ct
            end do
            do j = 1, np14(ii)
               uscale_ct(ip14(j,ii)) = u4scale_ct
            end do
            do k = i, npole
               kk = ipole(k)
               dukx = uind_ct(1,k)
               duky = uind_ct(2,k)
               dukz = uind_ct(3,k)
               pukx = uinp_ct(1,k)
               puky = uinp_ct(2,k)
               pukz = uinp_ct(3,k)
               do m = 1, ncell
                  xr = x(kk) - x(ii)
                  yr = y(kk) - y(ii)
                  zr = z(kk) - z(ii)
                  call imager (xr,yr,zr,m)
                  r2 = xr*xr + yr* yr + zr*zr
c
c     calculate the error function damping factors
c
                  if (r2 .le. cut2) then
                     r = sqrt(r2)
                     rr1 = 1.0d0 / r
                     rr2 = rr1 * rr1
                     rr3 = rr2 * rr1
                     rr5 = rr2 * rr3
                     ralpha = aewald * r
                     bn(0) = erfc(ralpha) * rr1
                     exp2a = exp(-ralpha**2)
                     aefac = aesq2n
                     do j = 1, 2
                        bfac = dble(j+j-1)
                        aefac = aesq2 * aefac
                        bn(j) = (bfac*bn(j-1)+aefac*exp2a) * rr2
                     end do
c
c     compute the polarization damping scale factors
c
                     scale3_ct = 1.0d0
                     scale5_ct = 1.0d0
                     damp_ct = pdi_ct * pdamp_ct(k)
                     damp_beta_ct = pti_beta_ct * beta_damp_ct(k)
                     if (damp_ct .ne. 0.0d0) then
                        pgamma_ct = min(pti_ct,thole_ct(k))
                        damp_ct = -pgamma_ct * (r/damp_ct)**3
                        if (damp_ct .gt. -50.0d0) then
                           expdamp_ct = exp(damp_ct)
                        scale3_ct = 1.0d0 - expdamp_ct*damp_beta_ct
                        scale5_ct = 1.0d0 - expdamp_ct*damp_beta_ct
     &                              *(1.0d0-damp_ct)
                        end if
                     end if
                     if (use_polymer) then
                        if (r2 .le. polycut2) then
                           scale3_ct = scale3_ct * uscale_ct(kk)
                           scale5_ct = scale5_ct * uscale_ct(kk)
                        end if
                     end if
                     bn(1) = -bn(1) + (1.0d0-scale3_ct)*rr3
                     bn(2) = bn(2) - 3.0d0*(1.0d0-scale5_ct)*rr5
c
c     find the field terms for the current interaction
c
                     duir = xr*duix + yr*duiy + zr*duiz
                     dukr = xr*dukx + yr*duky + zr*dukz
                     puir = xr*puix + yr*puiy + zr*puiz
                     pukr = xr*pukx + yr*puky + zr*pukz
                     fimd_ct(1) = bn(1)*dukx + bn(2)*dukr*xr
                     fimd_ct(2) = bn(1)*duky + bn(2)*dukr*yr
                     fimd_ct(3) = bn(1)*dukz + bn(2)*dukr*zr
                     fkmd_ct(1) = bn(1)*duix + bn(2)*duir*xr
                     fkmd_ct(2) = bn(1)*duiy + bn(2)*duir*yr
                     fkmd_ct(3) = bn(1)*duiz + bn(2)*duir*zr
                     fimp_ct(1) = bn(1)*pukx + bn(2)*pukr*xr
                     fimp_ct(2) = bn(1)*puky + bn(2)*pukr*yr
                     fimp_ct(3) = bn(1)*pukz + bn(2)*pukr*zr
                     fkmp_ct(1) = bn(1)*puix + bn(2)*puir*xr
                     fkmp_ct(2) = bn(1)*puiy + bn(2)*puir*yr
                     fkmp_ct(3) = bn(1)*puiz + bn(2)*puir*zr
c
c     increment the field at each site due to this interaction
c
                     do j = 1, 3
                        field_ct(j,i) = field_ct(j,i) + fimd_ct(j)
                        fieldp_ct(j,i) = fieldp_ct(j,i) + fimp_ct(j)
                        if (ii .ne. kk) then
                           field_ct(j,k) = field_ct(j,k) + fkmd_ct(j)
                           fieldp_ct(j,k) = fieldp_ct(j,k) + fkmp_ct(j)
                        end if
                     end do
                  end if
               end do
            end do
c
c     reset interaction scaling coefficients for connected atoms
c
            do j = 1, np11(ii)
               uscale_ct(ip11(j,ii)) = 1.0d0
            end do
            do j = 1, np12(ii)
               uscale_ct(ip12(j,ii)) = 1.0d0
            end do
            do j = 1, np13(ii)
               uscale_ct(ip13(j,ii)) = 1.0d0
            end do
            do j = 1, np14(ii)
               uscale_ct(ip14(j,ii)) = 1.0d0
            end do
         end do
      end if
c
c     perform deallocation of some local arrays
c
      deallocate (uscale_ct)
      return
      end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine umutual2b  --  Ewald real mutual field via list  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "umutual2b" computes the real space contribution of the induced
c     atomic dipole moments to the field via a neighbor list
c
c
      subroutine umutual2b_ct (field_ct,fieldp_ct)
      use sizes
      use mpole
      use polar
      use tarray
      use chargetransfer
      implicit none
      integer i,j,k,m
      real*8 fimd_ct(3),fkmd_ct(3)
      real*8 fimp_ct(3),fkmp_ct(3)
      real*8 field_ct(3,*)
      real*8 fieldp_ct(3,*)
      real*8, allocatable :: fieldt_ct(:,:)
      real*8, allocatable :: fieldtp_ct(:,:)
c
c
c     check for multipoles and set cutoff coefficients
c
      if (npole .eq. 0)  return
c
c     perform dynamic allocation of some local arrays
c
      allocate (fieldt_ct(3,npole))
      allocate (fieldtp_ct(3,npole))
c
c     set OpenMP directives for the major loop structure
c
c!$OMP PARALLEL default(private) shared(npole,uind_ct,uinp_ct,ntpair,tindex_ct,
c!$OMP& tdipdip,field_ct,fieldp_ct,fieldt_ct,fieldtp_ct)
c
c     initialize local variables for OpenMP calculation
c
c!$OMP DO collapse(2)
      do i = 1, npole
         do j = 1, 3
            fieldt_ct(j,i) = 0.0d0
            fieldtp_ct(j,i) = 0.0d0
         end do
      end do
c!$OMP END DO
c
c     find the field terms for each pairwise interaction
c
c!$OMP DO reduction(+:fieldt_ct,fieldtp_ct) schedule(guided)
      do m = 1, ntpair
         i = tindex(1,m)
         k = tindex(2,m)
       fimd_ct(1)= tdipdip(1,m)*uind_ct(1,k) + tdipdip(2,m)*uind_ct(2,k)
     &             + tdipdip(3,m)*uind_ct(3,k)
       fimd_ct(2)= tdipdip(2,m)*uind_ct(1,k) + tdipdip(4,m)*uind_ct(2,k)
     &             + tdipdip(5,m)*uind_ct(3,k)
       fimd_ct(3)= tdipdip(3,m)*uind_ct(1,k) + tdipdip(5,m)*uind_ct(2,k)
     &             + tdipdip(6,m)*uind_ct(3,k)
       fkmd_ct(1)= tdipdip(1,m)*uind_ct(1,i) + tdipdip(2,m)*uind_ct(2,i)
     &             + tdipdip(3,m)*uind_ct(3,i)
       fkmd_ct(2)= tdipdip(2,m)*uind_ct(1,i) + tdipdip(4,m)*uind_ct(2,i)
     &             + tdipdip(5,m)*uind_ct(3,i)
       fkmd_ct(3)= tdipdip(3,m)*uind_ct(1,i) + tdipdip(5,m)*uind_ct(2,i)
     &             + tdipdip(6,m)*uind_ct(3,i)
       fimp_ct(1)= tdipdip(1,m)*uinp_ct(1,k) + tdipdip(2,m)*uinp_ct(2,k)
     &             + tdipdip(3,m)*uinp_ct(3,k)
       fimp_ct(2)= tdipdip(2,m)*uinp_ct(1,k) + tdipdip(4,m)*uinp_ct(2,k)
     &             + tdipdip(5,m)*uinp_ct(3,k)
       fimp_ct(3)= tdipdip(3,m)*uinp_ct(1,k) + tdipdip(5,m)*uinp_ct(2,k)
     &             + tdipdip(6,m)*uinp_ct(3,k)
       fkmp_ct(1)= tdipdip(1,m)*uinp_ct(1,i) + tdipdip(2,m)*uinp_ct(2,i)
     &             + tdipdip(3,m)*uinp_ct(3,i)
       fkmp_ct(2)= tdipdip(2,m)*uinp_ct(1,i) + tdipdip(4,m)*uinp_ct(2,i)
     &             + tdipdip(5,m)*uinp_ct(3,i)
       fkmp_ct(3)= tdipdip(3,m)*uinp_ct(1,i) + tdipdip(5,m)*uinp_ct(2,i)
     &             + tdipdip(6,m)*uinp_ct(3,i)
c
c     increment the field at each site due to this interaction
c
         do j = 1, 3
            fieldt_ct(j,i) = fieldt_ct(j,i) + fimd_ct(j)
            fieldt_ct(j,k) = fieldt_ct(j,k) + fkmd_ct(j)
            fieldtp_ct(j,i) = fieldtp_ct(j,i) + fimp_ct(j)
            fieldtp_ct(j,k) = fieldtp_ct(j,k) + fkmp_ct(j)
         end do
      end do
c!$OMP END DO
c
c     end OpenMP directives for the major loop structure
c
c!$OMP DO
      do i = 1, npole
         do j = 1, 3
            field_ct(j,i) = fieldt_ct(j,i) + field_ct(j,i)
            fieldp_ct(j,i) = fieldtp_ct(j,i) + fieldp_ct(j,i)
         end do
      end do
c!$OMP END DO
c!$OMP END PARALLEL
c
c     perform deallocation of some local arrays
c
      deallocate (fieldt_ct)
      deallocate (fieldtp_ct)
      return
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine induce0d  --  Kirkwood SCRF induced dipoles  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "induce0d" computes the induced dipole moments at polarizable
c     sites for generalized Kirkwood SCRF and vacuum environments
c
c
      subroutine induce0d_ct
      use sizes
      use atoms
      use inform
      use iounit
      use mpole
      use polar
      use polpot
      use potent
      use units
      use uprior
      implicit none
      integer i,j,k,iter
      integer maxiter
      real*8 polmin
      real*8 eps,epsold
      real*8 epsd,epsp
      real*8 epsds,epsps
      real*8 udsum,upsum
      real*8 ussum,upssum
      real*8 a,ap,as,aps
      real*8 b,bp,bs,bps
      real*8 sum,sump
      real*8 sums,sumps
      real*8, allocatable :: poli(:)
      real*8, allocatable :: field(:,:)
      real*8, allocatable :: fieldp(:,:)
      real*8, allocatable :: fields(:,:)
      real*8, allocatable :: fieldps(:,:)
      real*8, allocatable :: udir(:,:)
      real*8, allocatable :: udirp(:,:)
      real*8, allocatable :: udirs(:,:)
      real*8, allocatable :: udirps(:,:)
      real*8, allocatable :: rsd(:,:)
      real*8, allocatable :: rsdp(:,:)
      real*8, allocatable :: rsds(:,:)
      real*8, allocatable :: rsdps(:,:)
      real*8, allocatable :: zrsd(:,:)
      real*8, allocatable :: zrsdp(:,:)
      real*8, allocatable :: zrsds(:,:)
      real*8, allocatable :: zrsdps(:,:)
      real*8, allocatable :: conj(:,:)
      real*8, allocatable :: conjp(:,:)
      real*8, allocatable :: conjs(:,:)
      real*8, allocatable :: conjps(:,:)
      real*8, allocatable :: vec(:,:)
      real*8, allocatable :: vecp(:,:)
      real*8, allocatable :: vecs(:,:)
      real*8, allocatable :: vecps(:,:)
      logical done
      character*6 mode
c
c
c     zero out the induced dipoles at each site; uind and uinp are
c     vacuum dipoles, uinds and uinps are SCRF dipoles
c
      do i = 1, npole
         do j = 1, 3
            uind(j,i) = 0.0d0
            uinp(j,i) = 0.0d0
            uinds(j,i) = 0.0d0
            uinps(j,i) = 0.0d0
         end do
      end do
      if (.not.use_polar .and. .not.use_solv)  return
c
c     set the switching function coefficients
c
      mode = 'MPOLE'
      call switch (mode)
c
c     perform dynamic allocation of some local arrays
c
      allocate (field(3,npole))
      allocate (fieldp(3,npole))
      allocate (fields(3,npole))
      allocate (fieldps(3,npole))
      allocate (udir(3,npole))
      allocate (udirp(3,npole))
      allocate (udirs(3,npole))
      allocate (udirps(3,npole))
c
c     compute the direct induced dipole moment at each atom, and
c     another set that also includes RF due to permanent multipoles
c
      call dfield0d_ct (field,fieldp,fields,fieldps)
c
c     set vacuum induced dipoles to polarizability times direct field;
c     set SCRF induced dipoles to polarizability times direct field
c     plus the GK reaction field due to permanent multipoles
c
      do i = 1, npole
         do j = 1, 3
            udir(j,i) = polarity(i) * field(j,i)
            udirp(j,i) = polarity(i) * fieldp(j,i)
            udirs(j,i) = polarity(i) * fields(j,i)
            udirps(j,i) = polarity(i) * fieldps(j,i)
            uind(j,i) = udir(j,i)
            uinp(j,i) = udirp(j,i)
            uinds(j,i) = udirs(j,i)
            uinps(j,i) = udirps(j,i)
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
            do i = 1, npole
               do j = 1, 3
                  udsum = 0.0d0
                  upsum = 0.0d0
                  ussum = 0.0d0
                  upssum = 0.0d0
                  do k = 1, nualt-1
                     udsum = udsum + bpred(k)*udalt(k,j,i)
                     upsum = upsum + bpredp(k)*upalt(k,j,i)
                     ussum = ussum + bpreds(k)*usalt(k,j,i)
                     upssum = upssum + bpredps(k)*upsalt(k,j,i)
                  end do
                  uind(j,i) = udsum
                  uinp(j,i) = upsum
                  uinds(j,i) = ussum
                  uinps(j,i) = upssum
               end do
            end do
         end if
c
c     perform dynamic allocation of some local arrays
c
         allocate (poli(npole))
         allocate (rsd(3,npole))
         allocate (rsdp(3,npole))
         allocate (rsds(3,npole))
         allocate (rsdps(3,npole))
         allocate (zrsd(3,npole))
         allocate (zrsdp(3,npole))
         allocate (zrsds(3,npole))
         allocate (zrsdps(3,npole))
         allocate (conj(3,npole))
         allocate (conjp(3,npole))
         allocate (conjs(3,npole))
         allocate (conjps(3,npole))
         allocate (vec(3,npole))
         allocate (vecp(3,npole))
         allocate (vecs(3,npole))
         allocate (vecps(3,npole))
c
c     set initial conjugate gradient residual and conjugate vector
c
         call ufield0d_ct (field,fieldp,fields,fieldps)
         do i = 1, npole
            poli(i) = max(polmin,polarity(i))
            do j = 1, 3
               rsd(j,i) = (udir(j,i)-uind(j,i))/poli(i)
     &                       + field(j,i)
               rsdp(j,i) = (udirp(j,i)-uinp(j,i))/poli(i)
     &                        + fieldp(j,i)
               rsds(j,i) = (udirs(j,i)-uinds(j,i))/poli(i)
     &                        + fields(j,i)
               rsdps(j,i) = (udirps(j,i)-uinps(j,i))/poli(i)
     &                         + fieldps(j,i)
               zrsd(j,i) = rsd(j,i) * poli(i)
               zrsdp(j,i) = rsdp(j,i) * poli(i)
               zrsds(j,i) = rsds(j,i) * poli(i)
               zrsdps(j,i) = rsdps(j,i) * poli(i)
               conj(j,i) = zrsd(j,i)
               conjp(j,i) = zrsdp(j,i)
               conjs(j,i) = zrsds(j,i)
               conjps(j,i) = zrsdps(j,i)
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
                  vecs(j,i) = uinds(j,i)
                  vecps(j,i) = uinps(j,i)
                  uind(j,i) = conj(j,i)
                  uinp(j,i) = conjp(j,i)
                  uinds(j,i) = conjs(j,i)
                  uinps(j,i) = conjps(j,i)
               end do
            end do
            call ufield0d_ct (field,fieldp,fields,fieldps)
            do i = 1, npole
               do j = 1, 3
                  uind(j,i) = vec(j,i)
                  uinp(j,i) = vecp(j,i)
                  uinds(j,i) = vecs(j,i)
                  uinps(j,i) = vecps(j,i)
                  vec(j,i) = conj(j,i)/poli(i) - field(j,i)
                  vecp(j,i) = conjp(j,i)/poli(i) - fieldp(j,i)
                  vecs(j,i) = conjs(j,i)/poli(i) - fields(j,i)
                  vecps(j,i) = conjps(j,i)/poli(i) - fieldps(j,i)
               end do
            end do
            a = 0.0d0
            ap = 0.0d0
            as = 0.0d0
            aps = 0.0d0
            sum = 0.0d0
            sump = 0.0d0
            sums = 0.0d0
            sumps = 0.0d0
            do i = 1, npole
               do j = 1, 3
                  a = a + conj(j,i)*vec(j,i)
                  ap = ap + conjp(j,i)*vecp(j,i)
                  as = as + conjs(j,i)*vecs(j,i)
                  aps = aps + conjps(j,i)*vecps(j,i)
                  sum = sum + rsd(j,i)*zrsd(j,i)
                  sump = sump + rsdp(j,i)*zrsdp(j,i)
                  sums = sums + rsds(j,i)*zrsds(j,i)
                  sumps = sumps + rsdps(j,i)*zrsdps(j,i)
               end do
            end do
            if (a .ne. 0.0d0)  a = sum / a
            if (ap .ne. 0.0d0)  ap = sump / ap
            if (as .ne. 0.0d0)  as = sums / as
            if (aps .ne. 0.0d0)  aps = sumps / aps
            do i = 1, npole
               do j = 1, 3
                  uind(j,i) = uind(j,i) + a*conj(j,i)
                  uinp(j,i) = uinp(j,i) + ap*conjp(j,i)
                  uinds(j,i) = uinds(j,i) + as*conjs(j,i)
                  uinps(j,i) = uinps(j,i) + aps*conjps(j,i)
                  rsd(j,i) = rsd(j,i) - a*vec(j,i)
                  rsdp(j,i) = rsdp(j,i) - ap*vecp(j,i)
                  rsds(j,i) = rsds(j,i) - as*vecs(j,i)
                  rsdps(j,i) = rsdps(j,i) - aps*vecps(j,i)
               end do
            end do
            b = 0.0d0
            bp = 0.0d0
            bs = 0.0d0
            bps = 0.0d0
            do i = 1, npole
               do j = 1, 3
                  zrsd(j,i) = rsd(j,i) * poli(i)
                  zrsdp(j,i) = rsdp(j,i) * poli(i)
                  zrsds(j,i) = rsds(j,i) * poli(i)
                  zrsdps(j,i) = rsdps(j,i) * poli(i)
                  b = b + rsd(j,i)*zrsd(j,i)
                  bp = bp + rsdp(j,i)*zrsdp(j,i)
                  bs = bs + rsds(j,i)*zrsds(j,i)
                  bps = bps + rsdps(j,i)*zrsdps(j,i)
               end do
            end do
            if (sum .ne. 0.0d0)  b = b / sum
            if (sump .ne. 0.0d0)  bp = bp / sump
            if (sums .ne. 0.0d0)  bs = bs / sums
            if (sumps .ne. 0.0d0)  bps = bps / sumps
            epsd = 0.0d0
            epsp = 0.0d0
            epsds = 0.0d0
            epsps = 0.0d0
            do i = 1, npole
               do j = 1, 3
                  conj(j,i) = zrsd(j,i) + b*conj(j,i)
                  conjp(j,i) = zrsdp(j,i) + bp*conjp(j,i)
                  conjs(j,i) = zrsds(j,i) + bs*conjs(j,i)
                  conjps(j,i) = zrsdps(j,i) + bps*conjps(j,i)
                  epsd = epsd + rsd(j,i)*rsd(j,i)
                  epsp = epsp + rsdp(j,i)*rsdp(j,i)
                  epsds = epsds + rsds(j,i)*rsds(j,i)
                  epsps = epsps + rsdps(j,i)*rsdps(j,i)
               end do
            end do
c
c     check the convergence of the mutual induced dipoles
c
            epsold = eps
            eps = max(epsd,epsp,epsds,epsps)
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
         deallocate (rsds)
         deallocate (rsdps)
         deallocate (zrsd)
         deallocate (zrsdp)
         deallocate (zrsds)
         deallocate (zrsdps)
         deallocate (conj)
         deallocate (conjp)
         deallocate (conjs)
         deallocate (conjps)
         deallocate (vec)
         deallocate (vecp)
         deallocate (vecs)
         deallocate (vecps)
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
   40       format (/,' INDUCE  --  Warning, Induced Charge Transfer',
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
      deallocate (fields)
      deallocate (fieldps)
      deallocate (udir)
      deallocate (udirp)
      deallocate (udirs)
      deallocate (udirps)
      return
      end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine dfield0d  --  generalized Kirkwood direct field  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "dfield0d" computes the direct electrostatic field due to
c     permanent multipole moments for use with with generalized
c     Kirkwood implicit solvation
c
c
      subroutine dfield0d_ct (field,fieldp,fields,fieldps)
      use sizes
      use atoms
      use couple
      use gkstuf
      use group
      use mpole
      use polar
      use polgrp
      use polpot
      use shunt
      use solute
      implicit none
      integer i,j,k,ii,kk
      real*8 xr,yr,zr
      real*8 xr2,yr2,zr2
      real*8 fgrp,r,r2
      real*8 rr3,rr5,rr7
      real*8 ci,uxi,uyi,uzi
      real*8 qxxi,qxyi,qxzi
      real*8 qyyi,qyzi,qzzi
      real*8 ck,uxk,uyk,uzk
      real*8 qxxk,qxyk,qxzk
      real*8 qyyk,qyzk,qzzk
      real*8 dir,dkr
      real*8 qix,qiy,qiz,qir
      real*8 qkx,qky,qkz,qkr
      real*8 damp,expdamp
      real*8 scale3,scale5
      real*8 scale7
      real*8 pdi,pti,pgamma
      real*8 rb2,rbi,rbk
      real*8 dwater,fc,fd,fq
      real*8 gf,gf2,gf3,gf5,gf7
      real*8 expterm,expc,expc1
      real*8 dexpc,expcdexpc
      real*8 a(0:3,0:2)
      real*8 gc(4),gux(10)
      real*8 guy(10),guz(10)
      real*8 gqxx(4),gqxy(4)
      real*8 gqxz(4),gqyy(4)
      real*8 gqyz(4),gqzz(4)
      real*8 fid(3),fkd(3)
      real*8, allocatable :: dscale(:)
      real*8, allocatable :: pscale(:)
      real*8 field(3,*)
      real*8 fieldp(3,*)
      real*8 fields(3,*)
      real*8 fieldps(3,*)
      real*8, allocatable :: fieldt(:,:)
      real*8, allocatable :: fieldtp(:,:)
      real*8, allocatable :: fieldts(:,:)
      real*8, allocatable :: fieldtps(:,:)
      logical proceed
c
c
c     zero out the value of the field at each site
c
      do i = 1, npole
         do j = 1, 3
            field(j,i) = 0.0d0
            fieldp(j,i) = 0.0d0
            fields(j,i) = 0.0d0
            fieldps(j,i) = 0.0d0
         end do
      end do
c
c     set dielectric constant and scaling factors for water
c
      dwater = 78.3d0
      fc = 1.0d0 * (1.0d0-dwater) / (1.0d0*dwater)
      fd = 2.0d0 * (1.0d0-dwater) / (1.0d0+2.0d0*dwater)
      fq = 3.0d0 * (1.0d0-dwater) / (2.0d0+3.0d0*dwater)
c
c     perform dynamic allocation of some local arrays
c
      allocate (dscale(n))
      allocate (pscale(n))
c
c     set arrays needed to scale connected atom interactions
c
      do i = 1, n
         dscale(i) = 1.0d0
         pscale(i) = 1.0d0
      end do
c
c     perform dynamic allocation of some local arrays
c
      allocate (fieldt(3,npole))
      allocate (fieldtp(3,npole))
      allocate (fieldts(3,npole))
      allocate (fieldtps(3,npole))
c
c     set OpenMP directives for the major loop structure
c
!$OMP PARALLEL default(private) shared(npole,ipole,pdamp,thole,rborn,
!$OMP& rpole,n12,n13,n14,n15,np11,np12,np13,np14,i12,i13,i14,i15,
!$OMP% ip11,ip12,ip13,ip14,p2scale,p3scale,p4scale,p41scale,p5scale,
!$OMP& d1scale,d2scale,d3scale,d4scale,use_intra,x,y,z,off2,fc,fd,fq,
!$OMP& gkc,field,fieldp,fields,fieldps)
!$OMP& firstprivate(dscale,pscale)
!$OMP% shared(fieldt,fieldtp,fieldts,fieldtps)
c
c     initialize local variables for OpenMP calculation
c
!$OMP DO collapse(2)
      do i = 1, npole
         do j = 1, 3
            fieldt(j,i) = 0.0d0
            fieldtp(j,i) = 0.0d0
            fieldts(j,i) = 0.0d0
            fieldtps(j,i) = 0.0d0
         end do
      end do
!$OMP END DO
c
c     find the field terms for each pairwise interaction
c
!$OMP DO reduction(+:fieldt,fieldtp,fieldts,fieldtps)
!$OMP& schedule(guided)
c
c     compute the direct induced dipole moment at each atom, and
c     another set that also includes RF due to permanent multipoles
c
      do i = 1, npole
         ii = ipole(i)
         pdi = pdamp(i)
         pti = thole(i)
         rbi = rborn(ii)
         ci = rpole(1,i)
         uxi = rpole(2,i)
         uyi = rpole(3,i)
         uzi = rpole(4,i)
         qxxi = rpole(5,i)
         qxyi = rpole(6,i)
         qxzi = rpole(7,i)
         qyyi = rpole(9,i)
         qyzi = rpole(10,i)
         qzzi = rpole(13,i)
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
         do k = i, npole
            kk = ipole(k)
            rbk = rborn(kk)
            proceed = .true.
            if (use_intra)  call groups (proceed,fgrp,ii,kk,0,0,0,0)
            if (proceed) then
               xr = x(kk) - x(ii)
               yr = y(kk) - y(ii)
               zr = z(kk) - z(ii)
               xr2 = xr * xr
               yr2 = yr * yr
               zr2 = zr * zr
               r2 = xr2 + yr2 + zr2
               if (r2 .le. off2) then
                  r = sqrt(r2)
                  ck = rpole(1,k)
                  uxk = rpole(2,k)
                  uyk = rpole(3,k)
                  uzk = rpole(4,k)
                  qxxk = rpole(5,k)
                  qxyk = rpole(6,k)
                  qxzk = rpole(7,k)
                  qyyk = rpole(9,k)
                  qyzk = rpole(10,k)
                  qzzk = rpole(13,k)
c
c     self-interactions for the solute field are skipped
c
                  if (i .ne. k) then
                     scale3 = 1.0d0
                     scale5 = 1.0d0
                     scale7 = 1.0d0
                     damp = pdi * pdamp(k)
                     if (damp .ne. 0.0d0) then
                        pgamma = min(pti,thole(k))
                        damp = -pgamma * (r/damp)**3
                        if (damp .gt. -50.0d0) then
                           expdamp = exp(damp)
                           scale3 = 1.0d0 - expdamp
                           scale5 = 1.0d0 - expdamp*(1.0d0-damp)
                           scale7 = 1.0d0 - expdamp
     &                                 *(1.0d0-damp+0.6d0*damp**2)
                        end if
                     end if
                     rr3 = scale3 / (r*r2)
                     rr5 = 3.0d0 * scale5 / (r*r2*r2)
                     rr7 = 15.0d0 * scale7 / (r*r2*r2*r2)
                     dir = uxi*xr + uyi*yr + uzi*zr
                     qix = qxxi*xr + qxyi*yr + qxzi*zr
                     qiy = qxyi*xr + qyyi*yr + qyzi*zr
                     qiz = qxzi*xr + qyzi*yr + qzzi*zr
                     qir = qix*xr + qiy*yr + qiz*zr
                     dkr = uxk*xr + uyk*yr + uzk*zr
                     qkx = qxxk*xr + qxyk*yr + qxzk*zr
                     qky = qxyk*xr + qyyk*yr + qyzk*zr
                     qkz = qxzk*xr + qyzk*yr + qzzk*zr
                     qkr = qkx*xr + qky*yr + qkz*zr
                     fid(1) = -xr*(rr3*ck-rr5*dkr+rr7*qkr)
     &                           - rr3*uxk + 2.0d0*rr5*qkx
                     fid(2) = -yr*(rr3*ck-rr5*dkr+rr7*qkr)
     &                           - rr3*uyk + 2.0d0*rr5*qky
                     fid(3) = -zr*(rr3*ck-rr5*dkr+rr7*qkr)
     &                           - rr3*uzk + 2.0d0*rr5*qkz
                     fkd(1) = xr*(rr3*ci+rr5*dir+rr7*qir)
     &                           - rr3*uxi - 2.0d0*rr5*qix
                     fkd(2) = yr*(rr3*ci+rr5*dir+rr7*qir)
     &                           - rr3*uyi - 2.0d0*rr5*qiy
                     fkd(3) = zr*(rr3*ci+rr5*dir+rr7*qir)
     &                           - rr3*uzi - 2.0d0*rr5*qiz
                     do j = 1, 3
                        fieldt(j,i) = fieldt(j,i) + fid(j)*dscale(kk)
                        fieldt(j,k) = fieldt(j,k) + fkd(j)*dscale(kk)
                        fieldtp(j,i) = fieldtp(j,i) + fid(j)*pscale(kk)
                        fieldtp(j,k) = fieldtp(j,k) + fkd(j)*pscale(kk)
                     end do
                  end if
c
c     set the reaction potential auxiliary terms
c
                  rb2 = rbi * rbk
                  expterm = exp(-r2/(gkc*rb2))
                  expc = expterm / gkc
                  dexpc = -2.0d0 / (gkc*rb2)
                  gf2 = 1.0d0 / (r2+rb2*expterm)
                  gf = sqrt(gf2)
                  gf3 = gf2 * gf
                  gf5 = gf3 * gf2
                  gf7 = gf5 * gf2
                  a(0,0) = gf
                  a(1,0) = -gf3
                  a(2,0) = 3.0d0 * gf5
                  a(3,0) = -15.0d0 * gf7
c
c     set the reaction potential gradient auxiliary terms
c
                  expc1 = 1.0d0 - expc
                  a(0,1) = expc1 * a(1,0)
                  a(1,1) = expc1 * a(2,0)
                  a(2,1) = expc1 * a(3,0)
c
c     dipole second reaction potential gradient auxiliary term
c
                  expcdexpc = -expc * dexpc
                  a(1,2) = expc1*a(2,1) + expcdexpc*a(2,0)
c
c     multiply the auxiliary terms by dielectric functions
c
                  a(0,1) = fc * a(0,1)
                  a(1,0) = fd * a(1,0)
                  a(1,1) = fd * a(1,1)
                  a(1,2) = fd * a(1,2)
                  a(2,0) = fq * a(2,0)
                  a(2,1) = fq * a(2,1)
c
c     unweighted dipole reaction potential tensor
c
                  gux(1) = xr * a(1,0)
                  guy(1) = yr * a(1,0)
                  guz(1) = zr * a(1,0)
c
c     unweighted reaction potential gradient tensor
c
                  gc(2) = xr * a(0,1)
                  gc(3) = yr * a(0,1)
                  gc(4) = zr * a(0,1)
                  gux(2) = a(1,0) + xr2*a(1,1)
                  gux(3) = xr * yr * a(1,1)
                  gux(4) = xr * zr * a(1,1)
                  guy(2) = gux(3)
                  guy(3) = a(1,0) + yr2*a(1,1)
                  guy(4) = yr * zr * a(1,1)
                  guz(2) = gux(4)
                  guz(3) = guy(4)
                  guz(4) = a(1,0) + zr2*a(1,1)
                  gqxx(2) = xr * (2.0d0*a(2,0)+xr2*a(2,1))
                  gqxx(3) = yr * xr2*a(2,1)
                  gqxx(4) = zr * xr2*a(2,1)
                  gqyy(2) = xr * yr2*a(2,1)
                  gqyy(3) = yr * (2.0d0*a(2,0)+yr2*a(2,1))
                  gqyy(4) = zr * yr2 * a(2,1)
                  gqzz(2) = xr * zr2 * a(2,1)
                  gqzz(3) = yr * zr2 * a(2,1)
                  gqzz(4) = zr * (2.0d0*a(2,0)+zr2*a(2,1))
                  gqxy(2) = yr * (a(2,0)+xr2*a(2,1))
                  gqxy(3) = xr * (a(2,0)+yr2*a(2,1))
                  gqxy(4) = zr * xr * yr * a(2,1)
                  gqxz(2) = zr * (a(2,0)+xr2*a(2,1))
                  gqxz(3) = gqxy(4)
                  gqxz(4) = xr * (a(2,0)+zr2*a(2,1))
                  gqyz(2) = gqxy(4)
                  gqyz(3) = zr * (a(2,0)+yr2*a(2,1))
                  gqyz(4) = yr * (a(2,0)+zr2*a(2,1))
c
c     unweighted dipole second reaction potential gradient tensor
c
                  gux(5) = xr * (3.0d0*a(1,1)+xr2*a(1,2))
                  gux(6) = yr * (a(1,1)+xr2*a(1,2))
                  gux(7) = zr * (a(1,1)+xr2*a(1,2))
                  gux(8) = xr * (a(1,1)+yr2*a(1,2))
                  gux(9) = zr * xr * yr * a(1,2)
                  gux(10) = xr * (a(1,1)+zr2*a(1,2))
                  guy(5) = yr * (a(1,1)+xr2*a(1,2))
                  guy(6) = xr * (a(1,1)+yr2*a(1,2))
                  guy(7) = gux(9)
                  guy(8) = yr * (3.0d0*a(1,1)+yr2*a(1,2))
                  guy(9) = zr * (a(1,1)+yr2*a(1,2))
                  guy(10) = yr * (a(1,1)+zr2*a(1,2))
                  guz(5) = zr * (a(1,1)+xr2*a(1,2))
                  guz(6) = gux(9)
                  guz(7) = xr * (a(1,1)+zr2*a(1,2))
                  guz(8) = zr * (a(1,1)+yr2*a(1,2))
                  guz(9) = yr * (a(1,1)+zr2*a(1,2))
                  guz(10) = zr * (3.0d0*a(1,1)+zr2*a(1,2))
c
c     generalized Kirkwood permanent reaction field
c
                  fid(1) = uxk*gux(2) + uyk*gux(3) + uzk*gux(4)
     &                        + 0.5d0 * (ck*gux(1) + qxxk*gux(5)
     &                            + qyyk*gux(8) + qzzk*gux(10)
     &                            + 2.0d0*(qxyk*gux(6)+qxzk*gux(7)
     &                                         +qyzk*gux(9)))
     &                        + 0.5d0 * (ck*gc(2) + qxxk*gqxx(2)
     &                            + qyyk*gqyy(2) + qzzk*gqzz(2)
     &                            + 2.0d0*(qxyk*gqxy(2)+qxzk*gqxz(2)
     &                                         +qyzk*gqyz(2)))
                  fid(2) = uxk*guy(2) + uyk*guy(3) + uzk*guy(4)
     &                        + 0.5d0 * (ck*guy(1) + qxxk*guy(5)
     &                            + qyyk*guy(8) + qzzk*guy(10)
     &                            + 2.0d0*(qxyk*guy(6)+qxzk*guy(7)
     &                                         +qyzk*guy(9)))
     &                        + 0.5d0 * (ck*gc(3) + qxxk*gqxx(3)
     &                            + qyyk*gqyy(3) + qzzk*gqzz(3)
     &                            + 2.0d0*(qxyk*gqxy(3)+qxzk*gqxz(3)
     &                                         +qyzk*gqyz(3)))
                  fid(3) = uxk*guz(2) + uyk*guz(3) + uzk*guz(4)
     &                        + 0.5d0 * (ck*guz(1) + qxxk*guz(5)
     &                            + qyyk*guz(8) + qzzk*guz(10)
     &                            + 2.0d0*(qxyk*guz(6)+qxzk*guz(7)
     &                                         +qyzk*guz(9)))
     &                        + 0.5d0 * (ck*gc(4) + qxxk*gqxx(4)
     &                            + qyyk*gqyy(4) + qzzk*gqzz(4)
     &                            + 2.0d0*(qxyk*gqxy(4)+qxzk*gqxz(4)
     &                                         +qyzk*gqyz(4)))
                  fkd(1) = uxi*gux(2) + uyi*gux(3) + uzi*gux(4)
     &                        - 0.5d0 * (ci*gux(1) + qxxi*gux(5)
     &                            + qyyi*gux(8) + qzzi*gux(10)
     &                            + 2.0d0*(qxyi*gux(6)+qxzi*gux(7)
     &                                         +qyzi*gux(9)))
     &                        - 0.5d0 * (ci*gc(2) + qxxi*gqxx(2)
     &                            + qyyi*gqyy(2) + qzzi*gqzz(2)
     &                            + 2.0d0*(qxyi*gqxy(2)+qxzi*gqxz(2)
     &                                         +qyzi*gqyz(2)))
                  fkd(2) = uxi*guy(2) + uyi*guy(3) + uzi*guy(4)
     &                        - 0.5d0 * (ci*guy(1) + qxxi*guy(5)
     &                            + qyyi*guy(8) + qzzi*guy(10)
     &                            + 2.0d0*(qxyi*guy(6)+qxzi*guy(7)
     &                                         +qyzi*guy(9)))
     &                        - 0.5d0 * (ci*gc(3) + qxxi*gqxx(3)
     &                            + qyyi*gqyy(3) + qzzi*gqzz(3)
     &                            + 2.0d0*(qxyi*gqxy(3)+qxzi*gqxz(3)
     &                                         +qyzi*gqyz(3)))
                  fkd(3) = uxi*guz(2) + uyi*guz(3) + uzi*guz(4)
     &                        - 0.5d0 * (ci*guz(1) + qxxi*guz(5)
     &                            + qyyi*guz(8) + qzzi*guz(10)
     &                            + 2.0d0*(qxyi*guz(6)+qxzi*guz(7)
     &                                         +qyzi*guz(9)))
     &                        - 0.5d0 * (ci*gc(4) + qxxi*gqxx(4)
     &                            + qyyi*gqyy(4) + qzzi*gqzz(4)
     &                            + 2.0d0*(qxyi*gqxy(4)+qxzi*gqxz(4)
     &                                         +qyzi*gqyz(4)))
c
c     scale the self-field by half, such that it sums to one below
c
                  if (i .eq. k) then
                     do j = 1, 3
                        fid(j) = 0.5d0 * fid(j)
                        fkd(j) = 0.5d0 * fkd(j)
                     end do
                  end if
                  do j = 1, 3
                     fieldts(j,i) = fieldts(j,i) + fid(j)
                     fieldts(j,k) = fieldts(j,k) + fkd(j)
                     fieldtps(j,i) = fieldtps(j,i) + fid(j)
                     fieldtps(j,k) = fieldtps(j,k) + fkd(j)
                  end do
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
!$OMP END DO
c
c     add local copies to global variables for OpenMP calculation
c
!$OMP DO
      do i = 1, npole
         do j = 1, 3
            field(j,i) = field(j,i) + fieldt(j,i)
            fieldp(j,i) = fieldp(j,i) + fieldtp(j,i)
            fields(j,i) = fields(j,i) + fieldts(j,i)
            fieldps(j,i) = fieldps(j,i) + fieldtps(j,i)
         end do
      end do
!$OMP END DO
c
c     combine permanent multipole field and GK reaction field
c
!$OMP DO
      do i = 1, npole
         do j = 1, 3
            fields(j,i) = field(j,i) + fields(j,i)
            fieldps(j,i) = fieldp(j,i) + fieldps(j,i)
         end do
      end do
!$OMP END DO
!$OMP END PARALLEL
c
c     perform deallocation of some local arrays
c
      deallocate (dscale)
      deallocate (pscale)
      deallocate (fieldt)
      deallocate (fieldtp)
      deallocate (fieldts)
      deallocate (fieldtps)
      return
      end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine ufield0d  --  generalized Kirkwood mutual field  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "ufield0d" computes the mutual electrostatic field due to
c     induced dipole moments for use with with generalized Kirkwood
c     implicit solvation
c
c
      subroutine ufield0d_ct (field,fieldp,fields,fieldps)
      use sizes
      use atoms
      use gkstuf
      use group
      use mpole
      use polar
      use polgrp
      use polpot
      use shunt
      use solute
      implicit none
      integer i,j,k,ii,kk
      real*8 xr,yr,zr
      real*8 xr2,yr2,zr2
      real*8 fgrp,r,r2
      real*8 rr3,rr5
      real*8 duix,duiy,duiz
      real*8 puix,puiy,puiz
      real*8 dukx,duky,dukz
      real*8 pukx,puky,pukz
      real*8 duir,dukr
      real*8 puir,pukr
      real*8 duixs,duiys,duizs
      real*8 puixs,puiys,puizs
      real*8 dukxs,dukys,dukzs
      real*8 pukxs,pukys,pukzs
      real*8 duirs,puirs
      real*8 dukrs,pukrs
      real*8 damp,expdamp
      real*8 scale3,scale5
      real*8 pdi,pti,pgamma
      real*8 rb2,rbi,rbk
      real*8 dwater,fd
      real*8 gf,gf2,gf3,gf5
      real*8 expterm,expc
      real*8 expc1,dexpc
      real*8 a(0:3,0:2)
      real*8 gux(10),guy(10)
      real*8 guz(10)
      real*8 fid(3),fkd(3)
      real*8 fip(3),fkp(3)
      real*8 fids(3),fkds(3)
      real*8 fips(3),fkps(3)
      real*8, allocatable :: dscale(:)
      real*8 field(3,*)
      real*8 fieldp(3,*)
      real*8 fields(3,*)
      real*8 fieldps(3,*)
      real*8, allocatable :: fieldt(:,:)
      real*8, allocatable :: fieldtp(:,:)
      real*8, allocatable :: fieldts(:,:)
      real*8, allocatable :: fieldtps(:,:)
      logical proceed
c
c
c     zero out the value of the field at each site
c
      do i = 1, npole
         do j = 1, 3
            field(j,i) = 0.0d0
            fieldp(j,i) = 0.0d0
            fields(j,i) = 0.0d0
            fieldps(j,i) = 0.0d0
         end do
      end do
c
c     set dielectric constant and scaling factor for water
c
      dwater = 78.3d0
      fd = 2.0d0 * (1.0d0-dwater) / (1.0d0+2.0d0*dwater)
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
c     perform dynamic allocation of some local arrays
c
      allocate (fieldt(3,npole))
      allocate (fieldtp(3,npole))
      allocate (fieldts(3,npole))
      allocate (fieldtps(3,npole))
c
c     set OpenMP directives for the major loop structure
c
!$OMP PARALLEL default(private) shared(npole,ipole,pdamp,thole,rborn,
!$OMP& uind,uinp,uinds,uinps,np11,np12,np13,np14,ip11,ip12,ip13,ip14,
!$OMP& u1scale,u2scale,u3scale,u4scale,use_intra,x,y,z,off2,fd,gkc,
!$OMP& field,fieldp,fields,fieldps)
!$OMP& firstprivate(dscale) shared(fieldt,fieldtp,fieldts,fieldtps)
c
c     initialize local variables for OpenMP calculation
c
!$OMP DO collapse(2)
      do i = 1, npole
         do j = 1, 3
            fieldt(j,i) = 0.0d0
            fieldtp(j,i) = 0.0d0
            fieldts(j,i) = 0.0d0
            fieldtps(j,i) = 0.0d0
         end do
      end do
!$OMP END DO
c
c     find the field terms for each pairwise interaction
c
!$OMP DO reduction(+:fieldt,fieldtp,fieldts,fieldtps)
!$OMP& schedule(guided)
c
c     compute the mutual electrostatic field at each atom,
c     and another field including RF due to induced dipoles
c
      do i = 1, npole
         ii = ipole(i)
         pdi = pdamp(i)
         pti = thole(i)
         rbi = rborn(ii)
         duix = uind(1,i)
         duiy = uind(2,i)
         duiz = uind(3,i)
         puix = uinp(1,i)
         puiy = uinp(2,i)
         puiz = uinp(3,i)
         duixs = uinds(1,i)
         duiys = uinds(2,i)
         duizs = uinds(3,i)
         puixs = uinps(1,i)
         puiys = uinps(2,i)
         puizs = uinps(3,i)
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
         do k = i, npole
            kk = ipole(k)
            rbk = rborn(kk)
            proceed = .true.
            if (use_intra)  call groups (proceed,fgrp,ii,kk,0,0,0,0)
            if (proceed) then
               xr = x(kk) - x(ii)
               yr = y(kk) - y(ii)
               zr = z(kk) - z(ii)
               xr2 = xr * xr
               yr2 = yr * yr
               zr2 = zr * zr
               r2 = xr2 + yr2 + zr2
               if (r2 .le. off2) then
                  r = sqrt(r2)
                  dukx = uind(1,k)
                  duky = uind(2,k)
                  dukz = uind(3,k)
                  pukx = uinp(1,k)
                  puky = uinp(2,k)
                  pukz = uinp(3,k)
                  dukxs = uinds(1,k)
                  dukys = uinds(2,k)
                  dukzs = uinds(3,k)
                  pukxs = uinps(1,k)
                  pukys = uinps(2,k)
                  pukzs = uinps(3,k)
                  if (i .ne. k) then
                     scale3 = dscale(kk)
                     scale5 = dscale(kk)
                     damp = pdi * pdamp(k)
                     if (damp .ne. 0.0d0) then
                        pgamma = min(pti,thole(k))
                        damp = -pgamma * (r/damp)**3
                        if (damp .gt. -50.0d0) then
                           expdamp = exp(damp)
                           scale3 = scale3 * (1.0d0-expdamp)
                           scale5 = scale5 * (1.0d0-(1.0d0-damp)
     &                                           *expdamp)
                        end if
                     end if
                     rr3 = scale3 / (r*r2)
                     rr5 = 3.0d0 * scale5 / (r*r2*r2)
                     duir = xr*duix + yr*duiy + zr*duiz
                     dukr = xr*dukx + yr*duky + zr*dukz
                     puir = xr*puix + yr*puiy + zr*puiz
                     pukr = xr*pukx + yr*puky + zr*pukz
                     duirs = xr*duixs + yr*duiys + zr*duizs
                     dukrs = xr*dukxs + yr*dukys + zr*dukzs
                     puirs = xr*puixs + yr*puiys + zr*puizs
                     pukrs = xr*pukxs + yr*pukys + zr*pukzs
                     fid(1) = -rr3*dukx + rr5*dukr*xr
                     fid(2) = -rr3*duky + rr5*dukr*yr
                     fid(3) = -rr3*dukz + rr5*dukr*zr
                     fkd(1) = -rr3*duix + rr5*duir*xr
                     fkd(2) = -rr3*duiy + rr5*duir*yr
                     fkd(3) = -rr3*duiz + rr5*duir*zr
                     fip(1) = -rr3*pukx + rr5*pukr*xr
                     fip(2) = -rr3*puky + rr5*pukr*yr
                     fip(3) = -rr3*pukz + rr5*pukr*zr
                     fkp(1) = -rr3*puix + rr5*puir*xr
                     fkp(2) = -rr3*puiy + rr5*puir*yr
                     fkp(3) = -rr3*puiz + rr5*puir*zr
                     fids(1) = -rr3*dukxs + rr5*dukrs*xr
                     fids(2) = -rr3*dukys + rr5*dukrs*yr
                     fids(3) = -rr3*dukzs + rr5*dukrs*zr
                     fkds(1) = -rr3*duixs + rr5*duirs*xr
                     fkds(2) = -rr3*duiys + rr5*duirs*yr
                     fkds(3) = -rr3*duizs + rr5*duirs*zr
                     fips(1) = -rr3*pukxs + rr5*pukrs*xr
                     fips(2) = -rr3*pukys + rr5*pukrs*yr
                     fips(3) = -rr3*pukzs + rr5*pukrs*zr
                     fkps(1) = -rr3*puixs + rr5*puirs*xr
                     fkps(2) = -rr3*puiys + rr5*puirs*yr
                     fkps(3) = -rr3*puizs + rr5*puirs*zr
                     do j = 1, 3
                        fieldt(j,i) = fieldt(j,i) + fid(j)
                        fieldt(j,k) = fieldt(j,k) + fkd(j)
                        fieldtp(j,i) = fieldtp(j,i) + fip(j)
                        fieldtp(j,k) = fieldtp(j,k) + fkp(j)
                        fieldts(j,i) = fieldts(j,i) + fids(j)
                        fieldts(j,k) = fieldts(j,k) + fkds(j)
                        fieldtps(j,i) = fieldtps(j,i) + fips(j)
                        fieldtps(j,k) = fieldtps(j,k) + fkps(j)
                     end do
                  end if
c
c     unweighted dipole reaction potential gradient tensor
c
                  rb2 = rbi * rbk
                  expterm = exp(-r2/(gkc*rb2))
                  expc = expterm / gkc
                  dexpc = -2.0d0 / (gkc*rbi*rbk)
                  gf2 = 1.0d0 / (r2+rb2*expterm)
                  gf = sqrt(gf2)
                  gf3 = gf2 * gf
                  gf5 = gf3 * gf2
                  a(1,0) = -gf3
                  a(2,0) = 3.0d0 * gf5
                  expc1 = 1.0d0 - expc
                  a(1,1) = expc1 * a(2,0)
                  gux(2) = fd * (a(1,0) + xr2*a(1,1))
                  gux(3) = fd * xr*yr*a(1,1)
                  gux(4) = fd * xr*zr*a(1,1)
                  guy(2) = gux(3)
                  guy(3) = fd * (a(1,0) + yr2*a(1,1))
                  guy(4) = fd * yr*zr*a(1,1)
                  guz(2) = gux(4)
                  guz(3) = guy(4)
                  guz(4) = fd * (a(1,0) + zr2*a(1,1))
                  fids(1) = dukxs*gux(2) + dukys*guy(2) + dukzs*guz(2)
                  fids(2) = dukxs*gux(3) + dukys*guy(3) + dukzs*guz(3)
                  fids(3) = dukxs*gux(4) + dukys*guy(4) + dukzs*guz(4)
                  fkds(1) = duixs*gux(2) + duiys*guy(2) + duizs*guz(2)
                  fkds(2) = duixs*gux(3) + duiys*guy(3) + duizs*guz(3)
                  fkds(3) = duixs*gux(4) + duiys*guy(4) + duizs*guz(4)
                  fips(1) = pukxs*gux(2) + pukys*guy(2) + pukzs*guz(2)
                  fips(2) = pukxs*gux(3) + pukys*guy(3) + pukzs*guz(3)
                  fips(3) = pukxs*gux(4) + pukys*guy(4) + pukzs*guz(4)
                  fkps(1) = puixs*gux(2) + puiys*guy(2) + puizs*guz(2)
                  fkps(2) = puixs*gux(3) + puiys*guy(3) + puizs*guz(3)
                  fkps(3) = puixs*gux(4) + puiys*guy(4) + puizs*guz(4)
                  if (i .eq. k) then
                     do j = 1, 3
                        fids(j) = 0.5d0 * fids(j)
                        fkds(j) = 0.5d0 * fkds(j)
                        fips(j) = 0.5d0 * fips(j)
                        fkps(j) = 0.5d0 * fkps(j)
                     end do
                  end if
                  do j = 1, 3
                     fieldts(j,i) = fieldts(j,i) + fids(j)
                     fieldts(j,k) = fieldts(j,k) + fkds(j)
                     fieldtps(j,i) = fieldtps(j,i) + fips(j)
                     fieldtps(j,k) = fieldtps(j,k) + fkps(j)
                  end do
               end if
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
!$OMP END DO
c
c     add local copies to global variables for OpenMP calculation
c
!$OMP DO
      do i = 1, npole
         do j = 1, 3
            field(j,i) = field(j,i) + fieldt(j,i)
            fieldp(j,i) = fieldp(j,i) + fieldtp(j,i)
            fields(j,i) = fields(j,i) + fieldts(j,i)
            fieldps(j,i) = fieldps(j,i) + fieldtps(j,i)
         end do
      end do
!$OMP END DO
!$OMP END PARALLEL
c
c     perform deallocation of some local arrays
c
      deallocate (dscale)
      deallocate (fieldt)
      deallocate (fieldtp)
      deallocate (fieldts)
      deallocate (fieldtps)
      return
      end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine induce0e  --  Poisson-Boltzmann induced dipoles  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "induce0e" computes the induced dipole moments at polarizable
c     sites for Poisson-Boltzmann SCRF and vacuum environments
c
c
      subroutine induce0e_ct
      use sizes
      use atoms
      use inform
      use iounit
      use mpole
      use polar
      use polpot
      use potent
      use units
      use uprior
      implicit none
      integer i,j,k,iter
      integer maxiter
      real*8 polmin
      real*8 eps,epsold
      real*8 epsd,epsp
      real*8 epsds,epsps
      real*8 udsum,upsum
      real*8 ussum,upssum
      real*8 a,ap,as,aps
      real*8 b,bp,bs,bps
      real*8 sum,sump
      real*8 sums,sumps
      real*8, allocatable :: poli(:)
      real*8, allocatable :: field(:,:)
      real*8, allocatable :: fieldp(:,:)
      real*8, allocatable :: fields(:,:)
      real*8, allocatable :: fieldps(:,:)
      real*8, allocatable :: udir(:,:)
      real*8, allocatable :: udirp(:,:)
      real*8, allocatable :: udirs(:,:)
      real*8, allocatable :: udirps(:,:)
      real*8, allocatable :: rsd(:,:)
      real*8, allocatable :: rsdp(:,:)
      real*8, allocatable :: rsds(:,:)
      real*8, allocatable :: rsdps(:,:)
      real*8, allocatable :: zrsd(:,:)
      real*8, allocatable :: zrsdp(:,:)
      real*8, allocatable :: zrsds(:,:)
      real*8, allocatable :: zrsdps(:,:)
      real*8, allocatable :: conj(:,:)
      real*8, allocatable :: conjp(:,:)
      real*8, allocatable :: conjs(:,:)
      real*8, allocatable :: conjps(:,:)
      real*8, allocatable :: vec(:,:)
      real*8, allocatable :: vecp(:,:)
      real*8, allocatable :: vecs(:,:)
      real*8, allocatable :: vecps(:,:)
      logical done
      character*6 mode
c
c
c     zero out the induced dipoles; uind and uinp are vacuum dipoles,
c     uinds and uinps are Poisson-Boltzmann SCRF dipoles
c
      do i = 1, npole
         do j = 1, 3
            uind(j,i) = 0.0d0
            uinp(j,i) = 0.0d0
            uinds(j,i) = 0.0d0
            uinps(j,i) = 0.0d0
         end do
      end do
      if (.not.use_polar .or. .not.use_solv)  return
c
c     set the switching function coefficients
c
      mode = 'MPOLE'
      call switch (mode)
c
c     perform dynamic allocation of some local arrays
c
      allocate (field(3,npole))
      allocate (fieldp(3,npole))
      allocate (fields(3,npole))
      allocate (fieldps(3,npole))
      allocate (udir(3,npole))
      allocate (udirp(3,npole))
      allocate (udirs(3,npole))
      allocate (udirps(3,npole))
c
c     compute the direct induced dipole moment at each atom, and
c     another set that also includes RF due to permanent multipoles
c
      call dfield0e_ct (field,fieldp,fields,fieldps)
c
c     set vacuum induced dipoles to polarizability times direct field;
c     SCRF induced dipoles are polarizability times direct field
c     plus the reaction field due to permanent multipoles
c
      do i = 1, npole
         do j = 1, 3
            udir(j,i) = polarity(i) * field(j,i)
            udirp(j,i) = polarity(i) * fieldp(j,i)
            udirs(j,i) = polarity(i) * fields(j,i)
            udirps(j,i) = polarity(i) * fieldps(j,i)
            uind(j,i) = udir(j,i)
            uinp(j,i) = udirp(j,i)
            uinds(j,i) = udirs(j,i)
            uinps(j,i) = udirps(j,i)
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
            do i = 1, npole
               do j = 1, 3
                  udsum = 0.0d0
                  upsum = 0.0d0
                  ussum = 0.0d0
                  upssum = 0.0d0
                  do k = 1, nualt-1
                     udsum = udsum + bpred(k)*udalt(k,j,i)
                     upsum = upsum + bpredp(k)*upalt(k,j,i)
                     ussum = ussum + bpreds(k)*usalt(k,j,i)
                     upssum = upssum + bpredps(k)*upsalt(k,j,i)
                  end do
                  uind(j,i) = udsum
                  uinp(j,i) = upsum
                  uinds(j,i) = ussum
                  uinps(j,i) = upssum
               end do
            end do
         end if
c
c     perform dynamic allocation of some local arrays
c
         allocate (poli(npole))
         allocate (rsd(3,npole))
         allocate (rsdp(3,npole))
         allocate (rsds(3,npole))
         allocate (rsdps(3,npole))
         allocate (zrsd(3,npole))
         allocate (zrsdp(3,npole))
         allocate (zrsds(3,npole))
         allocate (zrsdps(3,npole))
         allocate (conj(3,npole))
         allocate (conjp(3,npole))
         allocate (conjs(3,npole))
         allocate (conjps(3,npole))
         allocate (vec(3,npole))
         allocate (vecp(3,npole))
         allocate (vecs(3,npole))
         allocate (vecps(3,npole))
c
c     set initial conjugate gradient residual and conjugate vector
c
         call ufield0e_ct (field,fieldp,fields,fieldps)
         do i = 1, npole
            poli(i) = max(polmin,polarity(i))
            do j = 1, 3
               rsd(j,i) = (udir(j,i)-uind(j,i))/poli(i)
     &                       + field(j,i)
               rsdp(j,i) = (udirp(j,i)-uinp(j,i))/poli(i)
     &                        + fieldp(j,i)
               rsds(j,i) = (udirs(j,i)-uinds(j,i))/poli(i)
     &                        + fields(j,i)
               rsdps(j,i) = (udirps(j,i)-uinps(j,i))/poli(i)
     &                         + fieldps(j,i)
               zrsd(j,i) = rsd(j,i) * poli(i)
               zrsdp(j,i) = rsdp(j,i) * poli(i)
               zrsds(j,i) = rsds(j,i) * poli(i)
               zrsdps(j,i) = rsdps(j,i) * poli(i)
               conj(j,i) = zrsd(j,i)
               conjp(j,i) = zrsdp(j,i)
               conjs(j,i) = zrsds(j,i)
               conjps(j,i) = zrsdps(j,i)
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
                  vecs(j,i) = uinds(j,i)
                  vecps(j,i) = uinps(j,i)
                  uind(j,i) = conj(j,i)
                  uinp(j,i) = conjp(j,i)
                  uinds(j,i) = conjs(j,i)
                  uinps(j,i) = conjps(j,i)
               end do
            end do
            call ufield0e_ct (field,fieldp,fields,fieldps)
            do i = 1, npole
               do j = 1, 3
                  uind(j,i) = vec(j,i)
                  uinp(j,i) = vecp(j,i)
                  uinds(j,i) = vecs(j,i)
                  uinps(j,i) = vecps(j,i)
                  vec(j,i) = conj(j,i)/poli(i) - field(j,i)
                  vecp(j,i) = conjp(j,i)/poli(i) - fieldp(j,i)
                  vecs(j,i) = conjs(j,i)/poli(i) - fields(j,i)
                  vecps(j,i) = conjps(j,i)/poli(i) - fieldps(j,i)
               end do
            end do
            a = 0.0d0
            ap = 0.0d0
            as = 0.0d0
            aps = 0.0d0
            sum = 0.0d0
            sump = 0.0d0
            sums = 0.0d0
            sumps = 0.0d0
            do i = 1, npole
               do j = 1, 3
                  a = a + conj(j,i)*vec(j,i)
                  ap = ap + conjp(j,i)*vecp(j,i)
                  as = as + conjs(j,i)*vecs(j,i)
                  aps = aps + conjps(j,i)*vecps(j,i)
                  sum = sum + rsd(j,i)*zrsd(j,i)
                  sump = sump + rsdp(j,i)*zrsdp(j,i)
                  sums = sums + rsds(j,i)*zrsds(j,i)
                  sumps = sumps + rsdps(j,i)*zrsdps(j,i)
               end do
            end do
            if (a .ne. 0.0d0)  a = sum / a
            if (ap .ne. 0.0d0)  ap = sump / ap
            if (as .ne. 0.0d0)  as = sums / as
            if (aps .ne. 0.0d0)  aps = sumps / aps
            do i = 1, npole
               do j = 1, 3
                  uind(j,i) = uind(j,i) + a*conj(j,i)
                  uinp(j,i) = uinp(j,i) + ap*conjp(j,i)
                  uinds(j,i) = uinds(j,i) + as*conjs(j,i)
                  uinps(j,i) = uinps(j,i) + aps*conjps(j,i)
                  rsd(j,i) = rsd(j,i) - a*vec(j,i)
                  rsdp(j,i) = rsdp(j,i) - ap*vecp(j,i)
                  rsds(j,i) = rsds(j,i) - as*vecs(j,i)
                  rsdps(j,i) = rsdps(j,i) - aps*vecps(j,i)
               end do
            end do
            b = 0.0d0
            bp = 0.0d0
            bs = 0.0d0
            bps = 0.0d0
            do i = 1, npole
               do j = 1, 3
                  zrsd(j,i) = rsd(j,i) * poli(i)
                  zrsdp(j,i) = rsdp(j,i) * poli(i)
                  zrsds(j,i) = rsds(j,i) * poli(i)
                  zrsdps(j,i) = rsdps(j,i) * poli(i)
                  b = b + rsd(j,i)*zrsd(j,i)
                  bp = bp + rsdp(j,i)*zrsdp(j,i)
                  bs = bs + rsds(j,i)*zrsds(j,i)
                  bps = bps + rsdps(j,i)*zrsdps(j,i)
               end do
            end do
            if (sum .ne. 0.0d0)  b = b / sum
            if (sump .ne. 0.0d0)  bp = bp / sump
            if (sums .ne. 0.0d0)  bs = bs / sums
            if (sumps .ne. 0.0d0)  bps = bps / sumps
            epsd = 0.0d0
            epsp = 0.0d0
            epsds = 0.0d0
            epsps = 0.0d0
            do i = 1, npole
               do j = 1, 3
                  conj(j,i) = zrsd(j,i) + b*conj(j,i)
                  conjp(j,i) = zrsdp(j,i) + bp*conjp(j,i)
                  conjs(j,i) = zrsds(j,i) + bs*conjs(j,i)
                  conjps(j,i) = zrsdps(j,i) + bps*conjps(j,i)
                  epsd = epsd + rsd(j,i)*rsd(j,i)
                  epsp = epsp + rsdp(j,i)*rsdp(j,i)
                  epsds = epsds + rsds(j,i)*rsds(j,i)
                  epsps = epsps + rsdps(j,i)*rsdps(j,i)
               end do
            end do
c
c     check the convergence of the mutual induced dipoles
c
            epsold = eps
            eps = max(epsd,epsp,epsds,epsps)
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
         deallocate (rsds)
         deallocate (rsdps)
         deallocate (zrsd)
         deallocate (zrsdp)
         deallocate (zrsds)
         deallocate (zrsdps)
         deallocate (conj)
         deallocate (conjp)
         deallocate (conjs)
         deallocate (conjps)
         deallocate (vec)
         deallocate (vecp)
         deallocate (vecs)
         deallocate (vecps)
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
   40       format (/,' INDUCE  --  Warning, Induced Charge Transfer',
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
      deallocate (fields)
      deallocate (fieldps)
      deallocate (udir)
      deallocate (udirp)
      deallocate (udirs)
      deallocate (udirps)
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine dfield0e  --  Poisson-Boltzmann direct field  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "dfield0e" computes the direct electrostatic field due to
c     permanent multipole moments for use with in Poisson-Boltzmann
c
c
      subroutine dfield0e_ct (field,fieldp,fields,fieldps)
      use sizes
      use atoms
      use couple
      use group
      use mpole
      use pbstuf
      use polar
      use polgrp
      use polpot
      use shunt
      use solute
      implicit none
      integer i,j,k,ii,kk
      real*8 xr,yr,zr
      real*8 xr2,yr2,zr2
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
      real*8 field(3,*)
      real*8 fieldp(3,*)
      real*8 fields(3,*)
      real*8 fieldps(3,*)
      real*8, allocatable :: dscale(:)
      real*8, allocatable :: pscale(:)
      logical proceed
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
c     perform dynamic allocation of some local arrays
c
      allocate (dscale(n))
      allocate (pscale(n))
c
c     set arrays needed to scale connected atom interactions
c
      do i = 1, n
         pscale(i) = 1.0d0
         dscale(i) = 1.0d0
      end do
c
c     compute the direct electrostatic field at each atom, and
c     another field including RF due to permanent multipoles;
c     note self-interactions for the solute field are skipped
c
      do i = 1, npole
         ii = ipole(i)
         pdi = pdamp(i)
         pti = thole(i)
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
               xr2 = xr * xr
               yr2 = yr * yr
               zr2 = zr * zr
               r2 = xr2 + yr2 + zr2
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
                     pgamma = min(pti,thole(k))
                     damp = -pgamma * (r/damp)**3
                     if (damp .gt. -50.0d0) then
                        expdamp = exp(damp)
                        scale3 = 1.0d0 - expdamp
                        scale5 = 1.0d0 - expdamp*(1.0d0-damp)
                        scale7 = 1.0d0 - expdamp
     &                              *(1.0d0-damp+0.6d0*damp**2)
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
      deallocate (pscale)
c
c     find the Poisson-Boltzmann reaction field at each site
c
      call pbempole
c
c     combine permanent multipole field and PB reaction field
c
      do i = 1, npole
         ii = ipole(i)
         do j = 1, 3
            fields(j,i) = field(j,i) + pbep(j,ii)
            fieldps(j,i) = fieldp(j,i) + pbep(j,ii)
         end do
      end do
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine ufield0e  --  Poisson-Boltzmann mutual field  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "ufield0e" computes the mutual electrostatic field due to
c     induced dipole moments via a Poisson-Boltzmann solver
c
c
      subroutine ufield0e_ct (field,fieldp,fields,fieldps)
      use sizes
      use atoms
      use group
      use mpole
      use pbstuf
      use polar
      use polgrp
      use polpot
      use shunt
      use solute
      implicit none
      integer i,j,k,ii,kk
      real*8 xr,yr,zr
      real*8 xr2,yr2,zr2
      real*8 fgrp,r,r2
      real*8 rr3,rr5
      real*8 duix,duiy,duiz
      real*8 puix,puiy,puiz
      real*8 dukx,duky,dukz
      real*8 pukx,puky,pukz
      real*8 duir,puir
      real*8 dukr,pukr
      real*8 duixs,duiys,duizs
      real*8 puixs,puiys,puizs
      real*8 dukxs,dukys,dukzs
      real*8 pukxs,pukys,pukzs
      real*8 duirs,puirs
      real*8 dukrs,pukrs
      real*8 damp,expdamp
      real*8 scale3,scale5
      real*8 pdi,pti,pgamma
      real*8 fid(3),fkd(3)
      real*8 fip(3),fkp(3)
      real*8 fids(3),fkds(3)
      real*8 fips(3),fkps(3)
      real*8 field(3,*)
      real*8 fieldp(3,*)
      real*8 fields(3,*)
      real*8 fieldps(3,*)
      real*8, allocatable :: dscale(:)
      real*8, allocatable :: indpole(:,:)
      real*8, allocatable :: inppole(:,:)
      logical proceed
c
c
c     zero out the value of the field at each site
c
      do i = 1, npole
         do j = 1, 3
            field(j,i) = 0.0d0
            fieldp(j,i) = 0.0d0
            fields(j,i) = 0.0d0
            fieldps(j,i) = 0.0d0
         end do
      end do
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
c     compute the mutual electrostatic field at each atom,
c     and another field including RF due to induced dipoles
c
      do i = 1, npole
         ii = ipole(i)
         pdi = pdamp(i)
         pti = thole(i)
         duix = uind(1,i)
         duiy = uind(2,i)
         duiz = uind(3,i)
         puix = uinp(1,i)
         puiy = uinp(2,i)
         puiz = uinp(3,i)
         duixs = uinds(1,i)
         duiys = uinds(2,i)
         duizs = uinds(3,i)
         puixs = uinps(1,i)
         puiys = uinps(2,i)
         puizs = uinps(3,i)
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
            proceed = .true.
            if (use_intra)  call groups (proceed,fgrp,ii,kk,0,0,0,0)
            if (proceed) then
               xr = x(kk) - x(ii)
               yr = y(kk) - y(ii)
               zr = z(kk) - z(ii)
               xr2 = xr * xr
               yr2 = yr * yr
               zr2 = zr * zr
               r2 = xr2 + yr2 + zr2
               if (r2 .le. off2) then
                  r = sqrt(r2)
                  dukx = uind(1,k)
                  duky = uind(2,k)
                  dukz = uind(3,k)
                  pukx = uinp(1,k)
                  puky = uinp(2,k)
                  pukz = uinp(3,k)
                  dukxs = uinds(1,k)
                  dukys = uinds(2,k)
                  dukzs = uinds(3,k)
                  pukxs = uinps(1,k)
                  pukys = uinps(2,k)
                  pukzs = uinps(3,k)
                  scale3 = dscale(kk)
                  scale5 = dscale(kk)
                  damp = pdi * pdamp(k)
                  if (damp .ne. 0.0d0) then
                     pgamma = min(pti,thole(k))
                     damp = -pgamma * (r/damp)**3
                     if (damp .gt. -50.0d0) then
                        expdamp = exp(damp)
                        scale3 = scale3 * (1.0d0-expdamp)
                        scale5 = scale5 * (1.0d0-(1.0d0-damp)*expdamp)
                     end if
                  end if
                  rr3 = scale3 / (r*r2)
                  rr5 = 3.0d0 * scale5 / (r*r2*r2)
                  duir = xr*duix + yr*duiy + zr*duiz
                  dukr = xr*dukx + yr*duky + zr*dukz
                  puir = xr*puix + yr*puiy + zr*puiz
                  pukr = xr*pukx + yr*puky + zr*pukz
                  duirs = xr*duixs + yr*duiys + zr*duizs
                  dukrs = xr*dukxs + yr*dukys + zr*dukzs
                  puirs = xr*puixs + yr*puiys + zr*puizs
                  pukrs = xr*pukxs + yr*pukys + zr*pukzs
                  fid(1) = -rr3*dukx + rr5*dukr*xr
                  fid(2) = -rr3*duky + rr5*dukr*yr
                  fid(3) = -rr3*dukz + rr5*dukr*zr
                  fkd(1) = -rr3*duix + rr5*duir*xr
                  fkd(2) = -rr3*duiy + rr5*duir*yr
                  fkd(3) = -rr3*duiz + rr5*duir*zr
                  fip(1) = -rr3*pukx + rr5*pukr*xr
                  fip(2) = -rr3*puky + rr5*pukr*yr
                  fip(3) = -rr3*pukz + rr5*pukr*zr
                  fkp(1) = -rr3*puix + rr5*puir*xr
                  fkp(2) = -rr3*puiy + rr5*puir*yr
                  fkp(3) = -rr3*puiz + rr5*puir*zr
                  fids(1) = -rr3*dukxs + rr5*dukrs*xr
                  fids(2) = -rr3*dukys + rr5*dukrs*yr
                  fids(3) = -rr3*dukzs + rr5*dukrs*zr
                  fkds(1) = -rr3*duixs + rr5*duirs*xr
                  fkds(2) = -rr3*duiys + rr5*duirs*yr
                  fkds(3) = -rr3*duizs + rr5*duirs*zr
                  fips(1) = -rr3*pukxs + rr5*pukrs*xr
                  fips(2) = -rr3*pukys + rr5*pukrs*yr
                  fips(3) = -rr3*pukzs + rr5*pukrs*zr
                  fkps(1) = -rr3*puixs + rr5*puirs*xr
                  fkps(2) = -rr3*puiys + rr5*puirs*yr
                  fkps(3) = -rr3*puizs + rr5*puirs*zr
                  do j = 1, 3
                     field(j,i) = field(j,i) + fid(j)
                     field(j,k) = field(j,k) + fkd(j)
                     fieldp(j,i) = fieldp(j,i) + fip(j)
                     fieldp(j,k) = fieldp(j,k) + fkp(j)
                     fields(j,i) = fields(j,i) + fids(j)
                     fields(j,k) = fields(j,k) + fkds(j)
                     fieldps(j,i) = fieldps(j,i) + fips(j)
                     fieldps(j,k) = fieldps(j,k) + fkps(j)
                  end do
               end if
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
c
c     perform dynamic allocation of some global arrays
c
      if (.not. allocated(pbeuind))  allocate (pbeuind(3,n))
      if (.not. allocated(pbeuinp))  allocate (pbeuinp(3,n))
c
c     perform dynamic allocation of some local arrays
c
      allocate (indpole(3,n))
      allocate (inppole(3,n))
c
c     zero out the PB reaction field at each atomic site
c
      do i = 1, n
         do j = 1, 3
            indpole(j,i) = 0.0d0
            inppole(j,i) = 0.0d0
            pbeuind(j,i) = 0.0d0
            pbeuinp(j,i) = 0.0d0
         end do
      end do
c
c     find the Poisson-Boltzmann reaction field at each site
c
      do i = 1, npole
         ii = ipole(i)
         do j = 1, 3
            indpole(j,ii) = uinds(j,i)
            inppole(j,ii) = uinps(j,i)
         end do
      end do
      call apbsinduce (indpole,pbeuind)
      call apbsnlinduce (inppole,pbeuinp)
c
c     perform deallocation of some local arrays
c
      deallocate (indpole)
      deallocate (inppole)
c
c     combine mutual induced dipole field and PB reaction field
c
      do i = 1, npole
         ii = ipole(i)
         do j = 1, 3
            fields(j,i) = fields(j,i) + pbeuind(j,ii)
            fieldps(j,i) = fieldps(j,i) + pbeuinp(j,ii)
         end do
      end do
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine ulspred  --  induced dipole prediction coeffs  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "ulspred" uses standard extrapolation or a least squares fit
c     to set coefficients of an induced dipole predictor polynomial
c
c     literature references:
c
c     J. Kolafa, "Time-Reversible Always Stable Predictor-Corrector
c     Method for Molecular Dynamics of Polarizable Molecules", Journal
c     of Computational Chemistry, 25, 335-342 (2004)
c
c     W. Wang and R. D. Skeel, "Fast Evaluation of Polarizable Forces",
c     Journal of Chemical Physics, 123, 164107 (2005)
c
c
      subroutine ulspred_ct
      use sizes
      use mpole
      use uprior
      implicit none
      integer i,j,k,m
      real*8 coeff,udk,upk
      real*8 amax,apmax
      real*8 b_ct(maxualt_ct)
      real*8 bp_ct(maxualt_ct)
      real*8 a_ct(maxualt*(maxualt_ct+1)/2)
      real*8 ap_ct(maxualt*(maxualt_ct+1)/2)
      real*8 c_ct(maxualt_ct,maxualt_ct)
      real*8 cp_ct(maxualt_ct,maxualt_ct)
c
c
c     set the Gear predictor binomial coefficients
c
      if (polpred_ct .eq. 'GEAR') then
         do i = 1, nualt_ct
            coeff = gear_ct(i)
            bpred_ct(i) = coeff
            bpredp_ct(i) = coeff
            bpreds_ct(i) = coeff
            bpredps_ct(i) = coeff
         end do
c
c     set always stable predictor-corrector (ASPC) coefficients
c
      else if (polpred_ct .eq. 'ASPC') then
         do i = 1, nualt_ct
            coeff = aspc_ct(i)
            bpred_ct(i) = coeff
            bpredp_ct(i) = coeff
            bpreds_ct(i) = coeff
            bpredps_ct(i) = coeff
         end do
c
c     derive normal equations corresponding to least squares fit
c
      else
         do k = 1, nualt_ct
            b_ct(k) = 0.0d0
            bp_ct(k) = 0.0d0
            do m = k, nualt_ct
               c_ct(k,m) = 0.0d0
               cp_ct(k,m) = 0.0d0
            end do
         end do
         do i = 1, npole
            do j = 1, 3
               do k = 1, nualt_ct
                  udk = udalt_ct(k,j,i)
                  upk = upalt_ct(k,j,i)
                  do m = k, nualt
                     c_ct(k,m) = c_ct(k,m) + udk*udalt_ct(m,j,i)
                     cp_ct(k,m) = cp_ct(k,m) + upk*upalt_ct(m,j,i)
                  end do
               end do
            end do
         end do
         i = 0
         do k = 2, nualt_ct
            b_ct(k-1) = c_ct(1,k)
            bp_ct(k-1) = cp_ct(1,k)
            do m = k, nualt_ct
               i = i + 1
               a_ct(i) = c_ct(k,m)
               ap_ct(i) = cp_ct(k,m)
            end do
         end do
c
c     check for nonzero coefficients and solve normal equations
c
         k = nualt_ct - 1
         amax = 0.0d0
         apmax = 0.0d0
         do i = 1, k*(k+1)/2
            amax = max(amax,a_ct(i))
            apmax = max(apmax,ap_ct(i))
         end do
         if (amax .ne. 0.0d0)  call cholesky (k,a_ct,b_ct)
         if (apmax .ne. 0.0d0)  call cholesky (k,ap_ct,bp_ct)
c
c     transfer the final solution to the coefficient vector
c
         do k = 1, nualt_ct-1
            bpred_ct(k) = b_ct(k)
            bpredp_ct(k) = bp_ct(k)
            bpreds_ct(k) = b_ct(k)
            bpredps_ct(k) = bp_ct(k)
         end do
         bpred_ct(nualt_ct) = 0.0d0
         bpredp_ct(nualt_ct) = 0.0d0
         bpreds_ct(nualt_ct) = 0.0d0
         bpredps_ct(nualt_ct) = 0.0d0
      end if
      return
      end
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
      subroutine uscale0a_ct (mode,rsd_ct,rsdp_ct,zrsd_ct,zrsdp_ct)
      use sizes
      use atoms
      use limits
      use mpole
      use polar
      use polgrp
      use polpot
      use usolve
      use chargetransfer
      implicit none
      integer i,j,k,m
      integer ii,kk
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 r,r2,rr3,rr5
      real*8 pdi_ct,pti_ct,pti_beta_ct
      real*8 polmin_ct
      real*8 poli_ct,polik_ct
      real*8 damp_ct,expdamp_ct,damp_beta_ct
      real*8 pgamma_ct,off2
      real*8 scale3_ct,scale5_ct
      real*8 m1,m2,m3
      real*8 m4,m5,m6
      real*8, allocatable :: dscale_ct(:)
      real*8 rsd_ct(3,*)
      real*8 rsdp_ct(3,*)
      real*8 zrsd_ct(3,*)
      real*8 zrsdp_ct(3,*)
      character*6 mode
c
c
c     apply the preconditioning matrix to the current residual
c
      if (mode .eq. 'APPLY') then
c
c     use diagonal preconditioner elements as first approximation
c
         polmin_ct = 0.00000001d0
         do i = 1, npole
            poli_ct = udiag_ct * max(polmin_ct,charge_trans(i))
            do j = 1, 3
               zrsd_ct(j,i) = poli_ct * rsd_ct(j,i)
               zrsdp_ct(j,i) = poli_ct * rsdp_ct(j,i)
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
       zrsd_ct(1,i) = zrsd_ct(1,i) + m1*rsd_ct(1,k) + m2*rsd_ct(2,k)
     &                + m3*rsd_ct(3,k)
       zrsd_ct(2,i) = zrsd_ct(2,i) + m2*rsd_ct(1,k) + m4*rsd_ct(2,k)
     &                + m5*rsd_ct(3,k)
       zrsd_ct(3,i) = zrsd_ct(3,i) + m3*rsd_ct(1,k) + m5*rsd_ct(2,k)
     &                + m6*rsd_ct(3,k)
       zrsd_ct(1,k) = zrsd_ct(1,k) + m1*rsd_ct(1,i) + m2*rsd_ct(2,i)
     &                + m3*rsd_ct(3,i)
       zrsd_ct(2,k) = zrsd_ct(2,k) + m2*rsd_ct(1,i) + m4*rsd_ct(2,i)
     &                + m5*rsd_ct(3,i)
       zrsd_ct(3,k) = zrsd_ct(3,k) + m3*rsd_ct(1,i) + m5*rsd_ct(2,i)
     &                + m6*rsd_ct(3,i)
       zrsdp_ct(1,i) = zrsdp_ct(1,i) + m1*rsdp_ct(1,k) + m2*rsdp_ct(2,k)
     &                 + m3*rsdp_ct(3,k)
       zrsdp_ct(2,i) = zrsdp_ct(2,i) + m2*rsdp_ct(1,k) + m4*rsdp_ct(2,k)
     &                 + m5*rsdp_ct(3,k)
       zrsdp_ct(3,i) = zrsdp_ct(3,i) + m3*rsdp_ct(1,k) + m5*rsdp_ct(2,k)
     &                 + m6*rsdp_ct(3,k)
       zrsdp_ct(1,k) = zrsdp_ct(1,k) + m1*rsdp_ct(1,i) + m2*rsdp_ct(2,i)
     &                 + m3*rsdp_ct(3,i)
       zrsdp_ct(2,k) = zrsdp_ct(2,k) + m2*rsdp_ct(1,i) + m4*rsdp_ct(2,i)
     &                 + m5*rsdp_ct(3,i)
       zrsdp_ct(3,k) = zrsdp_ct(3,k) + m3*rsdp_ct(1,i) + m5*rsdp_ct(2,i)
     &                           + m6*rsdp_ct(3,i)
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
         allocate (dscale_ct(n))
c
c     set array needed to scale connected atom interactions
c
         do i = 1, n
            dscale_ct(i) = 1.0d0
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
            pdi_ct = pdamp_ct(i)
            pti_ct = thole_ct(i)
            poli_ct = charge_trans(i)
            pti_beta_ct = beta_damp_ct(i) 
            do j = i+1, npole
               dscale_ct(ipole(j)) = 1.0d0
            end do
            do j = 1, np11(ii)
               dscale_ct(ip11(j,ii)) = u1scale_ct
            end do
            do j = 1, np12(ii)
               dscale_ct(ip12(j,ii)) = u2scale_ct
            end do
            do j = 1, np13(ii)
               dscale_ct(ip13(j,ii)) = u3scale_ct
            end do
            do j = 1, np14(ii)
               dscale_ct(ip14(j,ii)) = u4scale_ct
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
                  scale3_ct = dscale_ct(kk)
                  scale5_ct = dscale_ct(kk)
                  damp_ct = pdi_ct * pdamp_ct(k)
                  damp_beta_ct = pti_beta_ct * beta_damp_ct(k)
                  if (damp_ct .ne. 0.0d0) then
                     pgamma_ct = min(pti_ct,thole_ct(k))
                     damp_ct = -pgamma_ct * (r/damp_ct)**3
                     if (damp_ct .gt. -50.0d0) then
                        expdamp_ct = exp(damp_ct)
                        scale3_ct = scale3_ct * expdamp_ct*damp_beta_ct
                        scale5_ct = scale5_ct * expdamp_ct*
     &                               (1.0d0-damp_ct)*damp_beta_ct
                     end if
                  end if
                  polik_ct = poli_ct * charge_trans(k)
                  rr3 = scale3_ct * polik_ct / (r*r2)
                  rr5 = 3.0d0 * scale5_ct * polik_ct / (r*r2*r2)
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
               dscale_ct(ip11(j,ii)) = 1.0d0
            end do
            do j = 1, np12(ii)
               dscale_ct(ip12(j,ii)) = 1.0d0
            end do
            do j = 1, np13(ii)
               dscale_ct(ip13(j,ii)) = 1.0d0
            end do
            do j = 1, np14(ii)
               dscale_ct(ip14(j,ii)) = 1.0d0
            end do
         end do
c
c     perform deallocation of some local arrays
c
         deallocate (dscale_ct)
      end if
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine uscale0b  --  dipole preconditioner via list  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "uscale0b" builds and applies a preconditioner for the conjugate
c     gradient induced dipole solver using a neighbor pair list
c
c
      subroutine uscale0b_ct (mode,rsd,rsdp,zrsd,zrsdp)
      use sizes
      use atoms
      use mpole
      use neigh
      use polar
      use polgrp
      use polpot
      use usolve
      use chargetransfer
      implicit none
      integer i,j,k,m
      integer ii,kk,kkk
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 r,r2,rr3,rr5
      real*8 pdi,pti,pti_beta_ct
      real*8 polmin
      real*8 poli,polik
      real*8 damp,expdamp,damp_beta
      real*8 pgamma
      real*8 scale3,scale5
      real*8 m1,m2,m3
      real*8 m4,m5,m6
      real*8, allocatable :: dscale(:)
      real*8 rsd(3,*)
      real*8 rsdp(3,*)
      real*8 zrsd(3,*)
      real*8 zrsdp(3,*)
      real*8, allocatable :: zrsdt(:,:)
      real*8, allocatable :: zrsdtp(:,:)
      character*6 mode
c
c
c     apply the preconditioning matrix to the current residual
c
      if (mode .eq. 'APPLY') then
c
c     perform dynamic allocation of some local arrays
c
         allocate (zrsdt(3,npole))
         allocate (zrsdtp(3,npole))
c
c     use diagonal preconditioner elements as first approximation
c
         polmin = 0.00000001d0
         do i = 1, npole
            poli = udiag * max(polmin,charge_trans(i))
            do j = 1, 3
               zrsd(j,i) = 0.0d0
               zrsdp(j,i) = 0.0d0
               zrsdt(j,i) = poli * rsd(j,i)
               zrsdtp(j,i) = poli * rsdp(j,i)
            end do
         end do
c
c     use the off-diagonal preconditioner elements in second phase
c
!$OMP PARALLEL default(private) shared(npole,mindex,minv,nulst,ulst,
!$OMP& rsd,rsdp,zrsd,zrsdp,zrsdt,zrsdtp)
!$OMP DO reduction(+:zrsdt,zrsdtp) schedule(guided)
         do i = 1, npole
            m = mindex(i)
            do kk = 1, nulst(i)
               k = ulst(kk,i)
               m1 = minv(m+1)
               m2 = minv(m+2)
               m3 = minv(m+3)
               m4 = minv(m+4)
               m5 = minv(m+5)
               m6 = minv(m+6)
               m = m + 6
               zrsdt(1,i) = zrsdt(1,i) + m1*rsd(1,k) + m2*rsd(2,k)
     &                        + m3*rsd(3,k)
               zrsdt(2,i) = zrsdt(2,i) + m2*rsd(1,k) + m4*rsd(2,k)
     &                        + m5*rsd(3,k)
               zrsdt(3,i) = zrsdt(3,i) + m3*rsd(1,k) + m5*rsd(2,k)
     &                        + m6*rsd(3,k)
               zrsdt(1,k) = zrsdt(1,k) + m1*rsd(1,i) + m2*rsd(2,i)
     &                        + m3*rsd(3,i)
               zrsdt(2,k) = zrsdt(2,k) + m2*rsd(1,i) + m4*rsd(2,i)
     &                        + m5*rsd(3,i)
               zrsdt(3,k) = zrsdt(3,k) + m3*rsd(1,i) + m5*rsd(2,i)
     &                        + m6*rsd(3,i)
               zrsdtp(1,i) = zrsdtp(1,i) + m1*rsdp(1,k) + m2*rsdp(2,k)
     &                         + m3*rsdp(3,k)
               zrsdtp(2,i) = zrsdtp(2,i) + m2*rsdp(1,k) + m4*rsdp(2,k)
     &                         + m5*rsdp(3,k)
               zrsdtp(3,i) = zrsdtp(3,i) + m3*rsdp(1,k) + m5*rsdp(2,k)
     &                         + m6*rsdp(3,k)
               zrsdtp(1,k) = zrsdtp(1,k) + m1*rsdp(1,i) + m2*rsdp(2,i)
     &                         + m3*rsdp(3,i)
               zrsdtp(2,k) = zrsdtp(2,k) + m2*rsdp(1,i) + m4*rsdp(2,i)
     &                         + m5*rsdp(3,i)
               zrsdtp(3,k) = zrsdtp(3,k) + m3*rsdp(1,i) + m5*rsdp(2,i)
     &                         + m6*rsdp(3,i)
            end do
         end do
!$OMP END DO
c
c     transfer the results from local to global arrays
c
!$OMP DO
         do i = 1, npole
            do j = 1, 3
               zrsd(j,i) = zrsdt(j,i) + zrsd(j,i)
               zrsdp(j,i) = zrsdtp(j,i) + zrsdp(j,i)
            end do
         end do
!$OMP END DO
!$OMP END PARALLEL
c
c     perform deallocation of some local arrays
c
         deallocate (zrsdt)
         deallocate (zrsdtp)
c
c     build the off-diagonal elements of preconditioning matrix
c
      else if (mode .eq. 'BUILD') then
         m = 0
         do i = 1, npole
            mindex(i) = m
            m = m + 6*nulst(i)
         end do
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
c     set OpenMP directives for the major loop structure
c
!$OMP PARALLEL default(private) shared(n,npole,ipole,x,y,z,pdamp,
!$OMP& thole,polarity,u1scale,u2scale,u3scale,u4scale,np11,ip11,
!$OMP& np12,ip12,np13,ip13,np14,ip14,nulst,ulst,mindex,minv)
!$OMP& firstprivate (dscale)
c
c     determine the off-diagonal elements of the preconditioner
c
!$OMP DO schedule(guided)
         do i = 1, npole
            ii = ipole(i)
            xi = x(ii)
            yi = y(ii)
            zi = z(ii)
            pdi = pdamp_ct(i)
            pti = thole_ct(i)
            poli = charge_trans(i)
            pti_beta_ct = beta_damp_ct(i) 
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
            m = mindex(i)
            do kkk = 1, nulst(i)
               k = ulst(kkk,i)
               kk = ipole(k)
               xr = x(kk) - xi
               yr = y(kk) - yi
               zr = z(kk) - zi
               call image (xr,yr,zr)
               r2 = xr*xr + yr* yr + zr*zr
               r = sqrt(r2)
               scale3 = dscale(kk)
               scale5 = dscale(kk)
               damp = pdi * pdamp_ct(k)
               damp_beta = pti_beta_ct * beta_damp_ct(k)
               if (damp .ne. 0.0d0) then
                  pgamma = min(pti,thole_ct(k))
                  damp = -pgamma * (r/damp)**3
                  if (damp .gt. -50.0d0) then
                     expdamp = exp(damp)
                     scale3 = scale3 * expdamp*damp_beta
                     scale5 = scale5 * expdamp*(1.0d0-damp)*damp_beta
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
!$OMP END DO
!$OMP END PARALLEL
c
c     perform deallocation of some local arrays
c
         deallocate (dscale)
      end if
      return
      end
