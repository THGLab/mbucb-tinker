
      subroutine dfield0a_mut_tensor(field,fieldp,T_mut)
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
      use iounit
      implicit none
      real*8, allocatable :: dscale_dir(:)
      real*8, allocatable :: dscale_mut(:)
      real*8, allocatable :: pscale(:)
      real*8 xr,yr,zr
      real*8 fgrp,r,r2
      real*8 rr3_dir,rr5_dir,rr7_dir
      real*8 rr3,rr5
      real*8 ci,dix,diy,diz
      real*8 duix,duiy,duiz
      real*8 puix,puiy,puiz
      real*8 qixx,qixy,qixz
      real*8 qiyy,qiyz,qizz
      real*8 ck,dkx,dky,dkz
      real*8 dukx,duky,dukz
      real*8 pukx,puky,pukz
      real*8 qkxx,qkxy,qkxz
      real*8 qkyy,qkyz,qkzz
      real*8 dir,duir,puir
      real*8 dkr,dukr,pukr
      real*8 qix,qiy,qiz,qir
      real*8 qkx,qky,qkz,qkr
      real*8 damp,expdamp
      real*8 scale3_dir,scale5_dir
      real*8 scale7_dir
      real*8 scale3_mut,scale5_mut
      real*8 pdi,pti,pgamma
      real*8 field(3,npole)
      real*8 fieldp(3,npole)
      real*8 fid(3),fkd(3),T_mut(3*npole,3*npole)
      real*8 fip(3),fkp(3)
      logical proceed
      character*6 mode
      integer i,ii,j,l1,k,m,kk
      integer l2,l3,k1,k2,i1,i2
      real*8 Txx,Txy,Txz,Tyx
      real*8 Tyy,Tyz,Tzx,Tzy,Tzz
c
c
c     zero out the value of the field at each site
c
      write (iout,11)
   11 format (/,' Note that the direct field field calculation',
     &   ' assumes that there is no cutoff nor periodic',
     &   ' boundary conditions.  use_replica is should be false,',
     &   ' and the portion of the code allowing for calculation',
     &   'between the central cell and image cells is omitted.', 
     &   'Refer to dfield0a in the original Tinker 7 distribution.')

      if(use_replica) then
             write (iout,21)
   21 format (/,' This code is for real-space only without',
     &          ' periodic boundary conditions.')
             call fatal
      end if

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

      allocate (dscale_dir(n))
      allocate (dscale_mut(n))
      allocate (pscale(n))
c
c     find the electrostatic field due to permanent multipoles,
c     as well the elements of the mutual interaction tensor

      do i = 1, npole-1
         i2 = 3*(i-1)
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
         do j = i+1, npole
            dscale_dir(ipole(j)) = 1.0d0
            pscale(ipole(j)) = 1.0d0
            dscale_mut(ipole(j)) = 1.0d0
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
            dscale_dir(ip11(j,ii)) = d1scale
            dscale_mut(ip11(j,ii)) = u1scale
         end do
         do j = 1, np12(ii)
            dscale_dir(ip12(j,ii)) = d2scale
            dscale_mut(ip12(j,ii)) = u2scale
         end do
         do j = 1, np13(ii)
            dscale_dir(ip13(j,ii)) = d3scale
            dscale_mut(ip13(j,ii)) = u3scale
         end do
         do j = 1, np14(ii)
            dscale_dir(ip14(j,ii)) = d4scale
            dscale_mut(ip14(j,ii)) = u4scale
         end do
         do k = i+1, npole
            k2=3*(k-1)
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
                  scale3_dir = 1.0d0
                  scale5_dir = 1.0d0
                  scale7_dir = 1.0d0
                  scale3_mut = dscale_mut(kk)
                  scale5_mut = dscale_mut(kk)

                  damp = pdi * pdamp(k)
                  if (damp .ne. 0.0d0) then
                     pgamma = min(pti,thole(k))
                     damp = -pgamma * (r/damp)**3
                     if (damp .gt. -50.0d0) then
                        expdamp = exp(damp)
                        scale3_dir = 1.0d0 - expdamp
                        scale5_dir = 1.0d0 - expdamp*(1.0d0-damp)
                        scale7_dir = 1.0d0 - expdamp
     &                              *(1.0d0-damp+0.6d0*damp**2)
                        scale3_mut = scale3_mut * (1.0d0-expdamp)
                        scale5_mut = scale5_mut * (1.0d0-expdamp
     &                                        *(1.0d0-damp))

                     end if
                  end if
                  rr3_dir = scale3_dir / (r*r2)
                  rr5_dir = 3.0d0 * scale5_dir / (r*r2*r2)
                  rr7_dir = 15.0d0 * scale7_dir / (r*r2*r2*r2)

                  rr3 = scale3_mut / (r*r2)
                  rr5 = 3.0d0 * scale5_mut / (r*r2*r2)

                  Txx = -(-rr3 + xr*xr*rr5)
                  Txy = -(xr*yr*rr5)
                  Txz = -(xr*zr*rr5)
                  Tyx = Txy
                  Tyy = -(-rr3 + yr*yr*rr5)
                  Tyz = -(yr*zr*rr5)
                  Tzx = Txz
                  Tzy = Tyz
                  Tzz = -(-rr3 + zr*zr*rr5)

                  T_mut(i2+1,k2+1) = Txx
                  T_mut(i2+1,k2+2) = Txy
                  T_mut(i2+1,k2+3) = Txz
                  T_mut(i2+2,k2+1) = Tyx
                  T_mut(i2+2,k2+2) = Tyy
                  T_mut(i2+2,k2+3) = Tyz
                  T_mut(i2+3,k2+1) = Tzx
                  T_mut(i2+3,k2+2) = Tzy
                  T_mut(i2+3,k2+3) = Tzz

                  T_mut(k2+1,i2+1) = Txx
                  T_mut(k2+1,i2+2) = Txy
                  T_mut(k2+1,i2+3) = Txz
                  T_mut(k2+2,i2+1) = Tyx
                  T_mut(k2+2,i2+2) = Tyy
                  T_mut(k2+2,i2+3) = Tyz
                  T_mut(k2+3,i2+1) = Tzx
                  T_mut(k2+3,i2+2) = Tzy
                  T_mut(k2+3,i2+3) = Tzz

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
                  fid(1) = -xr*(rr3_dir*ck-rr5_dir*dkr+rr7_dir*qkr)
     &                        - rr3_dir*dkx + 2.0d0*rr5_dir*qkx
                  fid(2) = -yr*(rr3_dir*ck-rr5_dir*dkr+rr7_dir*qkr)
     &                        - rr3_dir*dky + 2.0d0*rr5_dir*qky
                  fid(3) = -zr*(rr3_dir*ck-rr5_dir*dkr+rr7_dir*qkr)
     &                        - rr3_dir*dkz + 2.0d0*rr5_dir*qkz
                  fkd(1) = xr*(rr3_dir*ci+rr5_dir*dir+rr7_dir*qir)
     &                        - rr3_dir*dix - 2.0d0*rr5_dir*qix
                  fkd(2) = yr*(rr3_dir*ci+rr5_dir*dir+rr7_dir*qir)
     &                        - rr3_dir*diy - 2.0d0*rr5_dir*qiy
                  fkd(3) = zr*(rr3_dir*ci+rr5_dir*dir+rr7_dir*qir)
     &                        - rr3_dir*diz - 2.0d0*rr5_dir*qiz
                  do j = 1, 3
                     field(j,i) = field(j,i) + fid(j)*dscale_dir(kk)
                     field(j,k) = field(j,k) + fkd(j)*dscale_dir(kk)
                     fieldp(j,i) = fieldp(j,i) + fid(j)*pscale(kk)
                     fieldp(j,k) = fieldp(j,k) + fkd(j)*pscale(kk)
                  end do
               end if
            end if
         end do
      end do

      deallocate (dscale_dir)
      deallocate (dscale_mut)
      deallocate (pscale)
      return
      end

