      subroutine induce0a_aniso_invert
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
c      integer maxiter
c      real*8 polmin
c      real*8 eps,epsold
c      real*8 epsd,epsp
c      real*8 udsum,upsum
c      real*8 a,ap,b,bp
c      real*8 sum,sump
c      real*8, allocatable :: anisopoli(:,:)
      real*8, allocatable :: ranisopol_mat(:,:)
      real*8, allocatable :: field(:,:)
      real*8, allocatable :: fieldp(:,:)
      real*8, allocatable :: T_mut(:,:)
      real*8, allocatable :: Cmat(:,:)
c      real*8, allocatable :: udir(:,:)
c      real*8, allocatable :: udirp(:,:)
c      real*8, allocatable :: rsd(:,:)
c      real*8, allocatable :: rsdp(:,:)
c      real*8, allocatable :: zrsd(:,:)
c      real*8, allocatable :: zrsdp(:,:)
c      real*8, allocatable :: conj(:,:)
c      real*8, allocatable :: conjp(:,:)
c      real*8, allocatable :: vec(:,:)
c      real*8, allocatable :: vecp(:,:)
c      logical done
c      character*6 mode
      integer i2,k2,i1
c
c
c     zero out the induced dipoles at each site
c
      print*,"In induce_aniso_invert"
      do i = 1, npole
         if(i.eq.3) then
           do j=1,9
             print*,"ranisopol j=",j,ranisopolarity(j,i)
             print*,"anisopol j=",j,anisopolarity(j,i)
           end do
         end if 
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

      allocate (ranisopol_mat(3*npole,3*npole))
      allocate (T_mut(3*npole,3*npole))
      allocate (Cmat(3*npole,3*npole))      

      do i=1,3*npole
         do j=1,3*npole
            T_mut(j,i) = 0.0d0
            ranisopol_mat(j,i) = 0.0d0
            Cmat(j,i) = 0.0d0
         end do
      end do


c     Calculate direct field and mutual interaction tensor,T_mut
c     Note that the minus sign has already been applied
c     to T_mut in the subroutine below (-T_mut)
 
      call dfield0a_mut_tensor (field,fieldp,T_mut)

c     Place anisotropic polarizabilities in matrix
c      do i=1,npole-1
c         i2=3*(i-1)
c         do k=i+1,npole
c            k2 = 3*(k-1)
c  
c            ranisopol_mat(i2+1,k2+1) = ranisopolarity(1,i) !alpha_xx
c            ranisopol_mat(i2+1,k2+2) = ranisopolarity(2,i) !alpha_xy
c            ranisopol_mat(i2+1,k2+3) = ranisopolarity(3,i) !alpha_xz
c            ranisopol_mat(i2+2,k2+1) = ranisopolarity(4,i) !alpha_yx
c            ranisopol_mat(i2+2,k2+2) = ranisopolarity(5,i) !alpha_yy
c            ranisopol_mat(i2+2,k2+3) = ranisopolarity(6,i) !alpha_yz
c            ranisopol_mat(i2+3,k2+1) = ranisopolarity(7,i) !alpha_zx
c            ranisopol_mat(i2+3,k2+2) = ranisopolarity(8,i) !alpha_zy
c            ranisopol_mat(i2+3,k2+3) = ranisopolarity(9,i) !alpha_zz
c
c            ranisopol_mat(k2+1,i2+1) = ranisopolarity(1,i) !alpha_xx
c            ranisopol_mat(k2+1,i2+2) = ranisopolarity(2,i) !alpha_xy
c            ranisopol_mat(k2+1,i2+3) = ranisopolarity(3,i) !alpha_xz
c            ranisopol_mat(k2+2,i2+1) = ranisopolarity(4,i) !alpha_yx
c            ranisopol_mat(k2+2,i2+2) = ranisopolarity(5,i) !alpha_yy
c            ranisopol_mat(k2+2,i2+3) = ranisopolarity(6,i) !alpha_yz
c            ranisopol_mat(k2+3,i2+1) = ranisopolarity(7,i) !alpha_zx
c            ranisopol_mat(k2+3,i2+2) = ranisopolarity(8,i) !alpha_zy
c            ranisopol_mat(k2+3,i2+3) = ranisopolarity(9,i) !alpha_zz
c
c         end do
c      end do

      do i=1,npole
         i2=3*(i-1)

            ranisopol_mat(i2+1,i2+1) = ranisopolarity(1,i) !alpha_xx
            ranisopol_mat(i2+1,i2+2) = ranisopolarity(2,i) !alpha_xy
            ranisopol_mat(i2+1,i2+3) = ranisopolarity(3,i) !alpha_xz
            ranisopol_mat(i2+2,i2+1) = ranisopolarity(4,i) !alpha_yx
            ranisopol_mat(i2+2,i2+2) = ranisopolarity(5,i) !alpha_yy
            ranisopol_mat(i2+2,i2+3) = ranisopolarity(6,i) !alpha_yz
            ranisopol_mat(i2+3,i2+1) = ranisopolarity(7,i) !alpha_zx
            ranisopol_mat(i2+3,i2+2) = ranisopolarity(8,i) !alpha_zy
            ranisopol_mat(i2+3,i2+3) = ranisopolarity(9,i) !alpha_zz

      end do

c
c     Invert matrix of anisotropic polarizabilities
c
      call invert(3*npole,ranisopol_mat)

c
c     Construct 'C matrix' from inverted polarizability matrix and mutual
c     interaction tensor.  
c     Nomenclature and equations per Pengyu Ren et al., JCTC, 2011,
c     7:  3143-3161

      do i=1,3*npole
         do j=1,3*npole
            Cmat(j,i)=ranisopol_mat(j,i)+T_mut(j,i)
         end do
      end do


c
c     Invert C matrix    
c
      call invert(3*npole,Cmat)

         do i = 1, npole
            do i1 = 1, 3
               i2 = 3*(i-1) + i1
               do k = 1, npole
                  k2 = 3*(k-1)
                  uind(i1,i)=uind(i1,i)+ Cmat(i2,k2+1)*field(1,k)+
     &                                   Cmat(i2,k2+2)*field(2,k)+
     &                                   Cmat(i2,k2+3)*field(3,k)
                  uinp(i1,i)=uinp(i1,i)+ Cmat(i2,k2+1)*fieldp(1,k)+
     &                                   Cmat(i2,k2+2)*fieldp(2,k)+
     &                                   Cmat(i2,k2+3)*fieldp(3,k)

               end do
            end do
         end do

      deallocate (field)
      deallocate (fieldp)

      deallocate (ranisopol_mat)
      deallocate (T_mut)
      deallocate (Cmat) 

      return
      end
