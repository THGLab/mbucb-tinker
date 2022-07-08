c
c
c     ###################################################
c     ##  COPYRIGHT (C)  2012  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #################################################################
c     ##                                                             ##
c     ##  program testpol  --  check convergence of induced dipoles  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "testpol" computes the induced dipole moments for direct
c     polarization and increasing numbers of SCF iterations in
c     order to monitor the degree and rate of convergence
c
c
      program testpol
      use sizes
      use atoms
      use bound
      use inform
      use iounit
      use limits
      use polar
      use polpot
      use potent
      use rigid
      use units
      use kpolr_aniso
      implicit none
      integer i,j,k
      integer maxiter
      real*8 sum
      real*8, allocatable :: rms(:)
      real*8, allocatable :: drms(:)
      real*8, allocatable :: uexact(:,:)
      real*8, allocatable :: ustore(:,:,:)
      real*8 :: ksum
c
c
c     get the coordinates and required force field parameters
c
      call initial
      call getxyz
      call mechanic
c
c     check to make sure mutual polarization is being used
c
      if (.not.use_polar .or. poltyp.eq.'DIRECT') then
         write (iout,10)
   10    format (' TESTPOL  --  Mutual Polarization Model not in Use')
         call fatal
      end if
c
c     maintain any periodic boundary conditions
c
      if (use_bounds .and. .not.use_rigid) then
            if(use_anisopolz) then
             write (iout,11)
   11      format (/,' Anisotropic polarizability is not yet available',
     &              ' for Particle Mesh Ewald with Neighbor Lists')
             call fatal
           else
          call bounds
         end if
       end if 
c
c     update the pairwise interaction neighbor lists
c
      if (use_list)  then
            if(use_anisopolz) then
             write (iout,11)
   21      format (/,' Anisotropic polarizability is not yet available',
     &              ' for Particle Mesh Ewald with Neighbor Lists')
             call fatal
           else 
             call nblist
       end if
      end if
c
c     perform dynamic allocation of some local arrays
c
      maxiter = 500
      allocate (rms(0:maxiter))
      allocate (drms(maxiter))
      allocate (uexact(3,n))
      allocate (ustore(3,n,0:maxiter))
c
c     set tolerance and rotate multipoles to global frame
c
      debug = .false.
      poleps = 0.0000000001d0
      call chkpole
      call rotpole
c
c     get induced dipoles for noniterative direct polarization
c
      poltyp = 'DIRECT'
      if(use_anisopolz) then
      ! call induce_aniso
       call induce_aniso_akd
c       call induce0a_aniso_invert
c       call induce0a_aniso_invert_test
      else
       call induce
      end if
      do i = 1, n
         do j = 1, 3
            ustore(j,i,0) = debye * uind(j,i)
c             write (iout,*)'direct akd',uind(j,i)
         end do
      end do
c
c     compute induced dipoles for increasing iteration counts
c
      poltyp = 'MUTUAL'
      do k = 1, maxiter
         politer = k
      if(use_anisopolz) then
      ! call induce_aniso
c       call induce0a_aniso_invert
       call induce_aniso_akd
c       call induce0a_aniso_invert_test
      else
         call induce
      end if
         do i = 1, n
            do j = 1, 3
               ustore(j,i,k) = debye * uind(j,i)
c             write (iout,*)'mutual akd',uind(j,i),'ustore',ustore(j,i,k)
            end do
         end do
c
c     find the RMS difference per atom between iterations
c
         sum = 0.0d0
         do i = 1, n
            do j = 1, 3
               sum = sum + (ustore(j,i,k)-ustore(j,i,k-1))**2
            end do
         end do
         drms(k) = sqrt(sum/dble(npolar))
         if (drms(k) .lt. 0.5d0*poleps)  goto 20
      end do
   20 continue
      maxiter = politer
c
c     find the RMS of each iteration from converged dipoles
c
      do i = 1, n
         do j = 1, 3
            uexact(j,i) = ustore(j,i,maxiter)
         end do
      end do
      do k = 0, maxiter
         sum = 0.0d0
         do i = 1, n
            do j = 1, 3
               sum = sum + (ustore(j,i,k)-uexact(j,i))**2
c            write(iout,*)'akd diff',k, ustore(1,i,k)-uexact(1,i),
c     &                                 ustore(2,i,k)-uexact(2,i),
c     &                                 ustore(3,i,k)-uexact(3,i)
c            write(iout,*)'uexact diff',k,
c     &                   uexact(1,i),uexact(2,i),uexact(3,i)     
            end do
         end do
         rms(k) = sqrt(sum/dble(npolar))
      end do
c
c     print the RMS between iterations and to converged dipoles
c
      write (iout,30)
   30 format (/,' Convergence of Induced Dipole Moments :',
     &        //,2x,'Iteration',10x,'RMS Change',11x,'RMS vs Final')
      write (iout,40)  0,rms(0)
   40 format (/,i8,16x,'----',6x,f20.10)
      do k = 1, maxiter
         write (iout,50)  k,drms(k),rms(k)
   50    format (i8,3x,f20.10,3x,f20.10)
         if (rms(k) .lt. 0.5d0*poleps)  goto 60
      end do
   60 continue
      write (iout,*)'Total Iteration',k+1
c
c     perform deallocation of some local arrays
c
      deallocate (rms)
      deallocate (drms)
      deallocate (uexact)
      deallocate (ustore)
      end
