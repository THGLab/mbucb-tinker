c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1995  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine kpolar  --  assign polarizability parameters  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "kpolar" assigns atomic dipole polarizabilities to the atoms
c     within the structure and processes any new or changed values
c
c
      subroutine kpolar
      use sizes
      use atoms
      use inform
      use iounit
      use keys
      use kpolr
      use mpole
      use neigh
      use polar
      use polpot
      use potent
      use usolve
      use kpolr_aniso
      use polar_aniso
      use penetrate
      use chargetransfer
      use chgpen, only : ecp
      implicit none
      integer i,j,k
      integer npg,next
      integer pg(maxval)
      real*8 pol,thl
      real*8 sixth
      logical header
      character*20 keyword
      character*120 record
      character*120 string
      real*8 anisopol(9)
      real*8 pn
      real*8 cq_trans,thl_ct,beta_ct
      integer jj
! qtw >> in
      integer :: n_valence
      real*8  :: a_jp
      real*8  :: b_jp
      real*8  :: s_lo
      real*8  :: s_hi
      logical :: chgpen_header
! qtw << out

c
c
c     process keywords containing polarizability parameters
c
      header = .true.
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         if ((keyword(1:9) .eq. 'POLARIZE ').and.
     &       (.not.use_anisopolz)) then
            k = 0
            pol = 0.0d0
            thl = -1.0d0
            do j = 1, maxval
               pg(j) = 0
            end do
            call getnumb (record,k,next)
            string = record(next:120)
            read (string,*,err=10,end=10)  pol,thl,(pg(j),j=1,maxval)
   10       continue
            if (k .gt. 0) then
               if (header .and. .not.silent) then
                  header = .false.
                  write (iout,20)
   20             format (/,' Additional Atomic Dipole',
     &                       ' Polarizability Parameters :',
     &                    //,5x,'Atom Type',11x,'Alpha',8x,
     &                       'Damp',5x,'Group Atom Types'/)
               end if
               if (k .le. maxtyp) then
                  polr(k) = pol
                  athl(k) = thl
                  do j = 1, maxval
                     pgrp(j,k) = pg(j)
                     if (pg(j) .eq. 0) then
                        npg = j - 1
                        goto 30
                     end if
                  end do
   30             continue
                  if (.not. silent) then
                     write (iout,40)  k,pol,thl,(pg(j),j=1,npg)
   40                format (4x,i6,10x,f10.3,2x,f10.3,7x,20i5)
                  end if
               else
                  write (iout,50)
   50             format (/,' KPOLAR  --  Too many Dipole',
     &                       ' Polarizability Parameters')
                  abort = .true.
               end if
            end if
         else if( (keyword(1:14) .eq. 'ANISOPOLARIZE ').and.
     &              (use_anisopolz ) ) then
            k = 0
            do jj=1,9
             anisopol(jj) = 0.0d0
            end do
            thl = -1.0d0
            do j = 1, maxval
               pg(j) = 0
            end do
            call getnumb (record,k,next)
            string = record(next:120)
            read (string,*,err=11,end=11) anisopol(1),anisopol(2),
     &      anisopol(3),anisopol(4),anisopol(5),anisopol(6),
     &      anisopol(7),anisopol(8),anisopol(9),thl,
     &                                       (pg(j),j=1,maxval) 
   11       continue
            if (k .gt. 0) then
               if (header .and. .not.silent) then
                  header = .false.
                  write (iout,21)
   21             format (/,' Additional Atomic Dipole',
     &                       ' Polarizability Parameters :',
     &                    //,5x,'Atom Type',11x,'Alpha',8x,
     &                       'Damp',5x,'Group Atom Types'/)
               end if
               if (k .le. maxtyp) then
                  do jj=1,9
                   polr_aniso(jj,k) = anisopol(jj)
                  end do
                  athl(k) = thl
                  do j = 1, maxval
                     pgrp(j,k) = pg(j)
                     if (pg(j) .eq. 0) then
                        npg = j - 1
                        goto 31
                     end if
                  end do
   31             continue
                  if (.not. silent) then
                     write (iout,41) k,anisopol(1),anisopol(2),
     &                  anisopol(3),anisopol(4),anisopol(5),anisopol(6),
     &                  anisopol(7),anisopol(8),anisopol(9),
     &                  thl,(pg(j),j=1,npg)
   41                format (4x,i6,10x,f10.3,2x,f10.3,2x,f10.3,2x,
     &                       f10.3,2x,f10.3,2x,f10.3,2x,f10.3,2x,
     &                       f10.3,2x,f10.3,2x,f10.3,7x,20i5)
                  end if
               else
                  write (iout,51)
   51             format (/,' KPOLAR  --  Too many Dipole',
     &                       ' AnisoPolarizability Parameters')
                  abort = .true.
               end if
            end if

         else if( (keyword(1:9) .eq. 'POLARIZE ').and.
     &              (use_anisopolz ) ) then
            k = 0
            pol = 0.0d0
            do jj=1,9
             anisopol(jj) = 0.0d0
            end do
            thl = -1.0d0
            do j = 1, maxval
               pg(j) = 0
            end do
            call getnumb (record,k,next)
            string = record(next:120)
            read (string,*,err=12,end=12)  pol,thl,(pg(j),j=1,maxval)
   12       continue
            if (k .gt. 0) then
               if (header .and. .not.silent) then
                  header = .false.
                  write (iout,22)
   22             format (/,' Additional Atomic Dipole',
     &                       ' Polarizability Parameters :',
     &                    //,5x,'Atom Type',11x,'Alpha',8x,
     &                       'Damp',5x,'Group Atom Types'/)
               end if
               if (k .le. maxtyp) then
                  do jj=1,9
                   if((jj.eq.1).or.(jj.eq.5).or.(jj.eq.9)) then
                     polr_aniso(jj,k) = pol
                   else
                     polr_aniso(jj,k) = 0.0d0
                   end if
                  end do
                  athl(k) = thl
                  do j = 1, maxval
                     pgrp(j,k) = pg(j)
                     if (pg(j) .eq. 0) then
                        npg = j - 1
                        goto 32
                     end if
                  end do
   32             continue
                  if (.not. silent) then
                     write (iout,42)  k,pol,thl,(pg(j),j=1,npg)
   42                format (4x,i6,10x,f10.3,2x,f10.3,7x,20i5)
                  end if
               else
                  write (iout,52)
   52             format (/,' KPOLAR  --  Too many Dipole',
     &                       ' Polarizability Parameters')
                  abort = .true.
               end if
            end if

! qtw >> in
c
c     charge penetration parameters
c
         else if (keyword(1:10) .eq. 'CHGPENPRM ') then
            k = 0
            n_valence = 0
            a_jp = 0.0D0
            b_jp = 0.0D0

            call getnumb (record,k,next)
            string = record(next:120)
            ! read type N_val_JP alpha_JP beta_JP
            read (string,*,err=200,end=200) n_valence,a_jp, b_jp
  200       continue
            if (k .gt. 0) then
               if (chgpen_header) then
                  chgpen_header = .false.
                  write (iout,210)
  210             format (/,' Additional Charge',
     &                        ' Penetration Parameters:',
     &                    /,6x,'Atom',15x,'N_val',7x,'Alpha',7x,'Beta'/)
               end if

               if (k .le. maxtyp) then
                  ecp%n_val(k)  = n_valence
                  ecp%alp_jp(k) = a_jp
                  ecp%bet_jp(k) = b_jp

                  write (iout,220) k,n_valence,a_jp,b_jp
  220             format (4x,i6,10x,i10,2x,f10.3,2x,f10.3)
               else
                  write (iout,230)
  230             format (/,' KPOLAR  --  Too many',
     &                       ' Charge Penetration Parameters')
                  abort = .true.
               end if

            end if

         else if (keyword(1:14) .eq. 'CHGPEN-SWITCH ') then
            s_lo = 0.0D0
            s_hi = 0.0D0

            string = record(next:120)
            ! read switch_lo and switch_hi
            read (string,*,err=446,end=446) s_lo, s_hi
  446       continue
            ecp%switch_lo = s_lo
            ecp%switch_hi = s_hi

! qtw >> out

         else if(keyword(1:12) .eq. 'PENETRATION ') then
            k = 0
            pn = 0.0d0
            call getnumb (record,k,next)
            string = record(next:120)
            read (string,*,err=91,end=91)  pn
   91       continue
            if (k .gt. 0) then
               if (header .and. .not.silent) then
                  header = .false.
                  write (iout,92)
   92             format (/,' Additional',
     &               ' Penetration Parameters :',
     &               //,5x,'Atom Type',11x,'Pene'/)
               end if
               if (k .le. maxtyp) then
                 penr(k) = pn
                  if (.not. silent) then
                     write (iout,93)  k,pn
   93                format (4x,i6,10x,f10.3)
                  end if
               else
                  write (iout,94)
   94             format (/,' KPOLAR  --  Too many Penetration',
     &                       ' Parameters')
                  abort = .true.
               end if
            end if

         else if(keyword(1:15) .eq. 'CHARGETRANSFER ') then
            k = 0
            cq_trans = 0.0d0
            thl_ct = -1.0D0
            beta_ct = 1.0D0
            do j = 1, maxval
               pg(j) = 0
            end do
            call getnumb (record,k,next)
            string = record(next:120)
            read (string,*,err=191,end=191)cq_trans,thl_ct,beta_ct,
     &                    (pg(j),j=1,maxval)
  191       continue
            if (k .gt. 0) then
               if (header .and. .not.silent) then
                  header = .false.
                  write (iout,192)
  192             format (/,' Additional',
     &               ' Charge transfer Parameters :',
     &               //,5x,'Atom Type',11x,'charge_trans'/)
               end if
               if (k .le. maxtyp) then
                 cqr_trans(k) = cq_trans
                 rthole_ct(k) = thl_ct 
                 rbeta_charge(k) = beta_ct
                  do j = 1, maxval
                     pgrp(j,k) = pg(j)
                     if (pg(j) .eq. 0) then
                        npg = j - 1
                        goto 37
                     end if
                  end do
   37             continue
                  if (.not. silent) then
                     write (iout,193) k,cq_trans,thl_ct,beta_ct,
     &                             (pg(j),j=1,npg)
  193                format (4x,i6,10x,3f10.3,20i)
                  end if
               else
                  write (iout,194)
  194             format (/,' KPOLAR  --  Too many charge_trans',
     &                       ' Parameters')
                  abort = .true.
               end if
            end if
         end if
      end do
c
c     perform dynamic allocation of some global arrays
c
c      print*,"Before kpolar allocs"
      if (allocated(pen)) deallocate (pen)
      if (allocated(charge_trans)) deallocate (charge_trans)
      if (allocated(thole_ct))  deallocate (thole_ct)
      if (allocated(beta_charge))  deallocate (beta_charge)
      if (allocated(pdamp_ct))  deallocate (pdamp_ct)
      if (allocated(beta_damp_ct))  deallocate (beta_damp_ct)
      if (allocated(polarity))  deallocate (polarity)
      if (allocated(thole))  deallocate (thole)
      if (allocated(pdamp))  deallocate (pdamp)
      if (allocated(uind))  deallocate (uind)
      if (allocated(uinp))  deallocate (uinp)
      if (allocated(uinds))  deallocate (uinds)
      if (allocated(uinps))  deallocate (uinps)
      if (allocated(uind_ct))  deallocate (uind_ct)
      if (allocated(uinp_ct))  deallocate (uinp_ct)
      if (allocated(uinds_ct))  deallocate (uinds_ct)
      if (allocated(uinps_ct))  deallocate (uinps_ct)
      if (allocated(anisopolarity)) deallocate (anisopolarity)
      if (allocated(ranisopolarity)) deallocate (ranisopolarity)
      if(use_anisopolz) then
       allocate (anisopolarity(9,n))
       allocate (ranisopolarity(9,n))
      else
       allocate (polarity(n))
      end if
      allocate (thole(n))
      allocate (pdamp(n))
      allocate (uind(3,n))
      allocate (uinp(3,n))
      allocate (uinds(3,n))
      allocate (uinps(3,n))
      allocate (pen(n))
      allocate (charge_trans(n))
      allocate (thole_ct(n))
      allocate (beta_charge(n))
      allocate (pdamp_ct(n))
      allocate (beta_damp_ct(n))
      allocate (uind_ct(3,n))
      allocate (uinp_ct(3,n))
      allocate (uinds_ct(3,n))
      allocate (uinps_ct(3,n))
c      print*,"After kpolar allocs"
c
c     find and store the atomic dipole polarizability parameters
c
      do i = 1, n
        if(use_anisopolz) then
          do jj=1,9
           anisopolarity(jj,i)=polr_aniso(jj,type(i)) 
c       print*,"In kpolar i=",i,"anispolarityjj=",jj,anisopolarity(jj,i)
          end do
        else
           polarity(i) = polr(type(i))
        end if 
         thole(i) = athl(type(i))
         pen(i) = penr(type(i))
         charge_trans(i)= cqr_trans(type(i))
         thole_ct(i) = rthole_ct(type(i))
         beta_charge(i) = rbeta_charge(type(i)) 
      end do
c      print*,"After loading in anisopolarity array"
c
c     process keywords containing atom specific polarizabilities
c
      header = .true.
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         if ((keyword(1:9) .eq. 'POLARIZE ').and.
     &       (.not.use_anisopolz)) then
            k = 0
            pol = 0.0d0
            thl = 0.0d0
            call getnumb (record,k,next)
            if (k.lt.0 .and. k.ge.-n) then
               k = -k
               string = record(next:120)
               read (string,*,err=60,end=60)  pol,thl
   60          continue
               if (header) then
                  header = .false.
                  write (iout,70)
   70             format (/,' Additional Dipole Polarizabilities',
     &                       ' for Specific Atoms :',
     &                    //,6x,'Atom',15x,'Alpha',8x,'Damp',/)
               end if
               if (.not. silent) then
                  write (iout,80)  k,pol,thl
   80             format (4x,i6,10x,f10.3,2x,f10.3)
               end if
               polarity(k) = pol
               thole(k) = thl
            end if
         else if ( (keyword(1:14) .eq. 'ANISOPOLARIZE ') .and.
     &              (use_anisopolz ) ) then
            k = 0
            do jj=1,9
               anisopol(jj) = 0.0d0
            end do            
            thl = 0.0d0
            call getnumb (record,k,next)
            if (k.lt.0 .and. k.ge.-n) then
               k = -k
               string = record(next:120)
               read (string,*,err=61,end=61) anisopol(1),anisopol(2),
     &      anisopol(3),anisopol(4),anisopol(5),anisopol(6),
     &      anisopol(7),anisopol(8),anisopol(9),thl
   61          continue
               if (header) then
                  header = .false.
                  write (iout,71)
   71             format (/,' Additional Dipole Polarizabilities',
     &                       ' for Specific Atoms :',
     &                    //,6x,'Atom',15x,'Alpha',8x,'Damp',/)
               end if
               if (.not. silent) then
                  write (iout,81)  k,anisopol(1),anisopol(2),
     &      anisopol(3),anisopol(4),anisopol(5),anisopol(6),
     &      anisopol(7),anisopol(8),anisopol(9),thl
   81             format (4x,i6,10x,f10.3,2x,f10.3,2x,
     &                   f10.3,2x,f10.3,2x,f10.3,2x,f10.3,2x,
     &                   f10.3,2x,f10.3,2x,f10.3,2x,f10.3)
               end if
               do jj=1,9
                anisopolarity(jj,k) = anisopol(jj)
               end do
               thole(k) = thl
            end if
         else if ( (keyword(1:9) .eq. 'POLARIZE ') .and.
     &              (use_anisopolz ) ) then
            k = 0
            pol = 0.0d0
            thl = 0.0d0
            call getnumb (record,k,next)
            if (k.lt.0 .and. k.ge.-n) then
               k = -k
               string = record(next:120)
               read (string,*,err=62,end=62)  pol,thl
   62          continue
               if (header) then
                  header = .false.
                  write (iout,72)
   72             format (/,' Additional Dipole Polarizabilities',
     &                       ' for Specific Atoms :',
     &                    //,6x,'Atom',15x,'Alpha',8x,'Damp',/)
               end if
               if (.not. silent) then
                  write (iout,82)  k,pol,thl
   82             format (4x,i6,10x,f10.3,2x,f10.3)
               end if
               do jj=1,9
                if((jj.eq.1).or.(jj.eq.5).or.(jj.eq.9)) then
                 anisopolarity(jj,k) = pol
                else
                 anisopolarity(jj,k) = 0.0d0
                end if
               end do
               thole(k) = thl
            end if
! qtw >> in
c
c     charge penetration parameters
c
         else if (keyword(1:10) .eq. 'CHGPENPRM ') then
            k = 0
            n_valence = 0
            a_jp = 0.0D0
            b_jp = 0.0D0

            call getnumb (record,k,next)
            if (k.lt.0 .and. k.ge.-n) then
               k = -k
               string = record(next:120)
            ! read type N_val_JP alpha_JP beta_JP
               read (string,*,err=240,end=240) n_valence,a_jp,b_jp
  240          continue
               if (chgpen_header) then
                  chgpen_header = .false.
                  write (iout,250)
  250             format (/,' Additional Charge Penetration Parameters',
     &                       ' for Specific Atoms :',
     &                    /,6x,'Atom',15x,'N_val',7x,'Alpha',7x,'Beta'/)
               end if
               write (iout,260) k,n_valence,a_jp,b_jp
  260          format (4x,i6,10x,i10,2x,f10.3,2x,f10.3)

               ecp%n_val(k)  = n_valence
               ecp%alp_jp(k) = a_jp
               ecp%bet_jp(k) = b_jp
            end if
         else if (keyword(1:14) .eq. 'CHGPEN-SWITCH ') then
            s_lo = 0.0D0
            s_hi = 0.0D0

            string = record(next:120)
            ! read switch_lo and switch_hi
            read (string,*,err=447,end=447) s_lo, s_hi
  447       continue
            ecp%switch_lo = s_lo
            ecp%switch_hi = s_hi

! qtw >> out

         else if(keyword(1:12) .eq. 'PENETRATION ') then
            k = 0
            pn = 0.0d0
            call getnumb (record,k,next)
            if (k.lt.0 .and. k.ge.-n) then
               k = -k
               string = record(next:120)
               read (string,*,err=97,end=97)  pn
   97          continue
               if (header) then
                  header = .false.
                  write (iout,98)
   98             format (/,' Additional',
     &               ' Penetration Parameters :',
     &               //,5x,'Atom Type',11x,'Pene'/)
               end if
               if (.not. silent) then
                  write (iout,99)  k,pn
   99             format (4x,i6,10x,f10.3)
               end if
               pen(k) = pn
            end if

         else if(keyword(1:15) .eq. 'CHARGETRANSFER ') then
            k = 0
            cq_trans = 0.0d0
            thl_ct = 0.0D0
            beta_ct = 1.0D0 
            call getnumb (record,k,next)
            if (k.lt.0 .and. k.ge.-n) then
               k = -k
               string = record(next:120)
               read (string,*,err=197,end=197)  cq_trans,thl_ct,beta_ct
  197          continue
               if (header) then
                  header = .false.
                  write (iout,198)
  198             format (/,' Additional',
     &               ' Charge transfer Parameters :',
     &               //,5x,'Atom Type',11x,'Charge transfer'/)
               end if
               if (.not. silent) then
                  write (iout,199)  k,cq_trans,thl_ct,beta_ct
  199             format (4x,i6,10x,f10.3)
               end if
               charge_trans(k) = cq_trans
               thole_ct(k) = thl_ct
               beta_charge(k)=beta_ct
            end if
         end if
      end do
c
c     remove zero and undefined polarizable sites from the list
c
      npole = 0
      npolar = 0
      nctrans = 0
c      print*,"Before removal of zero and undefined polz"

      if (use_anisopolz) then
      do i = 1, n
         if (polsiz(i).ne.0 .or. anisopolarity(1,i).ne.0.0d0) then
            npole = npole + 1
            ipole(npole) = i
            pollist(i) = npole
            zaxis(npole) = zaxis(i)
            xaxis(npole) = xaxis(i)
            yaxis(npole) = yaxis(i)
            polaxe(npole) = polaxe(i)
            do k = 1, maxpole
               pole(k,npole) = pole(k,i)
            end do
            if (anisopolarity(1,i) .ne. 0.0d0)  npolar = npolar + 1
            do jj=1,9
             anisopolarity(jj,npole) = anisopolarity(jj,i)
            end do
            thole(npole) = thole(i)
         end if
      end do         
      else
      do i = 1, n
         if (polsiz(i).ne.0 .or. polarity(i).ne.0.0d0) then
            npole = npole + 1
            ipole(npole) = i
            pollist(i) = npole
            zaxis(npole) = zaxis(i)
            xaxis(npole) = xaxis(i)
            yaxis(npole) = yaxis(i)
            polaxe(npole) = polaxe(i)
            do k = 1, maxpole
               pole(k,npole) = pole(k,i)
            end do
            if (polarity(i) .ne. 0.0d0)  npolar = npolar + 1
            polarity(npole) = polarity(i)
            thole(npole) = thole(i)
         end if
      end do
      end if
cc akd
c      do i = 1, n
c         if (polsiz(i).ne.0 .or. charge_trans(i).ne.0.0d0) then
c            npole = npole + 1
c            ipole(npole) = i
c            pollist(i) = npole
c            zaxis(npole) = zaxis(i)
c            xaxis(npole) = xaxis(i)
c            yaxis(npole) = yaxis(i)
c            polaxe(npole) = polaxe(i)
c            do k = 1, maxpole
c               pole(k,npole) = pole(k,i)
c            end do
c            if (charge_trans(i) .ne. 0.0d0)  nctrans = nctrans + 1
c            charge_trans(npole) = charge_trans(i)
c            thole_ct(npole) = thole_ct(i)
c            beta_charge(npole) = beta_charge(i)
c            write(*,*)'dbug akd'
c         end if
c      end do
c akd
c      print*,"After removal of zero and undefined polz"
c
c     set the values used in the scaling of the polarizability
c
      sixth = 1.0d0 / 6.0d0
c      print*,"Before pdamp assignment"
      if(use_anisopolz) then
        do i=1, npole
           if (thole(i) .eq. 0.0d0) then
              pdamp(i) = 0.0d0
           else
              pdamp(i) = ( (anisopolarity(1,i)+ 
     &          anisopolarity(5,i)+ anisopolarity(9,i))/3.0d0)**sixth
           end if
        end do     
      else
        do i = 1, npole
           if (thole(i) .eq. 0.0d0) then
              pdamp(i) = 0.0d0
           else
              pdamp(i) = polarity(i)**sixth
           end if
        end do
      end if

      if(use_chargetranfer) then
        do i = 1, npole
           if (thole_ct(i) .eq. 0.0d0) then
              pdamp_ct(i) = 0.0d0
           else
              pdamp_ct(i) = charge_trans(i)**sixth
           end if
           if (beta_charge(i) .eq. 0.0d0) then
              beta_charge(i) = 0.0d0
           else
              beta_damp_ct(i) = beta_charge(i)
           end if
        end do
      end if
      
c      print*,"After pdamp assignment"
c
c     assign polarization group connectivity of each atom
c
c      print*,"Before call polargrp"

      call polargrp
c      print*,"After call polargrp"

c
c     test multipoles at chiral sites and invert if necessary
c
c      print*,"Before call chkpole"
      call chkpole
c      print*,"After call chkpole"
c
c     turn off polarizable multipole potential if it is not used
c
      if (npole .eq. 0)  use_mpole = .false.
      if (npolar .eq. 0)  use_polar = .false.
c
c     perform dynamic allocation of some global arrays
c
      if (use_polar) then
         if (allocated(mindex))  deallocate (mindex)
         if (allocated(minv))  deallocate (minv)
         allocate (mindex(npole))
         allocate (minv(3*maxulst*npole))
      end if
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine polargrp  --  polarization group connectivity  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "polargrp" generates members of the polarization group of
c     each atom and separate lists of the 1-2, 1-3 and 1-4 group
c     connectivities
c
c
      subroutine polargrp
      use sizes
      use atoms
      use couple
      use inform
      use iounit
      use kpolr
      use mpole
      use polgrp
      implicit none
      integer maxlist,maxkeep
      parameter (maxkeep=100)
      parameter (maxlist=1000)
      integer i,j,k
      integer it,jt
      integer jj,kk
      integer start,stop
      integer nlist,nkeep
      integer maxp11,maxp12
      integer maxp13,maxp14
      integer keep(maxkeep)
      integer list(maxlist)
      integer, allocatable :: mask(:)
      logical done
c
c
c     perform dynamic allocation of some global arrays
c
      maxp11 = 150
      maxp12 = 50
      maxp13 = 50
      maxp14 = 50
      if (allocated(np11))  deallocate (np11)
      if (allocated(np12))  deallocate (np12)
      if (allocated(np13))  deallocate (np13)
      if (allocated(np14))  deallocate (np14)
      if (allocated(ip11))  deallocate (ip11)
      if (allocated(ip12))  deallocate (ip12)
      if (allocated(ip13))  deallocate (ip13)
      if (allocated(ip14))  deallocate (ip14)
      allocate (np11(n))
      allocate (np12(n))
      allocate (np13(n))
      allocate (np14(n))
      allocate (ip11(maxp11,n))
      allocate (ip12(maxp12,n))
      allocate (ip13(maxp13,n))
      allocate (ip14(maxp14,n))
c
c     find the directly connected group members for each atom
c
      do i = 1, n
         np11(i) = 1
         ip11(1,i) = i
         it = type(i)
         do j = 1, n12(i)
            jj = i12(j,i)
            jt = type(jj)
            do k = 1, maxval
               kk = pgrp(k,it)
               if (kk .eq. 0)  goto 20
               if (pgrp(k,it) .eq. jt) then
                  np11(i) = np11(i) + 1
                  if (np11(i) .le. maxp11) then
                     ip11(np11(i),i) = jj
                  else
                     write (iout,10)
   10                format (/,' POLARGRP  --  Too many Atoms',
     &                          ' in Polarization Group')
                     abort = .true.
                     goto 30
                  end if
               end if
            end do
   20       continue
         end do
      end do
   30 continue
c
c     perform dynamic allocation of some local arrays
c
      allocate (mask(n))
c
c     find any other group members for each atom in turn
c
      do i = 1, n
         mask(i) = 0
      end do
      do i = 1, n
         done = .false.
         start = 1
         stop = np11(i)
         do j = start, stop
            jj = ip11(j,i)
            if (jj .lt. i) then
               done = .true.
               np11(i) = np11(jj)
               do k = 1, np11(i)
                  ip11(k,i) = ip11(k,jj)
               end do
            else
               mask(jj) = i
            end if
         end do
         do while (.not. done)
            done = .true.
            do j = start, stop
               jj = ip11(j,i)
               do k = 1, np11(jj)
                  kk = ip11(k,jj)
                  if (mask(kk) .ne. i) then
                     np11(i) = np11(i) + 1
                     if (np11(i) .le. maxp11) then
                        ip11(np11(i),i) = kk
                     else
                        write (iout,40)
   40                   format (/,' POLARGRP  --  Too many Atoms',
     &                             ' in Polarization Group')
                        abort = .true.
                        goto 50
                     end if
                     mask(kk) = i
                  end if
               end do
            end do
            if (np11(i) .ne. stop) then
               done = .false.
               start = stop + 1
               stop = np11(i)
            end if
         end do
         call sort (np11(i),ip11(1,i))
      end do
   50 continue
c
c     loop over atoms finding all the 1-2 group relationships
c
      do i = 1, n
         mask(i) = 0
      end do
      do i = 1, n
         do j = 1, np11(i)
            jj = ip11(j,i)
            mask(jj) = i
         end do
         nkeep = 0
         do j = 1, np11(i)
            jj = ip11(j,i)
            do k = 1, n12(jj)
               kk = i12(k,jj)
               if (mask(kk) .ne. i) then
                  nkeep = nkeep + 1
                  keep(nkeep) = kk
               end if
            end do
         end do
         nlist = 0
         do j = 1, nkeep
            jj = keep(j)
            do k = 1, np11(jj)
               kk = ip11(k,jj)
               nlist = nlist + 1
               list(nlist) = kk
            end do
         end do
         call sort8 (nlist,list)
         if (nlist .le. maxp12) then
            np12(i) = nlist
            do j = 1, nlist
               ip12(j,i) = list(j)
            end do
         else
            write (iout,60)
   60       format (/,' POLARGRP  --  Too many Atoms',
     &                 ' in 1-2 Polarization Group')
            abort = .true.
            goto 70
         end if
      end do
   70 continue
c
c     loop over atoms finding all the 1-3 group relationships
c
      do i = 1, n
         mask(i) = 0
      end do
      do i = 1, n
         do j = 1, np11(i)
            jj = ip11(j,i)
            mask(jj) = i
         end do
         do j = 1, np12(i)
            jj = ip12(j,i)
            mask(jj) = i
         end do
         nlist = 0
         do j = 1, np12(i)
            jj = ip12(j,i)
            do k = 1, np12(jj)
               kk = ip12(k,jj)
               if (mask(kk) .ne. i) then
                  nlist = nlist + 1
                  list(nlist) = kk
               end if
            end do
         end do
         call sort8 (nlist,list)
         if (nlist .le. maxp13) then
            np13(i) = nlist
            do j = 1, nlist
               ip13(j,i) = list(j)
            end do
         else
            write (iout,80)
   80       format (/,' POLARGRP  --  Too many Atoms',
     &                 ' in 1-3 Polarization Group')
            abort = .true.
            goto 90
         end if
      end do
   90 continue
c
c     loop over atoms finding all the 1-4 group relationships
c
      do i = 1, n
         mask(i) = 0
      end do
      do i = 1, n
         do j = 1, np11(i)
            jj = ip11(j,i)
            mask(jj) = i
         end do
         do j = 1, np12(i)
            jj = ip12(j,i)
            mask(jj) = i
         end do
         do j = 1, np13(i)
            jj = ip13(j,i)
            mask(jj) = i
         end do
         nlist = 0
         do j = 1, np13(i)
            jj = ip13(j,i)
            do k = 1, np12(jj)
               kk = ip12(k,jj)
               if (mask(kk) .ne. i) then
                  nlist = nlist + 1
                  list(nlist) = kk
               end if
            end do
         end do
         call sort8 (nlist,list)
         if (nlist .le. maxp14) then
            np14(i) = nlist
            do j = 1, nlist
               ip14(j,i) = list(j)
            end do
         else
            write (iout,100)
  100       format (/,' POLARGRP  --  Too many Atoms',
     &                 ' in 1-4 Polarization Group')
            abort = .true.
            goto 110
         end if
      end do
  110 continue
c
c     perform deallocation of some local arrays
c
      deallocate (mask)
      return
      end
