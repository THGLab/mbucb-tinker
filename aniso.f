      subroutine aniso
      use sizes
      use atoms
      use bound
      use hescut
      use keys
      use limits
      use neigh
      use polpot
      use tarray
      use kpolr_aniso
      implicit none
      integer i,next
      integer limit
      real*8 big,value
      logical truncate
      character*20 keyword
      character*120 record
      character*120 string

      use_anisopolz = .false.

      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         string = record(next:120)
c
c     get values related to use of Ewald summation
c
         if (keyword(1:10) .eq. 'ANISOPOLZ ') then
            use_anisopolz = .true.
   10    continue
         end if
      end do
      print*,"use_anisopolz=",use_anisopolz
      return
      end
