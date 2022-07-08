      subroutine cttrans
      use sizes
      use atoms
      use bound
      use hescut
      use keys
      use limits
      use neigh
      use polpot
      use tarray
      use chargetransfer
      implicit none
      integer i,next
      integer limit
      real*8 big,value
      logical truncate
      character*20 keyword
      character*120 record
      character*120 string

      use_chargetranfer = .false.

      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         string = record(next:120)
c
c     get values related to use of Ewald summation
c
         if (keyword(1:13) .eq. 'CHARGETRANSF ') then
            use_chargetranfer = .true.
   10    continue
         end if
      end do
c      print*,"use_chargetranfer=",use_chargetranfer
      return
      end
