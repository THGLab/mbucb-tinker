c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  module penetrate  
c     ##                                                            ##
c     ################################################################
c
c
c
c
      module penetrate 
      use sizes
      implicit none
      real*8, allocatable :: pen(:)
      real*8 penr(maxtyp)
      logical usepiquemalcdq
      save
      end
