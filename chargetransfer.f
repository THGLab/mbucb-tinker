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
      module chargetransfer
      use sizes
      implicit none
      integer nctrans
      real*8, allocatable :: charge_trans(:)
      real*8 cqr_trans(maxtyp)
      real*8, allocatable :: beta_charge(:)
      real*8 rbeta_charge(maxtyp)
      real*8, allocatable :: thole_ct(:)
      real*8 rthole_ct(maxtyp)
      real*8, allocatable :: pdamp_ct(:)
      real*8, allocatable :: beta_damp_ct(:)
      integer pgrp_ct(maxval,maxtyp)
      real*8, allocatable :: uind_ct(:,:)
      real*8, allocatable :: uinp_ct(:,:)
      real*8, allocatable :: uinds_ct(:,:)
      real*8, allocatable :: uinps_ct(:,:)
      logical use_chargetranfer
      save
      end
