c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##############################################################
c     ##                                                          ##
c     ##  module kpolr  --  polarizability forcefield parameters  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     polr_aniso   anisotropic dipole polarizability parameters for each atom type
c     use_anisopolz logical governing whether or not to use anisotropic polarizabilities 
c
c
      module kpolr_aniso
      use sizes
      implicit none
      real*8 polr_aniso(9,maxtyp)
      logical use_anisopolz
      save
      end
