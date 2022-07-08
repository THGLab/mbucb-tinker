c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #################################################################
c     ##                                                             ##
c     ##  module polar_aniso  --  anisotripoic polarizability  ##        
c     ##                                                             ##
c     #################################################################
c
c
c     anisopolarity  anisotropic dipole polarizability for each multipole site in the local coor. frame (Ang**3)
c     ranisopolarity anisotropic dipole polarizabilities rotated to the global coordinate system 
c
c
      module polar_aniso
      implicit none
      real*8, allocatable :: anisopolarity(:,:)
      real*8, allocatable :: ranisopolarity(:,:) 
      save
      end
