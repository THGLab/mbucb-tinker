c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  mutant.i  --  hybrid atoms for free energy perturbation  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     lambda     generic weighting between initial and final states
c     vlambda    state weighting value for electrostatic potentials
c     elambda    state weighting value for van der Waals potentials
c     scexp      scale factor for soft core buffered 14-7 potential
c     scalpha    scale factor for soft core buffered 14-7 potential
c     nmut       number of atoms mutated from initial to final state
c     imut       atomic sites differing in initial and final state
c     type0      atom type of each atom in the initial state system
c     class0     atom class of each atom in the initial state system
c     type1      atom type of each atom in the final state system
c     class1     atom class of each atom in the final state system
c     mut        true if an atom is to be mutated, false otherwise
c     mutintra   true if intramutant vdw are scaled, false otherwise
c
c
c      integer nmut,imut
c      integer type0,class0
c      integer type1,class1
c      real*8 lambda
c      real*8 vlambda,elambda
c      real*8 scexp,scalpha
c      logical mut
c      logical mutintra
c      common /mutant/ lambda,vlambda,elambda,scexp,scalpha,nmut,
c     &                imut(maxatm),type0(maxatm),class0(maxatm),
c     &                type1(maxatm),class1(maxatm),mut(maxatm),mutintra

      module mutant
      implicit none
      integer nmut
      integer, allocatable :: imut(:)
      integer, allocatable :: type0(:)
      integer, allocatable :: class0(:)
      integer, allocatable :: type1(:)
      integer, allocatable :: class1(:)
      real*8 lambda
      real*8 vlambda
      real*8 elambda
      real*8 scexp
      real*8 scalpha
      logical, allocatable :: mut(:)
      logical mutintra
      save
      end

