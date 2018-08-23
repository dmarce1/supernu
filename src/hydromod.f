*This file is part of SuperNu.  SuperNu is released under the terms of the GNU GPLv3, see COPYING.
*Copyright (c) 2013-2017 Ryan T. Wollaeger and Daniel R. van Rossum.  All rights reserved.
      module hydromod
c     --------------
      use gasmod
      implicit none
c
      real*8, parameter :: hydro_gamma = 5.0d0 / 3.0d0
c
      integer, parameter :: rho_i = 1
      integer, parameter :: px_i = 2
      integer, parameter :: py_i = 3
      integer, parameter :: pz_i = 4
      integer, parameter :: tau_i = 5
      integer, parameter :: egas_i = 6
      integer, parameter :: natom_i = 7
      integer, parameter :: nelec_i = 8
      integer, parameter :: frac_i = 9
      integer, parameter :: nfracs = 2 * gas_nchain + gas_nelem
      integer, parameter :: hydro_nf = 8 + nfracs
      integer :: hydro_nx, hydro_ny, hydro_nz, hydro_bw
      real(8), parameter :: des1 = 1d-3
      real(8), parameter :: des2 = 1d-2


      real*8, allocatable :: hydro_state(:,:,:,:)

      save
c
      contains

      subroutine hydromod_init
      use gridmod
      implicit none

      hydro_bw = 2
      hydro_nx = grd_nx + 2 * hydro_bw
      hydro_ny = grd_ny + 2 * hydro_bw
      hydro_nz = grd_nz + 2 * hydro_bw


      allocate(hydro_state(hydro_nx,hydro_ny,hydro_nz,hydro_nf))


      end subroutine

      subroutine hydromod_dealloc
      implicit none

      deallocate(hydro_state)

      end subroutine


      end module hydromod
c vim: fdm=marker
