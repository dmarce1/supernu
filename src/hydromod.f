*This file is part of SuperNu.  SuperNu is released under the terms of the GNU GPLv3, see COPYING.
*Copyright (c) 2013-2017 Ryan T. Wollaeger and Daniel R. van Rossum.  All rights reserved.
      module hydromod
c     --------------
      use gasmod
      implicit none
c
      logical :: hydro_ison = .true.
c
      integer, parameter :: px_i = 0
      integer, parameter :: py_i = 1
      integer, parameter :: pz_i = 2
      integer, parameter :: tau_i = 3
      integer, parameter :: egas_i = 4
      integer, parameter :: stable_i = 5
      integer, parameter :: frac_i = 6
      integer, parameter :: nfracs = 2 * gas_nchain + gas_nelem
      integer, parameter :: hydro_nf = 5 + nfracs
      integer :: hydro_nx, hydro_ny, hydro_nz, hydro_bw
      integer :: this_nx, this_ny, this_nz
      integer :: this_xb, this_yb, this_zb
      integer :: this_xe, this_ye, this_ze

      real*8, allocatable :: hydro_state(:,:,:,:)

      save
c
      contains

      subroutine hydromod_init
      use gridmod
      implicit none
      INCLUDE 'mpif.h'
      use mpimod

      hydro_bw = 1
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
