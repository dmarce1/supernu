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
      integer, parameter :: frac_i = 7
      integer, parameter :: nfracs = 2 * gas_nchain + gas_nelem
      integer, parameter :: hydro_nf = 6 + nfracs
      integer :: hydro_nx, hydro_ny, hydro_nz, hydro_bw

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

      pure subroutine hydro_velocity_at( x, y, z, vx, vy, vz, xi, yi,zi)
      use gridmod
      use inputparmod
      implicit none

      real*8, intent(in) :: x, y, z
      real*8, intent(out) ::vx, vy, vz
      integer, intent(in) :: xi, yi, zi

      real*8 :: dx, dy, dz

      if( in_hydro_on ) then
        if( grd_igeom .ne. 11 ) then
          dx = 0.5d0 * (x - grd_xarr(xi)) /
     &                 (grd_xarr(xi+1) - grd_xarr(xi))
          dy = 0.5d0 * (y - grd_yarr(yi)) /
     &                 (grd_yarr(yi+1) - grd_yarr(yi))
          dz = 0.5d0 * (z - grd_zarr(zi)) /
     &                 (grd_zarr(zi+1) - grd_zarr(zi))
          vx = grd_v(xi,yi,zi,1)
          vy = grd_v(xi,yi,zi,2)
          vz = grd_v(xi,yi,zi,3)
          vx = vx + grd_dvdx(xi,yi,zi,1,1) * dx
          vy = vy + grd_dvdx(xi,yi,zi,2,1) * dx
          vz = vz + grd_dvdx(xi,yi,zi,3,1) * dx
          vx = vx + grd_dvdx(xi,yi,zi,1,2) * dy
          vy = vy + grd_dvdx(xi,yi,zi,2,2) * dy
          vz = vz + grd_dvdx(xi,yi,zi,3,2) * dy
          vx = vx + grd_dvdx(xi,yi,zi,1,3) * dz
          vy = vy + grd_dvdx(xi,yi,zi,2,3) * dz
          vz = vz + grd_dvdx(xi,yi,zi,3,3) * dz
        else
          dx = 0.5d0 * (x - grd_xarr(xi)) /
     &                 (grd_xarr(xi+1) - grd_xarr(xi))
          vx = grd_v(xi,yi,zi,1)
          vx = vx + grd_dvdx(xi,yi,zi,1,1) * dx
          vy = 0.0d0
          vz = 0.0d0
        endif
      else
        vx = 0.0d0
        vy = 0.0d0
        vz = 0.0d0
       endif


      end subroutine


      end module hydromod
c vim: fdm=marker
