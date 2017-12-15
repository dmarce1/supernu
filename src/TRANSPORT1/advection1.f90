!This file is part of SuperNu.  SuperNu is released under the terms of the GNU GPLv3, see COPYING.
!Copyright (c) 2013-2017 Ryan T. Wollaeger and Daniel R. van Rossum.  All rights reserved.
pure subroutine advection1(pretrans,ptcl,ptcl2)

  use miscmod
  use timestepmod
  use gridmod
  use particlemod
  use transportmod
  use inputparmod
  implicit none
  logical,intent(in) :: pretrans
  type(packet),target,intent(inout) :: ptcl
  type(packet2),target,intent(inout) :: ptcl2
!-----------------------------------------------------------------------
! This routine computes the advection of IMC particles through the
! velocity grid.  It is geometry dependent
!-----------------------------------------------------------------------
!-- advection split parameter
  real*8,parameter :: alph2 = .5d0
  logical,parameter :: partstopper = .true.
!
  integer :: zholder, zfdiff
  real*8 :: help
  integer :: i
!-- pointers
  integer,pointer :: ix, iy, iz, ic, ig
  real*8,pointer :: x
!-- statement function
  integer :: l
  real*8 :: dx,dz,xm,dyac,ym
  dx(l) = grd_xarr(l+1) - grd_xarr(l)
  dz(l) = grd_zarr(l+1) - grd_zarr(l)
  xm(l) = 0.5d0*(grd_xarr(l+1) + grd_xarr(l))
  dyac(l) = grd_yacos(l) - grd_yacos(l+1)
  ym(l) = sqrt(1d0-0.25d0*(grd_yarr(l+1)+grd_yarr(l))**2)

  ix => ptcl2%ix
  iy => ptcl2%iy
  iz => ptcl2%iz
  ic => ptcl2%ic
  ig => ptcl2%ig
  x => ptcl%x

!
!-- setting tentative new position
  if(pretrans) then
    x = x*tsp_t/(tsp_t+alph2*tsp_dt)
  else
    x = x*(tsp_t + alph2*tsp_dt)/(tsp_t+tsp_dt)
  endif

!
!-- nothing to do
  if (x>=grd_xarr(ix)) return

!
!-- finding tentative new index
  zholder = binsrch(x,grd_xarr,grd_nx+1,.false.)

!
!-- quick exit if DDMC is active
  if(in_puretran .or. .not.partstopper) then
    ix = zholder
    ic = grd_icell(ix,iy,iz)
    return
  endif

  zfdiff = -1
  if(grd_isvelocity) then
     help = tsp_t
  else
     help = 1d0
  endif
  do i=ix-1,zholder,-1
     l = grd_icell(i,iy,iz)
     if((grd_sig(l)+grd_cap(ig,l))* &
          min(dx(i),xm(i)*dyac(iy),xm(i)*ym(iy)*dz(iz))* &
          help>=trn_tauddmc) then
        zfdiff = i
        exit
     endif
  enddo

!--
  if(zfdiff.ne.-1) then
     ix = zfdiff+1
     x = grd_xarr(ix)
  else
     ix = zholder
  endif
  ic = grd_icell(ix,iy,iz)

end subroutine advection1
! vim: fdm=marker
