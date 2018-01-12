

subroutine particle_labfact( part1, part2, labfact )
  use gridmod
  use particlemod
  use physconstmod
  implicit none
  type(packet),intent(in) :: part1
  type(packet2),intent(in) :: part2
  real*8,intent(out) :: labfact
  real*8 :: eta, xi, udotn, mu, om
  real*8 :: vx, vy, vz, x0, y0, z0, x, y, z
  integer :: i, j, k

  mu = part1%mu
  om = part1%om
  eta = sqrt(1d0-mu**2)*cos(om)
  xi = sqrt(1d0-mu**2)*sin(om)
  i = part2%ix
  j = part2%iy
  k = part2%iz
  x = part1%x
  y = part1%y
  z = part1%z

  x0 = (grd_xarr(i+1) + grd_xarr(i))*0.5d0
  y0 = (grd_yarr(j+1) + grd_yarr(j))*0.5d0
  z0 = (grd_zarr(k+1) + grd_zarr(k))*0.5d0

  vx = grd_v(i,j,k,1) + grd_dvdx(i,j,k,1,1) * (x - x0) + grd_dvdx(i,j,k,1,2) * (y - y0) + grd_dvdx(i,j,k,1,3) * (z - z0)
  vy = grd_v(i,j,k,2) + grd_dvdx(i,j,k,2,1) * (x - x0) + grd_dvdx(i,j,k,2,2) * (y - y0) + grd_dvdx(i,j,k,2,3) * (z - z0)
  vz = grd_v(i,j,k,3) + grd_dvdx(i,j,k,3,1) * (x - x0) + grd_dvdx(i,j,k,3,2) * (y - y0) + grd_dvdx(i,j,k,3,3) * (z - z0)


  select case(grd_igeom)
    case(11)
      udotn = mu * vx
    case(1,2)
      udotn = mu * vx + eta * vy + xi * vz
    case(3)
      udotn = mu * vz + eta * vx + xi * vy
  end select

  labfact = 1.0d0 - udotn / pc_c

end subroutine
