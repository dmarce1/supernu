

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




subroutine particle_dopfact( part1, part2, dopfact )
  use gridmod
  use particlemod
  use physconstmod
  implicit none
  type(packet),intent(in) :: part1
  type(packet2),intent(in) :: part2
  real*8,intent(out) :: dopfact
  real*8 :: eta, xi, mu, om, n_dot_dvdx(3), n_dot_dvdx_dot_n
  real*8 :: dvdx(3,3)
  integer :: i, j, k

  mu = part1%mu
  om = part1%om
  eta = sqrt(1d0-mu**2)*cos(om)
  xi = sqrt(1d0-mu**2)*sin(om)
  i = part2%ix
  j = part2%iy
  k = part2%iz

  dvdx = grd_dvdx(i,j,k,:,:)

  select case(grd_igeom)
    case(11)
      n_dot_dvdx_dot_n = 1.0d0
    case(1,2)
      n_dot_dvdx(1) = mu * dvdx(1,1) + eta * dvdx(2,1) + xi * dvdx(3,1)
      n_dot_dvdx(2) = mu * dvdx(1,2) + eta * dvdx(2,2) + xi * dvdx(3,2)
      n_dot_dvdx(3) = mu * dvdx(1,3) + eta * dvdx(2,3) + xi * dvdx(3,3)
      n_dot_dvdx_dot_n = mu * n_dot_dvdx(1) + eta * n_dot_dvdx(2) + xi * n_dot_dvdx(3)
    case(3)
      n_dot_dvdx(1) = mu * dvdx(3,1) + eta * dvdx(1,1) + xi * dvdx(2,1)
      n_dot_dvdx(2) = mu * dvdx(3,2) + eta * dvdx(1,2) + xi * dvdx(2,2)
      n_dot_dvdx(3) = mu * dvdx(3,3) + eta * dvdx(1,3) + xi * dvdx(2,3)
      n_dot_dvdx_dot_n = mu * n_dot_dvdx(3) + eta * n_dot_dvdx(1) + xi * n_dot_dvdx(2)
  end select



end subroutine
