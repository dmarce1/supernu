
subroutine particle_momentum( ptcl, mom )
  use particlemod
  use gridmod
  use physconstmod
  implicit none

  type(packet), intent(in) :: ptcl
  real*8, intent(out) :: mom(3)
  real*8 :: sin0, mu1, mu2, mu, om, mtot

  mu = ptcl%mu
  om = ptcl%om
  mtot = ptcl%e / pc_c
  if( grd_igeom .ne. 11 ) then
    sin0 = sqrt(1d0 - mu * mu)
    mu1 = sin0 * cos(om)
    mu2 = sin0 * sin(om)
  endif
  select case(grd_igeom)
    case(11)
      mom(2:3) = 0d0
      mom(1) = mtot * mu
    case(1)
      mom(1) = mtot * mu
      mom(2) = mtot * mu1
      mom(3) = mtot * mu2
    case(2)
      mom(2) = mtot * mu
      mom(1) = mtot * mu1
      mom(3) = mtot * mu2
    case(3)
      mom(3) = mtot * mu
      mom(1) = mtot * mu1
      mom(2) = mtot * mu2
  end select

end subroutine

