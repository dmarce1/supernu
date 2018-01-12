      subroutine velinterp( v, dvdx, nx, ny, nz )
      use hydromod
      use gridmod
      use timestepmod
      implicit none

      integer, intent(in) :: nx, ny, nz
      real*8, intent(out) :: v(nx,ny,nz,3)
      real*8, intent(out) :: dvdx(nx,ny,nz,3,3)
      real*8 :: vc(nx+1,ny+1,nz+1)
      real*8 :: u(0:nx+1,0:ny+1,0:nz+1)

      integer :: xb, xe, yb, ye, zb, ze, dm, i
      real :: dxinv, tfact
      xb = hydro_bw
      yb = xb
      zb = xb
      xe = hydro_nx - hydro_bw + 1
      ye = hydro_ny - hydro_bw + 1
      ze = hydro_nz - hydro_bw + 1

      do dm = 1, 3
        u(0:nx+1,0:ny+1,0:nz+1) =
     &             hydro_state(xb:xe,yb:ye,zb:ze,px_i-1+dm) /
     &             hydro_state(xb:xe,yb:ye,zb:ze,rho_i)
        vc(1:nx+1,1:ny+1,1:nz+1) =
     &         (u(0:nx+0,0:ny+0,0:nz+0) +
     &          u(0:nx+0,0:ny+0,1:nz+1) +
     &          u(0:nx+0,1:ny+1,0:nz+0) +
     &          u(0:nx+0,1:ny+1,1:nz+1) +
     &          u(1:nx+1,0:ny+0,0:nz+0) +
     &          u(1:nx+1,0:ny+0,1:nz+1) +
     &          u(1:nx+1,1:ny+1,0:nz+0) +
     &          u(1:nx+1,1:ny+1,1:nz+1)) * 0.125d0
        v(:,:,:,dm) =
     &         (vc(1:nx+0,1:ny+0,1:nz+0) +
     &          vc(1:nx+0,1:ny+0,2:nz+1) +
     &          vc(1:nx+0,2:ny+1,1:nz+0) +
     &          vc(1:nx+0,2:ny+1,2:nz+1) +
     &          vc(2:nx+1,1:ny+0,1:nz+0) +
     &          vc(2:nx+1,1:ny+0,2:nz+1) +
     &          vc(2:nx+1,2:ny+1,1:nz+0) +
     &          vc(2:nx+1,2:ny+1,2:nz+1)) * 0.125d0
        dvdx(:,:,:,dm,1) =
     &         ((-vc(1:nx+0,1:ny+0,1:nz+0)) +
     &          (-vc(1:nx+0,1:ny+0,2:nz+1)) +
     &          (-vc(1:nx+0,2:ny+1,1:nz+0)) +
     &          (-vc(1:nx+0,2:ny+1,2:nz+1)) +
     &          (+vc(2:nx+1,1:ny+0,1:nz+0)) +
     &          (+vc(2:nx+1,1:ny+0,2:nz+1)) +
     &          (+vc(2:nx+1,2:ny+1,1:nz+0)) +
     &          (+vc(2:nx+1,2:ny+1,2:nz+1))) * 0.25d0
        dvdx(:,:,:,dm,2) =
     &         ((-vc(1:nx+0,1:ny+0,1:nz+0)) +
     &          (-vc(1:nx+0,1:ny+0,2:nz+1)) +
     &          (+vc(1:nx+0,2:ny+1,1:nz+0)) +
     &          (+vc(1:nx+0,2:ny+1,2:nz+1)) +
     &          (-vc(2:nx+1,1:ny+0,1:nz+0)) +
     &          (-vc(2:nx+1,1:ny+0,2:nz+1)) +
     &          (+vc(2:nx+1,2:ny+1,1:nz+0)) +
     &          (+vc(2:nx+1,2:ny+1,2:nz+1))) * 0.25d0
        dvdx(:,:,:,dm,3) =
     &         ((-vc(1:nx+0,1:ny+0,1:nz+0)) +
     &          (+vc(1:nx+0,1:ny+0,2:nz+1)) +
     &          (-vc(1:nx+0,2:ny+1,1:nz+0)) +
     &          (+vc(1:nx+0,2:ny+1,2:nz+1)) +
     &          (-vc(2:nx+1,1:ny+0,1:nz+0)) +
     &          (+vc(2:nx+1,1:ny+0,2:nz+1)) +
     &          (-vc(2:nx+1,2:ny+1,1:nz+0)) +
     &          (+vc(2:nx+1,2:ny+1,2:nz+1))) * 0.25d0
      enddo

      if( grd_isvelocity ) then
        tfact = tsp_t
      else
        tfact = 1.0d0
      endif

      do i = 1, nx
        dxinv = 1.0d0/(grd_xarr(i+1) - grd_xarr(i))*tfact
        dvdx(i,:,:,:,1) = dvdx(i,:,:,:,1) * dxinv
      enddo
      do i = 1, ny
        dxinv = 1.0d0/(grd_yarr(i+1) - grd_yarr(i))*tfact
        dvdx(:,i,:,:,2) = dvdx(:,i,:,:,2) * dxinv
      enddo
      do i = 1, nz
        dxinv = 1.0d0/(grd_zarr(i+1) - grd_zarr(i))*tfact
        dvdx(:,:,i,:,2) = dvdx(:,:,i,:,3) * dxinv
      enddo


      end subroutine velinterp
