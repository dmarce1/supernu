      subroutine velinterp
      use hydromod
      use gridmod
      use timestepmod
      implicit none

      integer :: nx, ny, nz
      real*8 :: vc(grd_nx+1,grd_ny+1,grd_nz+1)
      real*8 :: u(0:grd_nx+1,0:grd_ny+1,0:grd_nz+1)
      integer :: xb, xe, yb, ye, zb, ze, dm, i
      real :: dxinv, tfact
      nx = grd_nx
      ny = grd_ny
      nz = grd_nz
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
        grd_v(:,:,:,dm) =
     &         (vc(1:nx+0,1:ny+0,1:nz+0) +
     &          vc(1:nx+0,1:ny+0,2:nz+1) +
     &          vc(1:nx+0,2:ny+1,1:nz+0) +
     &          vc(1:nx+0,2:ny+1,2:nz+1) +
     &          vc(2:nx+1,1:ny+0,1:nz+0) +
     &          vc(2:nx+1,1:ny+0,2:nz+1) +
     &          vc(2:nx+1,2:ny+1,1:nz+0) +
     &          vc(2:nx+1,2:ny+1,2:nz+1)) * 0.125d0
        grd_dvdx(:,:,:,dm,1) =
     &         ((-vc(1:nx+0,1:ny+0,1:nz+0)) +
     &          (-vc(1:nx+0,1:ny+0,2:nz+1)) +
     &          (-vc(1:nx+0,2:ny+1,1:nz+0)) +
     &          (-vc(1:nx+0,2:ny+1,2:nz+1)) +
     &          (+vc(2:nx+1,1:ny+0,1:nz+0)) +
     &          (+vc(2:nx+1,1:ny+0,2:nz+1)) +
     &          (+vc(2:nx+1,2:ny+1,1:nz+0)) +
     &          (+vc(2:nx+1,2:ny+1,2:nz+1))) * 0.25d0
        grd_dvdx(:,:,:,dm,2) =
     &         ((-vc(1:nx+0,1:ny+0,1:nz+0)) +
     &          (-vc(1:nx+0,1:ny+0,2:nz+1)) +
     &          (+vc(1:nx+0,2:ny+1,1:nz+0)) +
     &          (+vc(1:nx+0,2:ny+1,2:nz+1)) +
     &          (-vc(2:nx+1,1:ny+0,1:nz+0)) +
     &          (-vc(2:nx+1,1:ny+0,2:nz+1)) +
     &          (+vc(2:nx+1,2:ny+1,1:nz+0)) +
     &          (+vc(2:nx+1,2:ny+1,2:nz+1))) * 0.25d0
        grd_dvdx(:,:,:,dm,3) =
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
        grd_dvdx(i,:,:,:,1) = grd_dvdx(i,:,:,:,1) * dxinv
      enddo
      do i = 1, ny
        dxinv = 1.0d0/(grd_yarr(i+1) - grd_yarr(i))*tfact
        grd_dvdx(:,i,:,:,2) = grd_dvdx(:,i,:,:,2) * dxinv
      enddo
      do i = 1, nz
        dxinv = 1.0d0/(grd_zarr(i+1) - grd_zarr(i))*tfact
        grd_dvdx(:,:,i,:,2) = grd_dvdx(:,:,i,:,3) * dxinv
      enddo


      end subroutine velinterp
