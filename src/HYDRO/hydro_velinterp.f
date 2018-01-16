
      subroutine hydro_velinterp
      use hydromod
      use gridmod
      use timestepmod
      implicit none

      integer :: nx, ny, nz
      real*8 :: vc(grd_nx+1,grd_ny+1,grd_nz+1)
      real*8 :: u(0:grd_nx+1,0:grd_ny+1,0:grd_nz+1)
      integer :: xb, xe, yb, ye, zb, ze, dm, i, j
      real :: dxinv, tfact, mu, r, dr, rsin0
      nx = grd_nx
      ny = grd_ny
      nz = grd_nz
      xb = hydro_bw
      yb = xb
      zb = xb
      xe = hydro_nx - hydro_bw + 1
      ye = hydro_ny - hydro_bw + 1
      ze = hydro_nz - hydro_bw + 1


      if( grd_isvelocity ) then
        tfact = tsp_t
      else
        tfact = 1.0d0
      endif

      if( grd_igeom .eq. 1 ) then
          u(0:nx+1,0:ny+1,0:nz+1) =
     &             hydro_state(xb:xe,yb:ye,zb:ze,px_i) /
     &             hydro_state(xb:xe,yb:ye,zb:ze,rho_i)
        grd_v(:,:,:,2:3) = 0.0d0
        grd_v(1:nx,:,:,1) =
     &           0.5d0*u(1:nx,1:ny,1:nz) +
     &          0.25d0*u(0:nx-1,1:ny,1:nz) +
     &          0.25d0*u(2:nx+1,1:ny,1:nz)
        do i = 1, nx
          r = (grd_xarr(i+1) + grd_xarr(i))/2.0d0 / tfact
          dr = (grd_xarr(i+1) - grd_xarr(i)) / tfact
          grd_dvdx(i,:,:,1,1) =
     &         (u(i+1,1:ny,1:nz) - u(i-1,1:ny,1:nz))*0.5d0 / dr
          grd_dvdx(i,:,:,2,2) = grd_v(i,:,:,1) / r
          grd_dvdx(i,:,:,3,3) = grd_v(i,:,:,1) / r
        enddo
        grd_dvdx(:,:,:,1,2) = 0.0d0
        grd_dvdx(:,:,:,1,3) = 0.0d0
        grd_dvdx(:,:,:,2,1) = 0.0d0
        grd_dvdx(:,:,:,2,3) = 0.0d0
        grd_dvdx(:,:,:,3,1) = 0.0d0
        grd_dvdx(:,:,:,3,2) = 0.0d0
      else
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

      endif

      end subroutine
