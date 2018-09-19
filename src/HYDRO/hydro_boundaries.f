      subroutine hydro_boundaries(U,X,veldim,nx,ny,nz,nf,bw,t)
      use hydromod
      use gridmod
      use mpimod
      implicit none

      logical, parameter :: allow_inflow = .false.

      real*8,dimension(nx,ny,nz,nf),intent(inout) :: U
      logical, intent(in) :: veldim(3)
      real*8,dimension(nx,ny,nz,3),intent(in) :: X
      integer,intent(in) :: nx, ny, nz, bw, nf

      integer :: i, j, k, dm
      real*8, intent(in) :: t
      real*8,dimension(nx,ny,nz,3) :: dummy2
      integer :: xb,xe,yb,ye,zb,ze

      real*8 :: dummy
      dummy = t
      dummy2 = X

      xb = hydro_bw+1
      yb = hydro_bw+1
      zb = hydro_bw+1
      xe = hydro_nx-bw
      ye = hydro_ny-bw
      ze = hydro_nz-bw

c     Boundaries

c       pre-bound
c          if(.not.allow_inflow) then
            do dm = 1, 3

            exit



              if( veldim(dm) ) then
                U(xb:xe,yb:ye,zb:ze,egas_i)
     &            = U(xb:xe,yb:ye,zb:ze,egas_i)-
     &            U(xb:xe,yb:ye,zb:ze,px_i+dm-1)**2 * 0.5d0 /
     &            U(xb:xe,yb:ye,zb:ze,rho_i)
                U(xb:xe,yb:ye,zb:ze,px_i+dm-1) =
     &            U(xb:xe,yb:ye,zb:ze,px_i+dm-1) -
     &            U(xb:xe,yb:ye,zb:ze,rho_i)*X(xb:xe,yb:ye,zb:ze,dm) / t
                U(xb:xe,yb:ye,zb:ze,egas_i)
     &              = U(xb:xe,yb:ye,zb:ze,egas_i) +
     &                  U(xb:xe,yb:ye,zb:ze,px_i+dm-1)**2 * 0.5d0 /
     &                U(xb:xe,yb:ye,zb:ze,rho_i)
              endif
            enddo
c          endif

          do i = 1, bw
            select case( grd_igeom )

c     1D Spherical
              case(11)
                U(i,:,:,:) = U(2*bw-i+1,:,:,:)
                U(nx - i + 1,:, :, :) = U(nx-bw,:,:,:)
                U(i,:,:,px_i) = -U(i,:,:,px_i)
                U(:,i,:,:) = U(:,bw+1,:,:)
                U(:,ny-i+1,:,:) = U(:,ny-bw,:,:)
                U(:,:,i,:) = U(:,:,bw+1,:)
                U(:,:,nz-i+1,:) = U(:,:,nz-bw,:)
                if(.not.allow_inflow) then
                  U(nx - i + 1, :, :, px_i) =
     &              max( U(nx - i + 1, :, :, px_i), 0.d0 )
                endif
c     Spherical
              case(1)
c     Radial singularity at center and outflow at edge
                do j = 1, ny
                  k = mod(j + (nz-2*bw) / 2 - 1, nz) + 1
                  U(i,j,:,:) = U(2 * bw - i + 1,k,:,:)
                  U(i,j,:,px_i) = -U(i,k,:,px_i)
                enddo
                U(nx - i + 1,:, :, :) = U(nx-bw,:,:,:)
                if(.not.allow_inflow) then
                  U(nx - i + 1, :, :, px_i) =
     &              max( U(nx - i + 1, :, :, px_i), 0.d0 )
                endif
c     Azimuthal periodic
                U(:, :, i, :) = U(:, :, nz - bw - 1 + i, :)
                U(:, :, nz - i + 1, :) = U(:, :, bw + i, :)
c     Theta direction
                do j = 1, nz
                  k = mod(j + (nz-2*bw)/2-1, nz) + 1
                  U(:, i, j, :) = U(:, 2 * bw - i + 1, k, :)
                  U(:, i, j, pz_i) = -U(:, i, k, pz_i)
                  U(:,ny-i+1,j,:) = U(:,ny+i-bw-1,k,:)
                  U(:,ny-i+1,j,py_i) = -U(:,ny+i-bw-1,k,py_i)
                enddo


c     Cylindrical
              case(2)
c     Radial singularity at center and outflow at edge
                do j = 1, nz
                  k = mod(j + (nz-2*bw) / 2 - 1, nz) + 1
                  U(i,:,j,:) = U(2*bw-i+1,:,k,:)
                  U(i,:,j,px_i) = -U(i,:,k,px_i)
                enddo
                U(nx-i+1,:,:,:) = U(nx-bw,:,:,:)
                if(.not.allow_inflow) then
                  U(nx-i+1,:,:,px_i) =
     &              max( U(nx - i + 1, :,:,px_i), 0.d0 )
                endif
c     Azimuthal periodic
                U(:,:,i,:) = U(:,:,nz-bw-1+i,:)
                U(:,:,nz-i+1,:) = U(:,:,bw+i,:)
c     Vertical outflow both directions
                U(:,i,:,:) = U(:,bw+1,:,:)
                U(:,ny-i+1,:,:) = U(:,ny-bw,:,:)
                if(.not.allow_inflow) then
                  U(:,i,:,py_i) = min( U(:,i,:,py_i), 0.0d0 )
                  U(:,ny-i+1,:,py_i) = max( U(:,ny-i+1,:,py_i), 0.0d0 )
                endif

c     Cartesian
              case(3)
c     All dims are outflow
                U(i,:,:,:) = U(bw + 1,:,:,:)
                U(nx-i+1,:,:,:) = U(nx-bw,:,:,:)
                U(:,i,:,:) = U(:,bw+1,:,:)
                U(:,ny-i+1,:,:) = U(:,ny-bw,:,:)
                U(:,:,i,:) = U(:,:,bw+1,:)
                U(:,:,nz-i+1,:) = U(:,:,nz-bw,:)
                if(.not.allow_inflow) then
                  U(i,:,:,px_i) = min( U(i,:,:,px_i), 0.0d0 )
                  U(:,i,:,py_i) = min( U(:,i,:,py_i), 0.0d0 )
                  U(:,:,i,pz_i) = min( U(:,:,i,pz_i), 0.0d0 )
                  U(nx-i+1,:,:,px_i) = max( U(nx-i+1,:,:,px_i), 0.0d0 )
                  U(:,ny-i+1,:,py_i) = max( U(:,ny-i+1,:,py_i), 0.0d0 )
                  U(:,:,nz-i+1,pz_i) = max( U(:,:,nz-i+1,pz_i), 0.0d0 )
                endif
            end select
          enddo

c       post-bound
c          if(.not.allow_inflow) then
            do dm = 1, 3

             exit

              if( veldim(dm) ) then
                U(:,:,:,egas_i) = U(:,:,:,egas_i) -
     &                   U(:,:,:,px_i+dm-1)**2 * 0.5d0 / U(:,:,:,rho_i)
                U(:,:,:,px_i+dm-1) = U(:,:,:,px_i+dm-1) +
     &                   U(:,:,:,rho_i)*X(:,:,:,dm) / t
                U(:,:,:,egas_i) = U(:,:,:,egas_i) +
     &                   U(:,:,:,px_i+dm-1)**2 * 0.5d0 / U(:,:,:,rho_i)
              endif
            enddo
c          endif



      end subroutine
