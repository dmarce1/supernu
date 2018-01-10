


      subroutine hydro_update( t0, t1 )
      use hydromod
      use gridmod
      use mpimod
      implicit none
      logical, parameter :: allow_inflow = .true.
      real(8), intent(in) :: t0, t1

      real(8) :: t, dt
      real(8) :: dtinv_max

      integer :: dm, rk
      integer :: i, j, k

      real(8) :: U(hydro_nx,hydro_ny,hydro_nz,hydro_nf)
      real(8) :: U0(hydro_nx,hydro_ny,hydro_nz,hydro_nf)
      real(8) :: dU(hydro_nx,hydro_ny,hydro_nz,hydro_nf)
      real(8) :: UR(hydro_nx,hydro_ny,hydro_nz,hydro_nf)
      real(8) :: UL(hydro_nx,hydro_ny,hydro_nz,hydro_nf)
      real(8) :: slp(hydro_nx,hydro_ny,hydro_nz,hydro_nf)
      real(8) :: slp_p(hydro_nx,hydro_ny,hydro_nz,hydro_nf)
      real(8) :: slp_m(hydro_nx,hydro_ny,hydro_nz,hydro_nf)
      real(8) :: Fv(hydro_nx,hydro_ny,hydro_nz,hydro_nf)
      real(8) :: Fs(hydro_nx,hydro_ny,hydro_nz,hydro_nf)
      real(8) :: A(hydro_nx,hydro_ny,hydro_nz)
      real(8) :: scle(hydro_nx,hydro_ny,hydro_nz,3)
      real(8) :: area(hydro_nx,hydro_ny,hydro_nz,3)
      real(8) :: volinv(hydro_nx,hydro_ny,hydro_nz)
      real(8) :: X(hydro_nx,hydro_ny,hydro_nz,3)
      real(8) :: Xf(hydro_nx+1,hydro_ny+1,hydro_nz+1,3)
      real(8) :: dX(hydro_nx,hydro_ny,hydro_nz,3)
      real(8) :: kin(hydro_nx,hydro_ny,hydro_nz)
      real(8) :: ein(hydro_nx,hydro_ny,hydro_nz)
      real(8) :: kinR(hydro_nx,hydro_ny,hydro_nz)
      real(8) :: kinL(hydro_nx,hydro_ny,hydro_nz)
      real(8) :: velR(hydro_nx,hydro_ny,hydro_nz)
      real(8) :: velL(hydro_nx,hydro_ny,hydro_nz)
      real(8) :: cR(hydro_nx,hydro_ny,hydro_nz)
      real(8) :: cL(hydro_nx,hydro_ny,hydro_nz)
      real(8) :: einR(hydro_nx,hydro_ny,hydro_nz)
      real(8) :: einL(hydro_nx,hydro_ny,hydro_nz)
      real(8) :: tmp(hydro_nx,hydro_ny,hydro_nz)
      logical :: veldim(3)
      logical :: dimused(3)
      integer :: xb, xe, yb, ye, zb, ze
      integer :: nx, ny, nz, bw, nf
      real(8) :: gamma
      logical :: done
      logical, save :: first_call = .true.

      call gather_hydro
      if( first_call ) then
        call hydro_output()
        first_call = .false.
      endif

      gamma = hydro_gamma
      nx = hydro_nx
      ny = hydro_ny
      nz = hydro_nz
      bw = hydro_bw
      nf = hydro_nf

      t = t0

      xb = bw + 1
      yb = bw + 1
      zb = bw + 1
      xe = nx - bw
      ye = ny - bw
      ze = nz - bw


c     Find out which dimensions are in use
      dimused(1) = grd_nx .gt. 1
      if( grd_igeom .eq. 11 ) then
        dimused(2) = .false.
        dimused(3) = .false.
      else
        dimused(2) = grd_ny .gt. 1
        dimused(3) = grd_nz .gt. 1
      endif

c     Get face X from grd_?arr
      do i = 1+bw, nx-bw+1
      do j = 1, ny
      do k = 1, nz
        Xf(i,j,k,1) = grd_xarr(i-bw)
      enddo
      enddo
      enddo
      do i = 1, nx
      do j = 1+bw, ny-bw+1
      do k = 1, nz
        Xf(i,j,k,2) = grd_yarr(j-bw)
      enddo
      enddo
      enddo
      do i = 1, nx
      do j = 1, ny
      do k = 1+bw, nz-bw+1
        Xf(i,j,k,3) = grd_zarr(k-bw)
      enddo
      enddo
      enddo

      do i = 0, bw-1
          Xf(bw-i,:,:,1) = 2.0d0*Xf(bw-i+1,:,:,1) - Xf(bw-i+2,:,:,1)
          Xf(:,bw-i,:,2) = 2.0d0*Xf(:,bw-i+1,:,2) - Xf(:,bw-i+2,:,2)
          Xf(:,:,bw-i,3) = 2.0d0*Xf(:,:,bw-i+1,3) - Xf(:,:,bw-i+2,3)
      enddo
      do i = 1, bw
        Xf(nx-bw+i+1,:,:,1) =
     &    2.0d0*Xf(nx-bw+i,:,:,1) - Xf(nx-bw+i-1,:,:,1)
        Xf(:,ny-bw+i+1,:,2) =
     &    2.0d0*Xf(:,ny-bw+i,:,2) - Xf(:,ny-bw+i-1,:,2)
        Xf(:,:,nz-bw+i+1,3) =
     &    2.0d0*Xf(:,:,nz-bw+i,3) - Xf(:,:,nz-bw+i-1,3)
      enddo

c     Select dimensions where grid is in velocity space
      select case( grd_igeom )
        case(11 , 1)
          veldim(1) = .true. .and. grd_isvelocity
          veldim(2:3) = .false.
          Xf(:,:,:,1) = Xf(:,:,:,1) * t
        case(2)
          veldim((/1,3/)) = .true. .and. grd_isvelocity
          veldim(2) = .false.
          Xf(:,:,:,(/1,3/)) = Xf(:,:,:,(/1,3/)) * t
        case(3)
          veldim = .true. .and. grd_isvelocity
          Xf = Xf * t
      end select

      veldim = veldim .and. dimused

c     Compute cell centered X and dX from face X

      X(:,:,:,1) = (Xf(2:nx+1,1:ny,1:nz,1)+Xf(1:nx,1:ny,1:nz,1))*0.5d0
      X(:,:,:,2) = (Xf(1:nx,2:ny+1,1:nz,2)+Xf(1:nx,1:ny,1:nz,2))*0.5d0
      X(:,:,:,3) = (Xf(1:nx,1:ny,2:nz+1,3)+Xf(1:nx,1:ny,1:nz,3))*0.5d0

      dX(:,:,:,1) = (Xf(2:nx+1,1:ny,1:nz,1)-Xf(1:nx,1:ny,1:nz,1))
      dX(:,:,:,2) = (Xf(1:nx,2:ny+1,1:nz,2)-Xf(1:nx,1:ny,1:nz,2))
      dX(:,:,:,3) = (Xf(1:nx,1:ny,2:nz+1,3)-Xf(1:nx,1:ny,1:nz,3))

c     Compute geometrical scale, face area, and inverse cell volumes
c     (in units of dX)
      select case( grd_igeom )
        case(1,11)
          tmp = acos(X(:,:,:,2))
          tmp = sin(tmp)
          scle(:,:,:,1) = 1.0d0
          scle(:,:,:,2) = abs(X(:,:,:,1)) * tmp
          scle(:,:,:,3) = abs(X(:,:,:,1))
          area(:,:,:,1) = Xf(1:nx,1:ny,1:nz,1)**2 * tmp
          area(:,:,:,2) = abs(Xf(1:nx,1:ny,1:nz,1))
          area(:,:,:,3) = abs(Xf(1:nx,1:ny,1:nz,1)) * tmp
          j = hydro_bw + 1
          k = j
        case(2)
          scle(:,:,:,(/1,3/)) = 1.0d0
          scle(:,:,:,2) = X(:,:,:,1)
          area(:,:,:,1) = Xf(1:nx,1:ny,1:nz,1)
          area(:,:,:,2) = 1.0d0
          area(:,:,:,3) = Xf(1:nx,1:ny,1:nz,1)
        case(3)
          scle = 1.0d0
          area = 1.0d0
      end select
      tmp = Xf(1:nx,1:ny,1:nz,2)

      volinv(1:nx-1,1:ny-1,1:nz-1) = 1.0d0 /
     &                               scle(1:nx-1,1:ny-1,1:nz-1,1) /
     &                               scle(1:nx-1,1:ny-1,1:nz-1,2) /
     &                               scle(1:nx-1,1:ny-1,1:nz-1,3)

      done = .false.
c     Main loop - loop until desired time is reached
      U = hydro_state

      do while (.not. done)
        U0 = U

        do rk = 1,bw
          dU = 0.0d0
          dtinv_max = 0.0d0

c     Boundaries

c       pre-bound
          if(.not.allow_inflow) then
            do dm = 1, 3
              if( veldim(dm) ) then
                U(:,:,:,egas_i) = U(:,:,:,egas_i) -
     &                  U(:,:,:,px_i+dm-1)**2 * 0.5d0 / U(:,:,:,rho_i)
                U(: ,:,:,px_i+dm-1) = U(:,:,:,px_i+dm-1) -
     &                   U(:,:,:,rho_i)*X(:,:,:,px_i+dm-1) / t
              endif
            enddo
          endif

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
                do j = 1, ny
                  k = mod(j + (ny-2*bw) / 2 - 1, ny) + 1
                  U(i,j,:,:) = U(2*bw-i+1,k,:,:)
                  U(i,j,:,px_i) = -U(i,k,:,px_i)
                enddo
                U(nx-i+1,:,:,:) = U(nx-bw,:,:,:)
                if(.not.allow_inflow) then
                  U(nx-i+1,:,:,px_i) =
     &              max( U(nx - i + 1, :,:,px_i), 0.d0 )
                endif
c     Azimuthal periodic
                U(:,i,:,:) = U(:,ny-bw-1+i,:,:)
                U(:,ny-i+1,:,:) = U(:,bw+i,:,:)
c     Vertical outflow both directions
                U(:,:,i,:) = U(:,:,bw+1,:)
                U(:,:,nz-i+1,:) = U(:,:,nz-bw,:)
                if(.not.allow_inflow) then
                  U(:,:,i,pz_i) = min( U(:,:,i,pz_i), 0.0d0 )
                  U(:,:,nz-i+1,pz_i) = max( U(:,:,nz-i+1,pz_i), 0.0d0 )
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
          if(.not.allow_inflow) then
            do dm = 1, 3
              if( veldim(dm) ) then
                U(:,:,:,px_i+dm-1) = U(:,:,:,px_i+dm-1) +
     &                   U(:,:,:,rho_i)*X(:,:,:,px_i+dm-1) / t
                U(:,:,:,egas_i) = U(:,:,:,egas_i) +
     &                   U(:,:,:,px_i+dm-1)**2 * 0.5d0 / U(:,:,:,rho_i)
              endif
            enddo
          endif

c     Compute contribution to dudt in each flux direction
          do dm = 1, 3
            if( dimused(dm) ) then

c     Reconstruct face values

c         pre-recon
              do i = 0, 2
                U(:,:,:,px_i+i) = U(:,:,:,px_i+i) / U(:,:,:,rho_i)
                if( veldim(i+1) ) then
                   U(:,:,:,px_i+i) = U(:,:,:,px_i+i) - X(:,:,:,i+1) / t
                endif
              end do

c     Piecewise constant
              if( bw .eq. 1 ) then
                UR = U
                select case( dm )
                  case(1)
                    UL(2:nx,:,:,:) = U(1:nx-1,:,:,:)
                  case(2)
                    UL(:,2:ny,:,:) = U(:,1:ny-1,:,:)
                  case(3)
                    UL(:,:,2:nz,:) = U(:,:,1:nz-1,:)
                end select
c     Piecwise linear
              else if( bw .eq. 2 ) then
                select case(dm)
                  case(1)
                    slp(2:nx,:,:,:) = U(2:nx,:,:,:) - U(1:nx-1,:,:,:)
                    slp_m = slp
                    slp_p(1:nx-1,:,:,:) = slp(2:nx,:,:,:)
                  case(2)
                    slp(:,2:ny,:,:) = U(:,2:ny,:,:) - U(:,1:ny-1,:,:)
                    slp_m = slp
                    slp_p(:,1:ny-1,:,:) = slp(:,2:ny,:,:)
                  case(3)
                    slp(:,:,2:nz,:) = U(:,:,2:nz,:) - U(:,:,1:nz-1,:)
                    slp_m = slp
                    slp_p(:,:,1:nz-1,:) = slp(:,:,2:nz,:)
                end select
                slp = (sign(0.25d0,slp_m)+sign(0.25d0,slp_p))*
     &                   min(abs(slp_p),abs(slp_m))
                select case(dm)
                  case(1)
                    UR(xb:xe+1,yb:ye,zb:ze,:) =
     &               U(xb:xe+1,yb:ye,zb:ze,:)-slp(xb:xe+1,yb:ye,zb:ze,:)
                    UL(xb:xe+1,yb:ye,zb:ze,:) =
     &               U(xb-1:xe,yb:ye,zb:ze,:)+slp(xb-1:xe,yb:ye,zb:ze,:)
                  case(2)
                    UR(xb:xe,yb:ye+1,zb:ze,:) =
     &               U(xb:xe,yb:ye+1,zb:ze,:)-slp(xb:xe,yb:ye+1,zb:ze,:)
                   UL(xb:xe,yb:ye+1,zb:ze,:) =
     &               U(xb:xe,yb-1:ye,zb:ze,:)+slp(xb:xe,yb-1:ye,zb:ze,:)
                  case(3)
                    UR(xb:xe,yb:ye,zb:ze+1,:) =
     &               U(xb:xe,yb:ye,zb:ze+1,:)-slp(xb:xe,yb:ye,zb:ze+1,:)
                    UL(xb:xe,yb:ye,zb:ze+1,:) =
     &               U(xb:xe,yb:ye,zb-1:ze,:)+slp(xb:xe,yb:ye,zb-1:ze,:)
                end select

                else
                  write(*,*) 'bw not supported by hydro'
                  call abort()
                endif

c         post-recon
              do i = 0, 2
                if( veldim(i+1) ) then
                   U (:,:,:,px_i+i) = U (:,:,:,px_i+i) + X (:,:,:,i+1)/t
                   UR(:,:,:,px_i+i) = UR(:,:,:,px_i+i) +
     &                                          Xf(1:nx,1:ny,1:nz,i+1)/t
                   UL(:,:,:,px_i+i) = UL(:,:,:,px_i+i) +
     &                                          Xf(1:nx,1:ny,1:nz,i+1)/t
                endif
                U (:,:,:,px_i+i) = U (:,:,:,px_i+i) * U (:,:,:,rho_i)
                UR(:,:,:,px_i+i) = UR(:,:,:,px_i+i) * UR(:,:,:,rho_i)
                UL(:,:,:,px_i+i) = UL(:,:,:,px_i+i) * UL(:,:,:,rho_i)
              enddo


c     Compute face kinetic and internal energies and velocities
              if( grd_igeom .eq. 11 ) then
                kinL = 0.5d0 * UL(:,:,:,px_i)**2 / UL(:,:,:,rho_i)
                kinR = 0.5d0 * UR(:,:,:,px_i)**2 / UR(:,:,:,rho_i)
              else
                kinL = 0.5d0*(UL(:,:,:,px_i)**2 + UL(:,:,:,py_i)**2
     &                    + UL(:,:,:,pz_i)**2) / UL(:,:,:,rho_i)
                kinR = 0.5d0*(UR(:,:,:,px_i)**2 + UR(:,:,:,py_i)**2
     &                    + UR(:,:,:,pz_i)**2) / UR(:,:,:,rho_i)
              endif
              einL = UL(:,:,:,egas_i) - kinL
              einR = UR(:,:,:,egas_i) - kinR
              velL = UL(:,:,:,px_i+dm-1) / UL(:,:,:,rho_i)
              velR = UR(:,:,:,px_i+dm-1) / UR(:,:,:,rho_i)
c     Remove grid velocity
              if( veldim(dm) ) then
                velL = velL - Xf(1:nx,1:ny,1:nz,dm) / t
                velR = velR - Xf(1:nx,1:ny,1:nz,dm) / t
              endif

c     Apply dual energy formalism to compute final value for internal
c     energy
              do i = 1, nx
              do j = 1, ny
              do k = 1, nz
                if( einL(i,j,k) .lt. UL(i,j,k,egas_i) * 0.001d0 ) then
                  einL(i,j,k) = UL(i,j,k,tau_i)**gamma
                endif
                if( einR(i,j,k) .lt. UR(i,j,k,egas_i) * 0.001d0 ) then
                  einR(i,j,k) = UR(i,j,k,tau_i)**gamma
                endif
              enddo
              enddo
              enddo
c     Compute signal speeds
              cL = sqrt((gamma-1.0d0)*gamma*einL/UL(:,:,:,rho_i))
              cR = sqrt((gamma-1.0d0)*gamma*einR/UR(:,:,:,rho_i))
              A = max(cL + abs( velL ),cR + abs( velR ))

c     Compute maximum dt inverse
              dtinv_max = max(dtinv_max,maxval(A(xb:xe,yb:ye,zb:ze) /
     &          scle(xb:xe,yb:ye,zb:ze, dm)/dX(xb:xe,yb:ye,zb:ze, dm)))
              select case(dm)
                case(1)
                  dtinv_max=max(dtinv_max,
     &                      maxval(A(xb+1:xe+1,yb:ye,zb:ze)/
     &           scle(xb:xe,yb:ye,zb:ze, dm)/dX(xb:xe,yb:ye,zb:ze, dm)))
                case(2)
                  dtinv_max=max(dtinv_max,
     &                      maxval(A(xb:xe,yb+1:ye+1,zb:ze)/
     &           scle(xb:xe,yb:ye,zb:ze, dm)/dX(xb:xe,yb:ye,zb:ze, dm)))
                case(3)
                  dtinv_max=max(dtinv_max,
     &                      maxval(A(xb:xe,yb:ye,zb+1:ze+1)/
     &           scle(xb:xe,yb:ye,zb:ze, dm)/dX(xb:xe,yb:ye,zb:ze, dm)))
              end select


c     Compute advection terms for vector flux
              do i = 1, nf
                Fv(:,:,:,i)=(velL * UL(:,:,:,i)+
     &                       velR * UR(:,:,:,i))*0.5d0
                Fv(:,:,:,i) = Fv(:,:,:,i) -
     &             A * (UR(:,:,:,i) - UL(:,:,:,i)) * 0.5d0
              enddo

c     Add work term for energy
              Fv(:,:,:,egas_i) = Fv(:,:,:,egas_i) +
     &         0.5d0*(gamma-1.0d0)*(einL * velL + einR * velR)

              if( veldim(dm) ) then
                Fv(:,:,:,egas_i) = Fv(:,:,:,egas_i) +
     &                 0.5d0*(gamma-1.0d0)*(einL*Xf(1:nx,1:ny,1:nz,dm) +
     &                                 einR * Xf(1:nx,1:ny,1:nz,dm)) / t
              endif


c     Compute scalar fluxes
              Fs = 0.0d0
              Fs(:,:,:,px_i+dm-1) =
     &                             0.5d0*(gamma-1.0d0)*(einL + einR)
c     Apply face areas to vector fluxes
              do i = 1, nf
                Fv(:,:,:,i) = Fv(:,:,:,i) * area(:,:,:,dm)
              enddo

c     Add flux contribution to dudt
              do i = 1, nf
                select case( dm )
                  case(1)
                    dU(1:nx-1,:,:,i) = dU(1:nx-1,:,:,i) -
     &                (Fv(2:nx,:,:,i) - Fv(1:nx-1,:,:,i)) /
     &                  dX(1:nx-1,:,:,1) * volinv(1:nx-1,:,:)
     &               -
     &                (Fs(2:nx,:,:,i) - Fs(1:nx-1,:,:,i)) /
     &                  (dX(1:nx-1,:,:,1) * scle(1:nx-1,:,:,1))
                  case(2)
                    dU(:,1:ny-1,:,i) = dU(:,1:ny-1,:,i) -
     &                (Fv(:,2:ny,:,i) - Fv(:,1:ny-1,:,i)) /
     &                  dX(:,1:ny-1,:,2) * volinv(:,1:ny-1,:)
     &               -
     &                (Fs(:,2:ny,:,i) - Fs(:,1:ny-1,:,i)) /
     &                  (dX(:,1:ny-1,:,2) * scle(:,1:ny-1,:,2))
                  case(3)
                    dU(:,:,1:nz-1,i) = dU(:,:,1:nz-1,i) -
     &                (Fv(:,:,2:nz,i) - Fv(:,:,1:nz-1,i)) /
     &                  dX(:,:,1:nz-1,3) * volinv(:,:,1:nz-1)
     &               -
     &                (Fs(:,:,2:nz,i) - Fs(:,:,1:nz-1,i)) /
     &                  (dX(:,:,1:nz-1,3) * scle(:,:,1:nz-1,3))
                end select
              enddo

            endif

          enddo
c       Geometrical source terms
          if( grd_igeom .ne. 11 ) then
            if( grd_igeom .ne. 3 ) then
              tmp = 1.0d0 / (U(:,:,:,rho_i) * X(:,:,:,1))
              dU(:,:,:,px_i) = dU(:,:,:,px_i) +
     &          (U(:,:,:,py_i) * U(:,:,:,py_i)) * tmp
              dU(:,:,:,py_i) = dU(:,:,:,py_i) -
     &          (U(:,:,:,px_i) * U(:,:,:,py_i)) * tmp
            endif
            if( grd_igeom .eq. 1 ) then
              dU(:,:,:,px_i) = dU(:,:,:,px_i) +
     &          (U(:,:,:,pz_i) * U(:,:,:,pz_i)) * tmp
              dU(:,:,:,py_i) = dU(:,:,:,py_i) -
     &          (U(:,:,:,py_i) * U(:,:,:,pz_i)) / tan(X(:,:,:,3)) * tmp
              dU(:,:,:,pz_i) = dU(:,:,:,px_i) +
     &          (U(:,:,:,px_i) * U(:,:,:,pz_i)) * tmp
              dU(:,:,:,pz_i) = dU(:,:,:,py_i) -
     &          (U(:,:,:,py_i) * U(:,:,:,py_i)) / tan(X(:,:,:,3)) * tmp
            endif
          endif


c     Compute timestep
          if( rk .eq. 1 ) then
            dtinv_max = max(
     &              maxval(dU(xb:xe,yb:ye,zb:ze,(/rho_i,tau_i/)) /
     &                      (-U(xb:xe,yb:ye,zb:ze,(/rho_i,tau_i/)))),
     &                  dtinv_max)
            dt = 0.4d0 / dtinv_max
            if( dt .ge. t1 - t ) then
              dt = t1 - t
              done = .true.
            endif
          endif

c     If grid is moving, volume increases
          if( grd_isvelocity ) then
             dU = dU - 3.0d0 * U / t
          endif

c     Apply dudt
          if( rk .eq. 1 ) then
            U = U0 + dU * dt
          else
            U = U0 + (dU + U - U0) * 0.5d0 * dt
          endif

          if( rk .eq. 1 ) then

            if(any(U(xb:xe,yb:ye,zb:ze,rho_i).le.0.0d0)) then
               write(*,*) 'rho < 0.0'
               call abort()
            endif
            if(any(U(xb:xe,yb:ye,zb:ze,tau_i).le.0.0d0)) then
               write(*,*) 'tau < 0.0'
            endif

c     Upate X, Xf, and dX for moving grids
            do i = 1, 3
              if( veldim(i) ) then
                X(:,:,:,i) = X(:,:,:,i) * (1.0d0 + dt / t)
                Xf(:,:,:,i) = Xf(:,:,:,i) * (1.0d0 + dt / t)
                dX(:,:,:,i) = dX(:,:,:,i) * (1.0d0 + dt / t)
              endif
            enddo

            t = t + dt
          endif

        enddo


c     Apply dual energy formalism to update tau
        kin = U(:,:,:,px_i)**2+U(:,:,:,py_i)**2+U(:,:,:,pz_i)**2
        kin = kin * 0.5d0 / U(:,:,:,rho_i)
        ein = U(:,:,:,egas_i) - kin
        do i = 1, nx
        do j = 1, ny
        do k = 1, nz
          if( ein(i,j,k) .ge. U(i,j,k,egas_i) * 0.1d0 ) then
            U(i,j,k,tau_i) = ein(i,j,k)**(1.0d0/gamma)
          endif
        enddo
        enddo
        enddo
        write(*,*) t, t1, dt, dtinv_max


      enddo


      hydro_state = U

      call scatter_hydro

      call hydro_output()


      end subroutine


