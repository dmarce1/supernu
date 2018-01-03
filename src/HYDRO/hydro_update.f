




      subroutine hydro_update( t0, t1 )
      use hydromod
      use gridmod
      use mpimod
      implicit none
      real, intent(in) :: t0, t1

      real :: t, dt, avg, dif
      real :: tinv_max

      integer :: dm, dm_max
      integer :: i, j, k

      real :: U(hydro_nx,hydro_ny,hydro_nz,hydro_nf)
      real :: P(hydro_nx,hydro_ny,hydro_nz,hydro_nf)
      real :: dU(hydro_nx,hydro_ny,hydro_nz,hydro_nf)
      real :: UR(hydro_nx,hydro_ny,hydro_nz,hydro_nf)
      real :: UL(hydro_nx,hydro_ny,hydro_nz,hydro_nf)
      real :: Fv(hydro_nx,hydro_ny,hydro_nz,hydro_nf)
      real :: Fs(hydro_nx,hydro_ny,hydro_nz,hydro_nf)
      real :: A(hydro_nx,hydro_ny,hydro_nz)
      real :: scle(hydro_nx,hydro_ny,hydro_nz,3)
      real :: area(hydro_nx,hydro_ny,hydro_nz,3)
      real :: volinv(hydro_nx,hydro_ny,hydro_nz)
      real :: X(hydro_nx,hydro_ny,hydro_nz,3)
      real :: Xf(hydro_nx,hydro_ny,hydro_nz,3)
      real :: dX(hydro_nx,hydro_ny,hydro_nz,3)
      real :: kinR(hydro_nx,hydro_ny,hydro_nz)
      real :: kinL(hydro_nx,hydro_ny,hydro_nz)
      real :: velR(hydro_nx,hydro_ny,hydro_nz)
      real :: velL(hydro_nx,hydro_ny,hydro_nz)
      real :: einR(hydro_nx,hydro_ny,hydro_nz)
      real :: einL(hydro_nx,hydro_ny,hydro_nz)
      real :: tmp(hydro_nx,hydro_ny,hydro_nz)
      logical :: veldim(3)

      call gather_hydro
      t = t0

c     Get face X from grd_?arr
      do i = 1, hydro_nx
      do j = 1, hydro_ny
      do k = 1, hydro_nz
        Xf(i,j,k,1) = grd_xarr(i)
        Xf(i,j,k,2) = grd_yarr(j)
        Xf(i,j,k,3) = grd_zarr(k)
      enddo
      enddo
      enddo

c     Select dimensions where grid is in velocity space
      select case( grd_igeom )
        case(11 : 1)
          veldim(1) = .true. .and. grd_isvelocity
          veldim(2:3) = .false.
          Xf(:,:,:,1) = Xf(:,:,:,1)* t
        case(2)
          veldim((/1,3/)) = .true. .and. grd_isvelocity
          veldim(2) = .false.
          Xf(:,:,:,(/1,3/)) = Xf(:,:,:,(/1,3/)) * t
        case(3)
          veldim = .true. .and. grd_isvelocity
          Xf = Xf * t
      end select

c     Compute cell centered X and dX from face X
      X(1:hydro_nx-1,:,:,1) =
     &  (Xf(2:hydro_nx,:,:,1) + Xf(1:hydro_nx-1,:,:,1))*0.5d0
      X(:,1:hydro_ny-1,:,2) =
     &  (Xf(:,2:hydro_ny,:,2) + Xf(:,1:hydro_ny-1,:,2))*0.5d0
      X(:,:,1:hydro_nz-1,3) =
     &  (Xf(:,:,2:hydro_nz,3) + Xf(:,:,1:hydro_nz-1,3))*0.5d0
      dX(1:hydro_nx-1,:,:,1) =
     &  (Xf(2:hydro_nx,:,:,1) - Xf(1:hydro_nx-1,:,:,1))
      dX(:,1:hydro_ny-1,:,2) =
     &  (Xf(:,2:hydro_ny,:,2) - Xf(:,1:hydro_ny-1,:,2))
      dX(:,:,1:hydro_nz-1,3) =
     &  (Xf(:,:,2:hydro_nz,3) - Xf(:,:,1:hydro_nz-1,3))

c     Special case of 1D spherical
      if( grd_igeom .eq. 11 ) then
        dm_max = 1
      else
        dm_max = 3
      endif

c     Compute geometrical scale, face area, and inverse cell volumes
c     (in units of dX)
      select case( grd_igeom )
        case(1,11)
         scle(:,:,:,1) = 1.0d0
         scle(:,:,:,2) = X(:,:,:,1)
         scle(:,:,:,3) = X(:,:,:,1) * sin(X(:,:,:,2))
         area(:,:,:,1) = Xf(:,:,:,1)**2 * sin(Xf(:,:,:,2))
         area(:,:,:,2) = Xf(:,:,:,1)
         area(:,:,:,3) = Xf(:,:,:,1) * sin(Xf(:,:,:,2))
        case(2)
         scle(:,:,:,(/1,3/)) = 1.0d0
         scle(:,:,:,2) = X(:,:,:,1)
         area(:,:,:,1) = Xf(:,:,:,1)
         area(:,:,:,2) = 1.0d0
         area(:,:,:,3) = Xf(:,:,:,1)
        case(3)
         scle = 1.0d0
         area = 1.0d0
      end select
      volinv(:,:,:) = 1.0d0/(scle(:,:,:,1)*scle(:,:,:,2)*scle(:,:,:,3))

c     Main loop - loop until desired time is reached
      do while (t .lt. t1)

        dU = 0.0d0
        tinv_max = 0.0d0
        U = hydro_state

c     Boundaries
c     TODO: NEEDS COORD SYS SPECIALIZATION, THIS IS FOR CARTESIAN ONLY
        do i = 1, hydro_bw
          U(i, :, :, :) = U(hydro_bw + 1, :, :, :)
          U(:, i, :, :) = U(:, hydro_bw + 1, :, :)
          U(:, :, i, :) = U(:, :, hydro_bw + 1, :)
          U(hydro_nx - i + 1, :, :, :) = U(hydro_nx - hydro_bw, :, :, :)
          U(:, hydro_ny - i + 1, :, :) = U(:, hydro_ny - hydro_bw, :, :)
          U(:, :, hydro_nz - i + 1, :) = U(:, :, hydro_nz - hydro_bw, :)
          U(i, :, :, px_i) = min( U(i, :, :, px_i), 0.0d0 )
          U(:, i, :, py_i) = min( U(i, :, :, py_i), 0.0d0 )
          U(:, :, i, pz_i) = min( U(i, :, :, pz_i), 0.0d0 )
          U(hydro_nx - i + 1, :, :, px_i) =
     &      max( U(hydro_nx - i + 1, :, :, px_i), 0.d0 )
          U(:, hydro_ny - i + 1, :, py_i) =
     &      max( U(:, hydro_ny - i + 1, :, py_i), 0.0d0 )
          U(:, :, hydro_nz - i + 1, pz_i) =
     &      max( U(:, :, hydro_nz - i + 1, pz_i), 0.0d0 )
        enddo

c     Compute contribution to dudt in each flux direction
        do dm = 1, dm_max

c     Reconstruct face values (piecewise constant only)
c     TODO: HIGHER ORDER RECONSTRUCTIONS
          if( hydro_bw .eq. 1 ) then
            UR = hydro_state
            select case( dm )
            case(1)
              UL(2:hydro_nx,:,:,:) = U(1:hydro_nx-1,:,:,:)
            case(2)
              UL(:,2:hydro_ny,:,:) = U(:,1:hydro_ny-1,:,:)
            case(3)
              UL(:,:,2:hydro_nz,:) = U(:,:,1:hydro_nz-1,:)
            end select
          endif
c     Compute face kinetic and internal energies and velocities
          kinL = (UL(:,:,:,px_i)**2+UL(:,:,:,py_i)**2+UL(:,:,:,pz_i)**2)
     &                / UL(:,:,:,rho_i)
          kinR = (UR(:,:,:,px_i)**2+UR(:,:,:,py_i)**2+UR(:,:,:,pz_i)**2)
     &                / UR(:,:,:,rho_i)
          einL = UL(:,:,:,egas_i) - kinL
          einR = UR(:,:,:,egas_i) - kinR
          velL = UL(:,:,:,px_i+dm-1) / UL(:,:,:,rho_i)
          velR = UR(:,:,:,px_i+dm-1) / UR(:,:,:,rho_i)
c     Remove grid velocity
          if( veldim(dm) ) then
            velL = velL - Xf(:,:,:,dm) / t
            velR = velR - Xf(:,:,:,dm) / t
          endif
c     Apply dual energy formalism to compute final value for internal
c     energy
          do i = 1, hydro_nx
          do j = 1, hydro_ny
          do k = 1, hydro_nz
            if( einL(i,j,k) .lt. UL(i,j,k,egas_i) * 0.001d0 ) then
              einL(i,j,k) = UL(i,j,k,tau_i)**hydro_gamma
            endif
            if( einR(i,j,k) .lt. UR(i,j,k,egas_i) * 0.001d0 ) then
              einR(i,j,k) = UR(i,j,k,tau_i)**hydro_gamma
            endif
          enddo
          enddo
          enddo

c     Compute signal speeds
          A(:,:,:) = max(
     &            sqrt((hydro_gamma-1.0d0)*einL) + abs( velL ),
     &            sqrt((hydro_gamma-1.0d0)*einR) + abs( velR )
     &                     )

c     Compute maximum dt inverse
          tinv_max = max(tinv_max,maxval(A/dX(:,:,:,dm)))

c     Compute advection terms for vector flux
          do i = 1, hydro_nf
            Fv(:,:,:,i) = ((velL * UL(:,:,:,i) + velR * UR(:,:,:,i)) -
     &                         A * (UR(:,:,:,i) - UL(:,:,:,i))) * 0.5d0
          enddo

c     Add work term for energy
          Fv(:,:,:,egas_i) = Fv(:,:,:,egas_i) +
     &     0.5d0*(hydro_gamma-1.0d0)*(einL * velL + einR * velR)

          if( veldim(dm) ) then
            Fv(:,:,:,egas_i) = Fv(:,:,:,egas_i) +
     &               0.5d0*(hydro_gamma-1.0d0)*(einL * Xf(:,:,:,dm) +
     &                                          einR * Xf(:,:,:,dm)) / t
          endif

c     Compute scalar fluxes
          Fs = 0.0d0
          Fs(:,:,:,px_i+dm-1) = 0.5d0*(hydro_gamma-1.0d0)*(einL + einR)

c     Apply face areas to vector fluxes
          do i = 1, hydro_nf
            Fv(:,:,:,i) = Fv(:,:,:,i) * area(:,:,:,dm)
          enddo

c     Add flux contribution to dudt
          do i = 1, hydro_nf
            select case( dm )
            case(1)
              dU(:,:,:,i) = dU(:,:,:,i) -
     &          (Fv(2:hydro_nx,:,:,i) - Fv(1:hydro_nx-1,:,:,i)) /
     &            dX(1:hydro_nx-1,:,:,1) * volinv(1:hydro_nx-1,:,:)
     &         -
     &          (Fs(2:hydro_nx,:,:,i) - Fs(1:hydro_nx-1,:,:,i)) /
     &            (dX(1:hydro_nx-1,:,:,1) * scle(1:hydro_nx-1,:,:,1))
            case(2)
              dU(:,:,:,i) = dU(:,:,:,i) -
     &          (Fv(:,2:hydro_ny,:,i) - Fv(:,1:hydro_ny-1,:,i)) /
     &            dX(:,1:hydro_ny-1,:,2) * volinv(:,1:hydro_ny-1,:)
     &         -
     &          (Fs(:,2:hydro_ny,:,i) - Fs(:,1:hydro_ny-1,:,i)) /
     &            (dX(:,1:hydro_ny-1,:,2) * scle(:,1:hydro_ny-1,:,2))
            case(3)
              dU(:,:,:,i) = dU(:,:,:,i) -
     &          (Fv(:,:,2:hydro_nz,i) - Fv(:,:,1:hydro_nz-1,i)) /
     &            dX(:,:,1:hydro_nz-1,3) * volinv(:,:,1:hydro_nz-1)
     &         -
     &          (Fs(:,:,2:hydro_nz,i) - Fs(:,:,1:hydro_nz-1,i)) /
     &            (dX(:,:,1:hydro_nz-1,3) * scle(:,:,1:hydro_nz-1,3))
            end select
          enddo

        enddo

c       Geometrical source terms
        if( grd_igeom .ne. 3 ) then
            tmp = 1.0d0 / (U(:,:,:,rho_i) * X(:,:,:,1))
            dU(:,:,:,px_i) = dU(:,:,:,px_i) +
     &        (U(:,:,:,py_i) * U(:,:,:,py_i)) * tmp
            dU(:,:,:,py_i) = dU(:,:,:,py_i) -
     &        (U(:,:,:,px_i) * U(:,:,:,py_i)) * tmp
        endif
        if( (grd_igeom .eq. 1) .or. (grd_igeom .eq. 11) ) then
            dU(:,:,:,px_i) = dU(:,:,:,px_i) +
     &        (U(:,:,:,pz_i) * U(:,:,:,pz_i)) * tmp
            dU(:,:,:,py_i) = dU(:,:,:,py_i) -
     &        (U(:,:,:,py_i) * U(:,:,:,pz_i)) / tan(X(:,:,:,3)) * tmp
            dU(:,:,:,pz_i) = dU(:,:,:,px_i) +
     &        (U(:,:,:,px_i) * U(:,:,:,pz_i)) * tmp
            dU(:,:,:,pz_i) = dU(:,:,:,py_i) -
     &        (U(:,:,:,py_i) * U(:,:,:,py_i)) / tan(X(:,:,:,3)) * tmp
        endif

c     If grid is moving, volume increases
      if( grd_isvelocity ) then
        dU = dU - 3.0d0 * U / t
      endif

c     Compute timestep
        dt = 0.4d0 / tinv_max
        dt = min( dt, t1 - t )

c     Apply dudt
        U = U + dU * dt

c     Upate X, Xf, and dX for moving grids
        do i = 1, 3
          if( veldim(i) ) then
            X = X * (1.0d0 + dt / t)
            Xf = Xf * (1.0d0 + dt / t)
            dX = dX * (1.0d0 + dt / t)
          endif
        enddo

      enddo

      hydro_state = U

      call scatter_hydro

      end subroutine


