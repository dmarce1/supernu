




      subroutine hydro_update( t0, t1 )
      use hydromod
      use gridmod
      implicit none
      real, intent(in) :: t0, t1

      real :: t, dt, avg, dif
      real :: tinv_max

      integer :: dm
      integer :: i, j, k

      real :: U(hydro_nx,hydro_ny,hydro_nz,hydro_nf)
      real :: P(hydro_nx,hydro_ny,hydro_nz,hydro_nf)
      real :: dU(hydro_nx,hydro_ny,hydro_nz,hydro_nf)
      real :: UR(hydro_nx,hydro_ny,hydro_nz,hydro_nf)
      real :: UL(hydro_nx,hydro_ny,hydro_nz,hydro_nf)
      real :: F(hydro_nx,hydro_ny,hydro_nz,hydro_nf)
      real :: A(hydro_nx,hydro_ny,hydro_nz)
      real :: X(hydro_nx,hydro_ny,hydro_nz,3)
      real :: Xf(hydro_nx,hydro_ny,hydro_nz,3)
      real :: dX(hydro_nx,hydro_ny,hydro_nz,3)
      real :: kinR(hydro_nx,hydro_ny,hydro_nz)
      real :: kinL(hydro_nx,hydro_ny,hydro_nz)
      real :: velR(hydro_nx,hydro_ny,hydro_nz)
      real :: velL(hydro_nx,hydro_ny,hydro_nz)
      real :: einR(hydro_nx,hydro_ny,hydro_nz)
      real :: einL(hydro_nx,hydro_ny,hydro_nz)

      dU = 0.0d0
      tinv_max = 0.0d0

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
     &    max( U(hydro_nx - i + 1, :, :, px_i), 0.d0 )
        U(:, hydro_ny - i + 1, :, py_i) =
     &    max( U(:, hydro_ny - i + 1, :, py_i), 0.0d0 )
        U(:, :, hydro_nz - i + 1, pz_i) =
     &    max( U(:, :, hydro_nz - i + 1, pz_i), 0.0d0 )
      enddo

      do i = 1, hydro_nx
      do j = 1, hydro_ny
      do k = 1, hydro_nz

        Xf(i,j,k,1) = grd_xarr(i)
        Xf(i,j,k,2) = grd_yarr(j)
        Xf(i,j,k,3) = grd_zarr(k)

        avg = (grd_xarr(i+1)+grd_xarr(i))*0.5d0
        dif = (grd_xarr(i+1)-grd_xarr(i))
         X(i,j,k,1) = avg
        dX(i,j,k,1) = dif

        avg = (grd_yarr(j+1)+grd_yarr(j))*0.5d0
        dif = (grd_yarr(j+1)-grd_yarr(j))
         X(i,j,k,2) = avg
        dX(i,j,k,2) = dif

        avg = (grd_zarr(k+1)+grd_zarr(k))*0.5d0
        dif = (grd_zarr(k+1)-grd_zarr(k))
         X(i,j,k,3) = avg
        dX(i,j,k,3) = dif

      enddo
      enddo
      enddo

      if( grd_isvelocity ) then
        X = X * t
        dX = dX * t
      endif

      call gather_hydro

      U = hydro_state

      t = t0

      do while (t .lt. t1)

        do dm = 1, 3

          if( hydro_bw .eq. 1 ) then
            UR = hydro_state
            if( dm .eq. 1 ) then
              UL(2:hydro_nx,:,:,:) = U(1:hydro_nx-1,:,:,:)
            else if( dm .eq. 2 ) then
              UL(:,2:hydro_ny,:,:) = U(:,1:hydro_ny-1,:,:)
            else
              UL(:,:,2:hydro_nz,:) = U(:,:,1:hydro_nz-1,:)
            endif
          endif

          kinL = (UL(:,:,:,px_i)**2+UL(:,:,:,py_i)**2+UL(:,:,:,pz_i)**2)
     &                / UL(:,:,:,rho_i)
          kinR = (UR(:,:,:,px_i)**2+UR(:,:,:,py_i)**2+UR(:,:,:,pz_i)**2)
     &                / UR(:,:,:,rho_i)
          einL = UL(:,:,:,egas_i) - kinL
          einR = UR(:,:,:,egas_i) - kinR
          velL = UL(:,:,:,px_i+dm-1) / UL(:,:,:,rho_i)
          velR = UR(:,:,:,px_i+dm-1) / UR(:,:,:,rho_i)
          if( grd_isvelocity ) then
            velL = velL - Xf(:,:,:,dm) / t
            velR = velR - Xf(:,:,:,dm) / t
          endif
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

          A(:,:,:) = max(
     &            sqrt((hydro_gamma-1.0d0)*einL) + abs( velL ),
     &            sqrt((hydro_gamma-1.0d0)*einR) + abs( velR )
     &                     )

          tinv_max = max(tinv_max,maxval(A/dX(:,:,:,dm)))

          do i = 1, hydro_nf
            F(:,:,:,i) = ((velL * UL(:,:,:,i) + velR * UR(:,:,:,i)) -
     &                         A * (UR(:,:,:,i) - UL(:,:,:,i))) * 0.5d0
          enddo

          F(:,:,:,px_i+dm-1) = F(:,:,:,px_i+dm-1) +
     &     0.5d0*(hydro_gamma-1.0d0)*(einL + einR)

          F(:,:,:,egas_i) = F(:,:,:,egas_i) +
     &     0.5d0*(hydro_gamma-1.0d0)*(einL * velL + einR * velR)

          if( grd_isvelocity) then
            F(:,:,:,egas_i) = F(:,:,:,egas_i) +
     &               0.5d0*(hydro_gamma-1.0d0)*(einL * Xf(:,:,:,dm) +
     &                                          einR * Xf(:,:,:,dm)) / t
          endif

          do i = 1, hydro_nf
            if( dm .eq. 1 ) then
              dU(:,:,:,i) = dU(:,:,:,i) -
     &          (F(2:hydro_nx,:,:,i) - F(1:hydro_nx-1,:,:,i)) /
     &            dX(1:hydro_nx-1,:,:,1)
            else if( dm .eq. 2 ) then
              dU(:,:,:,i) = dU(:,:,:,i) -
     &          (F(:,2:hydro_ny,:,i) - F(:,1:hydro_ny-1,:,i)) /
     &            dX(:,1:hydro_ny-1,:,2)
            else
              dU(:,:,:,i) = dU(:,:,:,i) -
     &          (F(:,:,2:hydro_nz,i) - F(:,:,1:hydro_nz-1,i)) /
     &            dX(:,:,1:hydro_nz-1,3)
            endif
          enddo

        enddo

        if( grd_isvelocity ) then
          dU = dU + 3.0d0 * U / t
        endif

        dt = 0.4d0 / tinv_max
        dt = min( dt, t1 - t )

        U = U + dU * dt

        if( grd_isvelocity ) then
          X = X * (1.0d0 + dt / t)
          dX = dX * (1.0d0 + dt / t)
        endif

      enddo

      hydro_state = U

      call scatter_hydro

      end subroutine


