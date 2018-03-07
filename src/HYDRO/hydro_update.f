

      subroutine hydro_update( t0, t1, dt, dt_only )
      use hydromod
      use gridmod
      use mpimod
      implicit none
      logical, parameter :: allow_inflow = .true.
      real(8), intent(in) :: t0, t1
      real(8), intent(out) :: dt
      logical, intent(in) :: dt_only

      real(8) :: t
      real(8) :: dtinv_max

      integer :: dm, rk
      integer :: i, j, k, f, i0, j0, k0

      real(8), parameter :: des1 = 0.001d0
      real(8), parameter :: des2 = 0.1d0


      real(8) :: U(hydro_nx,hydro_ny,hydro_nz,hydro_nf)
      real(8) :: U0(hydro_nx,hydro_ny,hydro_nz,hydro_nf)
      real(8) :: dU(hydro_nx,hydro_ny,hydro_nz,hydro_nf)
      real(8) :: UR(hydro_nx,hydro_ny,hydro_nz,hydro_nf)
      real(8) :: UL(hydro_nx,hydro_ny,hydro_nz,hydro_nf)
      real(8) :: slp(hydro_nx,hydro_ny,hydro_nz)
      real(8) :: slp_p(hydro_nx,hydro_ny,hydro_nz)
      real(8) :: slp_m(hydro_nx,hydro_ny,hydro_nz)
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
      integer :: nx, ny, nz, bw, nf, l
      real(8) :: gamma, tfactor, tmp8
      logical :: done
      logical, save :: first_call = .true.

      if( first_call ) then
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
      if( grd_isvelocity ) then
        select case( grd_igeom )
          case(11 , 1)
            veldim(1) = .true.
            veldim(2:3) = .false.
            Xf(:,:,:,1) = Xf(:,:,:,1) * t
          case(2)
            veldim(1:2) = .true.
            veldim(3) = .false.
            Xf(:,:,:,1:2) = Xf(:,:,:,1:2) * t
          case(3)
            veldim = .true.
            Xf = Xf * t
        end select
      else
        veldim = .false.
      endif

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
        case(2)
          scle(:,:,:,1:2) = 1.0d0
          scle(:,:,:,3) = abs(X(:,:,:,1))
          area(:,:,:,1) = abs(Xf(1:nx,1:ny,1:nz,1))
          area(:,:,:,3) = 1.0d0
          area(:,:,:,2) = abs(Xf(1:nx,1:ny,1:nz,1))
        case(3)
          scle = 1.0d0
          area = 1.0d0
      end select

      volinv(1:nx-1,1:ny-1,1:nz-1) = 1.0d0 /
     &                               scle(1:nx-1,1:ny-1,1:nz-1,1) /
     &                               scle(1:nx-1,1:ny-1,1:nz-1,2) /
     &                               scle(1:nx-1,1:ny-1,1:nz-1,3)

      done = .false.
c     Main loop - loop until desired time is reached
      U = hydro_state

        write(*,*)
        write(*,*) t, t1
      do while (.not. done)

        U0 = U

        do rk = 1,bw

          dU = 0.0d0
          dtinv_max = 0.0d0

          call hydro_boundaries(U,X,veldim,nx,ny,nz,nf,bw,t)


c     Compute contribution to dudt in each flux direction
          do dm = 1, 3
            if( dimused(dm) ) then

c     Reconstruct face values

c         pre-recon
              do f = 1, nf
                if((f.ne.rho_i).and.(f.ne.tau_i).and.(f.ne.natom_i))then
                  if( f .eq. nelec_i .or. f .ge. frac_i ) then
                    l = natom_i
                  else
                    l = rho_i
                  endif
                  U(:,:,:,f) = U(:,:,:,f) / U(:,:,:,l)
                endif
              enddo
              U(:,:,:,egas_i) = U(:,:,:,egas_i) - 0.5*U(:,:,:,px_i)**2
              if( grd_igeom .ne. 11 ) then
               U(:,:,:,egas_i) = U(:,:,:,egas_i) - 0.5*U(:,:,:,py_i)**2
               U(:,:,:,egas_i) = U(:,:,:,egas_i) - 0.5*U(:,:,:,pz_i)**2
              endif
              do i = 0, 2
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
                do f = 1, nf
                 select case(dm)
                  case(1)
                    slp(2:nx,:,:) = (U(2:nx,:,:,f) - U(1:nx-1,:,:,f))/
     &                              (X(2:nx,:,:,1) - X(1:nx-1,:,:,1))
                    slp_m = slp
                    slp_p(1:nx-1,:,:) = slp(2:nx,:,:)
                  case(2)
                    slp(:,2:ny,:) = (U(:,2:ny,:,f) - U(:,1:ny-1,:,f))/
     &                              (X(:,2:ny,:,2) - X(:,1:ny-1,:,2))
                    slp_m = slp
                    slp_p(:,1:ny-1,:) = slp(:,2:ny,:)
                  case(3)
                    slp(:,:,2:nz) = (U(:,:,2:nz,f) - U(:,:,1:nz-1,f))/
     &                              (X(:,:,2:nz,3) - X(:,:,1:nz-1,3))
                    slp_m = slp
                    slp_p(:,:,1:nz-1) = slp(:,:,2:nz)
                 end select
                 slp = (sign(0.5d0,slp_m)+sign(0.5d0,slp_p))*
     &                   min(abs(slp_p),abs(slp_m))
                 select case(dm)
                  case(1)
                    UR(xb:xe+1,yb:ye,zb:ze,f) = U(xb:xe+1,yb:ye,zb:ze,f)
     &               - slp(xb:xe+1,yb:ye,zb:ze)
     &               *  dX(xb:xe+1,yb:ye,zb:ze,1)*0.5d0
                    UL(xb:xe+1,yb:ye,zb:ze,f) = U(xb-1:xe,yb:ye,zb:ze,f)
     &               + slp(xb-1:xe,yb:ye,zb:ze)
     &               *  dX(xb-1:xe,yb:ye,zb:ze,1)*0.5d0
                  case(2)
                    UR(xb:xe,yb:ye+1,zb:ze,f) = U(xb:xe,yb:ye+1,zb:ze,f)
     &               - slp(xb:xe,yb:ye+1,zb:ze)
     &               *  dX(xb:xe,yb:ye+1,zb:ze,2)*0.5d0
                    UL(xb:xe,yb:ye+1,zb:ze,f) = U(xb:xe,yb-1:ye,zb:ze,f)
     &               + slp(xb:xe,yb-1:ye,zb:ze)
     &               *  dX(xb:xe,yb-1:ye,zb:ze,2)*0.5d0
                  case(3)
                    UR(xb:xe,yb:ye,zb:ze+1,f) = U(xb:xe,yb:ye,zb:ze+1,f)
     &               - slp(xb:xe,yb:ye,zb:ze+1)
     &               *  dX(xb:xe,yb:ye,zb:ze+1,3)*0.5d0
                    UL(xb:xe,yb:ye,zb:ze+1,f) = U(xb:xe,yb:ye,zb-1:ze,f)
     &               + slp(xb:xe,yb:ye,zb-1:ze)
     &               *  dX(xb:xe,yb:ye,zb-1:ze,3)*0.5d0
                 end select
                enddo

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
              enddo
              U(:,:,:,egas_i) = U(:,:,:,egas_i) + 0.5*U(:,:,:,px_i)**2
              UR(:,:,:,egas_i) = UR(:,:,:,egas_i)+0.5*UR(:,:,:,px_i)**2
              UL(:,:,:,egas_i) = UL(:,:,:,egas_i)+0.5*UL(:,:,:,px_i)**2
              if( grd_igeom .ne. 11 ) then
               U(:,:,:,egas_i) = U(:,:,:,egas_i) + 0.5*U(:,:,:,py_i)**2
               U(:,:,:,egas_i) = U(:,:,:,egas_i) + 0.5*U(:,:,:,pz_i)**2
               UR(:,:,:,egas_i) = UR(:,:,:,egas_i)+0.5*UR(:,:,:,py_i)**2
               UR(:,:,:,egas_i) = UR(:,:,:,egas_i)+0.5*UR(:,:,:,pz_i)**2
               UL(:,:,:,egas_i) = UL(:,:,:,egas_i)+0.5*UL(:,:,:,py_i)**2
               UL(:,:,:,egas_i) = UL(:,:,:,egas_i)+0.5*UL(:,:,:,pz_i)**2
              endif
              do f = 1, nf
                if((f.ne.rho_i).and.(f.ne.tau_i).and.(f.ne.natom_i))then
                  if( f .eq. nelec_i .or. f .ge. frac_i ) then
                    l = natom_i
                  else
                    l = rho_i
                  endif
                  U(:,:,:,f) = U(:,:,:,f) * U(:,:,:,l)
                  UR(:,:,:,f) = UR(:,:,:,f) * UR(:,:,:,l)
                  UL(:,:,:,f) = UL(:,:,:,f) * UL(:,:,:,l)
                endif
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
              einL = max(UL(:,:,:,egas_i) - kinL,0d0)
              einR = max(UR(:,:,:,egas_i) - kinR,0d0)
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
                if( einL(i,j,k) .lt. UL(i,j,k,egas_i) * des1 ) then
                  einL(i,j,k) = UL(i,j,k,tau_i)**gamma
                endif
                if( einR(i,j,k) .lt. UR(i,j,k,egas_i) * des1 ) then
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
     &                 0.5d0*(gamma-1.0d0)*
     &      (einL + einR)*Xf(1:nx,1:ny,1:nz,dm)/t
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

c     If grid is moving, volume increases
          if( grd_isvelocity ) then
             dU = dU - 3.0d0 * U / t
          endif
          do i = xb, xe
          do j = yb, ye
          do k = zb, ze
            i0 = i - hydro_bw
            j0 = j - hydro_bw
            k0 = k - hydro_bw
            tmp8 = volinv(i,j,k) / product(dX(i,j,k,:)) / (t1-t0)
            dU(i,j,k,px_i) = dU(i,j,k,px_i)+grd_momdep(i0,j0,k0,1)*tmp8
            if( grd_igeom .ne. 11 ) then
              dU(i,j,k,py_i)=dU(i,j,k,py_i)+grd_momdep(i0,j0,k0,2)*tmp8
              dU(i,j,k,pz_i)=dU(i,j,k,pz_i)+grd_momdep(i0,j0,k0,3)*tmp8
            endif
          enddo
          enddo
          enddo


c     Compute timestep
          if( rk .eq. 1 ) then
            dtinv_max = max(
     &              maxval(dU(xb:xe,yb:ye,zb:ze,(/rho_i,tau_i/)) /
     &                      (-U(xb:xe,yb:ye,zb:ze,(/rho_i,tau_i/)))),
     &                  dtinv_max)
            dt = 0.1d0 / dtinv_max
            if( dt .ge. t1 - t .and. (.not. dt_only) ) then
              dt = t1 - t
              done = .true.
            endif
          endif

          if( .not. dt_only ) then

c     Apply dudt
            if( rk .eq. 1 ) then
              U = U0 + dU * dt
            else
              U = (U0 + U + dU*dt) * 0.5d0
            endif


            if( rk .eq. 1 ) then

              if(any(U(xb:xe,yb:ye,zb:ze,rho_i).le.0.0d0)) then
                write(*,*) 'rho < 0.0'
                call abort()
              endif
              if(any(U(xb:xe,yb:ye,zb:ze,tau_i).le.0.0d0)) then
                write(*,*) 'tau < 0.0'
              endif

              tfactor = 1.0d0 + dt / t

c     Upate X, Xf, and dX for moving grids
              do i = 1, 3
                if( veldim(i) ) then
                  X(:,:,:,i) = X(:,:,:,i) * tfactor
                  Xf(:,:,:,i) = Xf(:,:,:,i) * tfactor
                  dX(:,:,:,i) = dX(:,:,:,i) * tfactor
                endif
              enddo

c     Update geometrical quantities for moving grid
              if( grd_isvelocity ) then
                select case( grd_igeom )
                  case(1,11)
                    scle(:,:,:,2) = scle(:,:,:,2) * tfactor
                    scle(:,:,:,3) = scle(:,:,:,3) * tfactor
                    area(:,:,:,1) = area(:,:,:,1) * tfactor**2
                    area(:,:,:,2) = area(:,:,:,2) * tfactor
                    area(:,:,:,3) = area(:,:,:,3) * tfactor
                    volinv = volinv / (tfactor**2)
                  case(2)
                    scle(:,:,:,3) = scle(:,:,:,3) * tfactor
                    area(:,:,:,1) = area(:,:,:,1) * tfactor
                    area(:,:,:,2) = area(:,:,:,2) * tfactor
                    volinv = volinv / tfactor
                  case(3)
                end select
              endif
              t = t + dt

            endif

            call hydro_boundaries(U,X,veldim,nx,ny,nz,nf,bw,t)

c     Apply dual energy formalism to update tau
            kin = U(:,:,:,px_i)**2+U(:,:,:,py_i)**2+U(:,:,:,pz_i)**2
            kin = kin * 0.5d0 / U(:,:,:,rho_i)
            ein = U(:,:,:,egas_i) - kin
            do i = 2, nx-1
            do j = 2, ny-1
            do k = 2, nz-1
              tmp8 = U(i,j,k,egas_i)
              tmp8 = max(U(i-1,j,k,egas_i), tmp8)
              tmp8 = max(U(i+1,j,k,egas_i), tmp8)
              tmp8 = max(U(i,j-1,k,egas_i), tmp8)
              tmp8 = max(U(i,j+1,k,egas_i), tmp8)
              tmp8 = max(U(i,j,k-1,egas_i), tmp8)
              tmp8 = max(U(i,j,k+1,egas_i), tmp8)
              if( ein(i,j,k) .ge. tmp8 * des2 ) then
                U(i,j,k,tau_i) = ein(i,j,k)**(1.0d0/gamma)
              endif
            enddo
            enddo
            enddo

          endif

        enddo


        write(*,*) t, t1, dt, dtinv_max

        done = done .or. dt_only

      enddo

c      if( .not. dt_only ) then



        hydro_state = U


        call hydro_boundaries(hydro_state,X,veldim,nx,ny,nz,nf,bw,t)


        call scatter_hydro


        call hydro_velinterp

c     endif

      grd_momdep=0d0

      end subroutine


