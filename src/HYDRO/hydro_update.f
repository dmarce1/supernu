
      subroutine hydro_eos_from_temp
     &                           (ndens,nelec,temp,frac,energy,pressure)
      use ionsmod
      use physconstmod
      use gasmod
      use hydromod
      implicit none
      real(8), intent(inout) :: nelec
      real(8), intent(in) :: ndens, temp, frac(gas_nelem)
      real(8), intent(out) :: pressure, energy
      real(8) ::  pot
      integer :: iconv

        nelec = max( 0.00, nelec / ndens )
        call
     &   ions_solve_eos_energy(frac,temp,ndens,nelec,iconv,pot)
        nelec = nelec * ndens
        pressure = pc_kb * (ndens + nelec) * temp
        energy = pot + pressure / (hydro_gamma - 1d0)

      end subroutine

      subroutine hydro_eos_to_temp
     &                           (ndens,nelec,temp,frac,ein,pressure)
      use ionsmod
      use physconstmod
      use gasmod
      use hydromod
      implicit none
      real(8), intent(in) :: ndens, ein, frac(gas_nelem)
      real(8), intent(inout) :: nelec, temp
      integer :: iconv
      real(8) ::   pressure, test1, test2
      real(8) ::  tmp_elec
      real(8) ::  energy

      real(8) :: f, temp1, temp2, dfdt, err, c0, min_temp, min_ene

      min_temp = 1
      min_ene = 1.5d0 * pc_kb * ndens * min_temp
      temp = max( min_temp, temp)

      energy = ein
      if( energy .le. min_ene ) then
        energy = min_ene
        temp = min_temp
        pressure = (hydro_gamma - 1d0) * min_ene
        return
      endif

c      toler = 1d-9
      iconv = 0

      c0 = 1d0
      err = 1d0
c           write(*,*) ndens, energy, min_ene
      do while( err .gt. 1d-6 )
c        write(*,*) iconv, f, dfdt, temp, err, energy, c0
        temp1 = temp
        temp2 = temp1 * 1.0001d0
        call hydro_eos_from_temp
     &                     (ndens,nelec,temp2,frac,test2,pressure)
        tmp_elec = nelec
        call hydro_eos_from_temp
     &                     (ndens,tmp_elec,temp1,frac,test1,pressure)
        f = test1 - energy
        dfdt = (test2 - test1) / (temp2 - temp1)
        temp = temp - f / dfdt * c0
        temp = max( temp, temp1 / 2.)
c        temp = min( temp, temp1 * 2. )
        err = abs( f / dfdt / temp )
        iconv = iconv + 1
        if( iconv .ge. 500 ) then
           write(*,*) iconv, err, test2, test1, temp2,temp1,energy,nelec
        endif
        if( iconv .gt. 1000 ) then
           stop 'hydro_eos_to_temp: failed to converge'
        endif
        c0 = c0 * (0.1**(1d-3))
      enddo


      end subroutine

      subroutine hydro_eos_pressure(nfrac, nelec, iene, pressure, temp)
      use ionsmod
      use physconstmod
      use gasmod
      use hydromod
      implicit none
      real(8), intent(in) :: iene, nfrac(gas_nelem)
      real(8), intent(out) :: pressure
      real(8), intent(inout) :: nelec, temp

      real(8) :: ndens, frac(gas_nelem)

      ndens = sum(nfrac)
      if( ndens .gt. 0d0 .and. iene .gt. 0d0 ) then
         frac = nfrac / ndens
         call hydro_eos_to_temp(ndens,nelec,temp,frac,iene,pressure)
      else
         pressure = 0d0
      endif

      end subroutine


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
      real(8) :: dtinv_max, qin, qout

      integer :: dm, rk
      integer :: i, j, k, f, i0, j0, k0



      real(8) :: U(hydro_nx,hydro_ny,hydro_nz,hydro_nf)
      real(8), save, allocatable :: V(:,:,:,:)
      real(8) :: U0(hydro_nx,hydro_ny,hydro_nz,hydro_nf)
      real(8) :: dU(hydro_nx,hydro_ny,hydro_nz,hydro_nf)
      real(8) :: UR(hydro_nx,hydro_ny,hydro_nz,hydro_nf)
      real(8) :: UL(hydro_nx,hydro_ny,hydro_nz,hydro_nf)
      real(8) :: VR(hydro_nx,hydro_ny,hydro_nz,hydro_nf)
      real(8) :: VL(hydro_nx,hydro_ny,hydro_nz,hydro_nf)
      real(8) :: PR(hydro_nx,hydro_ny,hydro_nz)
      real(8) :: PL(hydro_nx,hydro_ny,hydro_nz)
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
      integer :: xb, xe, yb, ye, zb, ze, xe2, ye2, ze2
      integer :: nx, ny, nz, bw, nf, l
      real(8) :: gamma, tfactor, tmp8, tmp2, tau_floor, rho_floor
      logical :: done
      logical, save :: first_call = .true.

      if( first_call ) then
        allocate(V(hydro_nx,hydro_ny,hydro_nz,hydro_nf))
        V = 1d4
        first_call = .false.
      endif

      tau_floor = 1d-20
      rho_floor = 1d-20

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
        case(11)
          scle(:,:,:,1) = 1.0d0
          scle(:,:,:,2) = abs(X(:,:,:,1))
          scle(:,:,:,3) = abs(X(:,:,:,1))
          area(:,:,:,1) = Xf(1:nx,1:ny,1:nz,1)**2
          area(:,:,:,2) = abs(Xf(1:nx,1:ny,1:nz,1))
          area(:,:,:,3) = abs(Xf(1:nx,1:ny,1:nz,1))
        case(1)
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
      UR = 0d0
      UL = 0d0

        write(*,*)
        write(*,*) t, t1




      qin = 0d0
      do i = bw+1, nx-bw
      do j = bw+1, ny-bw
      do k = bw+1, nz-bw
        qin = qin +
     &U(i,j,k,tau_i)**hydro_gamma * product(dX(i,j,k,:)) / volinv(i,j,k)
      enddo
      enddo
      enddo

      do while (.not. done)

        U0 = U

        do rk = 1,bw

          dU = 0.0d0
          dtinv_max = 0.0d0

          call hydro_boundaries(U,X,veldim,nx,ny,nz,nf,bw,t)


c         pre-recon
          do f = 1, nf
            if( f .ne. tau_i ) then
              if((f.ne.rho_i).and.(f.ne.natom_i))then
                if( f .eq. nelec_i .or. f .ge. frac_i ) then
                  l = natom_i
                else
                  l = rho_i
                endif
                V(:,:,:,f) = U(:,:,:,f) / U(:,:,:,l)
              else
                V(:,:,:,f) = U(:,:,:,f)
              endif
            endif
          enddo
          V(:,:,:,egas_i) = V(:,:,:,egas_i) - 0.5*V(:,:,:,px_i)**2
          if( grd_igeom .ne. 11 ) then
            V(:,:,:,egas_i) = V(:,:,:,egas_i) - 0.5*V(:,:,:,py_i)**2
            V(:,:,:,egas_i) = V(:,:,:,egas_i) - 0.5*V(:,:,:,pz_i)**2
          endif
          do i = 0, 2
            if( veldim(i+1) ) then
               V(:,:,:,px_i+i) = V(:,:,:,px_i+i) - X(:,:,:,i+1) / t
            endif
          end do


          do i = hydro_bw+1, hydro_nx-hydro_bw
          do j = hydro_bw+1, hydro_ny-hydro_bw
          do k = hydro_bw+1, hydro_nz-hydro_bw
c            write(*,*) i, j, k
c            call hydro_eos_to_temp
c     &        (V(i,j,k,natom_i),U(i,j,k,nelec_i),
c     &         V(i,j,k,tau_i),V(i,j,k,frac_i:gas_nelem+frac_i-1),
c     &         V(i,j,k,egas_i)*V(i,j,k,rho_i),tmp2)
          enddo
          enddo
          enddo

          call hydro_boundaries(V,X,veldim,nx,ny,nz,nf,bw,t)

c     Compute contribution to dudt in each flux direction
          do dm = 1, 3
            if( dimused(dm) ) then
              select case( dm )
                case(1)
                  xe2 = xe + 1
                  ye2 = ye
                  ze2 = ze
                case(2)
                  xe2 = xe
                  ye2 = ye + 1
                  ze2 = ze
                case(3)
                  xe2 = xe
                  ye2 = ye
                  ze2 = ze + 1
              end select

c     Reconstruct face values


c     Piecewise constant
              if( bw .eq. 1 ) then
                VR = V
                select case( dm )
                  case(1)
                    VL(2:nx,:,:,:) = V(1:nx-1,:,:,:)
                  case(2)
                    VL(:,2:ny,:,:) = V(:,1:ny-1,:,:)
                  case(3)
                    VL(:,:,2:nz,:) = V(:,:,1:nz-1,:)
                end select
c     Piecwise linear
              else if( bw .eq. 2 ) then
                do f = 1, nf
                 select case(dm)
                  case(1)
                    slp(2:nx,:,:) = (V(2:nx,:,:,f) - V(1:nx-1,:,:,f))/
     &                              (X(2:nx,:,:,1) - X(1:nx-1,:,:,1))
                    slp_m = slp
                    slp_p(1:nx-1,:,:) = slp(2:nx,:,:)
                  case(2)
                    slp(:,2:ny,:) = (V(:,2:ny,:,f) - V(:,1:ny-1,:,f))/
     &                              (X(:,2:ny,:,2) - X(:,1:ny-1,:,2))
                    slp_m = slp
                    slp_p(:,1:ny-1,:) = slp(:,2:ny,:)
                  case(3)
                    slp(:,:,2:nz) = (V(:,:,2:nz,f) - V(:,:,1:nz-1,f))/
     &                              (X(:,:,2:nz,3) - X(:,:,1:nz-1,3))
                    slp_m = slp
                    slp_p(:,:,1:nz-1) = slp(:,:,2:nz)
                 end select
                 slp = (sign(0.5d0,slp_m)+sign(0.5d0,slp_p))*
     &                   min(abs(slp_p),abs(slp_m))
                 select case(dm)
                  case(1)
                    VR(xb:xe+1,yb:ye,zb:ze,f) = V(xb:xe+1,yb:ye,zb:ze,f)
     &               - slp(xb:xe+1,yb:ye,zb:ze)
     &               *  dX(xb:xe+1,yb:ye,zb:ze,1)*0.5d0
                    VL(xb:xe+1,yb:ye,zb:ze,f) = V(xb-1:xe,yb:ye,zb:ze,f)
     &               + slp(xb-1:xe,yb:ye,zb:ze)
     &               *  dX(xb-1:xe,yb:ye,zb:ze,1)*0.5d0
                  case(2)
                    VR(xb:xe,yb:ye+1,zb:ze,f) = V(xb:xe,yb:ye+1,zb:ze,f)
     &               - slp(xb:xe,yb:ye+1,zb:ze)
     &               *  dX(xb:xe,yb:ye+1,zb:ze,2)*0.5d0
                    VL(xb:xe,yb:ye+1,zb:ze,f) = V(xb:xe,yb-1:ye,zb:ze,f)
     &               + slp(xb:xe,yb-1:ye,zb:ze)
     &               *  dX(xb:xe,yb-1:ye,zb:ze,2)*0.5d0
                  case(3)
                    VR(xb:xe,yb:ye,zb:ze+1,f) = V(xb:xe,yb:ye,zb:ze+1,f)
     &               - slp(xb:xe,yb:ye,zb:ze+1)
     &               *  dX(xb:xe,yb:ye,zb:ze+1,3)*0.5d0
                    VL(xb:xe,yb:ye,zb:ze+1,f) = V(xb:xe,yb:ye,zb-1:ze,f)
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
                   VR(:,:,:,px_i+i) = VR(:,:,:,px_i+i) +
     &                                          Xf(1:nx,1:ny,1:nz,i+1)/t
                   VL(:,:,:,px_i+i) = VL(:,:,:,px_i+i) +
     &                                          Xf(1:nx,1:ny,1:nz,i+1)/t
                endif
              enddo
              VR(:,:,:,egas_i) = VR(:,:,:,egas_i)+0.5*VR(:,:,:,px_i)**2
              VL(:,:,:,egas_i) = VL(:,:,:,egas_i)+0.5*VL(:,:,:,px_i)**2
              if( grd_igeom .ne. 11 ) then
               VR(:,:,:,egas_i) = VR(:,:,:,egas_i)+0.5*VR(:,:,:,py_i)**2
               VR(:,:,:,egas_i) = VR(:,:,:,egas_i)+0.5*VR(:,:,:,pz_i)**2
               VL(:,:,:,egas_i) = VL(:,:,:,egas_i)+0.5*VL(:,:,:,py_i)**2
               VL(:,:,:,egas_i) = VL(:,:,:,egas_i)+0.5*VL(:,:,:,pz_i)**2
              endif
              do f = 1, nf
                if((f.ne.rho_i).and.(f.ne.natom_i))then
                  if( f .eq. nelec_i .or. f .ge. frac_i ) then
                    l = natom_i
                  else
                    l = rho_i
                  endif
                  UR(:,:,:,f) = VR(:,:,:,f) * VR(:,:,:,l)
                  UL(:,:,:,f) = VL(:,:,:,f) * VL(:,:,:,l)
                else
                  UR(:,:,:,f) = VR(:,:,:,f)
                  UL(:,:,:,f) = VL(:,:,:,f)
                endif
              enddo


c     Compute face kinetic and internal energies and velocities
              if( grd_igeom .eq. 11 ) then
                kinL(xb:xe2,yb:ye2,zb:ze2)
     &               = 0.5d0 * UL(xb:xe2,yb:ye2,zb:ze2,px_i)**2 /
     &                         UL(xb:xe2,yb:ye2,zb:ze2,rho_i)
                kinR(xb:xe2,yb:ye2,zb:ze2)
     &               = 0.5d0 * UR(xb:xe2,yb:ye2,zb:ze2,px_i)**2 /
     &                         UR(xb:xe2,yb:ye2,zb:ze2,rho_i)
              else
                kinL = 0.5d0*(UL(:,:,:,px_i)**2 + UL(:,:,:,py_i)**2
     &                    + UL(:,:,:,pz_i)**2) / UL(:,:,:,rho_i)
                kinR = 0.5d0*(UR(:,:,:,px_i)**2 + UR(:,:,:,py_i)**2
     &                    + UR(:,:,:,pz_i)**2) / UR(:,:,:,rho_i)
              endif
              einL = max(UL(:,:,:,egas_i) - kinL,0d0)
              einR = max(UR(:,:,:,egas_i) - kinR,0d0)
              velL(xb:xe2,yb:ye2,zb:ze2) =
     &               UL(xb:xe2,yb:ye2,zb:ze2,px_i+dm-1) /
     &               UL(xb:xe2,yb:ye2,zb:ze2,rho_i)
              velR(xb:xe2,yb:ye2,zb:ze2) =
     &               UR(xb:xe2,yb:ye2,zb:ze2,px_i+dm-1) /
     &               UR(xb:xe2,yb:ye2,zb:ze2,rho_i)
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

c     Compute pressures
              do i = xb, xe2
              do j = yb, ye2
              do k = zb, ze2
c                write(*,*) 'R', i,j,k
                 PL(i,j,k) = (2./3.)*einL(i,j,k)
                 PR(i,j,k) = (2./3.)*einR(i,j,k)
c                call hydro_eos_from_temp
c     &            (VR(i,j,k,natom_i),UR(i,j,k,nelec_i),
c     &             VR(i,j,k,tau_i),VR(i,j,k,frac_i:gas_nelem+frac_i-1),
c     &             VR(i,j,k,egas_i),PR(i,j,k))
c                call hydro_eos_from_temp
c     &            (VL(i,j,k,natom_i),UL(i,j,k,nelec_i),
c     &             VL(i,j,k,tau_i),VL(i,j,k,frac_i:gas_nelem+frac_i-1),
c     &             VL(i,j,k,egas_i),PL(i,j,k))
              enddo
              enddo
              enddo


c     Compute signal speeds

              cL(xb:xe2,yb:ye2,zb:ze2) =
     & sqrt(gamma*PL(xb:xe2,yb:ye2,zb:ze2) /
     &        UL(xb:xe2,yb:ye2,zb:ze2,rho_i))

              cR(xb:xe2,yb:ye2,zb:ze2) =
     & sqrt(gamma*PR(xb:xe2,yb:ye2,zb:ze2) /
     &            UR(xb:xe2,yb:ye2,zb:ze2,rho_i))
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
     &         0.5d0*(PL * velL + PR * velR)

              if( veldim(dm) ) then
                Fv(:,:,:,egas_i) = Fv(:,:,:,egas_i) +
     &                 0.5d0*(gamma-1.0d0)*
     &      (einL + einR)*Xf(1:nx,1:ny,1:nz,dm)/t
              endif


c     Compute scalar fluxes
              Fs = 0.0d0
              Fs(:,:,:,px_i+dm-1) = 0.5d0*(PL + PR)

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
             dU(:,:,:,tau_i) = dU(:,:,:,tau_i)
     &        + (3.0-3.0/gamma)*U(:,:,:,tau_i) / t
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
c           dtinv_max = max(
c     &              maxval(dU(xb:xe,yb:ye,zb:ze,rho_i) /
c     &                      (-U(xb:xe,yb:ye,zb:ze,rho_i))),
c     &              maxval(dU(xb:xe,yb:ye,zb:ze,tau_i) /
c     &                      (-U(xb:xe,yb:ye,zb:ze,tau_i))),
c     &                  dtinv_max)
            dt = 2d0 / dtinv_max / 15d0
            if( dt .ge. t1 - t .and. (.not. dt_only) ) then
              dt = t1 - t
              done = .true.
            else
              dt = (t1 - t) / ceiling((t1 - t) / dt)
            endif
c            write(*,*) 'dt=', dt, 0.4d0 / dtinv_max
          endif

          if( .not. dt_only ) then

c     Apply dudt
            if( rk .eq. 1 ) then
              U = U0 + dU * dt
            else
              U = (U0 + U + dU*dt) * 0.5d0
            endif


            if(any(U(xb:xe,yb:ye,zb:ze,rho_i).le.0.0d0)) then
              write(*,*) 'rho < 0.0'
              call abort()
            endif
            if(any(U(xb:xe,yb:ye,zb:ze,rho_i)/=
     &                   U(xb:xe,yb:ye,zb:ze,rho_i))) then
              write(*,*) 'rho nan'
              call abort()
            endif
            if(any(U(xb:xe,yb:ye,zb:ze,natom_i).le.0.0d0)) then
              write(*,*) 'natom < 0.0'
              call abort()
            endif
            if(any(U(xb:xe,yb:ye,zb:ze,natom_i)/=
     &         U(xb:xe,yb:ye,zb:ze,natom_i))) then
              write(*,*) 'natom nan'
              call abort()
            endif
            if(any(U(xb:xe,yb:ye,zb:ze,frac_i:gas_nelem+frac_i-1).
     &                    lt.0.0d0)) then
              write(*,*) 'frac < 0.0'
              call abort()
             endif

c              if(any(U(xb:xe,yb:ye,zb:ze,tau_i).le.0.0d0)) then
c                write(*,*) 'tau < 0.0'
c              endif

             if( rk .eq. 1 ) then
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
              if( grd_igeom .ne. 11 ) then
                tmp8 = max(U(i,j-1,k,egas_i), tmp8)
                tmp8 = max(U(i,j+1,k,egas_i), tmp8)
                tmp8 = max(U(i,j,k-1,egas_i), tmp8)
                tmp8 = max(U(i,j,k+1,egas_i), tmp8)
              endif
              if(ein(i,j,k).ge.tmp8*des2)then
                U(i,j,k,tau_i) = ein(i,j,k)**(1.0d0/gamma)
              endif
              U(i,j,k,tau_i) = max(tau_floor,U(i,j,k,tau_i))
              U(i,j,k,rho_i) = max(rho_floor,U(i,j,k,rho_i))
            enddo
            enddo
            enddo

          endif

c        pause

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





      qout = 0d0
      do i = bw+1, nx-bw
      do j = bw+1, ny-bw
      do k = bw+1, nz-bw
        qout = qout +
     &U(i,j,k,tau_i)**hydro_gamma * product(dX(i,j,k,:)) / volinv(i,j,k)
      enddo
      enddo
      enddo

      open(unit=1,file='sums.dat',status='unknown', position='append')
      write(1,*) (qout - qin) / (t1 - t0)
      close(1)

      end subroutine


