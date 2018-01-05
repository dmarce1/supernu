      subroutine hydro_setup
      use hydromod
      use gridmod
      implicit none
      hydro_state = 0.0d0
      end subroutine


      subroutine sedov_setup
      use hydromod
      use inputstrmod
      use gridmod
      use mpimod
      implicit none
c     http://cococubed.asu.edu/research_pages/sedov.shtml


c..declare
      integer          i,nstep
      real*8          time,zpos(grd_nx-1),
     1                 eblast,rho0,omega,vel0,ener0,pres0,cs0,gamma,
     2                 xgeom,
     3                 den(grd_nx-1),ener(grd_nx-1),pres(grd_nx-1),
     4                 vel(grd_nx-1),cs(grd_nx-1),
     5                 zlo,zhi,zstep,value
      integer :: j, k, l
      real*8 :: h_time

      write(*,*) '~~~~~~~~~~~~~~~~~~~~', grd_nx, grd_ny, grd_nz
      allocate(str_xleft(grd_nx))
      allocate(str_yleft(grd_ny))
      allocate(str_zleft(grd_nz))


      nstep = 120
      eblast = 0.851072d0
      xgeom  = 3.0d0
      omega  = 0.0d0


c..input parameters in cgs
      time   = 1.0d0
      rho0   = 1.0d0
      vel0   = 0.0d0
      ener0  = 0.0d0
      pres0  = 0.0d0
      cs0    = 0.0d0
      gamma  = 5.0d0/3.0d0




c..number of grid points, spatial domain, spatial step size.
c..to match hydrocode output, use the mid-cell points.
      zlo   = 0.0d0
      zhi   = 1.2d0
      zstep = (zhi - zlo)/float(nstep)
      do i=1,nstep
       zpos(i)   = zlo + 0.5d0*zstep + float(i-1)*zstep
      enddo


c..get the solution for all spatial points at once

       call sed_1d(time,nstep,zpos,
     1             eblast,omega,xgeom,
     2             rho0,vel0,ener0,pres0,cs0,gamma,
     3             den,ener,pres,vel,cs)



      h_time = 1.0d0 / maxval(vel/zpos)
      str_xleft(1) = 0.0d0
      do i = 2, grd_nx
        str_xleft(i) = 0.5d0 * (zpos(i) + zpos(i-1)) / h_time
      enddo
      str_yleft = 0.0d0
      str_zleft = 1.0d0

      do i = hydro_bw+1, hydro_nx - hydro_bw
      do j = hydro_bw+1, hydro_ny - hydro_bw
      do k = hydro_bw+1, hydro_nz - hydro_bw
        l = i - hydro_bw
        hydro_state(i,j,k,rho_i) = den(l)
        hydro_state(i,j,k,px_i) = den(l) * vel(l)
        hydro_state(i,j,k,py_i) = 0.0d0
        hydro_state(i,j,k,pz_i) = 0.0d0
        hydro_state(i,j,k,egas_i) = (ener(l) + 0.5d0*vel(l)**2) * den(l)
        hydro_state(i,j,k,tau_i) = (ener(l)*den(l))**(1.0d0/hydro_gamma)
        hydro_state(i,j,k,frac_i) = den(l)
        hydro_state(i,j,k,frac_i+1:hydro_nf) = 0.0d0
      enddo
      enddo
      enddo

      call scatter_hydro

      end subroutine
