      subroutine hydro_setup
      use hydromod
      use gridmod
      implicit none
      hydro_state = 0.0d0
      end subroutine


      subroutine sedov_setup
      use hydromod
      use inputstrmod
      use inputparmod
      use timestepmod
      use gridmod
      use mpimod
      use physconstmod
      implicit none
c     http://cococubed.asu.edu/research_pages/sedov.shtml


c..declare
      integer          i,nstep
      real*16          time,zpos(in_ndim(1)+1),
     1                 eblast,rho0,omega,vel0,ener0,pres0,cs0,gamma,
     2                 xgeom,
     3                 den(in_ndim(1)+1),ener(in_ndim(1)+1),
     4                 pres(in_ndim(1)+1),vel(in_ndim(1)+1),
     5                 cs(in_ndim(1)+1),zlo,zhi,zstep,value
      integer :: j, k
      real*16 :: h_time, vol, mint



      nstep = in_ndim(1)+1
      eblast = 1.0d0
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
      do i = 2, in_ndim(1)+1
        str_xleft(i) = 0.5d0 * (zpos(i) + zpos(i-1)) / in_tsp_tfirst
      enddo
      str_yleft(1) = -1.0d0
      str_yleft(2) = +1.0d0
      str_zleft(1) = 0.0d0
      str_zleft(2) = 2.0d0 * pc_pi


      mint = 1.0d+100
      do i = 1, in_ndim(1)
      do j = 1, in_ndim(2)
      do k = 1, in_ndim(3)
        vol = 4.0*3.14159/3.0*(str_xleft(i+1)**3 - str_xleft(i)**3)
        vol = vol * in_tsp_tfirst**3
        str_mass(i,j,k) = den(i)*vol
        str_vx(i,j,k) = vel(i)
        str_vy(i,j,k) = 0.0d0
        str_vz(i,j,k) = 0.0d0
        write(*,*) i,j,k,ener(i)
        str_temp(i,j,k) = ener(i) / (1.5d0*pc_kb/pc_mh)
        str_massfr(:,i,j,k) = 0.0d0
        str_massfr(1,i,j,k) = 1.0d0
        str_ye(i,j,k) = 1.0d0
        if( ener(i) .gt. 0.0d0 ) then
          mint = min(mint,ener(i))
        endif
      enddo
      enddo
      enddo

      mint = mint / (1.5d0*pc_kb/pc_mh)
      str_temp = max(str_temp, mint/10.0d0)
      str_mass = max(str_mass,1.0e-10)

      call scatter_hydro

      end subroutine
