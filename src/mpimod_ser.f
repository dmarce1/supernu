*This file is part of SuperNu.  SuperNu is released under the terms of the GNU GPLv3, see COPYING.
*Copyright (c) 2013-2017 Ryan T. Wollaeger and Daniel R. van Rossum.  All rights reserved.
      module mpimod
c     -------------
      implicit none
      integer :: MPI_COMM_WORLD=0
      integer,parameter :: MPI_MAX_PROCESSOR_NAME=13
      integer,private :: ierr=0
      integer :: impi=0  !mpi rank
      integer,parameter :: impi0=0 !master mpi rank
      logical :: lmpi0   !true on the master mpi rank
      integer :: nmpi=1  !number of mpi tasks
c
      save
c
      contains
c
      subroutine bcast_permanent
      end subroutine bcast_permanent
c
c
      subroutine scatter_inputstruct(ndim,icell1,ncell)
      use inputstrmod!{{{
      use gasmod
      implicit none
      integer,intent(in) :: ndim(3)
      integer,intent(out) :: icell1,ncell
************************************************************************
* mpi_scatter the input structure to all ranks in the worker comm.
************************************************************************
      integer :: dmy
c
      icell1 = 1
      ncell = str_nc
c
      dmy = ndim(1) !use the intent(in) variable
      allocate(str_massdd(ncell))
      str_massdd = reshape(str_massdc,[ncell]) !}}}
      if(str_nabund>0) then
       allocate(str_massfrdd(str_nabund,ncell))
       str_massfrdd = reshape(str_massfrdc,[str_nabund,ncell])
      endif
      if(str_ltemp) then
       allocate(str_tempdd(ncell))
       str_tempdd = reshape(str_tempdc,[ncell])
      endif
      if(str_lye) then
       allocate(str_yedd(ncell))
       str_yedd = reshape(str_yedc,[ncell])
      endif
      end subroutine scatter_inputstruct
c
c
      subroutine allgather_initialrad
c     -------------------------------
      use gridmod
      use gasmod
      grd_evolinit = reshape(gas_eraddens,[grd_ncell])
      end subroutine allgather_initialrad
c
c
      subroutine allgather_gammacap
c     ---------------------------
      use gridmod
      use gasmod
      grd_capgam = reshape(gas_capgam,[grd_ncell])
      grd_emitex = reshape(gas_emitex,[grd_ncell])
      end subroutine allgather_gammacap
c
c
      subroutine allreduce_gammaenergy
      end subroutine allreduce_gammaenergy
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      HYDRO LSU
      subroutine gather_hydro
      use gridmod
      use hydromod
      use gasmod
      use physconstmod
      implicit none

      integer :: i, j, k, l, f
      integer :: i0, j0, k0, f0
      real*8 :: eint

      do i = 1, grd_nx
      do j = 1, grd_ny
      do k = 1, grd_nz
        l = grd_icell(i,j,k)
        grd_vx(l) = grd_vx(l)+grd_momdep(i,j,k,1)/gas_mass(l)
        if( grd_igeom .ne. 11 ) then
        grd_vy(l) = grd_vy(l)+grd_momdep(i,j,k,2)/gas_mass(l)
        grd_vz(l) = grd_vz(l)+grd_momdep(i,j,k,3)/gas_mass(l)
        endif
      enddo
      enddo
      enddo
      grd_momdep=0d0

      do i = hydro_bw+1, hydro_nx - hydro_bw
      do j = hydro_bw+1, hydro_ny - hydro_bw
      do k = hydro_bw+1, hydro_nz - hydro_bw
         i0 = i - hydro_bw
         j0 = j - hydro_bw
         k0 = k - hydro_bw
         l = grd_icell(i0, j0, k0)
         hydro_state(i,j,k,rho_i) = gas_mass(l) / gas_vol(l)
         hydro_state(i,j,k,px_i) = grd_vx(l) * gas_rho(l)
         hydro_state(i,j,k,py_i) = grd_vy(l) * gas_rho(l)
         hydro_state(i,j,k,pz_i) = grd_vz(l) * gas_rho(l)
         eint =  1.5d0*pc_kb*(1.0d0+gas_nelec(l))
     &              * gas_natom(l) / gas_vol(l) * gas_temp(l)
         hydro_state(i,j,k,tau_i) = eint**(1.0d0 / hydro_gamma)
         hydro_state(i,j,k,natom_i) = gas_natom(l) / gas_vol(l)
         if( grd_igeom .eq. 11 ) then
           hydro_state(i,j,k,egas_i) = eint +
     &                               0.5d0*grd_vx(l)**2 * gas_rho(l)
         else
           hydro_state(i,j,k,egas_i) = eint + 0.5d0*
     &     (grd_vx(l)**2 + grd_vy(l)**2 + grd_vz(l)**2) * gas_rho(l)
         endif
         f = frac_i
         do f0 = 1, gas_nelem
           hydro_state(i,j,k,f) = gas_natom1fr(f0,l) * gas_natom(l)
     &                                               / gas_vol(l)
           f = f + 1
         enddo
         do f0 = -2*gas_nchain, -1
           hydro_state(i,j,k,f) = gas_natom1fr(f0,l) * gas_natom(l) /
     &                                                 gas_vol(l)
           f = f + 1
         enddo
c         write(*,*) i,j,k,gas_nelem,gas_natom1fr(1:gas_nelem,l)
         hydro_state(i,j,k,nelec_i) = gas_nelec(l) * gas_natom(l)
     &                                             / gas_vol(l)
      enddo
      enddo
      enddo
      end subroutine gather_hydro


      subroutine scatter_hydro
      use gridmod
      use hydromod
      use gasmod
      use physconstmod
      use timestepmod
      use elemdatamod
      implicit none

      integer :: i, j, k, l, f
      integer :: i0, j0, k0, f0
      real*8 :: natom
      real*8 :: eint, nnuc
      do i = hydro_bw+1, hydro_nx - hydro_bw
      do j = hydro_bw+1, hydro_ny - hydro_bw
      do k = hydro_bw+1, hydro_nz - hydro_bw
         i0 = i - hydro_bw
         j0 = j - hydro_bw
         k0 = k - hydro_bw
         l = grd_icell(i0, j0, k0)
         if( grd_isvelocity) then
           gas_vol(l) = gas_vol(l) * (1.0d0 + tsp_dt / tsp_t )**3
         endif
         gas_rho(l) = hydro_state(i,j,k,rho_i)
         gas_mass(l) = gas_rho(l) * gas_vol(l)
         grd_vx(l) = hydro_state(i,j,k,px_i) / gas_rho(l)
         grd_vy(l) = hydro_state(i,j,k,py_i) / gas_rho(l)
         grd_vz(l) = hydro_state(i,j,k,pz_i) / gas_rho(l)
         eint = hydro_state(i,j,k,egas_i) -
     &          (grd_vx(l)**2+grd_vy(l)**2+grd_vz(l)**2)*0.50d0*
     &              hydro_state(i,j,k,rho_i)
         if( eint .le. hydro_state(i,j,k,egas_i) * 0.1d0 ) then
           eint = hydro_state(i,j,k,tau_i)**(hydro_gamma)
         endif


         f = frac_i
         do f0 = 1, gas_nelem
           gas_natom1fr(f0,l) = hydro_state(i,j,k,f) * gas_vol(l)
           f = f + 1
         enddo
         do f0 = -2*gas_nchain, -1
           gas_natom1fr(f0,l) = hydro_state(i,j,k,f) * gas_vol(l)
           f = f + 1
         enddo
         gas_natom(l) = 0.0d0
         nnuc = 0.0d0
         do f = 1, gas_nelem
           natom = gas_natom1fr(f,l)
           nnuc = nnuc + natom * elem_data(f)%m
           gas_natom(l) = gas_natom(l) + natom
         enddo
         gas_natom(l) = hydro_state(i,j,k,natom_i) * gas_vol(l)
         gas_nelec(l) = hydro_state(i,j,k,nelec_i) * gas_vol(l)
     &                                             / gas_natom(l)
         gas_ye(l) = gas_nelec(l) * gas_natom(l) / nnuc
c         gas_nelec(l) = gas_nelec(l) / gas_natom(l)
         gas_bcoef(l) = 1.5d0*pc_kb*(1d0+gas_nelec(l))
     &              * gas_natom(l) / gas_vol(l)
         gas_temp(l) =    eint / gas_bcoef(l)
         gas_temp(l) = max(gas_temp(l),3d3)
c         write(*,*) l, gas_temp(l), eint, gas_bcoef(l)
         if( grd_isvelocity) then
           gas_vol(l) = gas_vol(l) / (1.0d0 + tsp_dt / tsp_t )**3
         endif
      enddo
      enddo
      enddo


      do i = 1, gas_ncell


        gas_natom1fr(28,i) = gas_natom1fr(28,i) -
     &   gas_natom1fr(gas_ini56,i)
        gas_natom1fr(27,i) = gas_natom1fr(27,i) -
     &   gas_natom1fr(gas_ico56,i)
        gas_natom1fr(26,i) = gas_natom1fr(26,i) -
     &   gas_natom1fr(gas_ife52,i)
        gas_natom1fr(25,i) = gas_natom1fr(25,i) -
     &   gas_natom1fr(gas_imn52,i)
        gas_natom1fr(24,i) = gas_natom1fr(24,i) -
     &   gas_natom1fr(gas_icr48,i)
        gas_natom1fr(23,i) = gas_natom1fr(23,i) -
     &   gas_natom1fr(gas_iv48,i)

       do l=1,gas_nelem
        gas_ye(i) = gas_ye(i) + gas_natom1fr(l,i)*l/elem_data(l)%m
       enddo
       if(gas_nchain/=3) stop 'massfr2natomfr: gas_nchain updated'
       gas_ye(i) = gas_ye(i) + gas_natom1fr(gas_ini56,i)*(28/56d0)
       gas_ye(i) = gas_ye(i) + gas_natom1fr(gas_ico56,i)*(27/56d0)
       gas_ye(i) = gas_ye(i) + gas_natom1fr(gas_ife52,i)*(26/52d0)
       gas_ye(i) = gas_ye(i) + gas_natom1fr(gas_imn52,i)*(25/52d0)
       gas_ye(i) = gas_ye(i) + gas_natom1fr(gas_icr48,i)*(24/48d0)
       gas_ye(i) = gas_ye(i) + gas_natom1fr(gas_iv48,i)*(23/48d0)
       gas_ye(i) = gas_ye(i)/sum(gas_natom1fr(:,i))

        gas_natom0fr(-2,i,1) = gas_natom1fr(gas_ini56,i)!unstable
        gas_natom0fr(-1,i,1) = gas_natom1fr(gas_ico56,i)!unstable
        gas_natom0fr(0:2,i,1) = gas_natom1fr(26:28,i)!stable
c-- fe/mn/cr
        gas_natom0fr(-2,i,2) = gas_natom1fr(gas_ife52,i)!unstable
        gas_natom0fr(-1,i,2) = gas_natom1fr(gas_imn52,i)!unstable
        gas_natom0fr(0:2,i,2) = gas_natom1fr(24:26,i)!stable
c-- cr/v/ti
        gas_natom0fr(-2,i,3) = gas_natom1fr(gas_icr48,i)!unstable
        gas_natom0fr(-1,i,3) = gas_natom1fr(gas_iv48,i)!unstable
        gas_natom0fr(0:2,i,3) = gas_natom1fr(22:24,i)!stable
c
        gas_natom1fr(28,i) = gas_natom1fr(28,i) +
     &   gas_natom1fr(gas_ini56,i)
        gas_natom1fr(27,i) = gas_natom1fr(27,i) +
     &   gas_natom1fr(gas_ico56,i)
        gas_natom1fr(26,i) = gas_natom1fr(26,i) +
     &   gas_natom1fr(gas_ife52,i)
        gas_natom1fr(25,i) = gas_natom1fr(25,i) +
     &   gas_natom1fr(gas_imn52,i)
        gas_natom1fr(24,i) = gas_natom1fr(24,i) +
     &   gas_natom1fr(gas_icr48,i)
        gas_natom1fr(23,i) = gas_natom1fr(23,i) +
     &   gas_natom1fr(gas_iv48,i)
c
        gas_natom(i) = sum(gas_natom1fr(1:,i))
        gas_natom1fr(:,i) = gas_natom1fr(:,i)/gas_natom(i)
        gas_natom0fr(:,i,:) = gas_natom0fr(:,i,:)/gas_natom(i)
      enddo


c      call abort()

      end subroutine scatter_hydro
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine bcast_nonpermanent
      use gridmod
      use gasmod
      use groupmod
      use sourcemod
      use particlemod
************************************************************************
* Broadcast the data that changes with time.
* - stub
************************************************************************
c-- domain decomposition
      grd_tempinv = reshape(1d0/gas_temp,[grd_ncell])
      grd_emit = reshape(gas_emit,[grd_ncell])


      grd_emitex = reshape(gas_emitex,[grd_ncell])
      grd_evolinit = reshape(gas_evolinit,[grd_ncell])
c
      grd_cap = reshape(gas_cap,[grp_ng,grd_ncell])
      grd_sig = reshape(gas_sig,[grd_ncell])
      grd_capgrey = reshape(gas_capgrey,[grd_ncell])
      grd_fcoef = reshape(gas_fcoef,[grd_ncell])
c
      src_nvacantall(1) = count(prt_isvacant)
      end subroutine bcast_nonpermanent
c
c
      subroutine allgather_leakage
      end subroutine allgather_leakage
c
c
      subroutine reduce_gridtally
************************************************************************
* Reduce the results from the packet transport that are needed for the
* temperature correction.
* - stub
************************************************************************
      use gridmod
      use gasmod
      gas_edep = reshape(grd_tally(1,:),[grd_ncell])
      gas_eraddens = reshape(grd_tally(2,:),[grd_ncell])
      end subroutine reduce_gridtally
c
c
      subroutine reduce_fluxtally
      end subroutine reduce_fluxtally
c
c
      subroutine reduce_fluxes
      end subroutine reduce_fluxes
c
      subroutine reduce_gastemp
c     -------------------------
      use gridmod
      use gasmod
      grd_tempinv = reshape(1d0/gas_temp,[grd_ncell])
      end subroutine reduce_gastemp
c
      subroutine scatter_restart_data
      end subroutine scatter_restart_data
c
c
      subroutine collect_restart_data
      end subroutine collect_restart_data
c
c
      subroutine mpimod_dealloc
      end subroutine mpimod_dealloc
c
c-- MPI intrinsics
c-----------------
      subroutine mpi_init(ierr_)
      implicit none
      integer :: ierr_
      ierr_ = ierr
      end subroutine mpi_init
c
      subroutine mpi_comm_rank(mpi_comm,impi_,ierr_)
      implicit none
      integer :: mpi_comm
      integer :: impi_,ierr_
      ierr_ = ierr
      impi_ = impi
      mpi_comm = MPI_COMM_WORLD
      end subroutine mpi_comm_rank
c
      subroutine mpi_comm_size(mpi_comm,nmpi_,ierr_)
      implicit none
      integer :: mpi_comm
      integer :: nmpi_,ierr_
      ierr_ = ierr
      nmpi_ = nmpi
      mpi_comm = MPI_COMM_WORLD
      end subroutine mpi_comm_size
c
      subroutine mpi_get_processor_name(pname,ilen_,ierr_)
      implicit none
      character*(MPI_MAX_PROCESSOR_NAME) :: pname
      integer :: ilen_,ierr_
      pname = 'NOT AVAILABLE'
      ierr_ = ierr
      ilen_ = 1
      end subroutine mpi_get_processor_name
c
      subroutine mpi_barrier(mpi_comm,ierr_)
      implicit none
      integer :: mpi_comm,ierr_
      ierr_ = ierr
      mpi_comm = MPI_COMM_WORLD
      end subroutine mpi_barrier
c
      subroutine mpi_finalize(ierr_)
      implicit none
      integer :: ierr_
      ierr_ = ierr
      end subroutine mpi_finalize
c
      end module mpimod
c vim: fdm=marker
