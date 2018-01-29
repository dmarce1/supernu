


      subroutine write_sedov_analytic(onum)
      use physconstmod
      use gridmod
      use hydromod
      use timestepmod
      implicit none

      integer          i,nstep
      real*16          time,zpos(grd_nx),
     1                 eblast,rho0,omega,vel0,ener0,pres0,cs0,gamma,
     2                 xgeom,
     3                 den(grd_nx),ener(grd_nx),
     4                 pres(grd_nx),vel(grd_nx),
     5                 cs(grd_nx)
      integer, intent(in) :: onum
      character*17 :: filename

      write(filename,'("hydro.a.",I0.4,".dat")') onum
      open(unit=2,file=filename,status='unknown')


      nstep = grd_nx
      eblast = 1.0d0
      xgeom  = 3.0d0
      omega  = 0.0d0
      zpos = (grd_xarr(2:nstep+1)+grd_xarr(1:nstep))/2.0d0 * tsp_t


c..input parameters in cgs
      time   = 1.0d0 + tsp_t - tsp_tfirst
      rho0   = 1.0d0
      vel0   = 0.0d0
      ener0  = 0.0d0
      pres0  = 0.0d0
      cs0    = 0.0d0
      gamma  = hydro_gamma




c..get the solution for all spatial points at once

       call sed_1d(time,nstep,zpos,
     1             eblast,omega,xgeom,
     2             rho0,vel0,ener0,pres0,cs0,gamma,
     3             den,ener,pres,vel,cs)


      do i = 1, grd_nx


        write(2,'(E14.6)', advance="no") zpos(i)
        write(2,'(E14.6)', advance="no") den(i)
        write(2,'(E14.6)', advance="no") den(i) * vel(i)
        write(2,'(E14.6)', advance="no") 0.0d0
        write(2,'(E14.6)', advance="no") 0.0d0
        write(2,'(E14.6)', advance="no")
     &            (ener(i)*den(i))**(1.0/gamma)
        write(2,'(E14.6)')
     &            ener(i)*den(i)+0.5d0*den(i)*vel(i)*vel(i)

      enddo

      close(2)

      end subroutine


      subroutine hydro_output


      use hydromod
      use gridmod
      use inputparmod
      use timestepmod
      implicit none

      integer, save :: onum = 0
      character*15 :: filename
      integer :: i, j
      real*8 :: x


      if( grd_igeom .eq. 11 ) then

        write(filename,'("hydro.",I0.4,".dat")') onum
        open(unit=1,file=filename,status='unknown')
        do i = hydro_bw+1, hydro_nx - hydro_bw
          x = (grd_xarr(i-hydro_bw) + grd_xarr(i+1-hydro_bw)) * 0.5d0
          if( grd_isvelocity ) then
            x = x * tsp_t
          endif
          write(1,'(E14.6)', advance="no") x
          do j = 1, hydro_nf
            write(1,'(E14.6)', advance="no")
     &                   hydro_state(i,hydro_bw+1,hydro_bw+1,j)
          enddo
          write(1,*)
        enddo
        close(1)

        select case( in_test_problem )
          case(1)
            call write_sedov_analytic(onum)
        end select

      else

        write(*,*)
     &        'Hydro output in this coordinate system not implemented'

      endif


      onum = onum + 1


      end subroutine


