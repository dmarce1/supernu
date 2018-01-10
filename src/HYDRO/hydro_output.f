


      subroutine hydro_output


      use hydromod
      use gridmod
      use timestepmod
      implicit none

      integer, save :: onum = 0
      character*15 :: filename
      integer :: i, j
      real :: x


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

      else

        write(*,*)
     &        'Hydro output in this coordinate system not implemented'

      endif


      onum = onum + 1


      end subroutine


