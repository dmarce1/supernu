      subroutine hydro_setup
      use hydromod
      use gridmod
      implicit none
      integer :: i, j, k, l, i0, j0, k0
      real*8 :: tmp
      hydro_state = 0.0d0
      end subroutine

      subroutine hydro_load_state
      use hydromod
      use gridmod
      implicit none

      integer :: i, j, k, l, i0, j0, k0

      do i0 = 1, grd_nx
      do j0 = 1, grd_ny
      do k0 = 1, grd_nz
         i = i0 + hydro_bw
         j = j0 + hydro_bw
         k = k0 + hydro_bw
         l = grd_icell(i0,j0,k0)
c         hydro_state(i,j,k,

      enddo
      enddo
      enddo


      end subroutine


      subroutine hydro_unload_state
      implicit none

      end subroutine
