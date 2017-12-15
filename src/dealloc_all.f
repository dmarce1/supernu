*This file is part of SuperNu.  SuperNu is released under the terms of the GNU GPLv3, see COPYING.
*Copyright (c) 2013-2017 Ryan T. Wollaeger and Daniel R. van Rossum.  All rights reserved.
      subroutine dealloc_all
c     ----------------------
      use mpimod
      use ionsmod
      use bbxsmod
      use gridmod
      use groupmod
      use gasmod
      use particlemod
      use fluxmod
      use sourcemod
      use randommod
      use timestepmod
      implicit none
************************************************************************
* deallocate all that was used till the end of the program. Any
* allocatable arrays remaining allocated after this had to be dealt
* with earlier.  This helps to catch memory leaks! (drr)
************************************************************************
c-- ionsmod
      call ions_dealloc
      call gas_dealloc
      call grid_dealloc
      call flux_dealloc
      deallocate(prt_particles,prt_isvacant)
      call mpimod_dealloc
      deallocate(grp_wl,grp_wlinv)
      deallocate(src_nvacantall)
      deallocate(rnd_states)
      if(allocated(bb_xs)) deallocate(bb_xs) !only impi==impi0, but only if nobbopac==f
      deallocate(tsp_tarr)

      end subroutine dealloc_all
c vim: fdm=marker
