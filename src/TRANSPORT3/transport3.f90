!This file is part of SuperNu.  SuperNu is released under the terms of the GNU GPLv3, see COPYING.
!Copyright (c) 2013-2017 Ryan T. Wollaeger and Daniel R. van Rossum.  All rights reserved.
pure subroutine transport3(ptcl,ptcl2,rndstate,edep,eraddens,eamp,totevelo,ierr)

  use randommod
  use miscmod
  use gridmod
  use groupmod
  use timestepmod
  use physconstmod
  use particlemod
  use transportmod
  use inputparmod
  implicit none
!
  type(packet),target,intent(inout) :: ptcl
  type(packet2),target,intent(inout) :: ptcl2
  type(rnd_t),intent(inout) :: rndstate
  real*8,intent(out) :: edep, eraddens, eamp
  real*8,intent(inout) :: totevelo
  integer,intent(out) :: ierr
!##################################################
  !This subroutine passes particle parameters as input and modifies
  !them through one IMC transport event.  If
  !the puretran boolean is set to false, this routine couples to the
  !corresponding DDMC diffusion routine.
!##################################################
  real*8,parameter :: cinv = 1d0/pc_c

  logical :: loutx,louty,loutz
  integer :: iznext,iynext,ixnext
  real*8 :: elabfact, eta, xi
  real*8 :: thelp, thelpinv, help
  real*8 :: alb, eps, beta, pp
  real*8 :: dcen,dcol,dthm,dbx,dby,dbz,ddop
  real*8 :: darr(7)
  real*8 :: r1, r2
!-- distance out of physical reach
  real*8 :: far
  real*8 :: emitlump

  integer,pointer :: ix, iy, iz, ic, ig
  real*8,pointer :: x,y,z,mu,om,e,e0,wl,d
!-- statement functions
  integer :: l,j,k
  real*8 :: dx,dy,dz,xm,ym,zm,rsq
  dx(l) = grd_xarr(l+1) - grd_xarr(l)
  dy(l) = grd_yarr(l+1) - grd_yarr(l)
  dz(l) = grd_zarr(l+1) - grd_zarr(l)
  xm(l) = .5*(grd_xarr(l+1) + grd_xarr(l))
  ym(l) = .5*(grd_yarr(l+1) + grd_yarr(l))
  zm(l) = .5*(grd_zarr(l+1) + grd_zarr(l))
  rsq(l,j,k) = xm(l)**2 + ym(j)**2 + zm(k)**2

  ix => ptcl2%ix
  iy => ptcl2%iy
  iz => ptcl2%iz
  ic => ptcl2%ic
  ig => ptcl2%ig
  d => ptcl2%dist
  x => ptcl%x
  y => ptcl%y
  z => ptcl%z
  mu => ptcl%mu
  om => ptcl%om
  e => ptcl%e
  e0 => ptcl%e0
  wl => ptcl%wl

!-- no error by default
  ierr = 0
!-- init
  edep = 0d0
  eraddens = 0d0
  eamp = 0d0
!
!-- projections
  eta = sqrt(1d0-mu**2)*sin(om)
  xi = sqrt(1d0-mu**2)*cos(om)
!
!-- setting vel-grid helper variables
  if(grd_isvelocity) then
!-- calculating initial transformation factors
     elabfact=1d0-(mu*z+eta*y+xi*x)*cinv
     thelp = tsp_t
  else
     elabfact = 1d0
     thelp = 1d0
  endif
!
!-- inverting vel-grid factor
  thelpinv = 1d0/thelp

!-- distance longer than distance to census
  far = 2d0*abs(pc_c*tsp_dt*thelpinv) !> dcen

!-- census distance
  dcen = abs(pc_c*(tsp_t1-ptcl%t)*thelpinv)
!
!-- boundary distances
  if(xi==0d0) then
     dbx = far
  else
     dbx = max((grd_xarr(ix)-x)/xi,(grd_xarr(ix+1)-x)/xi)
  endif
  if(eta==0d0) then
     dby = far
  else
     dby = max((grd_yarr(iy)-y)/eta,(grd_yarr(iy+1)-y)/eta)
  endif
  if(mu==0d0) then
     dbz = far
  else
     dbz = max((grd_zarr(iz)-z)/mu,(grd_zarr(iz+1)-z)/mu)
  endif
!
!-- Thomson scattering distance
  if(grd_sig(ic)>0d0) then
     call rnd_r(r1,rndstate)
     dthm = -log(r1)*thelpinv/(elabfact*grd_sig(ic))
  else
     dthm = far
  endif
!
!-- effective collision distance
  if(grd_cap(ig,ic)<=0d0) then
     dcol = far
  elseif(trn_isimcanlog) then
!-- calculating dcol for analog MC
     call rnd_r(r1,rndstate)
     dcol = -log(r1)*thelpinv/(elabfact*grd_cap(ig,ic))
  elseif(grd_fcoef(ic)<1d0.and.grd_fcoef(ic)>=0d0) then
     call rnd_r(r1,rndstate)
     dcol = -log(r1)*thelpinv/&
          (elabfact*(1d0-grd_fcoef(ic))*grd_cap(ig,ic))
  else
     dcol = far
  endif
!
!-- Doppler shift distance
  if(grd_isvelocity.and.ig<grp_ng) then
     ddop = pc_c*(elabfact-wl*grp_wlinv(ig+1))
     if(ddop<0d0) then
        ddop = far
     endif
  else
     ddop = far
  endif
!
!-- finding minimum distance
  darr = [dcen,dcol,dthm,ddop,dbx,dby,dbz]
  ptcl2%idist = minloc(darr,dim=1)
  d = darr(ptcl2%idist)
! if(any(darr/=darr) .or. d<0d0) then
!    write(0,*) darr
!    write(*,*) ix,iy,iz,x,y,z,mu,eta,xi,om
!    stop 'transport3: invalid distance'
! endif
  if(any(darr/=darr)) then
     ierr = 1
     return
  endif
  if(d<0d0) then
     ierr = 2
     return
  endif

!-- updating position
  if(d==dbx) then
     if(xi>=0d0) then
        ixnext = ix+1
        x = grd_xarr(ix+1)
     else
        ixnext = ix-1
        x = grd_xarr(ix)
     endif
  else
     x = x + xi*d
  endif

  if(d==dby) then
     if(eta>=0d0) then
        iynext = iy+1
        y = grd_yarr(iy+1)
     else
        iynext = iy-1
        y = grd_yarr(iy)
     endif
  else
     y = y + eta*d
  endif

  if(d==dbz) then
     if(mu>=0d0) then
        iznext = iz + 1
        z = grd_zarr(iz+1)
     else
        iznext = iz - 1
        z = grd_zarr(iz)
     endif
  else
     z = z + mu*d
  endif

!
!-- updating time
  ptcl%t = ptcl%t + thelp*cinv*d

!-- tallying energy densities
  if(trn_isimcanlog) then
!-- analog energy density
     eraddens = e*elabfact**2 * &
          d*thelp*cinv*tsp_dtinv
  else
!-- nonanalog energy density
     if(grd_fcoef(ic)*grd_cap(ig,ic)* &
          min(dx(ix),dy(iy),dz(iz))*thelp>1d-6) then
        eraddens = e* &
             (1.0d0-exp(-grd_fcoef(ic)*elabfact* &
             grd_cap(ig,ic)*d*thelp))* &
             elabfact/(grd_fcoef(ic)*elabfact * &
             grd_cap(ig,ic)*pc_c*tsp_dt)
     else
!-- analog energy density
        eraddens = e*elabfact**2 * &
             d*thelp*cinv*tsp_dtinv
     endif
!-- depositing nonanalog absorbed energy
     edep = e* &
          (1d0-exp(-grd_fcoef(ic)*grd_cap(ig,ic)* &
          elabfact*d*thelp))*elabfact
     if(edep/=edep) then
!       write(0,*) e,grd_fcoef(ic),grd_cap(ig,ic),elabfact,d,thelp
!       stop 'transport3: invalid energy deposition'
        ierr = 3
        return
     endif
!-- reducing particle energy
     e = e*exp(-grd_fcoef(ic)*grd_cap(ig,ic) * &
          elabfact*d*thelp)
  endif

!
!-- updating transformation factors
  if(grd_isvelocity) then
     elabfact=1d0-(xi*x+eta*y+mu*z)*cinv
  endif

!
!-- checking which event occurs from min distance

!
!-- census
  if(d==dcen) then
     ptcl2%stat = 'cens'
     return
  endif

!-- common manipulations for collisions
  if(d==dthm.or.d==dcol) then
!-- resampling direction
     call rnd_r(r1,rndstate)
     mu = 1d0 - 2d0*r1
     call rnd_r(r1,rndstate)
     om = pc_pi2*r1
!-- checking velocity dependence
     if(grd_isvelocity) then
        eta = sqrt(1d0-mu**2)*sin(om)
        xi = sqrt(1d0-mu**2)*cos(om)
!-- transforming mu
        mu = (mu+z*cinv)/(1d0+(mu*z+eta*y+xi*x)*cinv)
        if(mu>1d0) then
           mu = 1d0
        elseif(mu<-1d0) then
           mu = -1d0
        endif
!-- transforming om
        om = atan2(eta+y*cinv,xi+x*cinv)
        if(om<0d0) om=om+pc_pi2
!-- x,y lab direction cosines
        eta = sqrt(1d0-mu**2)*sin(om)
        xi = sqrt(1d0-mu**2)*cos(om)
     endif
  endif

!
!-- checking if escaped domain
  loutx = .false.
  if(d==dbx) then
     if(ixnext==grd_nx+1 .or. ixnext==0) then
        loutx = .true. !domain edge is reached
     elseif(grd_icell(ixnext,iy,iz)==grd_ivoid) then
        if(ixnext>ix.eqv.grd_xarr(ixnext)>0d0) then !away from center
           loutx = rsq(ixnext,iy,iz)>grd_rvoid**2 !enter void corners
        endif
     endif
  endif
  louty = .false.
  if(d==dby) then
     if(iynext==grd_ny+1 .or. iynext==0) then
        louty = .true. !domain edge is reached
     elseif(grd_icell(ix,iynext,iz)==grd_ivoid) then
        if(iynext>iy.eqv.grd_yarr(iynext)>0d0) then !away from center
           louty = rsq(ix,iynext,iz)>grd_rvoid**2 !enter void corners
        endif
     endif
  endif
  loutz = .false.
  if(d==dbz) then
     if(iznext==grd_nz+1 .or. iznext==0) then
        loutz = .true. !domain edge is reached
     elseif(grd_icell(ix,iy,iznext)==grd_ivoid) then
        if(iznext>iz.eqv.grd_zarr(iznext)>0d0) then !away from center
           loutz = rsq(ix,iy,iznext)>grd_rvoid**2 !enter void corners
        endif
     endif
  endif
  if(loutx.or.louty.or.loutz) then
!-- observer time correction
     eta = sqrt(1d0-mu**2)*sin(om)
     xi = sqrt(1d0-mu**2)*cos(om)
     ptcl%t=ptcl%t-(mu*z+eta*y+xi*x)*thelp*cinv
!-- ending particle
     ptcl2%stat = 'flux'
     return
  endif

!
!-- Thomson scatter
  if(d==dthm) then
!-- checking velocity dependence
     if(grd_isvelocity) then
!-- lab wavelength
        wl = wl*(1d0-(mu*z+eta*y+xi*x)*cinv)/elabfact
        help = elabfact/(1d0-(mu*z+eta*y+xi*x)*cinv)
!-- velocity effects accounting
        totevelo=totevelo+e*(1d0-help)
!
!-- energy weight
        e = e*help
        e0 = e0*help
     endif

!
!-- effective collision
  elseif(d==dcol) then
     call rnd_r(r1,rndstate)
!-- checking if analog
     if(trn_isimcanlog.and.r1<=grd_fcoef(ic)) then
!-- effective absorption
        ptcl2%stat = 'dead'
!-- adding comoving energy to deposition energy
        edep = e*elabfact
!-- velocity effects accounting
        totevelo = totevelo+e*(1d0-elabfact)
        return
     else
!-- effective scattering
!
!-- transforming to lab
        if(grd_isvelocity) then
           help = elabfact/(1d0-(mu*z+eta*y+xi*x)*cinv)
!-- velocity effects accounting
           totevelo=totevelo+e*(1d0-help)
!-- energy weight
           e = e*help
           e0 = e0*help
        endif
!
!-- sample group
        emitlump = grd_opaclump(8,ic)/grd_capgrey(ic)
        if(grp_ng==1) then
        elseif(emitlump<.99d0 .or. trn_nolumpshortcut .or. in_puretran) then
           call rnd_r(r1,rndstate)
           ig = emitgroup(r1,ic)
!-- checking if DDMC in new group
           if((grd_cap(ig,ic)+grd_sig(ic)) * &
              min(dx(ix),dy(iy),dz(iz))*thelp >= trn_tauddmc &
              .and. .not.in_puretran) ptcl2%itype = 2
!-- don't sample, it will end up in the lump anyway
        else
           ptcl2%itype = 2
!-- always put this in the single most likely group
           ig = nint(grd_opaclump(9,ic))
        endif
!
!-- sample wavelength
        if(ptcl2%itype==2) then
!-- transforming to cmf
           if(grd_isvelocity) then
!-- velocity effects accounting
              totevelo=totevelo+e*(mu*z+eta*y+xi*x)*cinv
!-- energy weight
              e = e*(1d0-(mu*z+eta*y+xi*x)*cinv)
              e0 = e0*(1d0-(mu*z+eta*y+xi*x)*cinv)
           endif
           wl = 0d0 !workaround ifort 13.1.3 bug
        else
!-- uniformly in new group
           call rnd_r(r1,rndstate)
           wl = 1d0/((1d0-r1)*grp_wlinv(ig)+r1*grp_wlinv(ig+1))
!-- converting comoving wavelength to lab frame wavelength
           if(grd_isvelocity) wl = wl*(1d0-(mu*z+eta*y+xi*x)*cinv)
        endif
     endif

!
!-- x-bound
  elseif(d==dbx) then

     l = grd_icell(ixnext,iy,iz)
     if((grd_cap(ig,l)+grd_sig(l)) * &
          min(dx(ixnext),dy(iy),dz(iz))*thelp < trn_tauddmc &
          .or.in_puretran) then
!-- IMC in adjacent cell
        ix = ixnext
        ic = grd_icell(ix,iy,iz)
     else
!-- DDMC in adjacent cell
        if(grd_isvelocity) then
!-- transforming x-cosine to cmf
           xi = (xi-x*cinv)/elabfact
           if(xi>1d0) then
              xi = 1
           elseif(xi<-1d0) then
              xi = -1
           endif
!-- amplification factor
           if(.not.trn_noampfact .and. &
                 ((xi<0d0.and.x>0d0).or.(xi>0d0.and.x<0d0))) then
              help=1d0/abs(xi)
              help = min(100d0, help) !-- truncate singularity
!
!-- velocity effects accounting
              totevelo=totevelo-e*2d0 * &
                   (0.55d0*help-1.25d0*abs(xi))*abs(x)*cinv
!
!-- apply the excess (higher than factor 2d0) to the energy deposition
              eamp = e*2d0*0.55d0*max(0d0,help-2d0)*abs(x)*cinv
!-- apply limited correction to the particle
              help = min(2d0,help)
              e0=e0*(1d0 + 2d0*(0.55d0*help-1.25d0*abs(xi))*abs(x)*cinv)
              e=e*(1d0 + 2d0*(0.55d0*help-1.25d0*abs(xi))*abs(x)*cinv)
           endif
        endif
        help = (grd_cap(ig,l)+grd_sig(l)) * &
             dx(ixnext)*thelp
        alb = grd_fcoef(l)*grd_cap(ig,l)/ &
             (grd_cap(ig,l)+grd_sig(l))
        eps = (4d0/3d0)*sqrt(3d0*alb)/(1d0+pc_dext*sqrt(3d0*alb))
        beta = 1.5d0*alb*help**2+sqrt(3d0*alb*help**2 + &
             2.25d0*alb**2*help**4)
        pp = 0.5d0*eps*beta/(beta-0.75*eps*help)
        !pp = 4d0/(3d0*help+6d0*pc_dext)
!-- sampling
        call rnd_r(r1,rndstate)
        if (r1 < pp * (1d0+1.5d0*abs(xi))) then
           ptcl2%itype = 2
           if(grd_isvelocity) then
!-- velocity effects accounting
              totevelo=totevelo+e*(1d0-elabfact)
!
              e = e*elabfact
              e0 = e0*elabfact
              wl = wl/elabfact
           endif
           ix = ixnext
           ic = grd_icell(ix,iy,iz)
        else
           call rnd_r(r1,rndstate)
           call rnd_r(r2,rndstate)
           xi = -(ixnext-ix)*max(r1,r2)
           call rnd_r(r1,rndstate)
           eta = sqrt(1d0-xi**2)*cos(pc_pi2*r1)
!-- resampling z-cosine
           mu = sqrt(1d0-xi**2)*sin(pc_pi2*r1)
!-- resampling azimuthal
           om = atan2(eta,xi)
           if(om<0d0) om=om+pc_pi2
           if(grd_isvelocity) then
!-- transforming mu to lab
              mu=(mu+z*cinv)/(1d0+(x*xi+y*eta+z*mu)*cinv)
              if(mu>1d0) then
                 mu = 1d0
              elseif(mu<-1d0) then
                 mu = -1d0
              endif
!-- transforming om to lab
              om = atan2(eta+y*cinv,xi+x*cinv)
              if(om<0d0) om=om+pc_pi2
           endif
        endif
     endif

!
!-- y-bound
  elseif(d==dby) then

     l = grd_icell(ix,iynext,iz)
     if((grd_cap(ig,l)+grd_sig(l)) * &
          min(dx(ix),dy(iynext),dz(iz))*thelp < trn_tauddmc &
          .or.in_puretran) then
!-- IMC in adjacent cell
        iy = iynext
        ic = grd_icell(ix,iy,iz)
     else
!-- DDMC in adjacent cell
        if(grd_isvelocity) then
!-- transforming y-cosine to cmf
           eta = (eta-y*cinv)/elabfact
           if(eta>1d0) then
              eta = 1
           elseif(eta<-1d0) then
              eta = -1
           endif
!-- amplification factor
           if(.not.trn_noampfact .and. &
                 ((eta<0d0.and.y>0d0).or.(eta>0d0.and.y<0d0))) then
              help=1d0/abs(eta)
              help = min(100d0, help) !-- truncate singularity
!
!-- velocity effects accounting
              totevelo=totevelo-e*2d0 * &
                   (0.55d0*help-1.25d0*abs(eta))*abs(y)*cinv
!
!-- apply the excess (higher than factor 2d0) to the energy deposition
              eamp = e*2d0*0.55d0*max(0d0,help-2d0)*abs(y)*cinv
!-- apply limited correction to the particle
              help = min(2d0,help)
              e0=e0*(1d0 + 2d0*(0.55d0*help-1.25d0*abs(eta))*abs(y)*cinv)
              e=e*(1d0 + 2d0*(0.55d0*help-1.25d0*abs(eta))*abs(y)*cinv)
           endif
        endif
        help = (grd_cap(ig,l)+grd_sig(l)) * &
             dy(iynext)*thelp
        alb = grd_fcoef(l)*grd_cap(ig,l)/ &
             (grd_cap(ig,l)+grd_sig(l))
        eps = (4d0/3d0)*sqrt(3d0*alb)/(1d0+pc_dext*sqrt(3d0*alb))
        beta = 1.5d0*alb*help**2+sqrt(3d0*alb*help**2 + &
             2.25d0*alb**2*help**4)
        pp = 0.5d0*eps*beta/(beta-0.75*eps*help)
        !pp = 4d0/(3d0*help+6d0*pc_dext)
!-- sampling
        call rnd_r(r1,rndstate)
        if (r1 < pp * (1d0+1.5d0*abs(eta))) then
           ptcl2%itype = 2
           if(grd_isvelocity) then
!-- velocity effects accounting
              totevelo=totevelo+e*(1d0-elabfact)
!
              e = e*elabfact
              e0 = e0*elabfact
              wl = wl/elabfact
           endif
           iy = iynext
           ic = grd_icell(ix,iy,iz)
        else
           call rnd_r(r1,rndstate)
           call rnd_r(r2,rndstate)
           eta = -(iynext-iy)*max(r1,r2)
           call rnd_r(r1,rndstate)
           xi = sqrt(1d0-eta**2)*cos(pc_pi2*r1)
!-- resampling z-cosine
           mu = sqrt(1d0-eta**2)*sin(pc_pi2*r1)
!-- resampling azimuthal
           om = atan2(eta,xi)
           if(om<0d0) om=om+pc_pi2
           if(grd_isvelocity) then
!-- transforming mu to lab
              mu=(mu+z*cinv)/(1d0+(x*xi+y*eta+z*mu)*cinv)
              if(mu>1d0) then
                 mu = 1d0
              elseif(mu<-1d0) then
                 mu = -1d0
              endif
!-- transforming om to lab
              om = atan2(eta+y*cinv,xi+x*cinv)
              if(om<0d0) om=om+pc_pi2
           endif
        endif
     endif

!
!-- z-bound
  elseif(d==dbz) then

     l = grd_icell(ix,iy,iznext)
     if((grd_cap(ig,l)+grd_sig(l)) * &
          min(dx(ix),dy(iy),dz(iznext))*thelp < trn_tauddmc &
          .or.in_puretran) then
!-- IMC in adjacent cell
        iz = iznext
        ic = grd_icell(ix,iy,iz)
     else
!-- DDMC in adjacent cell
        if(grd_isvelocity) then
!-- transforming z-cosine to cmf
           mu = (mu-z*cinv)/elabfact
           if(mu>1d0) then
              mu = 1
           elseif(mu<-1d0) then
              mu = -1
           endif
!-- amplification factor
           if(.not.trn_noampfact .and.&
                 ((mu<0d0.and.z>0d0).or.(mu>0d0.and.z<0d0))) then
              help=1d0/abs(mu)
              help = min(100d0, help) !-- truncate singularity
!
!-- velocity effects accounting
              totevelo=totevelo-e*2d0 * &
                   (0.55d0*help-1.25d0*abs(mu))*abs(z)*cinv
!
!-- apply the excess (higher than factor 2d0) to the energy deposition
              eamp = e*2d0*0.55d0*max(0d0,help-2d0)*abs(x)*cinv
!-- apply limited correction to the particle
              help = min(2d0,help)
              e0=e0*(1d0 + 2d0*(0.55d0*help-1.25d0*abs(mu))*abs(z)*cinv)
              e=e*(1d0 + 2d0*(0.55d0*help-1.25d0*abs(mu))*abs(z)*cinv)
           endif
        endif
        help = (grd_cap(ig,l)+grd_sig(l)) * &
             dz(iznext)*thelp
        alb = grd_fcoef(l)*grd_cap(ig,l)/ &
             (grd_cap(ig,l)+grd_sig(l))
        eps = (4d0/3d0)*sqrt(3d0*alb)/(1d0+pc_dext*sqrt(3d0*alb))
        beta = 1.5d0*alb*help**2+sqrt(3d0*alb*help**2 + &
             2.25d0*alb**2*help**4)
        pp = 0.5d0*eps*beta/(beta-0.75*eps*help)
        !pp = 4d0/(3d0*help+6d0*pc_dext)
        !-- sampling
        call rnd_r(r1,rndstate)
        if (r1 < pp * (1d0+1.5d0*abs(mu))) then
           ptcl2%itype = 2
           if(grd_isvelocity) then
!-- velocity effects accounting
              totevelo=totevelo+e*(1d0-elabfact)
!
              e = e*elabfact
              e0 = e0*elabfact
              wl = wl/elabfact
           endif
           iz = iznext
           ic = grd_icell(ix,iy,iz)
        else
           call rnd_r(r1,rndstate)
           call rnd_r(r2,rndstate)
!-- resampling z-cosine
           mu = -(iznext-iz)*max(r1,r2)
           call rnd_r(r1,rndstate)
!-- resampling azimuthal
           om = pc_pi2*r1
           xi = sqrt(1d0-mu**2)*cos(om)
           eta = sqrt(1d0-mu**2)*sin(om)
           if(grd_isvelocity) then
!-- transforming mu to lab
              mu=(mu+z*cinv)/(1d0+(x*xi+y*eta+z*mu)*cinv)
              if(mu>1d0) then
                 mu = 1d0
              elseif(mu<-1d0) then
                 mu = -1d0
              endif
!-- transforming om to lab
              om = atan2(eta+y*cinv,xi+x*cinv)
              if(om<0d0) om=om+pc_pi2
           endif
        endif
     endif

!
!-- Doppler shift
  elseif(d==ddop) then
     if(.not.grd_isvelocity) then
!       stop 'transport3: ddop and no velocity'
        ierr = 5
        return
     endif
     if(ig<grp_ng) then
!-- shifting group
        ig = ig+1
        wl = (grp_wl(ig)+1d-6*(grp_wl(ig+1)-grp_wl(ig)))*elabfact
     else
!-- resampling wavelength in highest group
        call rnd_r(r1,rndstate)
        wl=1d0/(r1*grp_wlinv(grp_ng+1) + (1d0-r1)*grp_wlinv(grp_ng))
        wl = wl*elabfact
     endif
!-- check if ddmc region
     if((grd_sig(ic)+grd_cap(ig,ic)) * &
          min(dx(ix),dy(iy),dz(iz))*thelp >= trn_tauddmc &
          .and..not.in_puretran) then
        ptcl2%itype = 2
        if(grd_isvelocity) then
!-- velocity effects accounting
           totevelo=totevelo+e*(1d0-elabfact)
!
           e = e*elabfact
           e0 = e0*elabfact
           wl = wl/elabfact
        endif
     endif
  else
!    stop 'transport3: invalid distance'
     ierr = 6
     return
  endif

end subroutine transport3
! vim: fdm=marker
