!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Module containing minimal fortran "physics unit (phunit)"
! for ridge based orographic gravity wave drag
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module gw_rdg_phunit_mod

use shr_kind_mod,   only: r8=>shr_kind_r8   !, cl=>shr_kind_cl
use gw_common, only: pver , GWBand
use coords_1d,  only: Coords1D
use physics_types,  only: physics_ptend
use constituent, only: pcnst
!----------------------------------
! Explcit dependencies
!----------------------------------
!  pcnst and pver are array sizes
!
!  Coords1D, GWBand, physics_ptend are derived types



  implicit none
  private

  !!public :: pcols
  public :: gw_rdg_phunit


  ! anisotropic ridge fields
  integer, parameter :: prdg = 16

  !!integer :: pcols


contains

subroutine gw_rdg_phunit( &
   type, band, ncol, lchnk, n_rdg, dt, &
   u, v, t, p, piln, zm, zi, &
   nm, ni, rhoi, kvtt, q, dse, &
   effgw_rdg, effgw_rdg_max, &
   hwdth, clngt, gbxar, &
   mxdis, angll, anixy, &
   rdg_cd_llb, trpd_leewv, &
   ptend, flx_heat)

   !use coords_1d,  only: Coords1D
   use gw_rdg,     only: gw_rdg_src, gw_rdg_belowpeak, gw_rdg_break_trap
   use gw_common,  only: gw_drag_prof, energy_change

   character(len=5), intent(in) :: type         ! BETA or GAMMA
   type(GWBand),     intent(in) :: band         ! I hate this variable  ... it just hides information from view
   integer,          intent(in) :: ncol         ! number of atmospheric columns
   integer,          intent(in) :: lchnk        ! chunk identifier
   integer,          intent(in) :: n_rdg
   real(r8),         intent(in) :: dt           ! Time step.

   real(r8),         intent(in) :: u(ncol,pver)    ! Midpoint zonal winds. ( m s-1)
   real(r8),         intent(in) :: v(ncol,pver)    ! Midpoint meridional winds. ( m s-1)
   real(r8),         intent(in) :: t(ncol,pver)    ! Midpoint temperatures. (K)
   type(Coords1D),   intent(in) :: p               ! Pressure coordinates.
   real(r8),         intent(in) :: piln(ncol,pver+1)  ! Log of interface pressures.
   real(r8),         intent(in) :: zm(ncol,pver)   ! Midpoint altitudes above ground (m).
   real(r8),         intent(in) :: zi(ncol,pver+1) ! Interface altitudes above ground (m).
   real(r8),         intent(in) :: nm(ncol,pver)   ! Midpoint Brunt-Vaisalla frequencies (s-1).
   real(r8),         intent(in) :: ni(ncol,pver+1) ! Interface Brunt-Vaisalla frequencies (s-1).
   real(r8),         intent(in) :: rhoi(ncol,pver+1) ! Interface density (kg m-3).
   real(r8),         intent(in) :: kvtt(ncol,pver+1) ! Molecular thermal diffusivity.
   real(r8),         intent(in) :: q(:,:,:)        ! Constituent array.
   real(r8),         intent(in) :: dse(ncol,pver)  ! Dry static energy.


   real(r8),         intent(in) :: effgw_rdg       ! Tendency efficiency.
   real(r8),         intent(in) :: effgw_rdg_max
   real(r8),         intent(in) :: hwdth(ncol,prdg) ! width of ridges.
   real(r8),         intent(in) :: clngt(ncol,prdg) ! length of ridges.
   real(r8),         intent(in) :: gbxar(ncol)      ! gridbox area

   real(r8),         intent(in) :: mxdis(ncol,prdg) ! Height estimate for ridge (m).
   real(r8),         intent(in) :: angll(ncol,prdg) ! orientation of ridges.
   real(r8),         intent(in) :: anixy(ncol,prdg) ! Anisotropy parameter.

   real(r8),         intent(in) :: rdg_cd_llb      ! Drag coefficient for low-level flow
   logical,          intent(in) :: trpd_leewv

   type(physics_ptend), intent(inout):: ptend   ! Parameterization net tendencies.

          ! this was dimensioned pcols before: But who the fuck understands when to use ncol or pcols
   real(r8),        intent(out) :: flx_heat(ncol)

   !---------------------------Local storage-------------------------------

   integer :: k, m, nn

   real(r8), allocatable :: tau(:,:,:)  ! wave Reynolds stress
   ! gravity wave wind tendency for each wave
   real(r8), allocatable :: gwut(:,:,:)
   ! Wave phase speeds for each column
   real(r8), allocatable :: c(:,:)

   ! Isotropic source flag [anisotropic orography].
   integer  :: isoflag(ncol)

   ! horiz wavenumber [anisotropic orography].
   real(r8) :: kwvrdg(ncol)

   ! Efficiency for a gravity wave source.
   real(r8) :: effgw(ncol)

   ! Indices of top gravity wave source level and lowest level where wind
   ! tendencies are allowed.
   integer :: src_level(ncol)
   integer :: tend_level(ncol)
   integer :: bwv_level(ncol)
   integer :: tlb_level(ncol)

   ! Projection of wind at midpoints and interfaces.
   real(r8) :: ubm(ncol,pver)
   real(r8) :: ubi(ncol,pver+1)

   ! Unit vectors of source wind (zonal and meridional components).
   real(r8) :: xv(ncol)
   real(r8) :: yv(ncol)

   ! Averages over source region.
   real(r8) :: ubmsrc(ncol) ! On-ridge wind.
   real(r8) :: usrc(ncol)   ! Zonal wind.
   real(r8) :: vsrc(ncol)   ! Meridional wind.
   real(r8) :: nsrc(ncol)   ! B-V frequency.
   real(r8) :: rsrc(ncol)   ! Density.

   ! normalized wavenumber
   real(r8) :: m2src(ncol)

   ! Top of low-level flow layer.
   real(r8) :: tlb(ncol)

   ! Bottom of linear wave region.
   real(r8) :: bwv(ncol)

   ! Froude numbers for flow/drag regimes
   real(r8) :: Fr1(ncol)
   real(r8) :: Fr2(ncol)
   real(r8) :: Frx(ncol)

   ! Wave Reynolds stresses at source level
   real(r8) :: tauoro(ncol)
   real(r8) :: taudsw(ncol)

   ! Surface streamline displacement height for linear waves.
   real(r8) :: hdspwv(ncol)

   ! Surface streamline displacement height for downslope wind regime.
   real(r8) :: hdspdw(ncol)

   ! Wave breaking level
   real(r8) :: wbr(ncol)

   real(r8) :: utgw(ncol,pver)       ! zonal wind tendency
   real(r8) :: vtgw(ncol,pver)       ! meridional wind tendency
   real(r8) :: ttgw(ncol,pver)       ! temperature tendency
   real(r8) :: qtgw(ncol,pver,pcnst) ! constituents tendencies

   ! Effective gravity wave diffusivity at interfaces.
   real(r8) :: egwdffi(ncol,pver+1)

   ! Temperature tendencies from diffusion and kinetic energy.
   real(r8) :: dttdf(ncol,pver)
   real(r8) :: dttke(ncol,pver)

   ! Wave stress in zonal/meridional direction
   real(r8) :: taurx(ncol,pver+1)
   real(r8) :: taurx0(ncol,pver+1)
   real(r8) :: taury(ncol,pver+1)
   real(r8) :: taury0(ncol,pver+1)

   ! U,V tendency accumulators
   real(r8) :: utrdg(ncol,pver)
   real(r8) :: vtrdg(ncol,pver)

   ! Energy change used by fixer.
   real(r8) :: de(ncol)

   character(len=1) :: cn
   character(len=9) :: fname(4)
   !----------------------------------------------------------------------------

   ! Allocate wavenumber fields.
   allocate(tau(ncol,  band%ngwv:band%ngwv  , pver+1))
   allocate(gwut(ncol,pver,band%ngwv:band%ngwv  ))
   allocate(c(ncol,band%ngwv:band%ngwv))

   ! initialize accumulated momentum fluxes and tendencies
   taurx = 0._r8
   taury = 0._r8 
   utrdg = 0._r8
   vtrdg = 0._r8

   do nn = 1, n_rdg
  
      kwvrdg  = 0.001_r8 / ( hwdth(:,nn) + 0.001_r8 ) ! this cant be done every time step !!!
      isoflag = 0   
      effgw   = effgw_rdg * ( hwdth(1:ncol,nn)* clngt(1:ncol,nn) ) / gbxar(1:ncol)
      effgw   = min( effgw_rdg_max , effgw )

    call gw_rdg_src(ncol, band, p, &
         u, v, t, mxdis(:,nn), angll(:,nn), anixy(:,nn), kwvrdg, isoflag, zi, nm, &
         src_level, tend_level, bwv_level, tlb_level, tau, ubm, ubi, xv, yv,  & 
         ubmsrc, usrc, vsrc, nsrc, rsrc, m2src, tlb, bwv, Fr1, Fr2, Frx, c)

#ifndef THREEDIM
write(511) mxdis(:,nn),angll(:,nn),anixy(:,nn),hwdth(:,nn),clngt(:,nn),kwvrdg,gbxar
write(511)  ubmsrc,usrc,vsrc,nsrc, tlb, bwv, Fr1, Fr2, Frx, rsrc
write(511)  src_level, tend_level, bwv_level, tlb_level
write(511)  tau(:,0,:),ubi,ubm,nm,zm,zi,u,v,t,p%mid,p%ifc
#endif

    call gw_rdg_belowpeak(ncol, band, rdg_cd_llb, &
         t, mxdis(:,nn), anixy(:,nn), kwvrdg, & 
         zi, nm, ni, rhoi, &
         src_level, tau, & 
         ubmsrc, nsrc, rsrc, m2src, tlb, bwv, Fr1, Fr2, Frx, & 
         tauoro, taudsw, hdspwv, hdspdw)

#ifndef THREEDIM
write(511)  tauoro,taudsw, hdspwv, hdspdw,effgw,tau(:,0,:)
#endif

    call gw_rdg_break_trap(ncol, band, &
         zi, nm, ni, ubm, ubi, rhoi, kwvrdg , bwv, tlb, wbr, & 
         src_level, tlb_level, hdspwv, hdspdw,  mxdis(:,nn), & 
         tauoro, taudsw, tau, & 
         ldo_trapped_waves=trpd_leewv)
     
#ifndef THREEDIM
write(511)  tau(:,0,:), wbr 
#endif

    call gw_drag_prof(ncol, band, p, src_level, tend_level, dt, &
         t,    &
         piln, rhoi, nm, ni, ubm, ubi, xv, yv,   &
         effgw, c, kvtt, q, dse, tau, utgw, vtgw, &
         ttgw, qtgw, egwdffi,   gwut, dttdf, dttke, &
         kwvrdg=kwvrdg, & 
         satfac_in = 1._r8 )

#ifndef THREEDIM
write(511)  tau(:,0,:)
#endif

      ! Add the tendencies from each ridge to the totals.
      do k = 1, pver
         ! diagnostics
         utrdg(:,k) = utrdg(:,k) + utgw(:,k)
         vtrdg(:,k) = vtrdg(:,k) + vtgw(:,k)
         ! physics tendencies
         ptend%u(:ncol,k) = ptend%u(:ncol,k) + utgw(:,k)
         ptend%v(:ncol,k) = ptend%v(:ncol,k) + vtgw(:,k)
         ptend%s(:ncol,k) = ptend%s(:ncol,k) + ttgw(:,k)
      end do

      write(*,*) "MAX abs(UTGW) ",nn,maxval(abs(utgw))*86400.

      do m = 1, pcnst
         do k = 1, pver
            ptend%q(:ncol,k,m) = ptend%q(:ncol,k,m) + qtgw(:,k,m)
         end do
      end do

      do k = 1, pver+1
         taurx0(:,k) =  tau(:,0,k)*xv
         taury0(:,k) =  tau(:,0,k)*yv
         taurx(:,k)  =  taurx(:,k) + taurx0(:,k)
         taury(:,k)  =  taury(:,k) + taury0(:,k)
      end do

      if (nn == 1) then
      end if

      if (nn <= 6) then
         write(cn, '(i1)') nn
      end if

   end do ! end of loop over multiple ridges

#ifdef THREEDIM
write(511)  utrdg
write(511)  vtrdg
#endif

   ! Calculate energy change for output to CAM's energy checker.
   call energy_change(dt, p, u, v, ptend%u(:ncol,:), &
          ptend%v(:ncol,:), ptend%s(:ncol,:), de)
   flx_heat(:ncol) = de


   if (trim(type) == 'BETA') then
      fname(1) = 'TAUGWX'
      fname(2) = 'TAUGWY'
      fname(3) = 'UTGWORO'
      fname(4) = 'VTGWORO'
   else if (trim(type) == 'GAMMA') then
      fname(1) = 'TAURDGGMX'
      fname(2) = 'TAURDGGMY'
      fname(3) = 'UTRDGGM'
      fname(4) = 'VTRDGGM'
   else
      call endrun('gw_rdg_calc: FATAL: type must be either BETA or GAMMA'&
                  //' type= '//type)
   end if


   deallocate(tau, gwut, c)

 end subroutine gw_rdg_phunit


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine endrun(msg)

   integer :: iulog

   character(len=*), intent(in), optional :: msg    ! string to be printed

    iulog=6

   if (present (msg)) then
      write(iulog,*)'ENDRUN:', msg
   else
      write(iulog,*)'ENDRUN: called without a message string'
   end if

   stop

end subroutine endrun


end module gw_rdg_phunit_mod
