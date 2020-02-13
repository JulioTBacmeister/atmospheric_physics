module gw_oro_phunit_mod
  
use shr_kind_mod,   only: r8=>shr_kind_r8, cl=>shr_kind_cl
use gw_common, only: pver , GWBand
use coords_1d,  only: Coords1D
use physics_types,  only: physics_ptend
use constituent, only: pcnst

! very wierd and convoluted way to get pcols into module
!use gw_rdg_phunit_mod, only : pcols

  implicit none
  private

  public :: gw_oro_phunit
  public :: set_band_oro


  ! anisotropic ridge fields
  integer, parameter :: prdg = 16

  !! integer :: pcols

  ! A mid-scale "band" with only stationary waves (l = 0).
  type(GWBand) :: band_oro

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
! nml settings with old GW scheme
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! &gw_drag_nl
! effgw_beres_dp		= 0.40D0
! effgw_cm		= 1.D0
! effgw_oro		= 0.125D0
! fcrit2		= 1.0
! frontgfc		= 3.0D-15
! gw_apply_tndmax		= .false.
! gw_dc		= 2.5D0
! gw_dc_long		= 0.D0
! gw_drag_file		= '/glade/p/cesmdata/cseg/inputdata/atm/waccm/gw/newmfspectra40_dc25.nc'
! gw_limit_tau_without_eff		= .true. 
! gw_lndscl_sgh		= .false.
! gw_oro_south_fac		= 2.d0   
! gw_polar_taper		= .true.
! gw_prndl		= 0.5d0  
! gw_qbo_hdepth_scaling		= 0.25D0
! pgwv		= 32
! pgwv_long		= 0
! tau_0_ubc		= .false.
! taubgnd		= 2.5D-3
!/
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

logical  :: gw_lndscl_sgh   = .false.	
!logical  :: gw_apply_tndmax = .false.	!- default in CAM5
logical  :: gw_apply_tndmax = .true.	!- default for Anisotropic: "Sean" limiters
real(r8) :: gw_oro_south_fac = 1.d0     ! 2.d0   

contains

subroutine gw_oro_phunit( band, &
   ncol, lchnk, dt, effgw_oro,  &
   u, v, t, p, piln, zm, zi,    &
   nm, ni, rhoi, kvtt, q, dse,  &
   sgh, landfrac,lats, &
   ptend, flx_heat)

   !!!use coords_1d,  only: Coords1D
   use gw_oro,     only: gw_oro_src
   use gw_common,  only: gw_drag_prof, energy_change

   type(GWBand),     intent(in) :: band         ! I hate this variable  ... it just hides information from view
   integer,          intent(in) :: ncol         ! number of atmospheric columns
   integer,          intent(in) :: lchnk        ! chunk identifier
   real(r8),         intent(in) :: dt           ! Time step.
   real(r8),         intent(in) :: effgw_oro

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

   real(r8),         intent(in) :: sgh(ncol)       ! Subgrid topo rms
   real(r8),         intent(in) :: landfrac(ncol)  ! land fraction
   real(r8),         intent(in) :: lats(ncol)      ! latitudes


   type(physics_ptend), intent(inout):: ptend   ! Parameterization net tendencies.

   real(r8),        intent(out) :: flx_heat(ncol)

   !---------------------------Local storage-------------------------------

   integer :: k, m, nn

   real(r8), allocatable :: tau(:,:,:)  ! wave Reynolds stress
   ! gravity wave wind tendency for each wave
   real(r8), allocatable :: gwut(:,:,:)
   ! Wave phase speeds for each column
   real(r8), allocatable :: c(:,:)

   ! Efficiency for a gravity wave source.
   real(r8) :: effgw(ncol)

   ! Indices of top gravity wave source level and lowest level where wind
   ! tendencies are allowed.
   integer :: src_level(ncol)
   integer :: tend_level(ncol)

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


   ! Wave Reynolds stresses at source level
   real(r8) :: tauoro(ncol)
   real(r8) :: taudsw(ncol)

   ! 
   real(r8) :: sgh_scaled(ncol)

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


   ! Energy change used by fixer.
   real(r8) :: de(ncol)

   character(len=1) :: cn
   character(len=9) :: fname(4)

   integer :: i,j

   !----------------------------------------------------------------------------

   ! Allocate wavenumber fields.
   allocate(tau(ncol,band%ngwv:band%ngwv,pver+1))
   allocate(gwut(ncol,pver,band%ngwv:band%ngwv))
   allocate(c(ncol,band%ngwv:band%ngwv))


   ! initialize accumulated momentum fluxes and tendencies
     if (gw_lndscl_sgh) then
        !!where (landfrac(:ncol) >= epsilon(1._r8))
        where (landfrac(:ncol) >= (1._r8)*1e-4 )
           effgw = effgw_oro * landfrac(:ncol)
           sgh_scaled = sgh(:ncol) / sqrt(landfrac(:ncol))
        elsewhere
           effgw = 0._r8
           sgh_scaled = 0._r8
        end where

        ! Determine the orographic wave source
        call gw_oro_src(ncol, band, p, &
           u, v, t, sgh_scaled, zm, nm, &
           src_level, tend_level, tau, ubm, ubi, xv, yv, c)
     else
        effgw = effgw_oro

        ! Determine the orographic wave source
        call gw_oro_src(ncol, band, p, &
             u, v, t, sgh(:ncol), zm, nm, &
             src_level, tend_level, tau, ubm, ubi, xv, yv, c)
     endif
     do i = 1, ncol
        if (lats(i) < 0._r8) then
           tau(i,:,:) = tau(i,:,:) * gw_oro_south_fac
        end if
     end do

#ifndef THREEDIM
write(511)  sgh,sgh_scaled
write(511)  src_level, tend_level
write(511)  tau(:,0,:),ubi,ubm,nm,zm,zi,u,v,t,p%mid,p%ifc
#endif

     ! satfac_in is 2 by default for CAM5

     ! Solve for the drag profile with orographic sources.
     call gw_drag_prof(ncol, band, p, src_level, tend_level,   dt,   &
          t,    &
          piln, rhoi,       nm,   ni, ubm,  ubi,  xv,    yv,   &
          effgw,c,          kvtt, q,  dse,  tau,  utgw,  vtgw, &
          ttgw, qtgw, egwdffi,  gwut, dttdf, dttke,            &
          satfac_in = 1._r8,                                   &
          lapply_effgw_in=gw_apply_tndmax)

write(511)  tau(:,0,:)

     ! For orographic waves, don't bother with taucd, since there are no
     ! momentum conservation routines or directional diagnostics.

     !  add the diffusion coefficients
     !do k = 1, pver+1
     !   egwdffi_tot(:,k) = egwdffi_tot(:,k) + egwdffi(:,k)
     !end do

     ! Add the orographic tendencies to the spectrum tendencies.
     ! Don't calculate fixers, since we are too close to the ground to
     ! spread momentum/energy differences across low layers.
     do k = 1, pver
        ptend%u(:ncol,k) = ptend%u(:ncol,k) + utgw(:,k)
        ptend%v(:ncol,k) = ptend%v(:ncol,k) + vtgw(:,k)
        ptend%s(:ncol,k) = ptend%s(:ncol,k) + ttgw(:,k)
     end do

#ifdef THREEDIM
write(511)  utgw
write(511)  vtgw
#endif

     ! Calculate energy change for output to CAM's energy checker.
     ! This is sort of cheating; we don't have a good a priori idea of the
     ! energy coming from surface stress, so we just integrate what we and
     ! actually have so far and overwrite flx_heat with that.
     call energy_change(dt, p, u, v, ptend%u(:ncol,:), &
          ptend%v(:ncol,:), ptend%s(:ncol,:), de)
     flx_heat(:ncol) = de

     do m = 1, pcnst
        do k = 1, pver
           ptend%q(:ncol,k,m) = ptend%q(:ncol,k,m) + qtgw(:,k,m)
        end do
     end do

   deallocate(tau, gwut, c)

end subroutine gw_oro_phunit



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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine set_band_oro(band_oro_in)
  type(GWBand) , intent(in) :: band_oro_in
! Just to bring band_oro into module
band_oro=band_oro_in
write(*,*) " Band ORO ",band_oro%ngwv
end subroutine set_band_oro


end module gw_oro_phunit_mod
