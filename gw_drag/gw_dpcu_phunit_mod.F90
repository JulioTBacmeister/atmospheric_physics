module gw_dpcu_phunit_mod

  
use shr_const_mod,   only: pi=>shr_const_pi    !, cl=>shr_kind_cl
use shr_kind_mod,   only: r8=>shr_kind_r8    !, cl=>shr_kind_cl
use gw_common, only: pver , GWBand
use coords_1d,  only: Coords1D
use physics_types,  only: physics_ptend
use constituent, only: pcnst
use gw_convect,     only: BeresSourceDesc


! very wierd and convoluted way to get pcols into module
!use gw_rdg_phunit_mod, only : pcols

  implicit none
  private

  public :: gw_dpcu_phunit


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
real(r8) :: gw_oro_south_fac = 1.d0     ! 2.d0   

contains

subroutine gw_dpcu_phunit( band, &
   ncol, lchnk, dt, effgw_dp,  &
   u, v, t, p, piln, zm, zi,    &
   nm, ni, rhoi, kvtt, q, dse,  &
   netdt,desc,lats, &
   ptend, flx_heat)

   !!!use coords_1d,  only: Coords1D
   use gw_convect,     only: gw_beres_src
   use gw_common,  only: gw_drag_prof, energy_change

   type(BeresSourceDesc), intent(in) :: desc
   type(GWBand),     intent(in) :: band         ! I hate this variable  ... it just hides information from view
   integer,          intent(in) :: ncol         ! number of atmospheric columns
   integer,          intent(in) :: lchnk        ! chunk identifier
   real(r8),         intent(in) :: dt           ! Time step.
   real(r8),         intent(in) :: effgw_dp

   real(r8),         intent(in) :: u(ncol,pver)      ! Midpoint zonal winds. ( m s-1)
   real(r8),         intent(in) :: v(ncol,pver)      ! Midpoint meridional winds. ( m s-1)
   real(r8),         intent(in) :: t(ncol,pver)      ! Midpoint temperatures. (K)
   real(r8),         intent(in) :: netdt(ncol,pver)  ! Convective heating rate (K s-1)
   type(Coords1D),   intent(in) :: p                 ! Pressure coordinates.
   real(r8),         intent(in) :: piln(ncol,pver+1) ! Log of interface pressures.
   real(r8),         intent(in) :: zm(ncol,pver)     ! Midpoint altitudes above ground (m).
   real(r8),         intent(in) :: zi(ncol,pver+1)   ! Interface altitudes above ground (m).
   real(r8),         intent(in) :: nm(ncol,pver)     ! Midpoint Brunt-Vaisalla frequencies (s-1).
   real(r8),         intent(in) :: ni(ncol,pver+1)   ! Interface Brunt-Vaisalla frequencies (s-1).
   real(r8),         intent(in) :: rhoi(ncol,pver+1) ! Interface density (kg m-3).
   real(r8),         intent(in) :: kvtt(ncol,pver+1) ! Molecular thermal diffusivity.
   real(r8),         intent(in) :: q(:,:,:)          ! Constituent array.
   real(r8),         intent(in) :: dse(ncol,pver)    ! Dry static energy.

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

   ! Wave breaking level
   real(r8) :: wbr(ncol)

   real(r8) :: utgw(ncol,pver)       ! zonal wind tendency
   real(r8) :: vtgw(ncol,pver)       ! meridional wind tendency
   real(r8) :: ttgw(ncol,pver)       ! temperature tendency
   real(r8) :: qtgw(ncol,pver,pcnst) ! constituents tendencies

   ! Heating depth [m] and maximum heating in each column.
   real(r8) :: hdepth(ncol), maxq0(ncol)

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

   logical  :: gw_apply_tndmax  	!- default .TRUE. for Anisotropic: "Sean" limiters

   
   !----------------------------------------------------------------------------

   ! Allocate wavenumber fields.
   allocate(tau(ncol,-band%ngwv:band%ngwv,pver+1))
   allocate(gwut(ncol,pver,-band%ngwv:band%ngwv))
   allocate(c(ncol,-band%ngwv:band%ngwv))

     gw_apply_tndmax  = .FALSE.


     ! Efficiency of gravity wave momentum transfer.
     ! This is really only to remove the pole points.
     where (pi/2._r8 - abs(lats(:ncol)) >= 1.e-4 )  !-4*epsilon(1._r8))
        effgw = effgw_dp
     elsewhere
        effgw = 0._r8
     end where

     ! Determine wave sources for Beres deep scheme
     call gw_beres_src(ncol, band , desc, &
          u, v, netdt, zm, src_level, tend_level, tau, &
          ubm, ubi, xv, yv, c, hdepth, maxq0)


     ! satfac_in is 2 by default for CAM5

     ! Solve for the drag profile with orographic sources.
     call gw_drag_prof(ncol, band, p, src_level, tend_level,   dt,   &
          t,    &
          piln, rhoi,       nm,   ni, ubm,  ubi,  xv,    yv,   &
          effgw,c,          kvtt, q,  dse,  tau,  utgw,  vtgw, &
          ttgw, qtgw, egwdffi,  gwut, dttdf, dttke,            &
          satfac_in = 1._r8,                                   &
          lapply_effgw_in=gw_apply_tndmax)


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

end subroutine gw_dpcu_phunit



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


end module gw_dpcu_phunit_mod
