module gw_convect

!
! This module handles gravity waves from convection, and was extracted from
! gw_drag in May 2013.
!

use shr_kind_mod,   only: r8=>shr_kind_r8    !, cl=>shr_kind_cl
use shr_const_mod,  only: pi => shr_const_pi
use gw_common,      only: gw_drag_prof, gw_prof, GWBand, gw_rair, gw_cpair
use gw_utils,       only: get_unit_vector, dot_2d, midpoint_interp



  

  
implicit none
private
save

!public :: BeresSourceDesc
!public :: gw_beres_src

public :: gw_beres_run
public :: gw_beres_init

type :: BeresSourceDesc
   ! Whether wind speeds are shifted to be relative to storm cells.
   logical :: storm_shift
   ! Index for level where wind speed is used as the source speed.
   integer :: k
   ! Heating depths below this value [m] will be ignored.
   real(r8) :: min_hdepth
   ! Table bounds, for convenience. (Could be inferred from shape(mfcc).)
   integer :: maxh
   integer :: maxuh
   ! Heating depths [m].
   real(r8), allocatable :: hd(:)
   ! Table of source spectra.
   real(r8), allocatable :: mfcc(:,:,:)
end type BeresSourceDesc

! Some description of GW spectrum
type(GWBand)   :: band         ! I hate this variable  ... it just hides information from view

! Source dscriptor for Beres Scheme (Deep). Created in gw_beres_init.
! This is needed in every processor,
type(BeresSourceDesc) :: desc

! Run time parameters
real(r8) :: qbo_hdepth_scaling

! Define PI
! real(r8) :: PI


contains

!==========================================================================

!------------------------------------
!> \section arg_table_gw_beres_init  Argument Table
!! \htmlinclude gw_beres_init.html
subroutine gw_beres_init (file_name , pver, pverp, pref_edge )
#include <netcdf.inc>

  
  
  character(len=*), intent(in) :: file_name
  !!type(GWBand), intent(inout) :: band

  integer,  intent(in) :: pver,pverp
  real(r8), intent(in) :: pref_edge( pverp )
  
  !!type(BeresSourceDesc), intent(inout) :: desc


  ! Stuff for Beres convective gravity wave source.
  real(r8), allocatable :: mfcc(:,:,:), hdcc(:)
  integer  :: hd_mfcc , mw_mfcc, ps_mfcc, ngwv_file, ps_mfcc_mid

  real(r8) :: gw_dc, wavelength
  integer  :: pgwv , k


  ! Vars needed by NetCDF operators
  integer  :: ncid, dimid, varid, status

  ! Set values for runtime parameters
  qbo_hdepth_scaling = 0.25_r8





  
  status = nf_open(file_name , 0, ncid)


    status = NF_INQ_DIMID(ncid, 'PS', dimid)
    IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
    status = NF_INQ_DIMLEN(ncid, dimid, ps_mfcc )

    status = NF_INQ_DIMID(ncid, 'MW', dimid)
    IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
    status = NF_INQ_DIMLEN(ncid, dimid, mw_mfcc )

    status = NF_INQ_DIMID(ncid, 'HD', dimid)
    IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
    status = NF_INQ_DIMLEN(ncid, dimid, hd_mfcc )

 !  allocate( mfcc(ps_mfcc , mw_mfcc, hd_mfcc ) )
  allocate( mfcc(hd_mfcc , mw_mfcc, ps_mfcc ) )
  allocate( hdcc(hd_mfcc) )
   
  status = NF_INQ_VARID(ncid, 'HD', varid)
  IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
  status = NF_GET_VAR_DOUBLE(ncid, varid, hdcc )
  IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)

  status = NF_INQ_VARID(ncid, 'mfcc', varid)
  IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
  status = NF_GET_VAR_DOUBLE(ncid, varid, mfcc )
  IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)




  status = nf_close (ncid)


! Need to call GWBand for convective waves
!-----------------
! band_mid = GWBand(pgwv, gw_dc, 1.0_r8, wavelength_mid)
!------------------------------------------
!!! WHERE is gw_dc set? Where is pgwv set?
!!! ==> They are set in atm_in. Out of the box values seem to be
!!!  pgwv = 32
!!!  gw_dc = 2.5D0  

  ! Hardwire for now
  gw_dc = 2.5_r8
  pgwv  = 32
  wavelength = 1.e5_r8
  band  = GWBand(pgwv, gw_dc, 1.0_r8, wavelength )

  ! Get HD (heating depth) dimension.
  desc%maxh = HD_MFCC  ! get_pio_dimlen(gw_file_desc, "HD", file_path)

  ! Get MW (mean wind) dimension.
  desc%maxuh = MW_MFCC ! get_pio_dimlen(gw_file_desc, "MW", file_path)

  ! Get PS (phase speed) dimension.
  ngwv_file = ps_mfcc  ! get_pio_dimlen(gw_file_desc, "PS", file_path)

  ! Number in each direction is half of total (and minus phase speed of 0).
  desc%maxuh = (desc%maxuh-1)/2
  ngwv_file = (ngwv_file-1)/2


    allocate(desc%hd(desc%maxh) , stat=status )

    allocate(desc%mfcc(desc%maxh,-desc%maxuh:desc%maxuh,&
       -band%ngwv:band%ngwv), stat=status )

    ! Don't understand point of having uh dimension go from -maxuh to +maxuh

    ! midpoint of spectrum in netcdf file is ps_mfcc (odd number) divided by 2, plus 1
    ! E.g., ps_mfcc = 5. Integer divide; ps_mfcc/2=2. Add one. So, ps_mfcc_mid = 3
    !       1  2  3  4  5
    !      -2 -1  0 +1 +2
    ps_mfcc_mid = INT( ps_mfcc/2) + 1


    desc%mfcc( : , -desc%maxuh:desc%maxuh ,  -band%ngwv:band%ngwv )  &
           = mfcc( :, :,  -band%ngwv+ps_mfcc_mid  :  band%ngwv+ ps_mfcc_mid )
  
    ! While not currently documented in the file, it uses kilometers. Convert
    ! to meters.
    desc%hd = hdcc *1000._r8

    desc%storm_shift=.TRUE.

    do k = 0, pver
        ! 700 hPa index
        if (pref_edge(k+1) < 70000._r8) desc%k = k+1
    end do

end subroutine gw_beres_init



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!  Main Interface
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!------------------------------------
!> \section arg_table_gw_beres_run  Argument Table
!! \htmlinclude gw_beres_run.html
subroutine gw_beres_run(  &
   ncol, pver, pverp, pcnst, dt, effgw_dp,  &
   u, v, t, pint, pmid, delp, piln, zm, zi,    &
   kvtt, q, dse,  &
   netdt, lats, &
   flx_heat, utgw, vtgw , ttgw, qtgw )


   !type(BeresSourceDesc), intent(in) :: desc
   !type(GWBand),     intent(in) :: band         ! I hate this variable  ... it just hides information from view
   integer,          intent(in) :: ncol         ! number of atmospheric columns
   integer,          intent(in) :: pver         ! Layer Vertical dimension.
   integer,          intent(in) :: pverp        ! Intfc Vertical dimension.
   integer,          intent(in) :: pcnst        ! consituent dimension
   real(r8),         intent(in) :: dt           ! Time step.
   real(r8),         intent(in) :: effgw_dp

   real(r8),         intent(in) :: u(ncol,pver)      ! Midpoint zonal winds. ( m s-1)
   real(r8),         intent(in) :: v(ncol,pver)      ! Midpoint meridional winds. ( m s-1)
   real(r8),         intent(in) :: t(ncol,pver)      ! Midpoint temperatures. (K)
   real(r8),         intent(in) :: netdt(ncol,pver)  ! Convective heating rate (K s-1)
   real(r8),         intent(in) :: delp(ncol,pver)   ! Delta interface pressures.
   real(r8),         intent(in) :: pint(ncol,pverp) ! Interface pressures.
   real(r8),         intent(in) :: pmid(ncol,pver)   ! Mid-point pressures.
   real(r8),         intent(in) :: piln(ncol,pverp) ! Log of interface pressures.
   real(r8),         intent(in) :: zm(ncol,pver)     ! Midpoint altitudes above ground (m).
   real(r8),         intent(in) :: zi(ncol,pverp)   ! Interface altitudes above ground (m).
   real(r8),         intent(in) :: kvtt(ncol,pverp) ! Molecular thermal diffusivity.
   real(r8),         intent(in) :: q(ncol,pver,pcnst) ! Constituent array.
   real(r8),         intent(in) :: dse(ncol,pver)     ! Dry static energy.

   real(r8),         intent(in) :: lats(ncol)      ! latitudes

   real(r8),        intent(out) :: flx_heat(ncol)
   real(r8),        intent(out) :: utgw(ncol,pver)       ! zonal wind tendency
   real(r8),        intent(out) :: vtgw(ncol,pver)       ! meridional wind tendency
   real(r8),        intent(out) :: ttgw(ncol,pver)       ! temperature tendency
   real(r8),        intent(out) :: qtgw(ncol,pver,pcnst) ! consituent tendency


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

   real(r8) :: nm(ncol,pver)   ! Midpoint Brunt-Vaisalla frequencies (s-1).
   real(r8) :: ni(ncol,pverp) ! Interface Brunt-Vaisalla frequencies (s-1).
   real(r8) :: rhoi(ncol,pverp) ! Interface density (kg m-3).

   ! Projection of wind at midpoints and interfaces.
   real(r8) :: ubm(ncol,pver)
   real(r8) :: ubi(ncol,pverp)

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


   ! Heating depth [m] and maximum heating in each column.
   real(r8) :: hdepth(ncol), maxq0(ncol)

   ! Effective gravity wave diffusivity at interfaces.
   real(r8) :: egwdffi(ncol,pverp)

   ! Temperature tendencies from diffusion and kinetic energy.
   real(r8) :: dttdf(ncol,pver)
   real(r8) :: dttke(ncol,pver)

   ! Wave stress in zonal/meridional direction
   real(r8) :: taurx(ncol,pverp)
   real(r8) :: taurx0(ncol,pverp)
   real(r8) :: taury(ncol,pverp)
   real(r8) :: taury0(ncol,pverp)


   ! Energy change used by fixer.
   real(r8) :: de(ncol)
   logical  :: gw_apply_tndmax  	!- default .TRUE. for Anisotropic: "Sean" limiters

   character(len=1) :: cn
   character(len=9) :: fname(4)

   integer :: i,j


   ! Calculate necessary thermodyanmic profiles for GW
   call gw_prof(ncol, pver,  pint, pmid, gw_cpair, t , rhoi, nm, ni)

   !----------------------------------------------------------------------------

   ! Allocate wavenumber fields.
   allocate(tau(ncol,-band%ngwv:band%ngwv,pverp))
   allocate(gwut(ncol,pver,-band%ngwv:band%ngwv))
   allocate(c(ncol,-band%ngwv:band%ngwv))

     gw_apply_tndmax  = .FALSE.

     flx_heat = 0._r8
     

     ! Efficiency of gravity wave momentum transfer.
     ! This is really only to remove the pole points.
     where (pi/2._r8 - abs(lats(:ncol)) >= 1.e-4 )  !-4*epsilon(1._r8))
        effgw = effgw_dp
     elsewhere
        effgw = 0._r8
     end where

     
     ! Determine wave sources for Beres deep scheme
     call gw_beres_src(ncol, pver ,   &
          u, v, netdt, zm, src_level, tend_level, tau, &
          ubm, ubi, xv, yv, c, hdepth, maxq0)


     ! satfac_in is 2 by default for CAM5

     ! Solve for the drag profile with orographic sources.
     call gw_drag_prof(ncol, pver, band, &
          pint, pmid, delp, src_level, tend_level,   dt,   &
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
     
   deallocate(tau, gwut, c)

 end subroutine gw_beres_run


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Non-interface routines
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!------------------------------------
subroutine gw_beres_src(ncol, pver , u, v, &
     netdt, zm, src_level, tend_level, tau, ubm, ubi, xv, yv, &
     c, hdepth, maxq0)
!-----------------------------------------------------------------------
! Driver for multiple gravity wave drag parameterization.
!
! The parameterization is assumed to operate only where water vapor
! concentrations are negligible in determining the density.
!
! Beres, J.H., M.J. Alexander, and J.R. Holton, 2004: "A method of
! specifying the gravity wave spectrum above convection based on latent
! heating properties and background wind". J. Atmos. Sci., Vol 61, No. 3,
! pp. 324-337.
!
!-----------------------------------------------------------------------
!------------------------------Arguments--------------------------------
  ! Column dimension.
  integer, intent(in) :: ncol
  ! Vertical dimension.
  integer, intent(in) :: pver

  ! Wavelengths triggered by convection.
  !!type(GWBand), intent(in) :: band

  ! Settings for convection type (e.g. deep vs shallow).
  !!type(BeresSourceDesc), intent(in) :: desc

  ! Midpoint zonal/meridional winds.
  real(r8), intent(in) :: u(ncol,pver), v(ncol,pver)
  ! Heating rate due to convection.
  real(r8), intent(in) :: netdt(:,:)
  ! Midpoint altitudes.
  real(r8), intent(in) :: zm(ncol,pver)

  ! Indices of top gravity wave source level and lowest level where wind
  ! tendencies are allowed.
  integer, intent(out) :: src_level(ncol)
  integer, intent(out) :: tend_level(ncol)

  ! Wave Reynolds stress.
  real(r8), intent(out) :: tau(ncol,-band%ngwv:band%ngwv,pver+1)
  ! Projection of wind at midpoints and interfaces.
  real(r8), intent(out) :: ubm(ncol,pver), ubi(ncol,pver+1)
  ! Unit vectors of source wind (zonal and meridional components).
  real(r8), intent(out) :: xv(ncol), yv(ncol)
  ! Phase speeds.
  real(r8), intent(out) :: c(ncol,-band%ngwv:band%ngwv)

  ! Heating depth [m] and maximum heating in each column.
  real(r8), intent(out) :: hdepth(ncol), maxq0(ncol)

!---------------------------Local Storage-------------------------------
  ! Column and level indices.
  integer :: i, k

  ! Zonal/meridional wind at roughly the level where the convection occurs.
  real(r8) :: uconv(ncol), vconv(ncol)

  ! Maximum heating rate.
  real(r8) :: q0(ncol)

  ! Bottom/top heating range index.
  integer  :: boti(ncol), topi(ncol)
  ! Index for looking up heating depth dimension in the table.
  integer  :: hd_idx(ncol)
  ! Mean wind in heating region.
  real(r8) :: uh(ncol)
  ! Min/max wavenumber for critical level filtering.
  integer :: Umini(ncol), Umaxi(ncol)
  ! Source level tau for a column.
  real(r8) :: tau0(-band%ngwv:band%ngwv)
  ! Speed of convective cells relative to storm.
  real(r8) :: CS(ncol)
  ! Index to shift spectra relative to ground.
  integer :: shift

  ! Heating rate conversion factor.
  real(r8), parameter :: CF = 20._r8
  ! Averaging length.
  real(r8), parameter :: AL = 1.0e5_r8

  !----------------------------------------------------------------------
  ! Initialize tau array
  !----------------------------------------------------------------------
  tau = 0.0_r8
  hdepth = 0.0_r8
  q0 = 0.0_r8
  tau0 = 0.0_r8

  !------------------------------------------------------------------------
  ! Determine wind and unit vectors approximately at the source level, then
  ! project winds.
  !------------------------------------------------------------------------

  ! Source wind speed and direction.
  uconv = u(:,desc%k)
  vconv = v(:,desc%k)

  ! Get the unit vector components and magnitude at the source level.
  call get_unit_vector(uconv, vconv, xv, yv, ubi(:,desc%k+1))

  ! Project the local wind at midpoints onto the source wind.
  do k = 1, pver
     ubm(:,k) = dot_2d(u(:,k), v(:,k), xv, yv)
  end do

  ! Compute the interface wind projection by averaging the midpoint winds.
  ! Use the top level wind at the top interface.
  ubi(:,1) = ubm(:,1)

  ubi(:,2:pver) = midpoint_interp(ubm)

  !-----------------------------------------------------------------------
  ! Calculate heating depth.
  !
  ! Heating depth is defined as the first height range from the bottom in
  ! which heating rate is continuously positive.
  !-----------------------------------------------------------------------

  ! First find the indices for the top and bottom of the heating range.
  boti = 0
  topi = 0
  do k = pver, 1, -1
     do i = 1, ncol
        if (boti(i) == 0) then
           ! Detect if we are outside the maximum range (where z = 20 km).
           if (zm(i,k) >= 20000._r8) then
              boti(i) = k
              topi(i) = k
           else
              ! First spot where heating rate is positive.
              if (netdt(i,k) > 0.0_r8) boti(i) = k
           end if
        else if (topi(i) == 0) then
           ! Detect if we are outside the maximum range (z = 20 km).
           if (zm(i,k) >= 20000._r8) then
              topi(i) = k
           else
              ! First spot where heating rate is no longer positive.
              if (.not. (netdt(i,k) > 0.0_r8)) topi(i) = k
           end if
        end if
     end do
     ! When all done, exit.
     if (all(topi /= 0)) exit
  end do

  ! Heating depth in m.
  hdepth = [ ( (zm(i,topi(i))-zm(i,boti(i))), i = 1, ncol ) ]

  ! J. Richter: this is an effective reduction of the GW phase speeds (needed to drive the QBO)
  hdepth = hdepth*qbo_hdepth_scaling

  hd_idx = index_of_nearest(hdepth, desc%hd)

  ! hd_idx=0 signals that a heating depth is too shallow, i.e. that it is
  ! either not big enough for the lowest table entry, or it is below the
  ! minimum allowed for this convection type.
  ! Values above the max in the table still get the highest value, though.
  where (hdepth < max(desc%min_hdepth, desc%hd(1))) hd_idx = 0

  ! Maximum heating rate.
  do k = minval(topi), maxval(boti)
     where (k >= topi .and. k <= boti)
        q0 = max(q0, netdt(:,k))
     end where
  end do

  !output max heating rate in K/day
  maxq0 = q0*24._r8*3600._r8

  ! Multipy by conversion factor
  q0 = q0 * CF

  if (desc%storm_shift) then

     ! Find the cell speed where the storm speed is > 10 m/s.
     ! Storm speed is taken to be the source wind speed.
     CS = sign(max(abs(ubm(:,desc%k))-10._r8, 0._r8), ubm(:,desc%k))

     ! Average wind in heating region, relative to storm cells.
     uh = 0._r8
     do k = minval(topi), maxval(boti)
        where (k >= topi .and. k <= boti)
           uh = uh + ubm(:,k)/(boti-topi+1)
        end where
     end do

     uh = uh - CS

  else

     ! For shallow convection, wind is relative to ground, and "heating
     ! region" wind is just the source level wind.
     uh = ubm(:,desc%k)

  end if

  ! Limit uh to table range.
  uh = min(uh, real(desc%maxuh, r8))
  uh = max(uh, -real(desc%maxuh, r8))

  ! Speeds for critical level filtering.
  Umini =  band%ngwv
  Umaxi = -band%ngwv
  do k = minval(topi), maxval(boti)
     where (k >= topi .and. k <= boti)
        Umini = min(Umini, nint(ubm(:,k)/band%dc))
        Umaxi = max(Umaxi, nint(ubm(:,k)/band%dc))
     end where
  end do

  Umini = max(Umini, -band%ngwv)
  Umaxi = min(Umaxi, band%ngwv)

  !-----------------------------------------------------------------------
  ! Gravity wave sources
  !-----------------------------------------------------------------------
  ! Start loop over all columns.
  !-----------------------------------------------------------------------
  do i=1,ncol

     !---------------------------------------------------------------------
     ! Look up spectrum only if the heating depth is large enough, else set
     ! tau0 = 0.
     !---------------------------------------------------------------------

     if (hd_idx(i) > 0) then

        !------------------------------------------------------------------
        ! Look up the spectrum using depth and uh.
        !------------------------------------------------------------------

        tau0 = desc%mfcc(hd_idx(i),nint(uh(i)),:)

        if (desc%storm_shift) then
           ! For deep convection, the wind was relative to storm cells, so
           ! shift the spectrum so that it is now relative to the ground.
           shift = -nint(CS(i)/band%dc)
           tau0 = eoshift(tau0, shift)
        end if

        ! Adjust magnitude.
        tau0 = tau0*q0(i)*q0(i)/AL

        ! Adjust for critical level filtering.
        tau0(Umini(i):Umaxi(i)) = 0.0_r8
 
        tau(i,:,topi(i)+1) = tau0

     end if ! heating depth above min and not at the pole

  enddo

  !-----------------------------------------------------------------------
  ! End loop over all columns.
  !-----------------------------------------------------------------------

  ! Output the source level.
  src_level = topi
  tend_level = topi

  ! Set phase speeds; just use reference speeds.
  c = spread(band%cref, 1, ncol)

end subroutine gw_beres_src



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

!=============================================
! Short routine to get the indices of a set of values rounded to their
! nearest points on a grid.
function index_of_nearest(x, grid) result(idx)
  real(r8), intent(in) :: x(:)
  real(r8), intent(in) :: grid(:)

  integer :: idx(size(x))

  real(r8) :: interfaces(size(grid)-1)
  integer :: i, n

  n = size(grid)
  interfaces = (grid(:n-1) + grid(2:))/2._r8

  idx = 1
  do i = 1, n-1
     where (x > interfaces(i)) idx = i + 1
  end do

end function index_of_nearest

end module gw_convect
