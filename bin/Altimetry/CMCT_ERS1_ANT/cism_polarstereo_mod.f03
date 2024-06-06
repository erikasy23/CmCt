!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!! NAME: cism_polarstereo_mod.f90
!!!
!!! PURPOSE: Sets up the standard polar stereographic coordinate system
!!! used by CISM.  Basically just initializes Jack Saba's
!!! polar_stereographic_mod.f90 for the standard CISM grid parameters and
!!! supplied grid spacing.  Currently only for Greenland, but plan to add
!!! Antarctica.  Call cism_ps_init_std_grn (which calls
!!! ps_SetAndCheckParms) first with the grid spacing (km) to initialize the
!!! projection, then call procedures from polar_stereographic_mod to make
!!! the transformations (see polar_stereographic_mod.f90 for details).
!!!
!!! Current Greenland parameters are taken from Jack Saba's
!!! bilinear_compare.pro, cism_psobj.pro, and ellipsoid.pro.
!!!
!!! USAGE:
!!!   use polar_stereographic_mod
!!!   use cism_polarstereo_mod
!!!   type( ps_params_def ) :: psparams            ! if wanted
!!!   ...
!!!   call cism_ps_init_std_grn( 4.0d0, status )   ! for 4.0 km grid spacing
!!!   ...
!!!   call latlon2pstxy( lat, lon, status, x, y )  ! x,y from lat,lon
!!!   call pstxy2latlon( x, y, lat, lon, status )  ! lat,lon from x,y
!!!   psparams = cism_ps_get_params()
!!!
!!!   X and Y are in units of the grid spacing (ie., if spacing=4.0km, then
!!!   x+1.0 is 4 km E of x); multiply by spacing to get km east and north.
!!!
!!! TBD: For Antarctica, add constant STD_ANT_CONSTS and subroutine
!!! cism_ps_init_std_ant, similar to the grn ones.  This really ought to be
!!! done with a configuration file rather than using hardcoded values,
!!! though.
!!!
!!! HISTORY: Author: Jeff Guerber, SigmaSpace/GSFC 615, June 2015.
!!!
!!! Author: Erika Simon, Innovim, 615, May 2018
!!! Added add constant STD_ANT_CONSTS
!!! Added subroutine cism_ps_init_std_ant
!!! 
!!!
!!!
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module cism_polarstereo_mod

  ! Use kinds_mod to ensure compatibility with polar_stereographic_mod
  use :: kinds_mod, only: I4b, R8b, L4b
  use :: polar_stereographic_mod

  implicit none
  private
  public :: ps_params_def
  public :: cism_ps_init_std_ant, cism_ps_get_params
  !! public :: cism_ps_init_std_grn, cism_ps_get_params


  !! PS_PARAMS_DEF: Define parameters needed for the polar stereographic
  !! projection. This is the set used by Jack Saba's
  !! polar_stereographic_mod.f90 (but with slightly different names).
  !! (Plus an "intialized" flag.)
  !!
  !! Jack describes orientation as "Longitude of the line extending from
  !! the pole downward for a N pole projection, and the longitude of the
  !! line extending from the pole upward for a S pole projection."
  !!
  !! Note that longitudes are in degrees EAST.
  !!
  !! spacing and a MUST be in the same units, which are km here, but we
  !! could as well use m (or even mi).
  !!
  !! Should probably add an indicator of whether we initialized for
  !! Greenland or Antarctica or whatever.
  type :: ps_params_def
     logical(kind=L4b) :: initialized = .false.  ! Has the PS system been init'd?
     real(kind=R8b)    :: a            ! semi-major axis, km
     real(kind=R8b)    :: e            ! eccentricity
     integer(kind=I4b) :: northsouth   ! 1=north, -1=south
     real(kind=R8b)    :: orientation  ! deg E
     real(kind=R8b)    :: spacing      ! grid spacing, km
     real(kind=R8b)    :: stdlat       ! lat of true scale, deg N
     real(kind=R8b)    :: xpole        ! X location of pole in grid coords
     real(kind=R8b)    :: ypole        ! Y location of pole in grid coords
  end type ps_params_def

  !! STD_ANT_PARAMS: Parameters for standard CISM PS grid for Antactica.
  !! Private.  Values are taken from Jack Saba's cism_psobj.pro and
  !! ellipsoid.pro.  (Structure constructor using keywords is an f2003
  !! feature).
  type( ps_params_def ), parameter :: STD_ANT_PARAMS = ps_params_def( &
       initialized = .false.,  &
       a = 6378.137000_R8b,    &   ! WGS-84, km
       e = 0.08181919104_R8b,  &   ! WGS-84
       northsouth = -1_I4b,     &   ! -1=north
       orientation = 0.0_R8b, &    ! lon 0 deg W
       spacing = -1.0_R8b,     &    ! DUMMY VALUE
       stdlat = -71.0_R8b,      &    ! deg N
       xpole = 0.0_R8b,        &    ! Pole at grid point (0,0)
       ypole = 0.0_R8b )

  !! CURRENT_PARAMS: Store the current set of parameters.
  type( ps_params_def ) :: current_params

contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! NAME: cism_ps_init_std_ant
  !!
  !! PURPOSE: Initialize polar_stereographic_mod with the standard CISM
  !! grid for Antarctica and the supplied grid spacing, by calling
  !! subroutine ps_SetAndCheckParms.  Returns the status from
  !! ps_SetAndCheckParms.
  !!
  !! DUMMY ARGUMENTS:
  !!   spacing: Grid cell spacing in km.  (Units same as a.)
  !!   status:  Status returned from ps_SetAndCheckParms.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine cism_ps_init_std_ant( spacing, status )

    real(kind=R8b), intent(in) :: spacing
    integer(kind=I4b), intent(out) :: status

    current_params = STD_ANT_PARAMS    ! NB: initialized = .false.
    current_params%spacing = spacing
    call ps_SetAndCheckParms ( &
         a_in                = current_params%a, &
         e_in                = current_params%e, &
         northsouth_in       = current_params%northsouth, &
         orientation_degE_in = current_params%orientation, &
         scale_in            = current_params%spacing, &
         stdLat_degN_in      = current_params%stdlat,  &
         xpole_in            = current_params%xpole, &
         ypole_in            = current_params%ypole, &
         status  = status )

    if (status /= 0) then
       ! This ought to go on stderr.  (Though polar_stereographic_mod.f90
       ! writes its errors to stdout.)
       print *, 'Error in cism_ps_init_std_ant: ps_SetAndCheckParams status = ', &
            status
    else
       ! We were successful
       current_params%initialized = .true.
    end if

    return
  end subroutine cism_ps_init_std_ant


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! NAME: cism_ps_get_params
  !!
  !! PURPOSE: Returns the current PS parameters, as a ps_params_def structure.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  type( ps_params_def ) function cism_ps_get_params()

    cism_ps_get_params = current_params

    return
  end function cism_ps_get_params


end module cism_polarstereo_mod
