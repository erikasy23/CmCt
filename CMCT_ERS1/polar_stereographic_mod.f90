! polar_stereographic_mod.f90: polar-stereographic coordinate transform module.
! 2000 Apr 24   Jack L. Saba
!
! SUBROUTINES INCLUDED:
!    ps_setAndCheckParms: sets parameters for a specific polar-stereographic
!                         coordinate system and makes sure they have acceptable
!                         values.
!    LatLon2PStXY: convert from lat/lon to p-s coordinates.
!    PStXY2LatLon: convert from p-s coordinates to lat/lon.
!
! NON-MODULE INTERFACE: Because module routines can only be called from Fortran
!      (90 or above) code (apparently), a set of interface routines was
!      written. These routines should NOT be called if accessing the P-S code
!      from Fortran because there is added overhead.
!
!      ps_parms: interface is identical to ps_setAndCheckParms. MUST BE CALLED
!                FIRST TO DEFINE THE PARAMETERS.
!      LL2PSt_Scalar: interface is identical to LatLon2PStXY. All arguments
!                must be scalars.
!      LL2PSt_Vector: interface is identical to LatLon2PStXY. All arguments
!                except nPts and status must be vectors.
!      LL2PSt_Tensor: interface is identical to LatLon2PStXY. All arguments
!                except nX, nY, and status must be 2-d arrays.
!      PSt2LL_Scalar: interface is identical to PStXY2LatLon. All arguments
!                must be scalars.
!      PSt2LL_Vector: interface is identical to PStXY2LatLon. All arguments
!                except nPts and status must be vectors.
!      PSt2LL_Tensor: interface is identical to PStXY2LatLon. All arguments
!                except nX, nY, and status must be 2-d arrays.
!
! USAGE:
!    CALL  ps_SetAndCheckParms ( a_in, e_in, northSouth_in, &
!                                Orientation_degE_in, Scale_in, &
!                                stdLat_degN_in, xPole_in, yPole_in, &
!                                status )
!    for scalar in/out:
!    CALL LatLon2PStXY ( Lat_degN, Lon_degE, status, x, y )
!    CALL PStXY2LatLon ( x, y, Lat_degN, Lon_degE, status )
!
!    for vector in/out:
!    CALL LatLon2PStXY ( Lat_degN, Lon_degE, nPts, status, x, y )
!    CALL PStXY2LatLon ( nPts, x, y, Lat_degN, Lon_degE, status )
!
!    for 2-d array in/out:
!    CALL LatLon2PStXY ( Lat_degN, Lon_degE, nX, nY, status, x, y )
!    CALL PStXY2LatLon ( nX, nY, x, y, Lat_degN, Lon_degE, status )
!
! NOTES:
!   1. ps_SetAndCheckParms MUST be called to set the transform paramters BEFORE
!      the transform subroutines are called.
!
!   2. The MODULE routines above CANNOT BE CALLED DIRECTLY from non-Fortran
!      programs. Instead, use the non-module interfaces described below.
!
!
! DESCRIPTION OF ARGUMENTS:
!   a_in                R8      Equatorial radius of the planetary body. Units
!                               must be same as units of Scale.
!   e_in                R8      Eccentricity of the planetary body.
!                               Dimensionless. Range [0,1).
!   northSouth_in       I4      1 ==> N pole map
!                               -1 ==> S pole map
!   Orientation_degE_in R8      Longitude of the line extending from the pole
!                               downward for a N pole projection, and the
!                               longitude of the line extending from the pole
!                               upward for a S pole projection.
!   Scale_in            R8      Scale of the map = length of one side of a cell
!                               at the standard latitude. Units must be the same
!                               as units of a. Range (0,Infinity).
!   Stdlat_Degn_in      R8      Latitude of true scale. Range -90 to +90.
!   xPole_in            R8      Location of the pole in grid coordinates
!   yPole_in            R8      Location of the pole in grid coordinates
!   Status              I4      0 = no problems
!                               -1 = bad input argument.
!   nPts                I4      Number of points in lat, lon, x and y vectors.
!   nX, nY              I4      Dimensions of lat, lon, x and y arrays, (nX,nY).
!   Lat_degN            R8      Latitude(s) of the point(s). This may be a
!                               scalar, vector, or 2-d array for calls to
!                               LatLon2PStXYand PStXY2LatLon. Use a variable of
!                               the appropriate shape for calls to the external
!                               subroutines.
!   Lon_degE            R8      Longitude(s) of the point(s). This may be a
!                               scalar, vector, or 2-d array for calls to
!                               LatLon2PStXYand PStXY2LatLon. Use a variable of
!                               the appropriate shape for calls to the external
!                               subroutines.
!   xcoord            R8        x coordinate(s) of the point(s). This may be a
!                               scalar, vector, or 2-d array for calls to
!                               LatLon2PStXYand PStXY2LatLon. Use a variable of
!                               the appropriate shape for calls to the external
!                               subroutines.
!   ycoord            R8        y coordinate(s) of the point(s). This may be a
!                               scalar, vector, or 2-d array for calls to
!                               LatLon2PStXYand PStXY2LatLon. Use a variable of
!                               the appropriate shape for calls to the external
!                               subroutines.
!                               subroutines.
!
! ERROR HANDLING:
!    If an error occurs, this module prints a message and returns a nonzero
!    status flag.
!
! VALUES OF PARAMETERS FOR PREVIOUS COORDINATE SYSTEMS
!
! Zwally radar altimetry grids
!                                       Grid Size
!                     ----------------------------------------------
!              Units      5 km            10 km           20 km
!              -----  -------------   -------------   --------------
!  Location           Gr        Ant   Gr        Ant   Gr         Ant
!  a           km        6378.273        6378.273        6378.273
!  e           none         0               0               0
!  northSouth  none   1         -1    1         -1     1         -1
!  orientation degE   315        0    315        0     315        0
!  scale       km       5.2286257       10.452171       20.955143
!  std lat     degN   90       -90    90        -90    90        -90
!  x pole      none        889             445             223
!  y pole      none        889             445             223
!  perim lat   degN*  50       -50    50        -50    50        -50
!  dscale      none*   2439.75121      1220.468499      608.754894
!  orientation degE*  45       270    45        270    45     270
!
! *: Parameters used in the old radar altimetry grids for completeness.
!    The orientation is defined differently but means the same thing (i.e.,
!    45 deg W is the downward vertical for Greenland, and 0 deg is the upward
!    vertical for Antarctica).
!
! SSMI 25-km grids
!              Units Greenland     Antarctica
!              ----- ------------  ----------
!  a           km    6378.273      6378.273
!  e           none  0.081816153   0.081816153
!  northSouth  none  1             -1
!  orientation degE  315           0
!  scale       km    25.0          25.0
!  std lat     degN  70            -70
!  x pole      none  154           158
!  y pole      none  234           174
!
!-------------------------------------------------------------------------------

MODULE Polar_Stereographic_MOD

   USE kinds_MOD, ONLY: I4b, R8b
   IMPLICIT NONE

   PRIVATE
   PUBLIC :: LatLon2PStXY, LatLon2PStXY_Scaler, LatLon2PStXY_Tensor, &
             LatLon2PStXY_Vector
   PUBLIC :: PStXY2LatLon, PStXY2LatLon_Scaler, PStXY2LatLon_Tensor, &
             PStXY2LatLon_Vector
   PUBLIC :: PS_SetAndCheckParms

!-------------------------------------------------------------------------------
!                              MODULE INTERFACES

!...Defined to allow input of vector arguments

   INTERFACE LatLon2PStXY
      MODULE PROCEDURE LatLon2PStXY_Scaler, LatLon2PStXY_Tensor, &
                       LatLon2PStXY_Vector
   END INTERFACE

   INTERFACE PStXY2LatLon
      MODULE PROCEDURE PStXY2LatLon_Scaler, PStXY2LatLon_Tensor, &
                       PStXY2LatLon_Vector
   END INTERFACE

!-------------------------------------------------------------------------------
!                              MODULE VARIABLES

!...Parameters

   REAL (KIND=R8b), PARAMETER  :: Epsilon = 100.0d0*TINY(1.0)
   REAL (KIND=R8b), PARAMETER  :: pi      = 3.1415926535897932d0
   REAL (KIND=R8b), PARAMETER  :: piOn2   = pi / 2.0d0
   REAL (KIND=R8b), PARAMETER  :: piOn4   = pi / 4.0d0
   REAL (KIND=R8b), PARAMETER  :: twoPi   = 2.0d0 * pi
   REAL (KIND=R8b), PARAMETER  :: DTOR    = pi / 180.0d0

!...Input arguments = parameters of the fit. These are initialized by a call to
!...ps_SetAndCheckParms and used in calls to LatLon2PStXY and PStXY2LatLon.
!...Where possible, they are initialized here to values that will cause trouble
!...if ps_SetAndCheckParms has not been called.

   ! orientation_degE: Longitude of the line extending from the pole downward
   !                   for a N pole projection, and the longitude of the line
   !                   extending from the pole upward for a S pole
   !                   projection. Range: anything, but will be converted to
   !                   [-180,180].
   ! Scale:            Scale of the map = length of one side of a cell at the
   !                   standard latitude. Units must be the same as units of
   !                   a. Range: Scale > 0.
   ! stdLat_degN:      Latitude of true scale. Range -90 to +90.

   ! NOTE: The three parameters above are not used directly except to set global
   !       parameters. Therefore they are not defined here.

   REAL (KIND=R8b), SAVE ::    a     = 0.0d0  ! Equatorial radius of the
                                                  !  planetary body. Units must
                                                  !  be same as units of Scale.
                                                  !  Range: a > 0.
   REAL (KIND=R8b), SAVE ::    e     = -1.0d0 ! Eccentricity of the
                                                  !  planetary body.
                                                  !  Dimensionless. Range [0,1).
   INTEGER (KIND=I4b), SAVE :: northSouth = 0 ! 1 = N pole map.
                                                  ! -1 = S pole map.
   REAL (KIND=R8b), SAVE ::    xPole = 0.0d0  ! Location of the pole in
   REAL (KIND=R8b), SAVE ::    yPole = 0.0d0  !  grid coordinates.

!...Module variables based on the input arguments, calculated in
!...ps_setAndCheckParms, used in LatLon2PStXY and/or PStXY2LatLon.

   REAL (KIND=R8b), SAVE :: D         = 0.0d0 ! Scaling factor used in
                                                  !  LatLon2PStXY when e=0.
   REAL (KIND=R8b), SAVE :: k0        = 0.0d0 ! Inverse of Scale. Used for
                                                  !  consistency with equations
                                                  !  in Snyder.
   REAL (KIND=R8b), SAVE :: Orientation_rad   ! Orientation of map converted
                                                  !  to radians in [-180,180].
   REAL (KIND=R8b), SAVE :: rhoFactor = 0.0d0 ! Constant factor used in
                                                  !  calculation of rho in
                                                  !  LatLon2PStXY
   REAL (KIND=R8b), SAVE :: rhoStd    = 0.0d0 ! Distance from pole to
                                                  !  circle at standard lat.
   REAL (KIND=R8b), SAVE :: term2Chi  = 0.0d0 ! Terms used in series
   REAL (KIND=R8b), SAVE :: term4Chi  = 0.0d0 !  expansion for latitude
   REAL (KIND=R8b), SAVE :: term6Chi  = 0.0d0 !  in PStXY2LatLon.
   REAL (KIND=R8b), SAVE :: tFactor   = 0.0d0 ! Factor used to calculate t.

CONTAINS

!===============================================================================

   SUBROUTINE ps_SetAndCheckParms ( a_in, e_in, northSouth_in, &
                                    Orientation_degE_in, Scale_in, &
                                    StdLat_degN_in, xPole_in, yPole_in, &
                                    status )

      IMPLICIT NONE
      ! Input arguments:
      REAL (KIND=R8b),    INTENT(IN) :: a_in    !...Equatorial radius
                                                    !...of the planetary body.
                                                    !...Units must be same as
                                                    !...units of Scale.
      REAL (KIND=R8b),    INTENT(IN) :: e_in    !...Eccentricity of
                                                    !...the planetary body.
                                                    !...Dimensionless. Range
                                                    !...[0,1).
      INTEGER (KIND=I4b), INTENT(IN) :: northSouth_in !...1 ==> N pole map
                                                          !...-1 ==> S pole map
      REAL (KIND=R8b),    INTENT(IN) :: Orientation_degE_in !...Longitude of
                                                    !...the line extending from
                                                    !...the pole downward for a
                                                    !...N pole projection, and
                                                    !...the longitude of the
                                                    !...line extending from the
                                                    !...pole upward for a S
                                                    !...pole projection.
      REAL (KIND=R8b),    INTENT(IN) :: Scale_in !...Scale of the
                                                    !...map = length of one
                                                    !...side of a cell at the
                                                    !...standard latitude.
                                                    !...Units must be the same
                                                    !...as units of a. Range
                                                    !...(0,Infinity).
      REAL (KIND=R8b),    INTENT(IN) :: Stdlat_Degn_in
                                                    !...Latitude of true scale.
                                                    !...Range -90 to +90.
      REAL (KIND=R8b),    INTENT(IN) :: xPole_in !...Location of the pole in
      REAL (KIND=R8b),    INTENT(IN) :: yPole_in !...grid coordinates
      !...Output arguments
      INTEGER (KIND=I4b), INTENT(OUT) :: Status !...0 = no problems
                                                    !...-1 = bad input argument.
      !...Local variables
      REAL (KIND=R8b) :: absStdLat_radN
      REAL (KIND=R8b) :: cosStd
      REAL (KIND=R8b) :: efactor      ! Factor in calculation of rho.
      REAL (KIND=R8b) :: esq          ! Square of eccentricity
      REAL (KIND=R8b) :: e4           ! Powers of
      REAL (KIND=R8b) :: e6           !    eccentricity
      REAL (KIND=R8b) :: mc           ! Factor in equation for rho when
                                          !  std lat = +/- 90 degs.
      REAL (KIND=R8b) :: sinStd
      REAL (KIND=R8b) :: tc           ! Factor in equation for rho when
                                          !  std lat = +/- 90 degs.

!------------------------------------------------------------------------------

!...Status = 0 ==> good data

      status = 0

!...Test values of parameters to make sure they are in bounds

      IF ( ABS(northSouth_in) /= 1 ) THEN
         PRINT *, 'PS_SETANDCHECKPARMS: northSouth must be +1 or -1'
         PRINT *, '   input northSouth: ', northSouth_in
         Status = -1
      ELSE
         northSouth = northSouth_in
      END IF

      IF ( a_in <= 0.0d0 ) THEN
         PRINT *, 'PS_SETANDCHECKPARMS: Planetary radius (a) <= 0'
         PRINT *, '   input a: ', a_in
         status = -1
      ELSE
         a = a_in
      END IF

      IF ( e_in < 0.0d0 .OR. e_in >= 1.0d0 ) THEN
         PRINT *, 'PS_SETANDCHECKPARMS: eccentricity must be in [0,1)'
         PRINT *, '   input e: ', e_in
         status = -1
      ELSE
         e = e_in
      END IF

      IF ( Scale_in <= 0.0d0 ) THEN
         PRINT *, 'PS_SETANDCHECKPARMS: map scale must be > 0'
         PRINT *, '   input Scale: ', Scale_in
         status = -1
      END IF

      IF ( northSouth*StdLat_degN_in < 0.0d0 &
      .OR.   northSouth*StdLat_degN_in > 90.0d0+Epsilon ) THEN
         PRINT *, 'PS_SETANDCHECKPARMS: standard latitude must be in proper' &
                  // ' hemisphere and <= 90 degs'
         PRINT *, '   input northSouth, StdLat_degN: ', northSouth, &
                  StdLat_degN_in
         status = -1
      END IF

      IF ( xPole_in < 0 ) THEN
         PRINT *, 'PS_SETANDCHECKPARMS: Pole position must be >= 0'
         PRINT *, '   input xPole: ', xPole_in
         status = -1
      ELSE
         xPole = xPole_in
      END IF

      IF ( yPole_in < 0 ) THEN
         PRINT *, 'PS_SETANDCHECKPARMS: Pole position must be >= 0'
         PRINT *, '   input yPole: ', yPole_in
         status = -1
      ELSE
         yPole = yPole_in
      END IF

!...Constants based on the input parameters. these are defined here so they
!...do not have to be calculated on each call to LatLon2PStXY or PStXY2LatLon.

      esq = e**2

      !...term... are factors used in expansion for latitude are constant.
      !...Calculate here for use in PStXY2LatLon

      e4       = esq**2
      e6       = esq**3
      term2chi = esq/2 + (5.0d0/24.0d0) * e4 + (1.0d0/12.0d0)*e6
      term4chi = (7.0d0/48.0d0)*e4 + (29.0d0/240.0d0)*e6
      term6chi = (7.0d0/120.0d0)*e6

      !...Define k0 for consistency with equations in Snyder
      k0 = 1.0d0 / Scale_in

      !...Convert angle to radians.  Make sure range of longitude is [-180,180].
      !...Compute trig fcns.

      IF ( Orientation_degE_in < 180.0d0 ) THEN
         Orientation_rad = DTOR * Orientation_degE_in
      ELSE
         Orientation_rad = DTOR * ( Orientation_degE_in - 360.0d0 )
      END IF
      absStdLat_radN = DTOR * DBLE(ABS(StdLat_degN_in))
      sinStd         = SIN ( absStdLat_radN )
      cosStd         = COS ( absStdLat_radN )

      !...Distance from pole to standard latitude in grid units.

      rhoStd = a * k0 * ( 1.0d0 + sinStd )

      !...The distance from pole to point is scaled differently for std lat not
      !...+/- 90 vs = +/- 90. Here we calculate the constant portion of the
      !...equation.

      IF ( ABS(StdLat_degN_in) < (90.0d0-Epsilon) ) THEN
         !...std lat is NOT at the pole.
         mc        = cosStd / SQRT ( 1.0d0 - esq * sinStd**2 )
         efactor   = ( (1.0d0+e*sinStd) / (1.0d0-e*sinStd) ) ** (e/2.0d0)
         tc        = TAN ( piOn4 - absStdLat_radN/2.0d0 ) * efactor
         rhoFactor = a * k0 * mc / tc
         tFactor   = tc / ( a * k0 * mc )
      ELSE
         !...Std lat at pole
         efactor   = SQRT ( (1.0d0+e)**(1.0d0+e) * (1.0d0-e)**(1.0d0-e) )
         rhoFactor = 2.0d0 * a * k0 / efactor
         tFactor   = efactor / ( 2.0d0 * a * k0 )
      END IF

      !...Scale factor is 2ak_o if std lat is 90 degs. Factor 2 is reduced
      !...otherwise. The general form from which the 2 is derived is
      !...1+sin(std lat). This is used only for spherical case (e=0).

      D = ( 1.0d0 + SIN(absStdLat_radN) ) * a * k0

!      PRINT *, 'PS_SETANDCHECKPARMS INPUT ARGUMENTS:'
!      PRINT *, '   a_in, e_in:          ', a_in, e_in
!      PRINT *, '   northSouth_in:       ', northSouth_in
!      PRINT *, '   Orientation_degE_in: ', Orientation_degE_in
!      PRINT *, '   Scale_in:            ', Scale_in
!      PRINT *, '   StdLat_degN_in:      ', StdLat_degN_in
!      PRINT *, '   xPole_in, yPole_in:  ', xPole_in, yPole_in
!      PRINT *, 'PS_SETANDCHECKPARMS LOCAL VARIABLES:'
!      PRINT *, '   absStdLat_radN:      ', absStdLat_radN
!      PRINT *, '   cosStd, sinStd:      ', cosStd, sinStd
!      PRINT *, '   efactor:             ', efactor
!      PRINT *, '   esq, e4, e6:         ', esq, e4, e6
!      PRINT *, '   mc, tc:              ', mc, tc
!      PRINT *, 'PS_SETANDCHECKPARMS MODULE VARIABLES:'
!      PRINT *, '   a, D, e:             ', a, D, e
!      PRINT *, '   k0:                  ', k0
!      PRINT *, '   northSouth:          ', northSouth
!      PRINT *, '   orientation_rad:     ', orientation_rad
!      PRINT *, '   rhoFactor, rhoStd:   ', rhoFactor, rhoStd
!      PRINT *, '   term (2,4,6) Chi:    ', term2Chi, term4Chi, term6Chi
!      PRINT *, '   tFactor:             ', tFactor
!      PRINT *, '   xPole, yPole:        ', xPole, yPole
!      PRINT *, 'PS_SETANDCHECKPARMS OUTPUT ARGUMENTS:'
!      PRINT *, '   status:              ', status

      RETURN

   END SUBROUTINE ps_SetAndCheckParms

!===============================================================================

   SUBROUTINE LatLon2PStXY_Scaler ( Lat_degN, Lon_degE, status, xcoord, ycoord )

      IMPLICIT NONE
   ! Input arguments
      REAL (KIND=R8b),    INTENT(IN)  :: Lat_degN
      REAL (KIND=R8b),    INTENT(IN)  :: Lon_degE
   ! Output arguments
      INTEGER (KIND=I4b), INTENT(OUT) :: status
      REAL (KIND=R8b),    INTENT(OUT) :: xcoord
      REAL (KIND=R8b),    INTENT(OUT) :: ycoord
   ! Local variabies
      REAL (KIND=R8b) :: absLat_radN
      REAL (KIND=R8b) :: dlon_rad     ! Angle between orientation angle and
                                          !  longitude of point.
      REAL (KIND=R8b) :: Lon          ! Longitude of current point.
      REAL (KIND=R8b) :: rho          ! Distance from pole to point in
                                          !  grid cell units.
      REAL (KIND=R8b) :: sinLat
      REAL (KIND=R8b) :: t            ! Factor in equation for rho.
      REAL (KIND=R8b) :: x            ! x coordinate relative to the pole.
      REAL (KIND=R8b) :: y            ! y coordinate relative to the pole.

!----------------------------------------------------------------------------

!...Check range of lat. Lon can have any value.

      IF ( ABS(Lat_degN) > 90.0d0+Epsilon ) THEN
         PRINT *, 'LatLon2PStXY_Scaler: Latitude out of bounds'
         PRINT *, '   Latitude: ', Lat_degN
         status = -1
         RETURN
      END IF

      status = 0

!...Special case 1: point at pole.

      IF ( ABS(Lat_degN) > (90.0d0-Epsilon) ) THEN
         xcoord = xPole
         ycoord = yPole
         GO TO 9999
      END IF

!...Convert angles to radians; make sure longitude is in [-180,180].

      absLat_radN = DTOR * ABS(Lat_degN)
      Lon = Lon_degE
      IF ( Lon < -180.0d0 ) THEN
         DO WHILE ( Lon < -180.0d0 )
            Lon = Lon + 360.0d0
         END DO
      END IF
      Lon = MOD ( Lon, 360.0d0 )

      dLon_rad = northSouth * ( DTOR * Lon - Orientation_rad )

!...Special case 2: spherical planet (e=0)

      IF ( e < Epsilon ) THEN
         rho    = D * TAN ( piOn4 - absLat_radN/2.0d0 )
         x      = rho * SIN(dLon_rad)
         y      = - rho * COS(dLon_rad)
         xcoord = xPole + northSouth * x
         ycoord = yPole - northSouth * y
         GO TO 9999
      END IF

!...Functions of angles needed for all further calculations

      sinLat = SIN ( absLat_radN )

!...Calculate distance from pole to point and break it into components

      t   = TAN ( piOn4 - absLat_radN/2.0d0 ) &
            * ( (1.0d0+e*sinLat) / (1.0d0-e*sinLat) ) ** (e/2.0d0)
      rho = rhoFactor * t
      x   =   rho * SIN ( dLon_rad )
      y   = - rho * COS ( dLon_rad )

!...Convert to upper left corner as origin.

      xcoord = xPole + northSouth * x
      ycoord = yPole - northSouth * y

 9999 CONTINUE  ! Jump here after completing calcs in different branches.

!      PRINT *, 'LATLON2PSTXY INPUT ARGUMENTS:'
!      PRINT *, '   Lat degN, Lon degE:  ', Lat_degN, Lon_degE
!      PRINT *, '   northSouth:          ', northSouth
!      PRINT *, '   a, e:                ', a, e
!      PRINT *, '   xPole, yPole:        ', xPole, yPole
!      PRINT *, 'LATLON2PSTXY LOCAL VARIABLES:'
!      PRINT *, '   abs lat, dlon (rad): ', absLat_radN, dlon_rad
!      PRINT *, '   sin lat:             ', sinlat
!      PRINT *, '   t, rho:              ', t, rho
!      PRINT *, '   x, y:                ', x, y
!      PRINT *, 'LATLON2PSTXY OUTPUT VARIABLES:'
!      PRINT *, '   status:              ', status
!      PRINT *, '   xcoord, ycoord:      ', xcoord, ycoord

     RETURN

   END SUBROUTINE LatLon2PStXY_Scaler

!===============================================================================

   SUBROUTINE LatLon2PStXY_Tensor ( Lat_degN, Lon_degE, nX, nY, status, &
                                    xcoord, ycoord )

      IMPLICIT NONE
   ! Input arguments
      INTEGER (KIND=I4b), INTENT(IN)  :: nX
      INTEGER (KIND=I4b), INTENT(IN)  :: nY
      REAL (KIND=R8b),    INTENT(IN)  :: Lat_degN(nX,nY)
      REAL (KIND=R8b),    INTENT(IN)  :: Lon_degE(nX,nY)
   ! Output arguments
      INTEGER (KIND=I4b), INTENT(OUT) :: status
      REAL (KIND=R8b),    INTENT(OUT) :: xcoord(nX,nY)
      REAL (KIND=R8b),    INTENT(OUT) :: ycoord(nX,nY)
   ! Local variabies
      INTEGER (KIND=I4b) :: ix
      INTEGER (KIND=I4b) :: iy

      DO iy=1,Ny
         DO ix=1,Nx
            CALL LatLon2PStXY_Scaler ( Lat_degN(ix,iy), Lon_degE(ix,iy), &
                                       status, xcoord(ix,iy), ycoord(ix,iy) )
!            IF ( Status /= 0 ) RETURN
         END DO
      END DO

      RETURN

   END SUBROUTINE LatLon2PStXY_Tensor

!===============================================================================

   SUBROUTINE LatLon2PStXY_Vector ( Lat_degN, Lon_degE, nPts, status, xcoord, &
                                    ycoord )

      IMPLICIT NONE
   ! Input arguments
      INTEGER (KIND=I4b), INTENT(IN) :: nPts
      REAL (KIND=R8b),    INTENT(IN)  :: Lat_degN(nPts)
      REAL (KIND=R8b),    INTENT(IN)  :: Lon_degE(nPts)
   ! Output arguments
      INTEGER (KIND=I4b), INTENT(OUT) :: status
      REAL (KIND=R8b),    INTENT(OUT) :: xcoord(nPts)
      REAL (KIND=R8b),    INTENT(OUT) :: ycoord(nPts)
   ! Local variabies
      INTEGER (KIND=I4b) :: i

      DO i=1,nPts
         CALL LatLon2PStXY_Scaler ( Lat_degN(i), Lon_degE(i), status, &
                                    xcoord(i), ycoord(i) )
         IF ( Status /= 0 ) RETURN
      END DO

      RETURN

   END SUBROUTINE LatLon2PStXY_Vector

!===============================================================================

   SUBROUTINE PStXY2LatLon_Scaler ( xcoord, ycoord, Lat_degN, Lon_degE, status )

      IMPLICIT  NONE
   !...Input arguments
      REAL (KIND=R8b),    INTENT(IN)  :: xcoord
      REAL (KIND=R8b),    INTENT(IN)  :: ycoord
   !...Output arguments
      REAL (KIND=R8b),    INTENT(OUT) :: Lat_degN
      REAL (KIND=R8b),    INTENT(OUT) :: Lon_degE
      INTEGER (KIND=I4b), INTENT(OUT) :: status

   !...Local variables
      REAL (KIND=R8b) :: chi
      REAL (KIND=R8b) :: Lat_radN
      REAL (KIND=R8b) :: Lon_radE
      REAL (KIND=R8b) :: lon
      REAL (KIND=R8b) :: phi
      REAL (KIND=R8b) :: rho
      REAL (KIND=R8b) :: t
      REAL (KIND=R8b) :: x
      REAL (KIND=R8b) :: y

!----------------------------------------------------------------------------

      status = 0

!...Special case 1: the pole.

      IF (  ( ABS(xcoord-xPole) < Epsilon ) &
      .AND. ( ABS(ycoord-yPole) < Epsilon ) ) THEN
         Lon_degE = 0.0d0
         Lat_degN = northSouth * 90.0d0
         GO TO 9999
      END IF

!...Coordinates relative to pole in unscaled units. Reverse signs of
!...appropriate input arguments if S pole to use same equations

      x = xcoord - xPole
      y = yPole - ycoord

!...Distance from pole to point in the projection plane in grid units.

      rho = SQRT ( x**2 + y**2 )

!...Special case 2: spherical planet.

      IF ( e < Epsilon ) THEN
         phi      = piOn2 - 2.0d0 * ATAN2 ( rho, rhoStd )
         Lat_degN = northSouth * phi / DTOR
         lon      = ATAN2 ( northSouth*x, -northSouth*y )
         Lon_radE = Orientation_rad + northSouth * Lon
         Lon_degE = Lon_radE / DTOR
         GO TO 9999
      END IF

!...The general case. To avoid iterations, the latitude is initially
!...defined using an approximation

      t   = rho * tFactor
      chi = piOn2 - 2.0d0 * ATAN ( t )

!...and a series expansion is used to calculate the final value

      Lat_radN = chi + term2chi * SIN(2.0d0*chi) &
                 + term4chi * SIN(4.0d0*chi) + term6chi * SIN(6.0d0*chi)

!...The longitude is simple

      Lon = ATAN2 ( northSouth*x, -northSouth*y )

      Lon_radE = Orientation_rad + northSouth * Lon

!...Convert angles from radians to degrees

      Lat_degN = northSouth * Lat_radN / DTOR
      Lon_degE = Lon_radE / DTOR

 9999 CONTINUE  ! Jump here after completing calcs in different branches.

!...Make sure longitudes are in bounds

      DO WHILE ( Lon_degE < 0.0d0 )
         Lon_degE = Lon_degE + 360.0d0
      END DO
      Lon_degE = MOD ( Lon_degE, 360.0d0 )

!      PRINT *, 'PSTXY2LATLON INPUT ARGUMENTS:'
!      PRINT *, '   xcoord, ycoord:     ', xcoord, ycoord
!      PRINT *, 'PSTXY2LATLON INTERNAL VARIABLES:'
!      PRINT *, '   chi:                ', chi
!      PRINT *, '   Lat_radN:           ', Lat_radN
!      PRINT *, '   Lon_radE:           ', Lon_radE
!      PRINT *, '   lon, phi:           ', lon, phi
!      PRINT *, '   rho, t:             ', rho, t
!      PRINT *, '   x, y:               ', x, y
!      PRINT *, 'PSTXY2LATLON OUTPUT ARGUMENTS:'
!      PRINT *, '   lat degN, lon degE: ', lat_degN, Lon_degE
!      PRINT *, '   status:             ', status

      RETURN

   END SUBROUTINE PStXY2LatLon_Scaler

!===============================================================================

   SUBROUTINE PStXY2LatLon_Tensor ( nX, nY, xcoord, ycoord, &
                                    Lat_degN, Lon_degE, status )

      IMPLICIT NONE
   ! Input arguments
      INTEGER (KIND=I4b), INTENT(IN) :: nX
      INTEGER (KIND=I4b), INTENT(IN) :: nY
      REAL (KIND=R8b),    INTENT(IN) :: xcoord(nX,nY)
      REAL (KIND=R8b),    INTENT(IN) :: ycoord(nX,nY)
   ! Output arguments
      REAL (KIND=R8b),    INTENT(OUT) :: Lat_degN(nX,nY)
      REAL (KIND=R8b),    INTENT(OUT) :: Lon_degE(nX,nY)
      INTEGER (KIND=I4b), INTENT(OUT) :: status
   ! Local variabies
      INTEGER (KIND=I4b) :: ix
      INTEGER (KIND=I4b) :: iy

      DO iy=1,Ny
         DO ix=1,Nx
            CALL PStXY2LatLon_Scaler ( xcoord(ix,iy), ycoord(ix,iy), &
                                       Lat_degN(ix,iy), Lon_degE(ix,iy), &
                                       status )
!            IF ( Status /= 0 ) RETURN
         END DO
      END DO

      RETURN

   END SUBROUTINE PStXY2LatLon_Tensor

!===============================================================================

   SUBROUTINE PStXY2LatLon_Vector ( nPts, xcoord, ycoord, &
                                    Lat_degN, Lon_degE, status )

      IMPLICIT NONE
   ! Input arguments
      INTEGER (KIND=I4b), INTENT(IN)  :: nPts
      REAL (KIND=R8b),    INTENT(IN)  :: xcoord(nPts)
            REAL (KIND=R8b),    INTENT(IN)  :: ycoord(nPts)
   ! Output arguments
      REAL (KIND=R8b),    INTENT(OUT) :: Lat_degN(nPts)
      REAL (KIND=R8b),    INTENT(OUT) :: Lon_degE(nPts)
      INTEGER (KIND=I4b), INTENT(OUT) :: status
   ! Local variabies
      INTEGER (KIND=I4b) :: i

      DO i=1,nPts
         CALL PStXY2LatLon_Scaler ( xcoord(i), ycoord(i), &
                                    Lat_degN(i), Lon_degE(i), status )
!         IF ( Status /= 0 ) RETURN
      END DO

      RETURN

   END SUBROUTINE PStXY2LatLon_Vector

!===============================================================================

END MODULE Polar_Stereographic_MOD

!===============================================================================
!===============================================================================

! The following code is needed for interfacing the module above to non-Fortran
! programs, since these programs can't access Fortran 90 modules directly as
! far as I can tell.

SUBROUTINE ps_parms ( a_in, e_in, northSouth_in, Orientation_degE_in, &
                      Scale_in, stdLat_degN_in, xPole_in, yPole_in, status )

   USE kinds_MOD, ONLY: I4b, R8b
   USE Polar_Stereographic_MOD
   IMPLICIT    NONE
   ! Input arguments: SEE ps_SetAndCheckParms for descriptions.
   REAL (KIND=R8b),    INTENT(IN) :: a_in
   REAL (KIND=R8b),    INTENT(IN) :: e_in
   INTEGER (KIND=I4b), INTENT(IN) :: northSouth_in
   REAL (KIND=R8b),    INTENT(IN) :: Orientation_degE_in
   REAL (KIND=R8b),    INTENT(IN) :: Scale_in
   REAL (KIND=R8b),    INTENT(IN) :: Stdlat_Degn_in
   REAL (KIND=R8b),    INTENT(IN) :: xPole_in
   REAL (KIND=R8b),    INTENT(IN) :: yPole_in
   !...Output arguments
   INTEGER (KIND=I4b), INTENT(OUT) :: Status !...0 = no problems
                                                    !...-1 = bad input argument.

   CALL  ps_SetAndCheckParms ( a_in, e_in, northSouth_in, Orientation_degE_in, &
                               Scale_in, stdLat_degN_in, xPole_in, yPole_in, &
                               status )
   RETURN
END SUBROUTINE ps_parms

!===============================================================================

SUBROUTINE LL2PSt_Scaler ( Lat_degN, Lon_degE, status, xcoord, ycoord )

   USE kinds_MOD, ONLY: I4b, R8b
   USE Polar_Stereographic_MOD
   IMPLICIT   NONE
!...Input arguments
   REAL (KIND=R8b),    INTENT(IN)  :: Lat_degN
   REAL (KIND=R8b),    INTENT(IN)  :: Lon_degE
!...Output arguments
   INTEGER (KIND=I4b), INTENT(OUT) :: status
   REAL (KIND=R8b),    INTENT(OUT) :: xcoord
   REAL (KIND=R8b),    INTENT(OUT) :: ycoord

   CALL LatLon2PStXY_Scaler ( Lat_degN, Lon_degE, status, xcoord, ycoord )

   RETURN
END SUBROUTINE LL2PSt_Scaler

!===============================================================================

SUBROUTINE LL2PSt_Vector ( Lat_degN, Lon_degE, nPts, status, xcoord, ycoord )

   USE kinds_MOD, ONLY: I4b, R8b
   USE Polar_Stereographic_MOD
   IMPLICIT   NONE
!...Input arguments
   INTEGER (KIND=I4b), INTENT(IN)  :: nPts
   REAL (KIND=R8b),    INTENT(IN)  :: Lat_degN(nPts)
   REAL (KIND=R8b),    INTENT(IN)  :: Lon_degE(nPts)
!...Output arguments
   INTEGER (KIND=I4b), INTENT(OUT) :: status
   REAL (KIND=R8b),    INTENT(OUT) :: xcoord(nPts)
   REAL (KIND=R8b),    INTENT(OUT) :: ycoord(nPts)

   CALL LatLon2PStXY_Vector ( Lat_degN, Lon_degE, nPts, status, xcoord, ycoord )

   RETURN
END SUBROUTINE LL2PSt_Vector

!===============================================================================

SUBROUTINE LL2PSt_Tensor ( Lat_degN, Lon_degE, nX, nY, status, xcoord, ycoord )

   USE kinds_MOD, ONLY: I4b, R8b
   USE Polar_Stereographic_MOD
   IMPLICIT   NONE
!...Input arguments
   INTEGER (KIND=I4b), INTENT(IN)  :: nX, nY
   REAL (KIND=R8b),    INTENT(IN)  :: Lat_degN(nX,nY)
   REAL (KIND=R8b),    INTENT(IN)  :: Lon_degE(nX,nY)
!...Output arguments
   INTEGER (KIND=I4b), INTENT(OUT) :: status
   REAL (KIND=R8b),    INTENT(OUT) :: xcoord(nX,nY)
   REAL (KIND=R8b),    INTENT(OUT) :: ycoord(nX,nY)

   CALL LatLon2PStXY_Tensor ( Lat_degN, Lon_degE, nX, nY, &
                              status, xcoord, ycoord )

   RETURN
END SUBROUTINE LL2PSt_Tensor

!===============================================================================

SUBROUTINE PSt2LL_Scaler ( xcoord, ycoord, Lat_degN, Lon_degE, status )

   USE kinds_MOD, ONLY: I4b, R8b
   USE Polar_Stereographic_MOD
   IMPLICIT   NONE
!...Input arguments
   REAL (KIND=R8b),    INTENT(IN)  :: xcoord
   REAL (KIND=R8b),    INTENT(IN)  :: ycoord
!...Output arguments
   REAL (KIND=R8b),    INTENT(OUT) :: Lat_degN
   REAL (KIND=R8b),    INTENT(OUT) :: Lon_degE
   INTEGER (KIND=I4b), INTENT(OUT) :: status

   CALL PStXY2LatLon_Scaler ( xcoord, ycoord, Lat_degN, Lon_degE, status )

   RETURN
END SUBROUTINE PSt2LL_Scaler

!===============================================================================

SUBROUTINE PSt2LL_Vector ( nPts, xcoord, ycoord, Lat_degN, Lon_degE, status )

   USE kinds_MOD, ONLY: I4b, R8b
   USE Polar_Stereographic_MOD
   IMPLICIT   NONE
!...Input arguments
   INTEGER (KIND=I4b), INTENT(IN)  :: nPts
   REAL (KIND=R8b),    INTENT(IN)  :: xcoord(nPts)
   REAL (KIND=R8b),    INTENT(IN)  :: ycoord(nPts)
!...Output arguments
   REAL (KIND=R8b),    INTENT(OUT) :: Lat_degN(nPts)
   REAL (KIND=R8b),    INTENT(OUT) :: Lon_degE(nPts)
   INTEGER (KIND=I4b), INTENT(OUT) :: status

   CALL PStXY2LatLon_Vector ( nPts, xcoord, ycoord, Lat_degN, Lon_degE, status )

   RETURN
END SUBROUTINE PSt2LL_Vector

!===============================================================================

SUBROUTINE PSt2LL_Tensor ( nX, nY, xcoord, ycoord, Lat_degN, Lon_degE, status )

   USE kinds_MOD, ONLY: I4b, R8b
   USE Polar_Stereographic_MOD
   IMPLICIT   NONE
!...Input arguments
   INTEGER (KIND=I4b), INTENT(IN)  :: nX, nY
   REAL (KIND=R8b),    INTENT(IN)  :: xcoord(nX,nY)
   REAL (KIND=R8b),    INTENT(IN)  :: ycoord(nX,nY)
!...Output arguments
   REAL (KIND=R8b),    INTENT(OUT) :: Lat_degN(nX,nY)
   REAL (KIND=R8b),    INTENT(OUT) :: Lon_degE(nX,nY)
   INTEGER (KIND=I4b), INTENT(OUT) :: status

   CALL PStXY2LatLon_Tensor ( nX, nY, xcoord, ycoord, &
                              Lat_degN, Lon_degE, status )

   RETURN
END SUBROUTINE PSt2LL_Tensor

!===============================================================================
