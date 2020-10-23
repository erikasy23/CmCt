!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! NAME: model_cism_txt_mod.f08
!!!
!!! PURPOSE: Module for using CISM model files in TEXT format, such as the
!!! "cism_usrf_yr_*.txt" files from Steve Price and which Jack Saba's
!!! bilinear_compare.pro reads.  usrf = upper surface, ie. elevation.
!!!
!!! These files have 6 columns, ASCII, separated by commas.  In Steve's
!!! cism_usrf_*.txt files, the single header line names the columns
!!! latitude, longitude, surf_elev_wgs84 (which Jack's program calls Elev),
!!! cism_x, cism_y, surf_elev_geoid (which Jack calls CISM_Elev).  cism_x
!!! and cism_y use a different grid spacing, apparently in meters with
!!! different origin, but this program doesn't use them, instead it
!!! converts lat and lon to our own PS grid.  We currently only use (or
!!! even read) the first 3 columns (lat, lon, surf_elev_wgs84).
!!!
!!! These files don't have a time dimension (and the value is not known,
!!! except perhaps from the filename), so this module returns a 1-element
!!! t, set to the invalid value.
!!!
!!! The grids are filled such that X and Y *increase*.
!!!
!!! The correct polar stereographic grid must be set up using
!!! cism_polarstereo_mod before calling this.
!!!
!!! TBD: 1) x(ix),y(iy) are rewritten for every row/column, should check
!!! that the new values agree with the earlier ones.  2) Should
!!! read_model_cism_txt return a status? Currently just stops.  3) When
!!! reading file, handle comment and blank lines.  4) Warn if xf,yf too far
!!! from integral values.  5) Convert to standard 1m grid.  6) Allow
!!! specification of which column to get the data from.
!!!
!!! AUTHOR: Jeff Guerber, GSFC 615/Sigma Space, June 2015
!!! 2015-10-28 JRG: Renamed from cism_txt_file_mod.f08.
!!!    Renamed read_cism_txt_file to read_model_cism_txt.
!!!    Changed argument "elevsWgs84" to more general "values".
!!!    Use NF90_FILL_DOUBLE as the invalid value and return as an argment.
!!! 2015-11-18 JRG: Output fixes.
!!! 2016-04-13 JRG: Values is now 3D, allocated to (Nx,Ny,1).
!!!    Added t argument.
!!! 2016-04-20 JRG: As a diagnostic, write out x,y,t,lats,lons,values at
!!!    two corners and the center.  Various minor changes mainly to printed
!!!    information.
!!! 2016-05-05 JRG: Handle I/O errors (more) gracefully.
!!! 2016-05-11 JRG: read_model_cism_txt: Deallocate the allocatable dummy
!!!    arguments before allocating them.  Calculate invflagval so we don't
!!!    have to require netcdf library just for that.
!!! 2016-05-20 JRG: Return time units (always "invalid").
!!!    Minor clarification on printout of model invflagval.
!!!
!!! Last SVN commit: $Id: model_cism_txt_mod.f08 79 2016-05-25 07:41:05Z jguerber $
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module model_cism_txt_mod

  use, intrinsic :: iso_fortran_env, only: REAL64, IOSTAT_END
  use :: polar_stereographic_mod
  use :: cism_polarstereo_mod

  implicit none
  private
  public read_model_cism_txt


contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! NAME: read_model_cism_txt
  !!
  !! PURPOSE: Read the text-format CISM model file.  Latitude and longitude
  !! are converted into CISM standard polar-stereographic grid coordinates,
  !! and the elevation data is placed in an array.  These files only
  !! contain one time step, so the returned elevation data array is
  !! dimensioned (Nx,Ny,1).  Points within the grid that the model file did
  !! not define are filled with NF90_FILL_DOUBLE (but we may change to IEEE
  !! NaN in the future once we have IEEE support in our gfortran).
  !!
  !! INPUT ARGUMENTS:
  !!   FILE: Name of the file to read.
  !!
  !! OUTPUT ARGUMENTS:
  !!   VALUES: The data values, read from the file, on the PS grid.
  !!       Eg., for Steve Price's cism_usrf_* files this is the surface
  !!       elevations relative to WGS84.  Allocatable 3D array of doubles
  !!       (1:Nx,1:Ny,1).
  !!   X: Allocatable 1D array of PS X values (1:Nx) of the grid points.
  !!   Y: Allocatable 1D array of PS Y values (1:Ny) of the grid points.
  !!        X,Y are in units of the grid spacing, E and N from the
  !!        origin. To convert to km, multiply by the grid spacing.
  !!   T: Allocatable 1D array of time values. These files don't have a
  !!        time dimension though, so returns 1 element, set to invflagval.
  !!   LATS, LONS: Latitudes and longitudes of each points on the grid.
  !!      Allocatable 2D arrays (1:Nx,1:Ny).  Since the grids are in (X,Y)
  !!      space, lat and lon vary.
  !!   INVFLAGVAL: Value used in the data arrays to flag invalid points.
  !!      Currently we use NF90_FILL_DOUBLE from NetCDF, but IEEE NaN
  !!      could be better.  Whichever way, it will be returned in this.
  !!   T_UNITS: Units of time.  Since time isn't defined for these files,
  !!      always set to "invalid".
  !!
  !!  It's OK if the allocatable actual arguments are already allocated.
  !!  They'll be reallocated here.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine read_model_cism_txt( file, values, x, y, t, lats, lons, &
       invflagval, t_units )

    !! Dummy arguments
    character(len=*), intent(in) :: file
    real(kind=REAL64), allocatable, intent(out) :: values(:,:,:)
    real(kind=REAL64), allocatable, intent(out) :: x(:), y(:), t(:)
    real(kind=REAL64), allocatable, intent(out) :: lats(:,:), lons(:,:)
    real(kind=REAL64), intent(out) :: invflagval
    character(len=*), intent(out) :: t_units

    !! Local variables

    ! the file: unit, i/o status, number of grid points, line
    integer :: lun, iostatus, nPnts
    character(len=256) :: line     ! a line read from the file
    character(len=256) :: iomsg    ! i/o status message

    ! lat, lon, elev read from the file (fl for "file"):
    real(kind=REAL64), allocatable :: flLats(:), flLons(:), flValues(:)

    ! CISM PS-grid x, y from the file, converted from flLats, flLons
    real(kind=REAL64), allocatable :: xf(:), yf(:)
    real(kind=REAL64) :: delta_x_intx, delta_y_inty  ! max x-int(x), y-int(y)

    ! Misc
    type( ps_params_def ) :: psparams   ! PS params from cism_polarstereo_mod
    integer :: psstatus    ! status from polar_stereographic_mod procedures
    integer :: ix, iy, Nx, Ny    ! indexes and sizes of grid arrays
    integer :: i
    real(kind=REAL64) :: xmin, ymin, xmax, ymax      ! limits of grid arrays

!!!------------------------------------------------------------------------

    print *
    print *, 'In model_cism_txt_mod%read_model_cism_txt'

    psparams = cism_ps_get_params()     ! From cism_polarstereo_mod
    if ( .not. psparams%initialized ) then
       print *
       print *, 'read_model_cism_txt: PS grid not initialized.'
       ! TBD: status?
       stop 1
    end if

    iostatus = 0
    iomsg    = ''

    !! Open the file
    open ( newunit=lun, file=file, form='formatted', access='sequential', &
         status='old', action='read', iostat=iostatus, iomsg=iomsg )
    if (iostatus /= 0) then
       print *
       print *, 'read_model_cism_text: Could not open cism-text model file: ', &
            trim(file)
       print '(a,i0,2a)', ' Status= ', iostatus, '  Msg= ', trim(iomsg)
       stop 1
    end if

    print *, '   Reading CISM-text format model file ', trim(file)

    ! Quick read through to find number of grid points, to allocate arrays
    ! (we can assume each line has a point, but probably good to make this
    ! more flexible).
    nPnts = 0
    iostatus = 0
    iomsg    = ''
    read (lun, *, iostat=iostatus, iomsg=iomsg)  ! ignore file's 1-line header
    do while (iostatus == 0)            ! read to end or error
       read (lun, *, iostat=iostatus, iomsg=iomsg)
       if (iostatus == 0) nPnts = nPnts + 1
    end do

    if (iostatus /= IOSTAT_END) then   ! error before reached end
       print *
       print *, 'read_model_cism_text: Error scanning cism-text model file:'
       print *, trim(file)
       print '(a,i0,2a)', ' Status= ', iostatus, '  Msg= ', trim(iomsg)
       stop 1
    end if

    if (nPnts == 0) then
       print *
       print *, 'read_model_cism_text: No points read from cism-text model file:'
       print *, trim(file)
       stop 1
    end if

    !! Allocate arrays for data read from the file.  They're grid points,
    !! but in a 1D list.  Then read them.  These files are CSV, so can
    !! simply use list-directed reads!  The files actually have their own X
    !! and Y fields (on a different PS grid spacing), but for now we ignore
    !! those and just use PS grids converted from the lat and lon.
    allocate ( flLats(nPnts), flLons(nPnts), flValues(nPnts), &
         xf(nPnts), yf(nPnts) )
    iostatus = 0
    iomsg    = ''
    rewind lun
    read (lun,*,iostat=iostatus,iomsg=iomsg)  ! ignore header line
    if (iostatus == 0) then
       do i = 1, nPnts
          read (lun,*,iostat=iostatus,iomsg=iomsg) &
               flLats(i), flLons(i), flValues(i)
          if (iostatus /= 0) exit
       end do
    end if
    close ( lun )

    if (iostatus /= 0) then
       print *
       print *, 'read_model_cism_text: Error reading data from cism-text model file:'
       print *, trim(file)
       print '(a,i0,2a)', ' Status= ', iostatus, '  Msg= ', trim(iomsg)
       stop 1
    end if

    print '(a,i0,a)', '    Read ', nPnts, ' records'

    !! Convert file lat,lon to CISM PS x,y.  LatLon2PstXY is in
    !! polar_stereographic_mod.
    call LatLon2PStXY( flLats, flLons, nPnts, psstatus, xf, yf )
    if (psstatus /= 0) then
       ! Only |lat|>90. triggers this.  LatLon2PStXY reports it itself.
       print *
       print *, 'read_model_cism_txt: Bad status from LatLon2PStXY'
       stop 1
    end if

    !! Real xf,yf should be (very close) to integral values in grid
    !! coordinates, since the model is (supposed to be!!) defined on the
    !! grid we're using, but have noise from the conversion.  Check to make
    !! sure that they are, so that we can use them to index the arrays.
    delta_x_intx = maxval( abs( xf - anint(xf) ) )
    delta_y_inty = maxval( abs( yf - anint(yf) ) )
    print '(a,2es12.3)', &
         '    Max errors of rounded x,y from lat/lon (should be small):', &
         delta_x_intx, delta_y_inty    ! Warn if too large?

    !! Snap xf,yf to grid points.  Allocate the grid arrays.
    xf = anint( xf )
    yf = anint( yf )

    xmin = minval( xf )
    ymin = minval( yf )
    xmax = maxval( xf )
    ymax = maxval( yf )
    Nx = nint(xmax) - nint(xmin) + 1
    Ny = nint(ymax) - nint(ymin) + 1

    ! Make sure the allocatable dummy arguments are not already allocated,
    ! then allocate them.
    if ( allocated( x ) )       deallocate ( x )
    if ( allocated( y ) )       deallocate ( y )
    if ( allocated( t ) )       deallocate ( t )
    if ( allocated( values ) )  deallocate ( values )
    if ( allocated( lats ) )    deallocate ( lats )
    if ( allocated( lons ) )    deallocate ( lons )
    allocate ( x(Nx), y(Ny), t(1), values(Nx,Ny,1), lats(Nx,Ny), lons(Nx,Ny) )

    ! invflagval is NF90_FILL_DOUBLE = 9.96920996838686905e+36 =
    ! 15.*(2**119).  Note that it converts to REAL32 without losing
    ! precision.
    invflagval = scale(15.0_REAL64, 119)
    x          = invflagval
    y          = invflagval
    t          = invflagval
    t_units    = 'invalid'
    values     = invflagval
    lats       = invflagval
    lons       = invflagval

    print '(a,f0.3)', '    Grid point spacing, km: ', psparams%spacing
    print '(a,i0,a,f0.2,a,f0.2)', '    Grid extent, X: ', Nx, &
         ' from ', xmin, ' to ', xmax
    print '(a,i0,a,f0.2,a,f0.2)', '    Grid extent, Y: ', Ny, &
         ' from ', ymin, ' to ', ymax
    print '(a,es15.7)', '    Model data invalid flag value: ', invflagval

    !! Now fill them up.  Note this fills the grids such that X and Y will
    !! be *increasing*.  (For decreasing X or Y, I think we would index by
    !! (xmax-xf(i)) or (ymax-yf(i)).)
    do i = 1, nPnts
       ix = nint( (xf(i) - xmin) ) + 1   ! Indices into the grid arrays
       iy = nint( (yf(i) - ymin) ) + 1
       ! TBD: SHOULD CHECK THAT XF,YF ARE "CLOSE" TO WHAT'S ALREADY AT THIS
       ! POINT, AND THIS POINT'S VALUE HASN'T ALREADY BEEN FILLED ALSO CHECK
       ! THAT THE VALUES ARE WITHIN THE GRID
       x(ix) = xf(i)
       y(iy) = yf(i)
       values(ix,iy,1) = flValues(i)
       lats(ix,iy) = flLats(i)
       lons(ix,iy) = modulo( flLons(i), 360.0_REAL64 )   ! force to (0,360] E
    end do

    print *, '   At (1,1), (Nx/2,Ny/2), (Nx,Ny):'
    print '(a,3g16.7)', '      x:     ', x(1), x(Nx/2), x(Nx)
    print '(a,3g16.7)', '      y:     ', y(1), y(Ny/2), y(Ny)
    print '(a,g16.7)',  '      t(1):  ', t(1)
    print '(a,3g16.7)', '      lats:  ', lats(1,1), lats(Nx/2,Ny/2), lats(Nx,Ny)
    print '(a,3g16.7)', '      lons:  ', lons(1,1), lons(Nx/2,Ny/2), lons(Nx,Ny)
    print '(a,3g16.7)', '      values:', values(1,1,1), values(Nx/2,Ny/2,1), &
         values(Nx,Ny,1)

    return
  end subroutine read_model_cism_txt

end module model_cism_txt_mod
