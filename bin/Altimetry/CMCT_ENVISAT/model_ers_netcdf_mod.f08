!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! NAME: model_netcdf_mod.f08
!!!
!!! PURPOSE: Module for reading model files in netCDF format.
!!!
!!! ASSUMPTIONS/NOTES:
!!!   * Variable to be read is 3D with dimensions (x,y,t). ((t,y,x) in C order.)
!!!   * The grid is rectangular in polar-stereographic projection space.
!!!   * Variables providing latitude and longitude are present.  They are
!!!     2D arrays (x,y) with standard_name attributes "latitude" and
!!!     "longitude".  (3D (x,y,1) lat/lon arrays have been observed and
!!!     found to work.)  X and Y are then calculated from them using our
!!!     own projection routines.  (The file's X and Y are not used, because
!!!     it's important that we have a consistent projection.)
!!!   * Data is not packed or scaled: Attributes "scale_factor" and
!!!     "add_offset" are currently NOT handled explicitly.  (Not clear to
!!!     me whether the library may do so automatically.  User's Guide and
!!!     lib docs don't say, "NetCDF Best Practices" says they don't.)
!!!   * _FillValue attribute is used if present, otherwise the netCDF
!!!     default NF90_FILL_DOUBLE is used.  (NF90_FILL_DOUBLE and
!!!     NF90_FILL_FLOAT are actually the same value, 15*(2**119), so they
!!!     should be interconvertable without issues.)  Attributes
!!!     "missing_value", "valid_range", "valid_min", and "valid_max" are
!!!     not handled.
!!!   * Output arrays are in order of *increasing* X and Y.  Input may be
!!!     either way, but will be reversed if needed.
!!!   * The "units" attribute of only the time variable is returned.  Units
!!!     of values, lat, and lon are assumed.
!!!
!!! The correct polar stereographic grid must be set up using
!!! cism_polarstereo_mod before calling this.
!!!
!!! AUTHOR: Jeff Guerber, GSFC 615/Sigma Space, April 2016
!!! 2016-04-20 JRG: Drop use of file's X/Y if lat/lon not present.  This
!!!    didn't work due to differences in the projected grids, and Sophie
!!!    says it's OK to require that the model file have lat/lon arrays.
!!!    Ensure that output arrays are in order of increasing X and Y.
!!!    As a diagnostic, write out x,y,t,lats,lons,values at two corners and
!!!    the center.
!!!    Many minor changes.
!!! 2016-05-11 JRG: read_model_netcdf: Deallocate the allocatable dummy
!!!    arguments before allocating them.
!!! 2016-05-20 JRG: Read time units.  Minor change to clarify double
!!!    printout of invflagval.
!!!
!!! Last SVN commit: $Id: model_netcdf_mod.f08 79 2016-05-25 07:41:05Z jguerber $
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module model_netcdf_mod

  use :: netcdf
  use :: polar_stereographic_mod
  use :: cism_polarstereo_mod
  use, intrinsic :: iso_fortran_env, only: REAL64

  implicit none
  private
  public read_model_netcdf

contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! NAME: read_model_netcdf
  !!
  !! PURPOSE: Read the netCDF-format model file.  Data is returned in an
  !! allocatable 3D array (Nx,Ny,Nt).  If latitude and longitude are in the
  !! file, they are converted into CISM polar-stereographic X and Y grid
  !! coordinates (even if X and Y are in the file, because we want the
  !! projection to be consistent).  If latitude and longitude are not in
  !! the file, but X and Y are, we derive lat and lon from the file's X and
  !! Y, trusting its projection.
  !!
  !! NOTE: Returned X, Y are on a grid in units of the grid spacing, not
  !! the ISMIP6 standard meters.  This is because upstream we use them to
  !! index the arrays.  Eventually this should be fixed.
  !!
  !! INPUT ARGUMENTS:
  !!   FILE: Name of the file to read.
  !!   VARNAME: Name of the netCDF variable to return.
  !!
  !! OUTPUT ARGUMENTS:
  !!   VALUES: Data values of file variable VARNAME, read from the file, on
  !!        the PS grid.  Allocatable 3D array of doubles (1:Nx,1:Ny,1:Nt).
  !!   X: Allocatable 1D array of PS X values (1:Nx) of the grid points.
  !!   Y: Allocatable 1D array of PS Y values (1:Ny) of the grid points.
  !!        NOTE: X,Y are in units of the grid spacing, E and N from the
  !!        origin. To convert to km, multiply by the grid spacing.
  !!   T: Allocatable 1D array of time values (1:Nt).  If time is not in
  !!      the file, value is set to invflagval.
  !!   LATS, LONS: Latitude and longitude of each point on the grid.
  !!      Allocatable 2D arrays (1:Nx,1:Ny).  Since the grids are in (X,Y)
  !!      space, lat and lon vary.
  !!   INVFLAGVAL: Value used in the data arrays to flag invalid points.
  !!      Attribute _FillValue if variable has it, else default to
  !!      NF90_FILL_DOUBLE.  (NB: Library automatically converts s.p. data
  !!      to double, including fill, plus I've verified that they really
  !!      are the same value.  Don't know if it converts integer data,
  !!      though.)
  !!   T_UNITS: Time variable units attribute.  'none' if no time has no
  !!      units attribute, 'invalid' if no time variable.
  !!
  !!  It's OK if the allocatable actual arguments are already allocated.
  !!  They'll be reallocated here.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine read_model_netcdf( file, varname, &
       values, x, y, t, lats, lons, invflagval, t_units )

    !! Dummy arguments
    character(len=*), intent(in) :: file
    character(len=*), intent(in) :: varname    ! name of variable to read
    real(kind=REAL64), allocatable, intent(out) :: values(:,:,:)  ! (x,y,t)
    real(kind=REAL64), allocatable, intent(out) :: x(:), y(:), t(:)
    real(kind=REAL64), allocatable, intent(out) :: lats(:,:), lons(:,:)
    real(kind=REAL64), intent(out) :: invflagval
    character(len=*), intent(out) :: t_units

    !! Local variables
    integer :: status, i
    character(len=500) :: msg    ! message for errors
    integer :: ncid      ! netcdf id for dataset
    integer :: nVars, nDims     ! number of variables, dimensions
    integer :: allDims(NF90_MAX_VAR_DIMS)     ! all dimensions in file
    character(len=50) :: stdname   ! standard_name attribute

    type info    ! information about each variable
       integer :: vid=-1    ! netcdf var id
       character(len=NF90_MAX_NAME) :: name    ! name of the variable
       integer :: ndims     ! number of dimensions
       integer :: dimids(NF90_MAX_VAR_DIMS)  !dim ids of this variable, in order
       integer :: xtype     ! type of the variable
    end type info
    type(info) ::  varInfo, timeInfo, latInfo, lonInfo

    integer :: Nx, Ny, Nt       ! dimensions of x, y, t
    real(kind=REAL64), allocatable :: xf(:,:), yf(:,:)   ! 2D x,y
    real(kind=REAL64) :: delta_x_intx, delta_y_inty  ! max x-int(x), y-int(y)
    real(kind=REAL64) :: fillValue, missingValue
    type( ps_params_def ) :: psparams   ! PS params from cism_polarstereo_mod
    real(kind=REAL64) :: xmin, ymin, xmax, ymax

!!!------------------------------------------------------------------------

    print *
    print *, 'In model_netcdf_mod%read_model_netcdf'

    psparams = cism_ps_get_params()
    if ( .not. psparams%initialized ) then
       print *, 'read_model_netcdf: PS grid not initialized.'
       ! TBD: status?
       stop 1
    end if

    status = nf90_open( path=file, mode=NF90_NOWRITE, ncid=ncid )
    call check_ncdf( status, 'Trying to open '// trim(file) )

    print *, '   Reading netCDF model file ', trim(file)
    print *, '   Looking for variable ', trim(varname)

    !! Get the numbers of vars and dims in this file, and each dimension
    status = nf90_inquire( ncid=ncid, nVariables=nVars, nDimensions=nDims )
    call check_ncdf( status, 'Inquiring for nVars and nDims')
    do i = 1, ndims
       status = nf90_inquire_dimension( ncid=ncid, dimid=i, len=allDims(i) )
       call check_ncdf( status, 'Looking for dimensions' )
    end do

    !! Is the requested variable in this file?  If so, get its varId,
    !! dimensions, and other information
    status = nf90_inq_varid( ncid=ncid, name=varName, varid=varInfo%vid )
    call check_ncdf( status, 'Requested variable not in the file! ' // &
         trim(varName) )
    varInfo%name = varName

    status = nf90_inquire_variable( ncid=ncid, varid=varInfo%vid, &
         ndims=varInfo%ndims, dimids=varInfo%dimids, xtype=varInfo%xtype )
    write (msg,*) 'Getting info about variable, varid=', varInfo%vid
    call check_ncdf( status, msg )
    Nx = allDims(varInfo%dimids(1))
    Ny = allDims(varInfo%dimids(2))
    Nt = allDims(varInfo%dimids(3))

    !! Get the invalid value.  _FillValue attribute if it exists, else use
    !! default.  Note that library converts to requested type.  (And I've
    !! verified that NF90_FILL_REAL and NF90_FILL_DOUBLE are the same.)
    !! Not currently handling integer types - don't know if those fills get
    !! converted properly.
    status = nf90_get_att( ncid=ncid, varid=varInfo%vid, &
         name='_FillValue', values=fillValue )
    if (status == NF90_NOERR ) then
       invflagval = fillValue
       print *, '   Using attribute _FillValue as invalid: ', invflagval
    else
       invflagval = NF90_FILL_DOUBLE
       print *, '   Using default double fill as invalid: ', invflagval
    end if

    !! Allocate the variable arrays.  First make sure the allocatable dummy
    !! arguments are not already allocated.
    if ( allocated( x ) )       deallocate ( x )
    if ( allocated( y ) )       deallocate ( y )
    if ( allocated( t ) )       deallocate ( t )
    if ( allocated( values ) )  deallocate ( values )
    if ( allocated( lats ) )    deallocate ( lats )
    if ( allocated( lons ) )    deallocate ( lons )
    allocate( x(Nx), y(Ny), t(Nt), values(Nx,Ny,Nt), lats(Nx,Ny), lons(Nx,Ny) )

    ! TBD: Should the others also be initialized?  model_cism_txt_mod does.
    values = invflagval

    !! Read the data variable
    status = nf90_get_var( ncid=ncid, varid=varInfo%vid, values=values )
    write (msg,*) 'Error reading Values, variable= ', &
         trim(varInfo%name), ', vid=', varInfo%vid
    call check_ncdf( status, msg )

    !! Find the varids of the other variables, that we need to get via
    !! their standard_name attributes.  Loop over all the variables and
    !! pick out the ones we may need.
    do i = 1, nVars

       status = nf90_get_att( ncid=ncid, varid=i, &
            name='standard_name', values=stdname )
       if (status == NF90_ENOTATT) cycle  ! just skip if this attr missing
       call check_ncdf( status, 'Reading standard_name attribute' )

       select case (stdname)
       case ('time')
          timeInfo%vid = i
       case('latitude')
          latInfo%vid = i
       case('longitude')
          lonInfo%vid = i
       end select

    end do

    !! Read each of the these variables.  For the coordinate
    !! variables, also check that they match the corresponding dimensions
    !! of value.  Coordinate variables also should all be 1D, and we should
    !! check that the ones we found are.

    if (timeInfo%vid /= -1) then          ! Time
       status = nf90_inquire_variable( ncid=ncid, varid=timeInfo%vid, &
            ndims=timeinfo%ndims, dimids=timeInfo%dimids )
       write (msg,*) 'Getting info about Time coord var, varid=', timeInfo%vid
       call check_ncdf( status, msg )

       if ( alldims(timeInfo%dimids(1)) /= Nt ) then
          print *, 'read_model_netcdf: T dimension mismatch! '
          print *, trim(timeInfo%name), ' = ', alldims(timeInfo%dimids(1))
          print *, trim(varInfo%name), ' = ', Nt
          stop 1
       end if

       status = nf90_get_var( ncid=ncid, varid=timeInfo%vid, values=t )
       write (msg,*) 'Error reading Time variable ', &
            trim(timeInfo%name), ' vid=', timeInfo%vid
       call check_ncdf( status, msg )

       status = nf90_get_att( ncid=ncid, varid=timeInfo%vid, &
            name='units', values=t_units )
       select case (status)
       case (NF90_NOERR)
          continue
       case (NF90_ENOTATT)
          t_units = 'none'
       case default
          write (msg,*) 'Getting time:units attribute'
          call check_ncdf( status, msg )
       end select

    else    ! no time variable
       t = invflagval
       t_units = 'invalid'
    end if

    if (latInfo%vid /= -1) then      ! Latitude
       status = nf90_inquire_variable( ncid=ncid, varid=latInfo%vid, &
            name=latInfo%name, ndims=latInfo%ndims, dimids=latInfo%dimids )
       write (msg,*) 'Getting info about latitude var, varid=', latInfo%vid
       call check_ncdf( status, msg )

       status = nf90_get_var( ncid=ncid, varid=latInfo%vid, values=lats )
       write (msg,*) 'Error reading latitude variable ', &
            trim(latInfo%name), ' vid=', latInfo%vid
       call check_ncdf( status, msg )
    end if

    if (lonInfo%vid /= -1) then    ! Longitude
       status = nf90_inquire_variable( ncid=ncid, varid=lonInfo%vid, &
            name=lonInfo%name, ndims=lonInfo%ndims, dimids=lonInfo%dimids )
       write (msg,*) 'Getting info about longitude var, varid=', lonInfo%vid
       call check_ncdf( status, msg )

       status = nf90_get_var( ncid=ncid, varid=lonInfo%vid, values=lons )
       write (msg,*) 'Error reading longitude variable ', &
            trim(lonInfo%name), ' vid=', lonInfo%vid
       call check_ncdf( status, msg )
    end if

    !!
    !! Calculate X, Y from Lats, Lons using our projection and grid.
    !!

    if ((latInfo%vid == -1) .or. (lonInfo%vid == -1) ) then
       print *, 'read_model_netcdf: File has no lat/lon arrays!'
       stop 1
    end if
    print *, '   Found lat & lon, deriving x,y'

    !! Convert lats, lons into 2D x, y arrays on our own projection.
    !! LatLon2PStXY is in polar_stereographic_mod, should have already
    !! been set up.
    allocate( xf(Nx,Ny), yf(Nx,Ny) )
    call LatLon2PStXY( lats, lons, Nx, Ny, status, xf, yf )
    if (status /= 0) then
       ! Only |lat|>90. triggers this.  LatLon2PStXY reports it itself.
       print *, 'read_model_netcdf: Bad status from LatLon2PStXY'
       stop 1
    end if
    lons =  modulo( lons, 360.0_REAL64 )   ! force to (0,360] E

    !! Real xf,yf should be (very close) to integral values in grid
    !! coordinates, since the model is (supposed to be!!) defined on the
    !! grid we're using, but have noise from the conversion.  Check to make
    !! sure that they are, so that we can use them to index the arrays.
    delta_x_intx = maxval( abs( xf - anint(xf) ) )
    delta_y_inty = maxval( abs( yf - anint(yf) ) )
    print '(a,2es12.3)', &
         '    Max errors of rounded x,y from lat/lon (should be small):', &
         delta_x_intx, delta_y_inty    ! Warn if too large?

    !! TBD: We should also check to see how well each row and column
    !! match the others.

    !! Return appropriate slices
    x = anint( xf(:,1) )
    y = anint( yf(1,:) )

    !!  Make sure the arrays are in order of increasing x and y, reverse
    !!  them if necessary.
    if (x(1) > x(Nx)) then
       x = x(Nx:1:-1)
       values = values(Nx:1:-1, :, :)
       lats = lats(Nx:1:-1, :)
       lons = lons(Nx:1:-1, :)
    end if

    if (y(1) > y(Ny)) then
       y = y(Ny:1:-1)
       values = values(:, Ny:1:-1, :)
       lats = lats(:, Ny:1:-1)
       lons = lons(:, Ny:1:-1)
    end if

    print '(a,f0.3)', '    Grid point spacing, km: ', psparams%spacing
    print '(a,i0,a,f0.2,a,f0.2)', '    Grid extent, X: ', Nx, &
         ' from ', x(1), ' to ', x(Nx)
    print '(a,i0,a,f0.2,a,f0.2)', '    Grid extent, Y: ', Ny, &
         ' from ', y(1), ' to ', y(Ny)
    print '(a,es15.7)', '    Model data invalid flag value: ', invflagval

    print *, '   At (1,1,1), (Nx/2,Ny/2,1), (Nx,Ny,1):'
    print '(a,3g16.7)', '      x:     ', x(1), x(Nx/2), x(Nx)
    print '(a,3g16.7)', '      y:     ', y(1), y(Ny/2), y(Ny)
    print '(a,g16.7)',  '      t(1):  ', t(1)
    print '(a,3g16.7)', '      lats:  ', lats(1,1), lats(Nx/2,Ny/2), lats(Nx,Ny)
    print '(a,3g16.7)', '      lons:  ', lons(1,1), lons(Nx/2,Ny/2), lons(Nx,Ny)
    print '(a,3g16.7)', '      values:', values(1,1,1), values(Nx/2,Ny/2,1), &
         values(Nx,Ny,1)

    return
  end subroutine read_model_netcdf


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Private subroutine check_ncdf
!!! Checks the status of netcdf calls.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine check_ncdf ( status, msg )
    integer, intent(in) :: status
    character(len=*), intent(in) :: msg

    if ( status /= nf90_NoErr ) then
       print *
       print *, 'model_netcdf_mod: NetCDF error ' // msg
       print *, '   nf90 status = ', status
       print *, '   ', trim( nf90_strerror(status) )
       stop 1
    end if
    return

  end subroutine check_ncdf

end module model_netcdf_mod
