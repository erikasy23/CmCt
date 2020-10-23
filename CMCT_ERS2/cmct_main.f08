!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!! NAME: cmct_main.f08
!!!
!!! PURPOSE: Cryospheric Model Comparison Tool, main program.
!!!
!!! NOTE: Currently assumes that the units of the values being
!!! differenced are meters.
!!!
!!! AUTHOR: Jeff Guerber (JRG), SigmaSpace/GSFC 615, July 2015
!!! 2015-10-28 JRG: Changed a bunch of variable names to be more readable.
!!!    cism_txt_file_mod is now model_cism_txt_mod and read_cism_txt_file
!!!    is now read_model_cism_txt, with arg list change and invalid.
!!! 2015-11-03 JRG: Added overall stats - diff_MomsTot
!!! 2015-11-18 JRG: Bug fixes. Better output. Debug level replaced by
!!!    verbosity, read from ctl file.
!!! 2015-11-23 JRG: Add _FillValue keywords, set invalid elements of
!!!    diff_Means, diff_SDs to default fill value.
!!! 2015-12-01 JRG: Added x_grid, y_grid coordinate variables to
!!!    .meansdgrid.nc file.  Minor output tweaks.
!!! 2015-12-04 JRG: Changed run config icesat-glas section to
!!!    mission-options. Added o/p of model information on stdout.
!!!    Added netCDF global attributes. baseOut as base of o/p filenames.
!!! 2015-12-31 JRG: Added histogram.
!!! 2016-01-05 JRG: Addded histogram cdf.
!!! 2016-03-03 JRG: Clean observation record.
!!! 2016-03-09 JRG: Misc small fixes based on recent test.
!!! 2016-03-31 JRG: Moved run-config file parsing out into cmct_runcfg_mod.
!!! 2016-04-09 JRG: Changed delta from (obs-model) to (model-obs), so +
!!!    now "model higher than observed". Sophie says this is typical.
!!!    Also netCDF output file is now netCDF4/HDF5 format.
!!! 2016-04-13 JRG: Read netCDF model files. Model time axis: model_T,
!!!    iTime, and time. model_Data now 3D (x,y,t).
!!! 2016-04-16 JRG: Added bug banner.
!!! 2016-04-21 JRG: Changed format of baseOut. Include model file and time
!!!    in reclist and histogram headers.  Several minor fixes.
!!! 2016-04-21 JRG: Log model time and index.
!!! 2016-04-29 JRG: Loop over comparisons in the run config file.
!!! 2016-05-12 JRG: Set moments invalid value. Various minor fixes.
!!! 2016-05-20 JRG: Check that time index is within bounds of time array.
!!!    Get time units from models, write to output files.
!!!    Write model_time_index and model_time global attributes in the
!!!    meansdgrids.nc file (probably should really be a 1-elem time axis).
!!! 2016-05-24 JRG: Treat iTime range for cism-text model files (whose iTime
!!!    should always be 1) same as for netcdf model files.
!!!    Added units and standard_name attributes to the netcdf MeanSDGrid file.
!!! 2016-06-08 JRG: Use new datasets_class to get observations. ds now a
!!!    datasets_class object.
!!! 2016-07-13 JRG: Put both dataset name and description in log and o/p files.
!!! 2016-07-14 JRG: Bumped version to BETA-1.
!!! 2016-07-15 JRG: Bug fix: Call datasets_info before ds%open.
!!!
!!! Last SVN commit: $Id: cmct_main.f08 116 2016-07-18 22:37:33Z jguerber $
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program cmct_main

  use, intrinsic :: iso_fortran_env, only: INT32, REAL32, REAL64
  use :: fson
  use :: fson_value_m
  use :: cism_polarstereo_mod
  use :: polar_stereographic_mod
  use :: model_cism_txt_mod
  use :: model_netcdf_mod
  use :: moments2d_class_mod
  use :: datasets_mod
  use :: netcdf
  use :: histogram_class_mod
  use :: cmct_runcfg_mod

  implicit none

  !!
  !! LOCAL VARIABLE DECLARATIONS
  !!

  character(len=*), parameter :: CMCT_VERSION="BETA-1"

  !
  ! Main Control file: fson parse tree pointer, and items read from it
  !
  type( fson_value ), pointer :: ctl_File=>null()    ! top level
  character(len=256) :: ctl_FileName  ! Command line arg: main control file name
  character(len=256) :: models_Dir    ! Directory where the models are
  character(len=256) :: run_CfgFileName     ! Run Configuration file name
  integer :: verbosity     ! verbosity of output, 0=normal, (>0)=more

  !
  ! Run configuration.  cmct_run_config is defined in cmc_runcfg_mod.f08.
  ! The run config file is designed to include more than one comparison in
  ! a run.
  !
  type( cmct_run_config ) :: runcfg
  integer :: comparison    ! Which comparison is this

  !
  ! Gridded model.  cism-text format model files don't have a time
  ! dimension, so just make it 1.
  !
  character(len=500) :: model_FullName        ! Full path to the model file
  real(kind=REAL64), allocatable :: model_Data(:,:,:)    ! data array (x,y,t)
  real(kind=REAL64), allocatable :: model_X(:), model_Y(:)   ! PS x, y
  real(kind=REAL64), allocatable :: model_T(:)       ! time array
  character(len=50) :: model_Tunits                  ! units of time
  real(kind=REAL64), allocatable :: model_Lats(:,:), model_Lons(:,:)
  real(kind=REAL64) :: model_Invalid  ! data invalid value
  real(kind=REAL64) :: xmin, xmax, ymin, ymax     ! bounds of the grid
  integer           :: iTime    ! index into model_T(:)
  real(kind=REAL64) :: time     ! value of model_T(iTime)

  !
  ! Observation dataset.  Managing the different datasets is now
  ! consolodated within datasets_class which is defined in
  ! datasets_mod.f08.  Which dataset to use is specified in the run config
  ! file.
  !
  type( datasets_class ) :: ds     ! observations datasets class
  real(kind=REAL64) :: obsLat, obsLon  ! observation loc from ds%next_obs
  real(kind=REAL32) :: obsValue    ! obs value at obsLat,obsLon from ds%next_obs
  real(kind=REAL64) :: obsX, obsY  ! obsLat,obsLon converted to grid coords
  character(len=DS_DESCLEN) :: ds_desc  ! dataset description string

  !
  ! Interpolation of the observation point into the model grid.
  !
  real(kind=REAL64) :: cell(0:1,0:1) ! grid data surounding the obs point
  integer :: xfloor, yfloor          ! grid points to lower left of obs pt
  real(kind=REAL64) :: p, q, f       ! bilinear interpolation, f=result
  real(kind=REAL64) :: delta         ! = f - obsValue

  !
  ! Results
  !
  integer :: szXout, szYout    ! OUTPUT grid sizes in X, Y directions
  ! Difference moments objects:  diff_Moms: on grid, diff_MomsTot: overall
  !   (NB: args must be 2D, allocatable, & deferred-shape, even if 1x1.)
  type(moments2d_class) :: diff_Moms, diff_MomsTot
  ! nPts: Number of data points in each grid cell.  nPtsTot: total data points
  integer(kind=INT32), allocatable :: nPts(:,:), nPtsTot(:,:)
  ! means, variances, and std deviations of (data - model), each grid cell:
  real(kind=REAL64), allocatable ::   &
       diff_Means(:,:), diff_Vars(:,:), diff_SDs(:,:)
  ! Overall mean, var, std dev of (data-model):
  real(kind=REAL64), allocatable ::  &
       diff_MeanTot(:,:), diff_VarTot(:,:), diff_SDTot(:,:)
  real(kind=REAL64) :: diff_InvFlag    ! invalid-value flag from moments2d_class
  real(kind=REAL64) :: minDelta, maxDelta     ! min and max delta
  ! numbers of obs records read, bad, outside grid, outside model, accepted:
  integer :: nRecsRead, nBad, nOutGrid, nOutModel, nGood
  ! Distance to nearest model node
  real(kind=REAL64) :: distNode, maxDistNode, distNode_max
  integer :: nDistMax
  ! Histogram
  type(histogram_class) :: hist1m      ! 1m histogram
  integer(kind=INT32) :: histnbins, histtotal, histbelow, histabove
  real(kind=REAL32) :: histmin, histmax, histbinsize
  integer(kind=INT32), allocatable :: histogram(:)
  real(kind=REAL32), allocatable :: histbins(:), histpdf(:), histcdf(:)

  !
  ! Output files
  !
  character(len=256) :: baseOut, listOut, meanOut, histOut  ! output file names
  integer :: listLun, histLun       ! LUNs for record list, histogram files
  character(len=80)  :: listFmt     ! format string for listLun records
  ! NetCDF variables:
  integer :: meanNCID, xdimID, ydimID    ! file and dimension IDs
  integer :: xVID, yVID, latsVID, lonsVID, nVID, meansVID, sdVID  !Variable IDs

  !
  ! Misc
  !
  integer :: status, i
  real    :: cputime_1, cputime_last   ! starting, ending CPU times

  call cpu_time( cputime_1 )

  print *, 'CRYOSPHERE MODEL COMPARISON TOOL'
  print *, 'NASA Goddard Space Flight Center, Cryospheric Sciences Lab'
  print *, 'by Erika Simon Innovim LLC./NASA GSFC 615'
  print *, 'cmct_main version ', CMCT_VERSION

  print *
  print *, '**********************************************'
  print *, '*      THIS IS BETA-LEVEL SOFTWARE!!!        *'
  print *, '*      THERE MAY WELL STILL BE BUGS!!!       *'
  print *, '*   Report them to erika.g.simon@nasa.gov     *'
  print *, '**********************************************'
  print *

  !!
  !! READ CONTROL/CONFIGURATION FILES
  !!

  !! Read main control file name from command line, then parse its JSON
  ctl_FileName = ''
  call get_command_argument( number=1, value=ctl_FileName, status=status )
  print *, 'cmct_main: Control file name = "', trim(ctl_FileName), '"'
  if (status /= 0) then
     print *, "cmct_main: get_command_argument status = ", status
     print *, "cmct_main: Usage: cmct_main control_file"
     stop 1
  end if
  ctl_File => fson_parse( ctl_FileName )
  call fson_get( ctl_File, 'cmct_main.verbosity', verbosity )
  call fson_get( ctl_File, 'cmct_main.models_dir', models_Dir  )

  !! Read info about known datasets from datasets.json
  call datasets_init( verbosity=verbosity )

  !! Run Config file holds the info returned by the web form.  Its name is
  !! in the control file.
  run_CfgFileName = ''
  call fson_get( ctl_File, 'cmct_main.cmct_run_config', run_CfgFileName )
  print *
  print *, 'cmct_main: Run config file name = "', trim(run_CfgFileName), '"'

  runcfg = parse_runcfg( run_CfgFileName )

  print *, 'cmct_main: Run ID       = ', trim(runcfg%run%runid)
  print *, 'cmct_main: Submitted on = ', trim(runcfg%run%date)
  print *, 'cmct_main: Submitted by = ', trim(runcfg%run%actualusername)
  print *, 'cmct_main: User Title   = "', trim(runcfg%run%user_run_title), '"'

  !!
  !! DO THE COMPARISONS
  !!
  !! Loop over all the entries in the run configuration file's
  !! "comparisons" array (runcfg%comps array).
  !!
  COMPLOOP:  do comparison = 1, runcfg%nComps

     print *
     print '( x, 65("=") )'
     print *
     print '(a,i0)', ' cmct_main: Comparison # ', comparison


     associate( modelInfo => runcfg%comps(comparison)%model, &
          obsInfo => runcfg%comps(comparison)%obs )

       ! baseOut is base of the output file names
       write (baseOut, '("CMCT_", a, ".", i2.2)' ) &
            trim(runcfg%run%runid), comparison

       print *
       print *, 'cmct_main: Model File Info:'
       print *, '   Name: ', trim(modelInfo%modelname)
       print *, '   File: ', trim(modelInfo%filename)
       print *, '   Format: ', trim(modelInfo%Format)
       print *, '   Region: ', trim(modelInfo%Region)
       print *, '   Variable: ', trim(modelInfo%variable)
       print *, '   Comments: ', trim(modelInfo%user_comments)

       !! Initialize the polar stereographic module for the grid and region used
       !! by this model.  (Region selection probably ought to be done inside
       !! cism_polarstereo_mod.)
       select case ( trim(modelInfo%Region) )
       case ( 'greenland' )
          call cism_ps_init_std_grn( modelInfo%spacing_km, status )
       case default
          print *,'cmct_main: Unknown region='//trim(modelInfo%Region)// &
               '. Stopping.'
          stop 1
       end select
       if (status /= 0) then
          print *,'cmct_main: cism_ps_init_std_grn error, stopping, status=', status
          stop 1
       end if

       !!
       !! Get the model data, PS x, PS y, t, lat, and lon arrays from the file
       !!
       model_FullName = trim(models_Dir) // '/' // trim(modelInfo%filename)
       select case ( trim(modelInfo%Format) )

       case ( 'cism-text' )
          ! read_model_cism_txt from model_cism_txt_mod.
          ! These files only contain 1 time, and model_time_index in the
          ! run config file should always be 1 for them.
          call read_model_cism_txt( file=model_FullName, values=model_Data, &
               x=model_X, y=model_Y, t=model_T, &
               lats=model_Lats, lons=model_Lons, &
               invflagval=model_Invalid, t_units=model_Tunits )

       case ( 'netcdf' )
          ! read_model_netcdf from model_netcdf_mod.
          ! Use time specified in run config file.
          call read_model_netcdf( file=model_FullName, varname=modelInfo%variable,&
               values=model_Data, &
               x=model_X, y=model_Y, t=model_T, &
               lats=model_Lats, lons=model_Lons, &
               invflagval=model_Invalid, t_units=model_Tunits )

       case default
          print *, 'cmct_main: Unknown model format='//trim(modelInfo%Format)// &
               '. Stopping.'
          stop 1
       end select

       ! Check that model_time_index (iTime) is within range:
       iTime = modelInfo%model_time_index
       if ( (lbound(model_T,1)<=iTime) .and. (iTime<=ubound(model_T,1)) ) then
          time = model_T(iTime)
       else
          print *
          print '(a,i0,a,i0,a,i0,a)', ' cmct_main: Time index (', iTime, &
               ') outside bounds of model times (', &
               lbound(model_T), ':', ubound(model_T), ')'
          print *, 'Comparison aborted.'
          cycle COMPLOOP
       end if

       ! Caution! This probably won't work if model_Invalid is NaN:
       xmin = minval( model_X, mask = model_X /= model_Invalid )
       xmax = maxval( model_X, mask = model_X /= model_Invalid )
       ymin = minval( model_Y, mask = model_Y /= model_Invalid )
       ymax = maxval( model_Y, mask = model_Y /= model_Invalid )

       print *
       print '(a,i0)',    ' cmct_main: Using model time index: ', iTime
       print '(a,g16.7)', '    Model time at index: ', time

       !!
       !! OPEN THE APPROPRIATE OBSERVATION FILE OBJECT.
       !!
       call datasets_info( dataset=obsInfo%dataset, desc=ds_desc, status=status )
       print *
       print *, 'cmct_main:  Comparison dataset: ' // trim(obsInfo%dataset)
       print *, '   ', trim(ds_desc)

       call ds%open( mission=obsInfo%mission, dataset=obsInfo%dataset, &
            ers2_campaign=obsInfo%campaign, &
            status=status, verbosity=verbosity )
       if (status /= DS_OPEN_OK) then
          print *
          print *, 'cmct_main:  Could not open dataset reader!'
          print *, 'Comparison aborted.'
          cycle COMPLOOP
       end if

       !!
       !! SET UP DIFFERENCE MOMENTS.  One for grid cells, one for everything.
       !!
       ! The grid points form the corners of the cells, so there's 1 less
       ! cell in each direction than model grid points.  We've standardized
       ! on using NetCDF default fill values for invalid values.
       szXout = size( model_X ) - 1
       szYout = size( model_Y ) - 1
       call diff_Moms%new( ubi=szXout, ubj=szYout, invflag=NF90_FILL_DOUBLE )
       ! overall moments obj only has 1 "cell":
       call diff_MomsTot%new( ubi=1, ubj=1, invflag=NF90_FILL_DOUBLE )

       !!
       !! INITIALIZE HISTOGRAM OBJECT
       !!
       ! 1m histogram.  -700 to +1000 is the range used in cism_results.pro (it
       ! overrides its histrange argument).  (cism_results.pro actually collects
       ! this range at 25m and only does 1m for the central 90%, but this could
       ! be rebinned to 25m easily in postprocessing.)  TBD: These values should
       ! be in configuration.
       call hist1m%new( min=-700., max=1000., binsize=1. )

       !!
       !! OPEN RECORD LIST OUTPUT FILE.  Lists every observation record we
       !! compare to the model.
       !!
       listOut = trim(baseOut) // '.reclist.txt'
       open ( newunit=listLun, file=listOut, status='new', form='formatted', &
            action='write' )    ! TBD: need to handle errors

       ! Header.  TBD: Some lines are currently icesat-specific
       write (listLun, '(a)')    '# CMCT RECORD LIST'
       write (listLun, '(2a)')   '# RUNID      = ', trim(runcfg%run%runid)
       write (listLun, '(2a)')   '# TITLE      = ', trim(runcfg%run%user_run_title)
       write (listLun, '(2a)')   '# DATE       = ', trim(runcfg%run%date)
       write (listLun, '(a,i0)') '# COMPARISON = ', comparison
       write (listLun, '(2a)')   '# MODEL      = ', trim(modelInfo%modelname)
       write (listLun, '(2a)')   '# MODEL FILE = ', trim(modelInfo%filename)
       write (listLun, '(2a)')   '# VARIABLE   = ', trim(modelInfo%variable)
       write (listLun, '(a,i0)') '# TIME INDEX = ', iTime
       write (listLun, '(a,g0)') '# TIME       = ', time
       write (listLun, '(2a)')   '# TIME UNITS = ', trim(model_Tunits)
       write (listLun, '(2a)')   '# OBS DATASET = ', trim(obsInfo%dataset)
       write (listLun, '(2a)')   '# DATASET DESCRIPTION = ', trim(ds_desc)
       write (listLun, '(2a)')   '# MISSION    = ', trim(obsInfo%mission)
       write (listLun, '(2a)')   '# CAMPAIGN   = ', trim(obsInfo%campaign)
       write (listLun, '(2a)')   '# REGION     = ', trim(modelInfo%Region)
       write (listLun, '(a,f0.3)') '# GRID UNITS, km = ', modelInfo%spacing_km
       write (listLun, '(2a)')   '# CREATED BY = cmct_main version ', &
            trim(CMCT_VERSION)

       ! TBD: We should get the column names from datasets and model file
       write (listLun, '(a)') '# NCOLS = 8'
       write (listLun, '(a)') '# COL 1 = latitude, degrees_north'
       write (listLun, '(a)') '# COL 2 = longitude, degrees_east'
       write (listLun, '(a,f0.1,a)') '# COL 3 = x_grid, ', modelInfo%spacing_km,' km'
       write (listLun, '(a,f0.1,a)') '# COL 4 = y_grid, ', modelInfo%spacing_km,' km'
       write (listLun, '(a)') '# COL 5 = observed elevation, m'
       write (listLun, '(a)') '# COL 6 = interpolated model elevation, m'
       write (listLun, '(a)') '# COL 7 = delta elevation (model - observed), m'
       write (listLun, '(a)') '# COL 8 = distance to nearest model node, km'

       write (listLun, '(a)') '# END'  ! Last line in header

       !!
       !! LOOP OVER ALL THE OBSERVATION DATASET RECORDS
       !!
       nRecsRead = 0   ! number of records read
       nBad      = 0   ! failed quality checks
       nOutGrid  = 0   ! outside grid
       nOutModel = 0   ! outside model cells
       nGood     = 0   ! number accepted
       ! distNode_max is expected maximum distance to the nearest node == half
       ! the diagonal size of a cell.  0.0001 fudge in case point is smack in
       ! the center.  maxDistNode is the maximum distance.  nDistMax is the
       ! number of points that exceeded distNodeMax.
       distNode_max = ( modelInfo%spacing_km * sqrt(2.0_REAL64) / 2.0_REAL64 ) &
            + 0.0001_REAL64
       maxDistNode  = 0.0_REAL64
       nDistMax     = 0

       OBSLOOP:  do

          !! Get the next observation.
          call ds%next_obs( lat=obsLat, lon=obsLon, value=obsValue, &
               status=status, verbosity=verbosity )
          select case (status)
          case (DS_NEXTOBS_OK)
             ! Got a good observation
             nRecsRead = nRecsRead + 1
          case (DS_NEXTOBS_DONE)
             ! End of the dataset
             exit OBSLOOP
          case (DS_NEXTOBS_UNCLEAN)
             ! Read a record but it failed the cleanliness check.  This
             ! also filters out ones with invalid elevation values.
             nRecsRead = nRecsRead + 1
             nBad = nBad + 1
             cycle OBSLOOP
          case (DS_NEXTOBS_ERROR)
             ! ds%next_obs couldn't get the next record
             ! Specifc error should have been reported by the reader object
             print *, 'cism_main: Error reading observation, skipping rest'
             exit OBSLOOP
          case default
             print *
             print *, 'cmct_main: Unhandled status from ds: ', status
             print *, 'Stopping.'
             stop 1
          end select

          !! We have a good observation. Process it:
          !! - convert to grid coordinates
          !! - interpolate grid at position of record
          !! - add difference from interpolant to gridded and overall
          !!   moments objects
          !! - write to listOut

          !! Convert observation to grid coordinates.  LatLon2PStXY is in
          !! polar_stereographic_mod.f90.
          call LatLon2PStXY( obsLat, obsLon, status, obsX, obsY )
          if (status /= 0) then
             print *, 'cmct_main: Error from LatLon2PStXY, status=', status
             stop 1
          end if

          !! Reject immediately if outside the grid (we need to calculate which
          !! grid cell it's in).
          if ( (obsX <= xmin) .or. (xmax <= obsX) .or. &
               (obsY <= ymin) .or. (ymax <= obsY) ) then
             nOutGrid = nOutGrid + 1
             cycle OBSLOOP
          end if

          !!  Calculate grid cell within the model the obs point falls in, and
          !!  reject ones that are not surrounded by 4 valid model points.
          !!  Yes, this means that the model is likely to extend a bit further
          !!  than the accepted observation points.  (This is the way Jack did
          !!  it in cism_getdata.pro.).  xfloor, yfloor are now the indices of
          !!  the grid point towards the (-x,-y) direction.  NOTE: Assumes
          !!  that X and Y are monotonically increasing.
          xfloor = floor( obsX - xmin ) + 1
          yfloor = floor( obsY - ymin ) + 1
          cell = model_Data( xfloor:xfloor+1, yfloor:yfloor+1, iTime )
          ! Caution! This probably won't work if model_Invalid is NaN:
          if ( (cell(0,0) == model_Invalid) .or. &
               (cell(1,0) == model_Invalid) .or. &
               (cell(0,1) == model_Invalid) .or. &
               (cell(1,1) == model_Invalid) ) then
             nOutModel = nOutModel + 1
             cycle OBSLOOP
          end if

          nGood = nGood + 1

          !! Bilinear Interpolate obs point into the model grid.  The
          !! formulation is from Abramowitz and Stegun, eqn. 25.2.66.
          p = ( obsX - model_X(xfloor) ) / ( model_X(xfloor+1) - model_X(xfloor) )
          q = ( obsY - model_Y(yfloor) ) / ( model_Y(yfloor+1) - model_Y(yfloor) )
          f = (1.0_REAL64 - p) * (1.0_REAL64 - q) * cell(0,0) + &
               p  * (1.0_REAL64 - q) * cell(1,0) + &
               (1.0_REAL64 - p) *               q  * cell(0,1) + &
               p  *               q  * cell(1,1)
          delta = f - obsValue      ! model - obs
          call diff_Moms%add_point( delta, xfloor, yfloor )
          call diff_MomsTot%add_point( delta, 1, 1 )
          call hist1m%accumulate( real(delta) )    ! histogram_class uses reals

          if (nGood == 1) then
             maxDelta = delta
             minDelta = delta
          else
             maxDelta = max( delta, maxDelta )
             minDelta = min( delta, minDelta )
          end if

          !! Calculate distance (in same units as modelInfo%spacing_km) from
          !! observation point to nearest grid point.  By now, we know the
          !! observation is surrounded by 4 good model points, so we can simply
          !! round obsx/obsy to their nearest integers and take the distance to
          !! there.
          distNode = sqrt( ((obsx-anint(obsx))**2) + ((obsy-anint(obsy))**2) ) &
               * modelInfo%spacing_km
          maxDistNode = max( maxDistNode, distNode )
          if ( distNode .gt. distNode_max ) then
             nDistMax = nDistMax + 1
          end if

          !! write per-record output
          listFmt = '(2(f13.8,2x), 2(f10.4,2x), x, 3(f9.3,2x), x, f6.3)'
          write (listLun, fmt=listFmt)  obsLat, obsLon, obsX, obsY, &
               obsValue, f, delta, distNode

       end do OBSLOOP

       close (listLun)

       print *
       print '(a,f9.3)', &
            ' cmct_main: Max distance of an obs to closest model node, km: ', &
            maxDistNode
       print '(a,f0.3,a,i0)', &
            ' cmct_main: Number farther from node than expected max (', &
            distNode_max, ' km):  ', nDistMax

       !! Get mean and std dev arrays
       call diff_Moms%get( n=nPts, mean=diff_Means, var=diff_Vars, &
            stddev=diff_SDs, invflag=diff_InvFlag )
       call diff_MomsTot%get( n=nPtsTot, mean=diff_MeanTot, var=diff_VarTot, &
            stddev=diff_SDTot )

       !! Statistics
       print *
       print *, 'cmct_main:  SUMMARY STATISTICS'
       print '(a,i10)', ' Observation records read:  ', nRecsRead
       print '(a,i10)', '    Records failed quality: ', nBad
       print '(a,i10)', '    Records outside grid:   ', nOutGrid
       print '(a,i10)', '    Records outside model:  ', nOutModel
       print '(a,i10)', '    Records acceptable:     ', nGood
       print '(a,i10)', ' Total observations used:   ', nPtsTot
       print *
       print '(a,f9.3)', &
            ' Overall Mean of Delta (Interp model - Obs), m: ', &
            diff_MeanTot
       print '(a,f9.3)', &
            ' Overall Std Dev of Delta (Interp model - Obs): ', &
            diff_SDTot
       print '(a,f0.3,a,f0.3)', &
            ' Min Delta: ', minDelta, '   Max Delta: ', maxDelta

       !!
       !! WRITE THE GRID OUTPUTS TO A NETCDF FILE.  (This might be better
       !! all separated off into a subroutine in a module.)  Sophie
       !! requested NetCDF-4 (HDF5) format.  Standardized on
       !! NF90_FILL_DOUBLE as flag for invalid values.  TBD: Make this more
       !! CF-compliant.
       !!
       print *
       meanOut = trim(baseOut) // '.meansdgrids.nc'
       print *, 'cmct_main: Writing netCDF MeanSDGrids file: '// trim(meanOut)
       print *, 'cmct_main: .meansdgrids.nc file invalid flag value: ', &
            NF90_FILL_DOUBLE

       status = nf90_create( path=meanOut, cmode=ior(NF90_NETCDF4,NF90_NOCLOBBER),&
            ncid=meanNCID )
       call check_ncdf( status, 'creating '//trim(meanOut) )

       ! X and Y dimensions.
       status = nf90_def_dim( meanNCID, 'x_grid', szXout, xdimID )
       call check_ncdf( status, 'defining X_grid dim' )

       status = nf90_def_dim( meanNCID, 'y_grid', szYout, ydimID )
       call check_ncdf( status, 'defining Y_grid dim' )

       ! Global attributes.

       status = nf90_put_att( meanNCID, NF90_GLOBAL, 'runid', &
            trim(runcfg%run%runid) )
       call check_ncdf( status, 'putting runid attribute' )

       status = nf90_put_att( meanNCID, NF90_GLOBAL, 'title', &
            trim(runcfg%run%user_run_title) )
       call check_ncdf( status, 'putting title attribute' )

       status = nf90_put_att( meanNCID, NF90_GLOBAL, 'date', &
            trim(runcfg%run%date) )
       call check_ncdf( status, 'putting date attribute' )

       status = nf90_put_att( meanNCID, NF90_GLOBAL, 'creator', &
            'cmct_main version ' // trim(CMCT_VERSION) )
       call check_ncdf( status, 'putting creator attribute' )


       status = nf90_put_att( meanNCID, NF90_GLOBAL, 'model_name', &
            trim(modelInfo%modelname) )
       call check_ncdf( status, 'putting model_name attribute' )

       status = nf90_put_att( meanNCID, NF90_GLOBAL, 'model_filename', &
            trim(modelInfo%filename) )
       call check_ncdf( status, 'putting model_filename attribute' )

       status = nf90_put_att( meanNCID, NF90_GLOBAL, 'model_variable', &
            trim(modelInfo%variable) )
       call check_ncdf( status, 'putting model_variable attribute' )

       status = nf90_put_att( meanNCID, NF90_GLOBAL, 'model_time_index', iTime )
       call check_ncdf( status, 'putting model_time_index attribute' )

       status = nf90_put_att( meanNCID, NF90_GLOBAL, 'model_time', time )
       call check_ncdf( status, 'putting model_time attribute' )

       status = nf90_put_att( meanNCID, NF90_GLOBAL, &
            'model_time_units', trim(model_Tunits) )
       call check_ncdf( status, 'putting model_time_units attribute' )

       status = nf90_put_att( meanNCID, NF90_GLOBAL, 'model_usercomments', &
            trim(modelInfo%user_comments) )
       call check_ncdf( status, 'putting model_usercomments attribute' )

       ! TBD: Some of these attributes are icesat-glas specific.
       status = nf90_put_att( meanNCID, NF90_GLOBAL, 'observation_dataset', &
            trim(obsInfo%dataset) )
       call check_ncdf( status, 'putting observation_dataset attribute' )

       status = nf90_put_att( meanNCID, NF90_GLOBAL, 'dataset_description', &
            trim(ds_desc) )
       call check_ncdf( status, 'putting dataset_description attribute' )

       status = nf90_put_att( meanNCID, NF90_GLOBAL, 'mission', &
            trim(obsInfo%mission) )
       call check_ncdf( status, 'putting mission attribute' )

       status = nf90_put_att( meanNCID, NF90_GLOBAL, 'campaign', &
            trim(obsInfo%campaign) )
       call check_ncdf( status, 'putting campaign attribute' )

       status = nf90_put_att( meanNCID, NF90_GLOBAL, 'region', &
            trim(modelInfo%Region) )
       call check_ncdf( status, 'putting region attribute' )


       ! Define the variables.  Set _FillValue attributes to NetCDF default.
       ! Note: nf90_def_var_fill doesn't work in this version of netcdf.

       ! x_grid, y_grid are the coordinate variables
       status = nf90_def_var( meanNCID, 'x_grid', NF90_DOUBLE, xdimID, xVID )
       call check_ncdf( status, 'defining var x_grid')
       status = nf90_put_att( meanNCID, xVID, 'standard_name', &
            'projection_x_coordinate' )
       call check_ncdf( status, 'setting standard_name for x_grid')

       status = nf90_def_var( meanNCID, 'y_grid', NF90_DOUBLE, ydimID, yVID )
       call check_ncdf( status, 'defining var y_grid')
       status = nf90_put_att( meanNCID, yVID, 'standard_name', &
            'projection_y_coordinate' )
       call check_ncdf( status, 'setting standard_name for y_grid')

       ! Latitude grid
       status = nf90_def_var( meanNCID, 'latitude', NF90_DOUBLE, [xdimID,ydimID], &
            latsVID )
       call check_ncdf( status, 'defining var latitude')
       status = nf90_put_att( meanNCID, latsVID, 'units', 'degrees_north' )
       call check_ncdf( status, 'setting units for latitude' )
       status = nf90_put_att( meanNCID, latsVID, 'standard_name', 'latitude' )
       call check_ncdf( status, 'setting standard_name for latitude' )
       status = nf90_put_att( meanNCID, latsVID, '_FillValue', NF90_FILL_DOUBLE )
       call check_ncdf( status, 'setting fill for latitude' )

       ! Longitude grid
       status = nf90_def_var( meanNCID, 'longitude', NF90_DOUBLE, [xdimID,ydimID], &
            lonsVID )
       call check_ncdf( status, 'defining var longitude')
       status = nf90_put_att( meanNCID, lonsVID, 'units', 'degrees_east' )
       call check_ncdf( status, 'setting units for longitude' )
       status = nf90_put_att( meanNCID, lonsVID, 'standard_name', 'longitude' )
       call check_ncdf( status, 'setting standard_name for longitude' )
       status = nf90_put_att( meanNCID, lonsVID, '_FillValue', NF90_FILL_DOUBLE )
       call check_ncdf( status, 'setting fill for longitude' )

       ! Number grid
       status = nf90_def_var( meanNCID, 'number', NF90_INT, [xdimID,ydimID], nVID )
       call check_ncdf( status, 'defining var number')

       ! Means grid
       status = nf90_def_var( meanNCID, 'means', NF90_DOUBLE, [xdimID,ydimID], &
            meansVID )
       call check_ncdf( status, 'defining var means')
       status = nf90_put_att( meanNCID, meansVID, 'units', 'meters' )
       call check_ncdf( status, 'setting units for means' )
       status = nf90_put_att( meanNCID, meansVID, '_FillValue', NF90_FILL_DOUBLE )
       call check_ncdf( status, 'setting fill for means' )

       ! Std_devs grid
       status = nf90_def_var( meanNCID, 'std_devs', NF90_DOUBLE, [xdimID,ydimID], &
            sdVID )
       call check_ncdf( status, 'defining var std_devs')
       status = nf90_put_att( meanNCID, sdVID, 'units', 'meters' )
       call check_ncdf( status, 'setting units for std_devs' )
       status = nf90_put_att( meanNCID, sdVID, '_FillValue', NF90_FILL_DOUBLE )
       call check_ncdf( status, 'setting fill for std_devs' )

       status = nf90_enddef( meanNCID )
       call check_ncdf( status, 'enddef' )

       ! Write the variables.  Invalid values of diff_Means and diff_SDs
       ! are set to netCDF default fill.  TBD: Should this be done for
       ! model_lats, model_lons too?  (Although currently
       ! model_cism_txt_mod uses that value anyway.)

       status = nf90_put_var( meanNCID, xVID, model_X(1:szXout), [1], [szXout] )
       call check_ncdf( status, 'putting x_grid' )

       status = nf90_put_var( meanNCID, yVID, model_Y(1:szYout), [1], [szYout] )
       call check_ncdf( status, 'putting y_grid' )


       status = nf90_put_var( meanNCID, latsVID, model_lats(1:szXout,1:szYout), &
            [1,1], [szXout,szYout] )
       call check_ncdf( status, 'putting latitude' )

       status = nf90_put_var( meanNCID, lonsVID, model_lons(1:szXout,1:szYout), &
            [1,1], [szXout,szYout] )
       call check_ncdf( status, 'putting longitude' )


       status = nf90_put_var( meanNCID, nVID, nPts, [1,1], [szXout,szYout] )
       call check_ncdf( status, 'putting number' )

       where (diff_Means == diff_InvFlag)  diff_Means = NF90_FILL_DOUBLE
       status = nf90_put_var( meanNCID, meansVID, diff_Means, &
            [1,1], [szXout,szYout] )
       call check_ncdf( status, 'putting diff_Means' )

       where (diff_SDs == diff_InvFlag)  diff_SDs = NF90_FILL_DOUBLE
       status = nf90_put_var( meanNCID, sdVID, diff_SDs, &
            [1,1], [szXout,szYout] )
       call check_ncdf( status, 'putting diff_SDs' )


       status = nf90_close( meanNCID )
       call check_ncdf( status, 'closing meanNCID' )

       !!
       !! WRITE HISTOGRAM FILE
       !!
       histOut = trim(baseOut) // '.histogram.txt'
       print *
       print *, 'cmct_main: Writing histogram file: '// trim(histOut)

       call hist1m%get( min=histmin, max=histmax, binsize=histbinsize, &
            nbins=histnbins, bins=histbins, histogram=histogram, &
            n=histtotal, nbelow=histbelow, nabove=histabove, &
            pdf=histpdf, cdf=histcdf )

       open ( newunit=histLun, file=histOut, status='new', form='formatted', &
            action='write' )    ! TBD: need to handle errors

       ! Header.  Yes, many of these are the same as the reclist file.  TBD:
       ! Should probably put these in a procedure.  TBD: Some lines are
       ! currently icesat-specific
       write (histLun, '(a)')    '# CMCT HISTOGRAM'
       write (histLun, '(2a)')   '# RUNID      = ', trim(runcfg%run%runid)
       write (histLun, '(2a)')   '# TITLE      = ', trim(runcfg%run%user_run_title)
       write (histLun, '(2a)')   '# DATE       = ', trim(runcfg%run%date)
       write (histLun, '(a,i0)') '# COMPARISON = ', comparison
       write (histLun, '(2a)')   '# MODEL      = ', trim(modelInfo%modelname)
       write (histLun, '(2a)')   '# MODEL FILE = ', trim(modelInfo%filename)
       write (histLun, '(2a)')   '# VARIABLE   = ', trim(modelInfo%variable)
       write (histLun, '(a,i0)') '# TIME INDEX = ', iTime
       write (histLun, '(a,g0)') '# TIME       = ', time
       write (histLun, '(2a)')   '# TIME UNITS = ', trim(model_Tunits)
       write (histLun, '(2a)')   '# OBS DATASET = ', trim(obsInfo%dataset)
       write (histLun, '(2a)')   '# DATASET DESCRIPTION = ', trim(ds_desc)
       write (histLun, '(2a)')   '# MISSION    = ', trim(obsInfo%mission)
       write (histLun, '(2a)')   '# CAMPAIGN   = ', trim(obsInfo%campaign)
       write (histLun, '(2a)')   '# REGION     = ', trim(modelInfo%Region)
       write (histLun, '(2a)')   '# CREATED BY = cmct_main version ', &
            trim(CMCT_VERSION)

       write (histLun, '(a,f0.2)') '# HISTOGRAM MIN, m = ', histmin
       write (histLun, '(a,f0.2)') '# HISTOGRAM MAX, m = ', histmax
       write (histLun, '(a,i0)')   '# HISTOGRAM NUMBER BINS = ', histnbins
       write (histLun, '(a,f0.2)') '# HISTOGRAM BINSIZE, m  = ', histbinsize
       write (histLun, '(a,i0)')   '# HISTOGRAM TOTAL IN RANGE = ', histtotal
       write (histLun, '(a,i0)')   '# HISTOGRAM BELOW RANGE    = ', histbelow
       write (histLun, '(a,i0)')   '# HISTOGRAM ABOVE RANGE    = ', histabove
       write (histLun, '(a)')      '# NOTE bins lowerbound <= x < upperbound'

       write (histLun, '(a)') '# NCOLS = 5'
       write (histLun, '(a)') '# COL 1 = lower bin boundary, m'
       write (histLun, '(a)') '# COL 2 = upper bin boundary, m'
       write (histLun, '(a)') '# COL 3 = counts'
       write (histLun, '(a)') '# COL 4 = probability density function'
       write (histLun, '(a)') '# COL 5 = empirical cumulative distribution function'

       write (histLun, '(a)') '# END'  ! Last line in header

       do i = lbound(histogram,1), ubound(histogram,1)
          write (histLun, '(2(f8.2,2x), i8, 3x, e12.5, 2x, e12.5)')  &
               histbins(i), histbins(i+1), histogram(i), histpdf(i), histcdf(i)
       end do

       close (histLun)

     end associate

  end do COMPLOOP

  !!
  !! WRAP UP
  !!
  print *
  print '( x, 65("=") )'
  print *
  call cpu_time( cputime_last )
  print '(a,g10.5,a)', ' cmct_main: Elapsed CPU: ', &
       cputime_last - cputime_1, ' sec'
  print *, 'cmct_main: Done'
  stop 0   ! normal exit

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Internal procedure check_ncdf
!!! Checks the status of netcdf calls.  Written as an internal
!!! procedure so it has access to all of the host variables.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine check_ncdf ( status, msg )
  integer, intent(in) :: status
  character(len=*), intent(in) :: msg

  if ( status /= nf90_NoErr ) then
     print *
     print *, 'cmct_main: NetCDF error ' // msg
     print *, '   nf90 status = ', status
     print *, '   ', trim( nf90_strerror(status) )
     stop 1
  end if
  return

end subroutine check_ncdf

end program
