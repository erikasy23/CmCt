!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! NAME: datasets_mod.f08
!!!
!!! PURPOSE: Module for accessing the observation datasets.  Reads
!!! datasets.json and maintains a list of available observation datasets
!!! and their specific information (currently just their specific
!!! configuration files).  Provides an object class (datasets_class) for
!!! accessing the datasets in a generic way, accessing the appropriate
!!! specific dataset class.  It should be OK to reuse objects of the
!!! datasets_class for new datasets by calling %open again.
!!!
!!! USAGE: Schematic of how this module may be used. See procedure prologs
!!! for details and options.
!!!   type(datasets_class) :: ds                          ! in declarations
!!!   call datasets_init                                  ! call once
!!!   call datasets_info( dataset, configfile, desc, status) ! call as needed
!!!   call ds%open( mission, dataset, campaign, status )  ! open new dataset
!!!   obsloop: do                                         ! loop on observations
!!!      call ds%next_obs(lat, lon, value, status)        ! get next obs in ds
!!!      if (status /= DS_NEXTOBS_OK) then exit obsloop   ! couldn't get obs
!!!      [process this observation]                       ! did get obs
!!!   end do obsloop
!!!
!!! AUTHOR: Jeff Guerber, GSFC 615/SigmaSpace, 2015-08-20
!!! 2015-11-27 JRG: Added verbosity.
!!! 2016-05-11 JRG: datasets_init: If sets is already allocated, deallocate.
!!! 2016-06-07 JRG: Extended to add datasets_class and its methods open(),
!!!    next_obs(), etc.
!!! 2016-06-20 JRG: Changed CMCT netCDF datasets to icesat-cmctnc-*,
!!!    separate from flat-binary icesat-cmct-* (not currently supported
!!!    here, but there's been talk...).  icesat_cmctnc_class%clean() now
!!!    needs arg specifying which elevation to check.
!!! 2016-07-02 JRG: Bug fixes: in open(), misspelled one of the dataset
!!!    names, and error message wrote mission instead of dataset.
!!! 2016-07-13 JRG: Read and store the dataset description strings (already
!!!    in datasets.json) and have datasets_info return them.
!!!    Open() returns status constants instead of halting on errors, so
!!!    cmct_main can go on to next comparison.
!!!    Provide length constants for the strings in type dsinfo.
!!!    Moved dsfilep inside datasets_init, don't need for rest of module.
!!!    General cleanup esp. in datasets_init and open.
!!!
!!! Last SVN commit: $Id: datasets_mod.f08 116 2016-07-18 22:37:33Z jguerber $
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE datasets_mod

  use icesat_cism_class_mod
  use icesat_cism_rec_mod
  use icesat_cmctnc_class_mod
  use icesat_cmct_rec_mod
  use fson
  use fson_value_m, only: fson_value_count, fson_value_get
  use iso_fortran_env, only: REAL32, REAL64

  implicit none
  private
  ! Note: The object methods (open, next_obs) must be private here but
  ! public in the type definition.
  public :: datasets_init, datasets_info, datasets_class
  public :: DS_OPEN_OK, DS_OPEN_ERROR
  public :: DS_NEXTOBS_OK, DS_NEXTOBS_UNCLEAN, DS_NEXTOBS_DONE, DS_NEXTOBS_ERROR
  public :: DS_NAMELEN, DS_DESCLEN, DS_CONFIGFILELEN   ! char string len consts

  !! List of information about the datasets, read from DSFILENAME
  character(len=*), parameter :: DSFILENAME = 'datasets.json'
  integer, parameter :: DS_NAMELEN = 50    ! max len of dataset name
  integer, parameter :: DS_DESCLEN = 150   ! max len of dataset description
  integer, parameter :: DS_CONFIGFILELEN = 150  ! max len of configfile name
  type dsinfo
     character(len=DS_NAMELEN)       :: name       ! name of the dataset
     character(len=DS_DESCLEN)       :: desc       ! dataset description string
     character(len=DS_CONFIGFILELEN) :: configfile ! name of its config file
  end type dsinfo
  type(dsinfo), dimension(:), allocatable :: sets  ! list of dataset info

  !!
  !! TYPE FOR THE DATASET-ACCESS CLASS.
  !!
  type datasets_class

     private
     character(len=25) :: mission       ! mission name
     character(len=DS_NAMELEN) :: dataset       ! name of the dataset
     character(len=4)  :: icesat_campaign !for mission=icesat-glas, the campaign

     !! A reader object for each *type* of dataset supported.  Will only
     !! use one of these, per object instance:
     type( icesat_cism_class )   :: icesat_cism
     type( icesat_cmctnc_class ) :: icesat_cmctnc

     contains
       ! The procedures that are type-bound to this type, which comprise
       ! the methods of the class:
       procedure :: open
       procedure :: next_obs

    end type datasets_class

    !! OPEN() STATUS CODES
    integer, parameter :: DS_OPEN_OK    = 0   ! Open successful
    integer, parameter :: DS_OPEN_ERROR = 1   ! Error opening dataset

    !! NEXT_OBS() STATUS CODES
    integer, parameter :: DS_NEXTOBS_OK    = 0    ! Observation is good
    integer, parameter :: DS_NEXTOBS_UNCLEAN = 1  ! Obs failed clean()
    integer, parameter :: DS_NEXTOBS_DONE  = 2    ! No more observations avail
    integer, parameter :: DS_NEXTOBS_ERROR = 3    ! Could not get an observation

CONTAINS

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! NAME: datasets_init
  !!
  !! PURPOSE: Read the datasets.json file and store its information in the
  !! module for use by the object and retrieval through datasets_info.
  !!
  !! ARGUMENT:
  !!    verbosity: Opt Input. Verbosity level, integer. Defaults to 0.
  !!
  !!  Note that if fson fails, it (usually) just halts with its own error
  !!  message.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE datasets_init( verbosity )

    ! Dummy argument
    integer, intent(in), optional :: verbosity

    ! Local variables
    type(fson_value), pointer :: dsfilep  ! fson parse of file
    type(fson_value), pointer :: knownp   ! fson ptr to file's list of
                                          !    known datasets
    type(fson_value), pointer :: fp       ! general fson ptr
    integer :: i, nSets       ! nSets = number of known datasets
    integer :: verb           ! local copy of verbosity dummy arg

    !!---------------------------------------------------------------

    ! Set verbosity locally rather than for object, so can use different
    ! values when calling other methods.
    verb = 0
    if (present(verbosity))  verb = verbosity

    if (verb >= 1) then
       print *
       print *, 'datasets_init: Reading known-datasets info from ', &
            trim(DSFILENAME)
    end if
    dsfilep => fson_parse( DSFILENAME )

    !! Get the names of the known datasets from the known_datsets section
    !! of the datasets.json file.  Fson can't extract string arrays or
    !! strings directly from an array, so have to march a pointer through
    !! the array and get them one at a time.
    call fson_get( dsfilep, 'known_datasets', knownp )
    nSets = fson_value_count( knownp )
    if (verb >= 1) print '(a,i0)', 'datasets_init: Known datasets: ', nSets

    if ( allocated( sets ) )  deallocate( sets )
    allocate ( sets(nSets) )
    do i = 1, nSets
       fp => fson_value_get( knownp, i )
       call fson_get( this=fp, value=sets(i)%name )
       if (verb >= 1)  print '(5x,i0,2a)', i, ': dataset name = ', &
            trim(sets(i)%name)
    end do

    !! Now look up the information for each dataset.  datasets.json has a
    !! separate section for each one.
    do i = 1, nSets
       call fson_get( dsfilep, trim(sets(i)%name)//'.desc', &
            sets(i)%desc )
       if (verb >= 1)  print '(5x,3a)', trim(sets(i)%name), &
            '.desc = ', trim(sets(i)%desc)
       call fson_get( dsfilep, trim(sets(i)%name)//'.configfile', &
            sets(i)%configfile )
       if (verb >= 1)  print '(5x,3a)', trim(sets(i)%name), &
            '.configfile = ', trim(sets(i)%configfile)
    end do

    call fson_destroy( dsfilep )

    return
  END SUBROUTINE datasets_init

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! NAME: datasets_info
  !!
  !! PURPOSE: Return the information about the named dataset.
  !!
  !! TBD: May need to add checks for which dataset, because they could
  !! have different items.
  !!
  !! ARGUMENTS:
  !!    dataset: Input. Name of the dataset to get information about.
  !!    configfile: Opt Output. Name of the configuration file.
  !!    desc: Opt Output. Description string for this dataset,
  !!       from datasets.json
  !!    status: Output. 0=OK, found dataset; 1=dataset not found in known list
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE datasets_info ( dataset, configfile, desc, status )

    ! Dummy arguments
    character(len=*), intent(in) :: dataset
    character(len=*), intent(out), optional :: configfile, desc
    integer, intent(out) :: status  ! 0=OK, 1=dataset not known

    ! Local variables
    integer :: i

    !!---------------------------------------------------------------

    do i = 1, size( sets )
       if ( trim( sets(i)%name ) .eq. trim( dataset ) ) then
          status = 0

          if ( present(configfile) ) then
             configfile = sets(i)%configfile
          end if

          if ( present(desc) ) then
             desc = sets(i)%desc
          end if

          !! Add future dataset items here (and in the dummy argument list).

          return
       end if
    end do

    !! If we got to here, we didn't find the dataset in our list
    print *
    print *, 'datasets_info: Dataset not found in ' // trim(DSFILENAME) &
         // ': ' // trim(dataset)
    status = 1
    return
  END SUBROUTINE datasets_info


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! NAME: open
  !!
  !! PURPOSE: Open a particular dataset in this object instance.  Object
  !!    method, type-bound to type(datasets_class).  Should be OK to call
  !!    this again to open the object instance on a different dataset.
  !!
  !! NOTE: Currently just stops on an error.
  !!
  !! ARGUMENTS:
  !!    mission:  Input. Mission name.
  !!    dataset:  Input. Dataset name.
  !!    icesat_campaign:  Input. ICEsat campaign. Only for mission icesat-glas.
  !!    status:  Output. Status of the operation, DS_OPEN_OK or DS_OPEN_ERROR
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE open( self, mission, dataset, icesat_campaign, status, verbosity )

    !! Dummy arguments
    ! self is the class variable, don't include as an actual argument
    class( datasets_class ), intent(inout) :: self
    character(len=*), intent(in) :: mission
    character(len=*), intent(in) :: dataset
    character(len=*),  intent(in), optional :: icesat_campaign
    integer, intent(out) :: status
    integer, intent(in), optional :: verbosity


    !! Local variables
    character(len=150) :: dsCfgFile   ! dataset-specific configuration file
    integer :: dsstatus       ! status from the datasets
    integer :: verb           ! local copy of verbosity dummy arg

    !!---------------------------------------------------------------

    verb = 0
    if (present(verbosity))  verb = verbosity

    if ( .not. allocated(sets) ) then
       print *, 'datasets_mod%open: Module has not been initialized, need ', &
            'to call datasets_init. Stopping.'
       stop 1      ! FATAL
    end if

    self%mission = trim(mission)

    !! Get configuration file name for this dataset.
    call datasets_info( dataset=dataset, configfile=dsCfgFile, status=dsstatus )
    if (dsstatus /= 0) then  ! currently Unknown is only error from datasets_info
       print *, 'datasets_mod%open: Unknown dataset requested: ', &
            trim(dataset)
       status = DS_OPEN_ERROR
       return
    end if
    self%dataset = dataset
    self%icesat_campaign = ''   ! will be reset below if needed


    SEL_MISSION: select case ( trim(self%mission) )

    !! MISSION ICESAT-GLAS
    case ( 'icesat-glas' ) SEL_MISSION

       if (present(icesat_campaign)) then
          self%icesat_campaign = icesat_campaign
       else
          print *, 'datasets_mod%open: Campaign needed for mission icesat-glas'
          status = DS_OPEN_ERROR
          return
       end if

       !!
       !! Initialize the appropriate dataset object.
       !! Note that the case selectors must be elements of known_datasets
       !! in datasets.json.
       !!
       SEL_DS: select case ( trim(self%dataset) )

       !!
       !! DATASET: ICESAT-CISM-ELEV-GRN
       !!
       case ( 'icesat-cism-elev-grn' ) SEL_DS
          call self%icesat_cism%init( configFile=dsCfgFile, &
               campaign=self%icesat_campaign, &
               status=dsstatus, verbosity=verb )
          if ( dsstatus /= ICM_INIT_OK ) then
             print *, 'datasets_mod%open: Error initializing icesat_cism_class, status=', &
                  dsstatus
             status = DS_OPEN_ERROR
             return
          end if

       !!
       !! DATASETS: ICESAT-CMCTNC-*-GRN
       !!
       case ( "icesat-cmctnc-tpelev-grn", &
            "icesat-cmctnc-wgselev-grn", &
            "icesat-cmctnc-egmmtelev-grn", &
            "icesat-cmctnc-egmtfelev-grn" ) SEL_DS
          call self%icesat_cmctnc%init( configFile=dsCfgFile, &
               campaign=self%icesat_campaign, &
               status=dsstatus, verbosity=verb )
          if ( dsstatus /= ICT_INIT_OK ) then
             print *, 'datasets_mod%open: Error initializing icesat_cmctnc_class, status=', &
                  dsstatus
             status = DS_OPEN_ERROR
             return
          end if

       !!
       !! Other datasets for this mission will go here, in new SEL_DS case
       !! blocks.
          !!
          !! This is where I am adding the datasets for Antarctica for the ICESat Glas mission
          !! SO the only thing I am changing here is the dataset files that are for the case
          !!

       !!
       !! DATASETS: ICESAT-CMCTNC-*-ANT
       !!
       case ( "icesat-cmctnc-tpelev-ant" ) SEL_DS
          call self%icesat_cmctnc%init( configFile=dsCfgFile, &
               campaign=self%icesat_campaign, &
               status=dsstatus, verbosity=verb )
          if ( dsstatus /= ICT_INIT_OK ) then
             print *, 'datasets_mod%open: Error initializing icesat_cmctnc_class, status=', &
                  dsstatus
             status = DS_OPEN_ERROR
             return
          end if

       case default SEL_DS
          print *, 'datasets_mod%open: Unknown dataset = ' // trim(self%dataset)
          status = DS_OPEN_ERROR
          return

       end select SEL_DS

    !!
    !! Other missions will go here, in new SEL_MISSION case blocks.
    !! This is where I will add the ERS 1 and ERS 2 datasets.
    !!   

    case default SEL_MISSION
       print *, 'datasets_mod%open: Unknown mission=' // trim(self%mission)
       status = DS_OPEN_ERROR
       return

    end select SEL_MISSION

    !!  Success!
    status = DS_OPEN_OK
    return
  END SUBROUTINE open

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! NAME: next_obs
  !!
  !! PURPOSE: Returns the next observation from the appropriate dataset for
  !!    this object: lat, lon, value (currently elevation).  Type-bound to
  !!    type(datasets_class).
  !!
  !! ARGUMENTS:
  !!    lat, lon: Output. Latitude and longitude
  !!    value: Output. Value of the observation, depends on dataset.
  !!    status: Output. See possible values in header.
  !!    verbosity: Opt Input.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE next_obs ( self, lat, lon, value, status, verbosity )

    !! Dummy arguments
    ! self is the class variable, don't include as an actual argument
    class( datasets_class ), intent(inout) :: self
    real(kind=REAL64), intent(out) :: lat, lon
    real(kind=REAL32), intent(out) :: value
    integer, intent(out) :: status
    integer, intent(in), optional :: verbosity

    !! Local variables
    ! Specific records for each dataset type
    type(icesat_cism_rec) :: icism_rec
    type(icesat_cmct_rec) :: icmct_rec
    integer :: nrstatus       ! Status from next_rec calls
    integer :: recQual        ! result of quality check
    integer :: verb           ! local copy of verbosity dummy arg

    !!---------------------------------------------------------------

    verb = 0
    if (present(verbosity))  verb = verbosity

    select case ( trim(self%dataset) )

    !!
    !! DATASET: ICESAT-CISM-ELEV-GRN
    !!
    case ( 'icesat-cism-elev-grn' )
       icism_rec = self%icesat_cism%next_rec( status=nrstatus, verbosity=verb )
       if (nrstatus == ICM_NEXT_REC_DONE) then
          ! No more records in the dataset
          status = DS_NEXTOBS_DONE
          return
       else if (nrstatus /= ICM_NEXT_REC_OK) then
          ! Some other error, we couldn't get a record
          status = DS_NEXTOBS_ERROR
          return
       end if

       ! Else, we got a valid observation (though it may yet fail clean())
       lat   = icism_rec%lat_degN
       lon   = icism_rec%lon_degE
       value = icism_rec%elev_m

       ! Clean data, eg. check if points have flag values or are out of
       ! reasonable elevation range.  (This dataset has no invalid values.)
       ! See icesat_cism%clean for specifics.  Just set status and let
       ! caller decide what to do with it.  TBD: Should we return recQual?
       recQual = self%icesat_cism%clean( rec=icism_rec, verbosity=verbosity )
       if (recQual == ICM_RECORD_GOOD) then
          status = DS_NEXTOBS_OK           ! passed the clean() checks
       else
          status = DS_NEXTOBS_UNCLEAN      ! no it didn't
       end if

       return

    !!
    !! DATASETS: ICESAT-CMCTNC-*-GRN
    !! From the ICESat-GLAS CMCT Greenland datasets in netCDF format
    !!
    case ( "icesat-cmctnc-tpelev-grn", &
         "icesat-cmctnc-wgselev-grn", &
         "icesat-cmctnc-egmmtelev-grn", &
         "icesat-cmctnc-egmtfelev-grn" )
       icmct_rec = self%icesat_cmctnc%next_rec( status=nrstatus, verbosity=verb )
       if (nrstatus == ICT_NEXT_REC_DONE) then
          ! No more records in the dataset
          status = DS_NEXTOBS_DONE
          return
       else if (nrstatus /= ICT_NEXT_REC_OK) then
          ! Some other error, we couldn't get a record
          status = DS_NEXTOBS_ERROR
          return
       end if

       ! We got a valid observation.  Now set the return value and clean
       ! the data, eg. check if point has flag values such as invalid, or
       ! is out of reasonable elevation range.  See icesat_cmctnc%clean for
       ! specifics.  Just set status and let caller decide what to do with
       ! it.  TBD: Should we return recQual?
       lat = icmct_rec%lat_degN
       lon = icmct_rec%lon_degE
       select case (self%dataset)
       case ( "icesat-cmctnc-tpelev-grn" )
          value = icmct_rec%TopexElev_m
          recQual = self%icesat_cmctnc%clean( rec=icmct_rec,  &
               fieldname='topexelev_m', verbosity=verb )

       case ( "icesat-cmctnc-wgselev-grn" )
          value = icmct_rec%WGS84Elev_m
          recQual = self%icesat_cmctnc%clean( rec=icmct_rec, &
               fieldname='wgs84elev_m', verbosity=verb )

       case ( "icesat-cmctnc-egmmtelev-grn" )
          recQual = self%icesat_cmctnc%clean( rec=icmct_rec, &
               fieldname='egm08meantideelev_m', verbosity=verb )
          value = icmct_rec%EGM08MeanTideElev_m

       case ( "icesat-cmctnc-egmtfelev-grn" )
          value = icmct_rec%EGM08TideFreeElev_m
          recQual = self%icesat_cmctnc%clean( rec=icmct_rec, &
               fieldname='egm08tidefreeelev_m', verbosity=verb )

       end select

       if (recQual == ICT_RECORD_GOOD) then
          status = DS_NEXTOBS_OK           ! passed the clean() checks
       else
          status = DS_NEXTOBS_UNCLEAN      ! no it didn't
       end if

       return

       !!
       !! New datasets go here, in new case blocks
       !! Below I put the setup for the ICESat Glas Antarctic data.
       !! Technically this should work easily just the same way as the Greenland data.
       !! 
       
    !!
    !! DATASETS: ICESAT-CMCTNC-*-ANT
    !! From the ICESat-GLAS CMCT Greenland datasets in netCDF format
    !!
    case ( "icesat-cmctnc-tpelev-ant" )
       icmct_rec = self%icesat_cmctnc%next_rec( status=nrstatus, verbosity=verb )
       if (nrstatus == ICT_NEXT_REC_DONE) then
          ! No more records in the dataset
          status = DS_NEXTOBS_DONE
          return
       else if (nrstatus /= ICT_NEXT_REC_OK) then
          ! Some other error, we couldn't get a record
          status = DS_NEXTOBS_ERROR
          return
       end if

       ! We got a valid observation.  Now set the return value and clean
       ! the data, eg. check if point has flag values such as invalid, or
       ! is out of reasonable elevation range.  See icesat_cmctnc%clean for
       ! specifics.  Just set status and let caller decide what to do with
       ! it.  TBD: Should we return recQual?
       lat = icmct_rec%lat_degN
       lon = icmct_rec%lon_degE
       select case (self%dataset)
       case ( "icesat-cmctnc-tpelev-ant" )
          value = icmct_rec%TopexElev_m
          recQual = self%icesat_cmctnc%clean( rec=icmct_rec,  &
               fieldname='topexelev_m', verbosity=verb )

       end select

       if (recQual == ICT_RECORD_GOOD) then
          status = DS_NEXTOBS_OK           ! passed the clean() checks
       else
          status = DS_NEXTOBS_UNCLEAN      ! no it didn't
       end if

       return

       !!
       !! New datasets go here, in new case blocks
       !! 

    end select

    !! Shouldn't get here!
    print *, 'datasets_mod%next_obs: Internal error, past select!  Stopping.'
    stop 1
  END SUBROUTINE next_obs

END MODULE datasets_mod
