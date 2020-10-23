!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! NAME: ers2_cmctnc_class_mod.f08
!!!
!!! PURPOSE: Defines a class for reading the ers2-CMCT files in netCDF
!!! format (GLA12 derivitives by Jack Saba, April 2016 processing).  Uses
!!! "icesat_config_mod" to parse the JSON configuration file and
!!! "libnetcdff" to read the netCDF files.
!!!
!!! TBD: Currently assumes that the config file's top-level key is
!!! "ers2_cmctnc_config".  Should this be changed to distinguish grn from
!!! ant?  Or maybe not because they have the same format.
!!!
!!! HISTORY: Author: Jeff Guerber, Sigma Space/GSFC 615, 2016-05-27
!!!    Adapted from ers2_cism_class_mod.f08.
!!! 2016-05-31 JRG: In clean(), check all 4 elevations.
!!! 2016-06-09 JRG: Class methods need to be private at module level.
!!! 2016-06-14 JRG: Check that parse_icesat_config returned a file list.
!!!    Bug fix in next_file: fix file open loop.
!!! 2016-06-15 JRG: Parameter FILLATTNM to work around problem in data files.
!!!    Bug fix in next_rec: make iRecStr longer.
!!! 2016-06-21 JRG: Rewrote clean() to only check a requested field.
!!! 2016-06-22 JRG: ~20x slower than flat files! Trying to improve
!!!    performance with nf90_open(...,cache_elems=).
!!! 2016-07-15 JRG: Bug fix: set self%iRec=0 in next_file() so is done always.
!!!    Improve verbose output in clean().
!!!
!!! Last SVN commit: $Id: icesat_cmctnc_class_mod.f08 116 2016-07-18 22:37:33Z jguerber $
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE ers2_cmctnc_class_mod

  use ers2_cmct_rec_mod
  use ers2_config_mod
  use netcdf

  implicit none
  private
  ! Note: The object methods (init, next_rec, next_file, clean, destroy,
  ! etc.) must be private here but public in the type definition.
  public :: ers2_cmctnc_class
  public :: ICT_INIT_OK, ICT_INIT_NOFILE, ICT_INIT_ERROR
  public :: ICT_NEXT_REC_OK, ICT_NEXT_REC_DONE, ICT_NEXT_REC_NOFILE
  public :: ICT_NEXT_FILE_OK, ICT_NEXT_FILE_END
  public :: ICT_RECORD_GOOD, ICT_RECORD_TOO_HIGH, ICT_RECORD_DEM_DIFFERENCE, &
       ICT_RECORD_VALUE_INVALID

  !! Constants
  integer, parameter :: lenstr = 256    ! Length of most strings, eg. filenames

  !! Status constants returned by the methods.  Begin with ICT_
  !! (for Icesat cmCT).  0 means successful, as usual.
  integer, parameter :: ICT_INIT_OK = 0
  integer, parameter :: ICT_INIT_NOFILE = 1
  integer, parameter :: ICT_INIT_ERROR = 2

  integer, parameter :: ICT_NEXT_REC_OK = 0
  integer, parameter :: ICT_NEXT_REC_DONE = 1
  integer, parameter :: ICT_NEXT_REC_NOFILE = 2
  integer, parameter :: ICT_NEXT_REC_ERROR = 3

  integer, parameter :: ICT_NEXT_FILE_OK = 0
  integer, parameter :: ICT_NEXT_FILE_END = 1

  !! Record quality codes from clean()
  integer, parameter :: ICT_RECORD_GOOD = 0       ! record passed
  integer, parameter :: ICT_RECORD_TOO_HIGH = 1   ! elevation is too high
  integer, parameter :: ICT_RECORD_DEM_DIFFERENCE = 2  ! elev too far from DEM
  integer, parameter :: ICT_RECORD_VALUE_INVALID = 3   ! value is invalid

  !! FILLATTNM: Fill value Attribute Name.  In some versions of Jack's
  !! netCDF GLAS data files, he used attribute _FILLVALUE instead of the
  !! correct _FillValue.  Once he fixes them, just swap this parameter.
  character(len=*), parameter :: FILLATTNM = '_FillValue'
  ! character(len=*), parameter :: FILLATTNM = '_FILLVALUE'

  !!
  !! DT for all the netCDF ids we need
  !!
  type :: netcdf_ids_t
     integer :: ncid
     ! Varids for each of the variables, correspond to the fields in
     ! icesat_cmct_rec (although abbreviated).
     integer :: time, lat, lon, tpx, wgs, egmmean, egmtidefree
     integer :: elevsd, slope, reflunc, wfsigma, gain, dem
     integer :: drainsys, icetype
  end type netcdf_ids_t

  !!
  !! Base type for the class
  !!
  type :: ers2_cmctnc_class
     private

     ! Note: gfortran 4.6.3 doesn't yet support deferred-len character
     ! components, ie. making these (len=:).
     character(len=lenstr) :: configFile  ! name of the config file
     character(len=lenstr) :: top         ! config file's top level JSON key
     character(len=4)      :: campaign    ! ers2 campaign ("L2A", etc.)
     character(len=lenstr) :: dir         ! dir for the data files for campaign
     character(len=lenstr), allocatable :: fileList(:)  ! list of campaign files
     integer  :: cFileIdx=0    ! index in fileList of file that's currently open
     type( netcdf_ids_t )  :: id         ! netcdf ids
     type( ers2_cmct_rec ) :: inv    ! invalid (fill) values read from file
     integer :: iRec, nRecs  ! record index & number of records in current file

     contains
       ! These procedures are type-bound to this type.  They comprise the
       ! methods of the class.  All are public.
       procedure :: init
       procedure :: next_rec
       procedure :: next_file
       procedure :: clean
       procedure :: destroy

  end type ers2_cmctnc_class


CONTAINS

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! init: Initializes the object, gets the list of data files from the
  !! config file, and calls next_file to open the first data file.
  !! Type-bound to type(icesat_cmctnc_class).  If already initialized,
  !! destroys self first and reinitializes.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE init( self, configFile, campaign, status, verbosity )

    !! The class variable.  Don't include as an actual argument.
    class( ers2_cmctnc_class ), intent(inout) :: self

    !! Other dummy args
    character(len=*), intent(in) :: configFile  ! JSON w/ icesat_cmctnc_config
    character(len=*), intent(in) :: campaign    ! Icesat campaign
    integer, intent(out) :: status
    integer, intent(in), optional :: verbosity   ! >0: print extra output

    !! Local variables
    integer :: nfStatus   ! status from next_file
    integer :: verb       ! local copy of verbosity

!!!---------------------------------------------------------------

    print *
    print *, 'In ers2_cmctnc_class_mod%init'

    ! Set verbosity locally rather than for object, so can use different
    ! values when calling other methods.
    verb = 0
    if (present(verbosity))  verb = verbosity

    ! Reset
    call self%destroy

    self%top        = 'ers2_cmctnc_config'
    self%campaign   = trim(adjustl(campaign))
    self%configFile = trim(adjustl(configFile))

    !! Parse the config file.  parse_ers2_config from icesat_config_mod.
    !! No status because fson just halts if there's a problem.
    call parse_ers2_config( configfile=self%configFile, top=self%top, &
         campaign=self%campaign, verbosity=verb, &
         dir=self%dir, fileList=self%fileList )
    if (.not. allocated(self%fileList)) then
       print *, 'ers2_cmctnc_class_mod%init: No file list found.'
       status = ICT_INIT_NOFILE
       return
    end if
    print *

    !! Call next_file() to open the first file.
    call self%next_file( status=nfStatus, verbosity=verb )
    if (nfStatus /= ICT_NEXT_FILE_OK) then
       print *, 'ers2_cmctnc_class_mod%init: next_file status =', nfStatus
       status = ICT_INIT_NOFILE
       return
    end if

    status = ICT_INIT_OK
    return

  END SUBROUTINE init

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! next_rec: returns the next record, of type(icesat_cmct_rec) from the
  !! current file.  If the end of the current file is reached, returns the
  !! first record from the next file.
  !!
  !! NOTE: STATUS only says whether a record was obtained, and if not why.
  !! Caller should call the CLEAN() method on the record and check its
  !! return code to see if the record should be used in calculations.
  !! (There may be situations where the record is desired even if it didn't
  !! pass, eg. for diagnostics on the file.)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  FUNCTION next_rec( self, status, invalid, verbosity ) result( retrec )

    !! The class variable.  Don't include as an actual argument.
    class( ers2_cmctnc_class ), intent(inout) :: self

    !! Dummy args & return value
    !! INVALID is an ers2_cmct_rec structure filled with the invalid
    !! values (_FillValue) for the variables in the current file.
    type( ers2_cmct_rec ) :: retrec
    integer, intent(out) :: status
    type( ers2_cmct_rec ), intent(out), optional :: invalid
    integer, intent(in), optional :: verbosity

    !! Local variables
    integer :: readstat    ! i/o status of the read statement
    integer :: nfStatus    ! status returned by next_file
    integer :: ncStatus    ! netCDF status
    character(len=lenstr) :: readmsg   ! error message from the read
    integer :: i
    integer :: verb     ! local copy of verbosity
    character(len=20) :: iRecStr    ! self%iRec written in a string
    integer :: iRecArr(1)   ! nf90_get_var needs start to be an array

!!!---------------------------------------------------------------

    verb = 0
    if (present(verbosity))  verb = verbosity

    !!
    !! We've reached the end of this file, so ask for the next one
    !!
    if ( self%iRec == self%nRecs ) then

       if (verb >= 1)  print *, "ers2_cmctnc_class_mod%next_rec: ", &
            "EOF reached, trying next file"

       call self%next_file( status=nfStatus, verbosity=verbosity )

       NEXT_FILE_STATUS: select case (nfStatus)
       case (ICT_NEXT_FILE_OK)
          ! OK: Got a file to read, so fall through to read the record
          continue
       case (ICT_NEXT_FILE_END)
          ! That was the last file, so we're all done!
          if (verb >= 1) print *, 'ers2_cmctnc_class_mod%next_rec: ', &
               'No more records, no more files'
          status = ICT_NEXT_REC_DONE
          return
       case default
          ! Some other status from next_file
          print '(a,a,i0)', 'ers2_cmctnc_class_mod%next_rec: ', &
               'Error obtaining next file, status = ', nfStatus
          status = ICT_NEXT_REC_NOFILE
          return
       end select NEXT_FILE_STATUS

    end if

    ! nf90_get_var requires start to be array even if only 1 elem
    self%iRec = self%iRec + 1
    iRecArr = [ self%iRec ]
    write (iRecStr, '("rec = ",i0)'), self%iRec

    retrec = self%inv     ! load with the invalid values
    if (present(invalid)) invalid = self%inv

    ncstatus = nf90_get_var( ncid=self%id%ncid, varid=self%id%time, &
         values=retrec%mjdtime_j2000s, start=iRecArr )
    !call check_ncdf( ncstatus, 'reading TIME variable, '//iRecStr )
    call check_ncdf( ncstatus, 'reading time variable, '//iRecStr )

    ncstatus = nf90_get_var( ncid=self%id%ncid, varid=self%id%lat, &
         values=retrec%Lat_degN, start=iRecArr )
    !call check_ncdf( ncstatus, 'reading LAT_DEGN variable, '//iRecStr )
    call check_ncdf( ncstatus, 'reading lat_degn variable, '//iRecStr )

    ncstatus = nf90_get_var( ncid=self%id%ncid, varid=self%id%lon, &
         values=retrec%Lon_degE, start=iRecArr )
    !call check_ncdf( ncstatus, 'reading LON_DEGE variable, '//iRecStr )
    call check_ncdf( ncstatus, 'reading lon_dege variable, '//iRecStr )

    !ncstatus = nf90_get_var( ncid=self%id%ncid, varid=self%id%tpx, &
    !     values=retrec%TopexElev_m, start=iRecArr )
    !call check_ncdf( ncstatus, 'reading TOPEXELEV_M variable, '//iRecStr )

    ncstatus = nf90_get_var( ncid=self%id%ncid, varid=self%id%wgs, &
         values=retrec%elev_wgs84_m, start=iRecArr )
    !call check_ncdf( ncstatus, 'reading WGS84ELEV_M variable, '//iRecStr )
    call check_ncdf( ncstatus, 'reading elev_wgs84_m variable, '//iRecStr )

    ncstatus = nf90_get_var( ncid=self%id%ncid, varid=self%id%egmmean, &
         values=retrec%elev_egm08meantide_m, start=iRecArr )
    !call check_ncdf(ncstatus, 'reading EGM08MEANTIDEELEV_M variable, '//iRecStr)
    call check_ncdf(ncstatus, 'reading elev_egm08meantide_m variable, '//iRecStr)

    ncstatus = nf90_get_var( ncid=self%id%ncid, varid=self%id%egmtidefree, &
         values=retrec%elev_egm08tidefree_m, start=iRecArr )
    !call check_ncdf(ncstatus, 'reading EGM08TIDEFREEELEV_M variable, '//iRecStr)
    call check_ncdf(ncstatus, 'reading elev_egm08tidefree_m variable, '//iRecStr)

    ncstatus = nf90_get_var( ncid=self%id%ncid, varid=self%id%elevsd, &
         values=retrec%elev_stddev_m, start=iRecArr )
    !call check_ncdf( ncstatus, 'reading ELEV_STDDEV_M variable, '//iRecStr )
    call check_ncdf( ncstatus, 'reading elev_stddev_m variable, '//iRecStr )

    ncstatus = nf90_get_var( ncid=self%id%ncid, varid=self%id%slope, &
         values=retrec%slope_deg, start=iRecArr )
    !call check_ncdf( ncstatus, 'reading SLOPE_DEG variable, '//iRecStr )
    call check_ncdf( ncstatus, 'reading slope_deg variable, '//iRecStr )

   ! ncstatus = nf90_get_var( ncid=self%id%ncid, varid=self%id%reflunc, &
   !      values=retrec%refl_uncorr, start=iRecArr )
   ! call check_ncdf( ncstatus, 'reading REFL_UNCORR variable, '//iRecStr )

    !ncstatus = nf90_get_var( ncid=self%id%ncid, varid=self%id%wfsigma, &
    !     values=retrec%wf_fit_sigma, start=iRecArr )
    !call check_ncdf( ncstatus, 'reading WF_FIT_SIGMA variable, '//iRecStr )

    !ncstatus = nf90_get_var( ncid=self%id%ncid, varid=self%id%gain, &
    !     values=retrec%gain, start=iRecArr )
    !call check_ncdf( ncstatus, 'reading GAIN variable, '//iRecStr )

    ncstatus = nf90_get_var( ncid=self%id%ncid, varid=self%id%dem, &
         values=retrec%dem_m, start=iRecArr )
    !call check_ncdf( ncstatus, 'reading DEM_M variable, '//iRecStr )
    call check_ncdf( ncstatus, 'reading dem_m variable, '//iRecStr )

    ncstatus = nf90_get_var( ncid=self%id%ncid, varid=self%id%time, &
         values=retrec%drainagesystem, start=iRecArr )
    !call check_ncdf( ncstatus, 'reading DRAINAGESYSTEM variable, '//iRecStr )
    call check_ncdf( ncstatus, 'reading drainagesystem variable, '//iRecStr )

    ncstatus = nf90_get_var( ncid=self%id%ncid, varid=self%id%drainsys, &
         values=retrec%icetype, start=iRecArr )
    !call check_ncdf( ncstatus, 'reading ICETYPE variable, '//iRecStr )
    call check_ncdf( ncstatus, 'reading icetype variable, '//iRecStr )

    status = ICT_NEXT_REC_OK
    return

  END FUNCTION next_rec


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! subroutine next_file: Keeps track of which file is currently open.
  !! Open the next one when requested.  Fill in all the netcdf varids and
  !! get the invalid (fill) values for each variable.
  !! These files have 1D variables whose only dimension is a record
  !! (unlimited) one.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE next_file( self, status, verbosity )

    !! The class variable.  Don't include as an actual argument.
    class( ers2_cmctnc_class ), intent(inout) :: self

    !! Dummy args
    integer, intent(out) :: status
    integer, intent(in), optional :: verbosity

    !! Local variables
    integer :: iFile
    character(len=2*lenstr) :: full   ! full name of the file
    integer :: ncstatus               ! netCDF status
    integer :: verb                   ! local verbosity flag
    character(len=lenstr) :: msg      ! error message
    integer :: recDimID         ! dimid of the record (unlimited) dimension

    !!---------------------------------------------------------------

    verb = 0
    if (present(verbosity))  verb = verbosity

    !! Close the currently open file.  Not a problem if it fails, just
    !! means there was no file open already.
    ncstatus = nf90_close( self%id%ncid )

    !! Starting at next file (which may be the first), loop up through the
    !! array of file names until we find one we can open.
    status = ICT_NEXT_FILE_END
    do iFile = self%cFileIdx + 1, size(self%fileList)

       !! Try to open the next file.
       !!
       !! It's not clear to me just how netCDF-4 chunk caching works,
       !! especially how cache_size and cache_nelems interact.  These files
       !! have 15 fields and ~9500 to ~35000 records.  For the set Jack
       !! created with 5000 element chunks, try making the cache big enough
       !! to hold several chunks from all the variables.
       full = trim( self%dir ) // trim( self%fileList(iFile) )
       ncstatus = nf90_open( path=full, mode=NF90_NOWRITE, ncid=self%id%ncid, &
            cache_nelems=15*5000*5 )

       ! !! Try to open the next file.
       ! !!
       ! !! It's not clear to me just how netCDF-4 chunk caching works,
       ! !! especially how cache_size and cache_nelems interact.  These files
       ! !! have 15 fields and ~9500 to ~35000 records.  In Jack's latest
       ! !! versions, he wrote them with the chunk size depending on the
       ! !! number of elements; processing eg. L3G still takes ~330 sec
       ! !! without chunk caching here.  Let's try making it big enough to
       ! !! hold at least 1 chunk from all the variables in the largest
       ! !! dataset:
       ! full = trim( self%dir ) // trim( self%fileList(iFile) )
       ! ncstatus = nf90_open( path=full, mode=NF90_NOWRITE, ncid=self%id%ncid, &
       !      cache_nelems=15*35000 )

       !! Success, so record the index and get out of loop
       if ( ncstatus == NF90_NOERR ) then
          print '(a,i0,a)',' ers2_cmctnc_class_mod%next_file: Reading file ', &
               iFile, ':'
          print *, trim(full)
          self%cFileIdx = iFile
          status = ICT_NEXT_FILE_OK
          exit
       end if

       !! If couldn't open it, complain and try the next one
       write (msg, '(a,i0,a,a)'), 'in next_file: opening file ', iFile, &
            ':', trim(full)
       call check_ncdf( ncstatus, msg )

    end do

    !! If we get here and status is still ICT_NEXT_FILE_END, we've
    !! exhausted the list of files, so report that and return.
    if (status == ICT_NEXT_FILE_END) then
       if (verb >= 1) print '(a,i0)', &
            'ers2_cmctnc_class_mod%next_file: File list exhausted. ifile=', &
            ifile
       return
    end if

    !! Get the number of records in this file (size of the unlimited
    !! dimension)

    ncstatus = nf90_inquire( ncid=self%id%ncid, unlimitedDimID=recDimID )
    call check_ncdf( ncstatus, 'getting ID of unlimited dimension' )
    ncstatus = nf90_inquire_dimension( ncid=self%id%ncid, dimid=recDimID, &
         len=self%nRecs )
    call check_ncdf( ncstatus, 'getting number of records' )
    print '(a,i0)', ' records = ', self%nRecs

    self%iRec = 0

    !!
    !! Now get all the variable ids and invalid (fill) values
    !!
    ncstatus = nf90_inq_varid( ncid=self%id%ncid, &
         name='time', varid=self%id%time)
    call check_ncdf( ncstatus, 'getting time varid')
    !ncstatus = nf90_get_att( ncid=self%id%ncid, varid=self%id%time, &
    !     name=FILLATTNM, values=self%inv%utctime_j2000s )
    ncstatus = nf90_get_att( ncid=self%id%ncid, varid=self%id%time, &
         name=FILLATTNM, values=self%inv%mjdtime_j2000s )
    call check_ncdf( ncstatus, 'getting time fill value')

    ncstatus = nf90_inq_varid( ncid=self%id%ncid, &
         name='lat_degn', varid=self%id%lat)
    call check_ncdf( ncstatus, 'getting lat_degn varid')
    ncstatus = nf90_get_att( ncid=self%id%ncid, varid=self%id%lat, &
         name=FILLATTNM, values=self%inv%lat_degn )
    call check_ncdf( ncstatus, 'getting lat_degn fill value')

    ncstatus = nf90_inq_varid( ncid=self%id%ncid, &
         name='lon_dege', varid=self%id%lon)
    call check_ncdf( ncstatus, 'getting lon_dege varid')
    ncstatus = nf90_get_att( ncid=self%id%ncid, varid=self%id%lon, &
         name=FILLATTNM, values=self%inv%lon_dege )
    call check_ncdf( ncstatus, 'getting lon_dege fill value')

   ! ncstatus = nf90_inq_varid( ncid=self%id%ncid, &
    !     name='TOPEXELEV_M', varid=self%id%tpx)
   ! call check_ncdf( ncstatus, 'getting TOPEXELEV_M varid')
   ! ncstatus = nf90_get_att( ncid=self%id%ncid, varid=self%id%tpx, &
   !      name=FILLATTNM, values=self%inv%topexelev_m )
   ! call check_ncdf( ncstatus, 'getting TOPEXELEV_M fill value')

    ncstatus = nf90_inq_varid( ncid=self%id%ncid, &
         name='elev_wgs84_m', varid=self%id%wgs)
    call check_ncdf( ncstatus, 'getting elev_wgs84_m varid')
    ncstatus = nf90_get_att( ncid=self%id%ncid, varid=self%id%wgs, &
         name=FILLATTNM, values=self%inv%elev_wgs84_m)
    call check_ncdf( ncstatus, 'getting elev_wgs84_m fill value')

    ncstatus = nf90_inq_varid( ncid=self%id%ncid, &
         name='elev_egm08meantide_m', varid=self%id%egmmean)
    call check_ncdf( ncstatus, 'getting elev_egm08meantide_m varid')
    ncstatus = nf90_get_att( ncid=self%id%ncid, varid=self%id%egmmean, &
         name=FILLATTNM, values=self%inv%elev_egm08meantide_m )
    call check_ncdf( ncstatus, 'getting elev_egm08meantide_m fill value')

    ncstatus = nf90_inq_varid( ncid=self%id%ncid, &
         name='elev_egm08tidefree_m', varid=self%id%egmtidefree)
    call check_ncdf( ncstatus, 'getting elev_egm08tidefree_m varid')
    ncstatus = nf90_get_att( ncid=self%id%ncid, varid=self%id%egmtidefree, &
         name=FILLATTNM, values=self%inv%elev_egm08tidefree_m )
    call check_ncdf( ncstatus, 'getting elev_egm08tidefree_m fill value')

    ncstatus = nf90_inq_varid( ncid=self%id%ncid, &
         name='elev_stddev_m', varid=self%id%elevsd)
    call check_ncdf( ncstatus, 'getting elev_stddev_m varid')
    ncstatus = nf90_get_att( ncid=self%id%ncid, varid=self%id%elevsd, &
         name=FILLATTNM, values=self%inv%elev_stddev_m )
    call check_ncdf( ncstatus, 'getting elev_stddev_m fill value')

    ncstatus = nf90_inq_varid( ncid=self%id%ncid, &
         name='slope_deg', varid=self%id%slope)
    call check_ncdf( ncstatus, 'getting slope_deg varid')
    ncstatus = nf90_get_att( ncid=self%id%ncid, varid=self%id%slope, &
         name=FILLATTNM, values=self%inv%slope_deg )
    call check_ncdf( ncstatus, 'getting slope_deg fill value')

   ! ncstatus = nf90_inq_varid( ncid=self%id%ncid, &
   !      name='REFL_UNCORR', varid=self%id%reflunc)
   ! call check_ncdf( ncstatus, 'getting REFL_UNCORR varid')
   ! ncstatus = nf90_get_att( ncid=self%id%ncid, varid=self%id%reflunc, &
   !      name=FILLATTNM, values=self%inv%refl_uncorr )
   ! call check_ncdf( ncstatus, 'getting REFL_UNCORR fill value')

    !ncstatus = nf90_inq_varid( ncid=self%id%ncid, &
    !     name='WF_FIT_SIGMA_M', varid=self%id%wfsigma)
    !call check_ncdf( ncstatus, 'getting WF_FIT_SIGMA_M varid')
    !ncstatus = nf90_get_att( ncid=self%id%ncid, varid=self%id%wfsigma, &
    !     name=FILLATTNM, values=self%inv%wf_fit_sigma )
    !call check_ncdf( ncstatus, 'getting WF_FIT_SIGMA_M fill value')

    !ncstatus = nf90_inq_varid( ncid=self%id%ncid, &
    !     name='GAIN', varid=self%id%gain)
    !call check_ncdf( ncstatus, 'getting GAIN varid')
    !ncstatus = nf90_get_att( ncid=self%id%ncid, varid=self%id%gain, &
    !     name=FILLATTNM, values=self%inv%gain )
    !call check_ncdf( ncstatus, 'getting GAIN fill value')

    ncstatus = nf90_inq_varid( ncid=self%id%ncid, &
         name='dem_m', varid=self%id%dem)
    call check_ncdf( ncstatus, 'getting dem_m varid')
    ncstatus = nf90_get_att( ncid=self%id%ncid, varid=self%id%dem, &
         name=FILLATTNM, values=self%inv%dem_m )
    call check_ncdf( ncstatus, 'getting dem_m fill value')

    ncstatus = nf90_inq_varid( ncid=self%id%ncid, &
         name='drainagesystem', varid=self%id%drainsys )
    call check_ncdf( ncstatus, 'getting drainagesystem varid')
    ncstatus = nf90_get_att( ncid=self%id%ncid, varid=self%id%drainsys, &
         name=FILLATTNM, values=self%inv%drainagesystem )
    call check_ncdf( ncstatus, 'getting drainagesystem fill value')

    ncstatus = nf90_inq_varid( ncid=self%id%ncid, &
         name='icetype', varid=self%id%icetype)
    call check_ncdf( ncstatus, 'getting icetype varid')
    ncstatus = nf90_get_att( ncid=self%id%ncid, varid=self%id%icetype, &
         name=FILLATTNM, values=self%inv%icetype )
    call check_ncdf( ncstatus, 'getting icetype fill value')

    return
  END SUBROUTINE next_file


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! clean: Apply quality filters to a field of the record, returns code
  !! ICT_RECORD_GOOD (== 0) or a reason the record should be rejected (/=
  !! 0).  The filters are based on those in Jack's cism_getdata.pro and
  !! read_glas_data.pro.  Currently checks any of the 4 elevations in the
  !! record.
  !!
  !! Currently implemented filters are:
  !!   - Field does not have its INVALID value.
  !!   - Elevation over Greenland < 3300 m
  !!   - Elevation over Antarctica < 4500 m
  !!   - Elevation is within 200 m of the GIMP DEM
  !!
  !! TBD: Jack may have included these when he created the data files.
  !! This may not be necessary at all!
  !!
  !! ARGUMENTS:
  !!   rec (icesat_cmct_rec, in): The record to check.
  !!   fieldname (char, in): Name of the specific field to verify.
  !!      Currently only does the four elevation fields. In LOWERCASE.
  !!   verbosity (int, in):
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  FUNCTION clean( self, rec, fieldname, verbosity )

    use iso_fortran_env, only: REAL32

    integer :: clean

    !! The class variable.  Don't include as an actual argument.
    class( ers2_cmctnc_class ), intent(inout) :: self

    !! Dummy args
    type( ers2_cmct_rec ), intent(in) :: rec
    character(len=*), intent(in)        :: fieldname
    integer, intent(in), optional       :: verbosity

    !! Local variables
    integer :: verb                 ! local verbosity flag
    real(kind=REAL32) :: elev, inv  ! specific elevation value and its invalid
    character(len=30) :: elevName   ! name of the elevation

    !! Constants.  These values are taken from Jack's cism_getdata.pro and
    !! read_glas_data.pro.  They should probably be read from configuration.
    real(kind=REAL32), parameter :: GRN_MAX_ELEV = 3300.0
    real(kind=REAL32), parameter :: ANT_MAX_ELEV = 4500.0
    real(kind=REAL32), parameter :: DEM_MAX_DIFF = 200.0

    !!---------------------------------------------------------------

    verb = 0
    if (present(verbosity))  verb = verbosity

    !! Pull out the requested field (for now, only doing elevations),
    !! appropriate invalied value (may need to change how this is done if
    !! extend to other fields, to handle other types), and set
    select case ( trim(fieldname) )
    !case ( 'topexelev_m' )
    !   elev = rec%topexelev_m
    !   inv  = self%inv%topexelev_m
    !   elevName = 'TopexElev_m'

    case ( 'elev_wgs84_m' )
       elev = rec%elev_wgs84_m
       inv  = self%inv%elev_wgs84_m
       elevName = 'elev_wgs84_m'

    case ( 'elev_egm08meantide_m' )
       elev = rec%elev_egm08meantide_m
       inv  = self%inv%elev_egm08meantide_m
       elevName = 'elev_egm08meantide_m'

    case ( 'elev_egm08tidefree_m' )
       elev = rec%elev_egm08tidefree_m
       inv  = self%inv%elev_egm08tidefree_m
       elevName = 'elev_egm08tidefree_m'

    case default
       print *, 'ers2_cmctnc_class_mod%clean: Unhandled field requested: ',&
            trim(fieldname)
       print *, 'Stopping.'
       stop 1
    end select

    !! Invalid value
    if ( elev == inv ) then
       clean = ICT_RECORD_VALUE_INVALID
       if (verb >= 1) then
          print *, 'ers2_cmctnc_class_mod%clean: ', &
               trim(elevName), ' has "Invalid" value'
          print *, elev, ' at ', rec%lat_degN, rec%lon_degE, ' rec ', self%iRec
       end if
       return
    end if

    !! Elevation is whacko - Greenland
    if ( (rec%lat_degN > 0.0d0) .and. (elev >= GRN_MAX_ELEV) ) then
       clean = ICT_RECORD_TOO_HIGH
       if (verb >= 1) then
          print *, 'ers2_cmctnc_class_mod%clean: ', &
               trim(elevName), ' exceeds Greenland max ', GRN_MAX_ELEV
          print *, elev, ' at ', rec%lat_degN, rec%lon_degE, ' rec ', self%iRec
       end if
       return
    end if

    !! Elevation is whacko - Antarctica
    if ( (rec%lat_degN <= 0.0d0) .and. (elev >= ANT_MAX_ELEV) ) then
       clean = ICT_RECORD_TOO_HIGH
       if (verb >= 1) then
          print *, 'ers2_cmctnc_class_mod%clean: ', &
               trim(elevName), ' exceeds Antarctica max ', ANT_MAX_ELEV
          print *, elev, ' at ', rec%lat_degN, rec%lon_degE, ' rec ', self%iRec
       end if
       return
    end if

    !! Elevation is too far from the GIMP DEM
    if ( abs(elev - rec%dem_m) > DEM_MAX_DIFF ) then
       clean = ICT_RECORD_DEM_DIFFERENCE
       if (verb >= 1) then
          print *, 'ers2_cmctnc_class_mod%clean: ', &
               trim(elevName), ' is too far from DEM ', DEM_MAX_DIFF
          print *, elev, rec%dem_m, elev - rec%dem_m, &
               ' at ', rec%lat_degN, rec%lon_degE, ' rec ', self%iRec
       end if
       return
    end if

    clean = ICT_RECORD_GOOD
    return

  END FUNCTION clean


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! destroy - Close file, free allocated arrays.
  !!
  !! This should be a finalizer, so the file is closed and the array
  !! deallocated if the object goes out of scope, but gfortran doesn't
  !! support them until 4.9.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE destroy ( self )

    ! self is the class variable, don't include as an actual argument
    class( ers2_cmctnc_class ), intent(inout) :: self

    ! local variables
    integer :: status

    ! OK if nf90_close fails, just means ncid isn't open which is fine.
    status = nf90_close( self%id%ncid )
    self%id%ncid = -1
    self%cFileIdx = 0

    if ( allocated( self%fileList ) )  deallocate ( self%fileList )

    return

  END SUBROUTINE destroy


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Private module subroutine check_ncdf
  !! Checks the status of netcdf calls.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine check_ncdf ( status, msg )
    integer, intent(in) :: status
    character(len=*), intent(in) :: msg

    if ( status /= nf90_NoErr ) then
       print *
       print *, 'ers2_cmctnc_class_mod: NetCDF error ' // trim(msg)
       print *, '   nf90 status = ', status
       print *, '   ', trim( nf90_strerror(status) )
       stop 1
    end if
    return

  end subroutine check_ncdf


END MODULE ers2_cmctnc_class_mod
