!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! NAME: icesat_cism_class_mod.f08
!!!
!!! PURPOSE: Defines a class for reading the ICESAT-CISM files (GLA12
!!! derivitives by Jack Saba, orignal 2014 processing, flat-binary format).
!!!
!!! HISTORY: Author: Jeff Guerber, Sigma Space/GSFC 615, 2015-07-01
!!! 2015-11-18 JRG: Changed debug to verbosity, now integer. Fixed up prints.
!!! 2016-03-03 JRG: Added method function CLEAN(), ICM_RECORD_* codes.
!!! 2016-05-11 JRG: init: Deallocate self%fileList if already allocated.
!!! 2016-05-13 JRG: Added method DESTROY. Call from INIT.
!!! 2016-05-25 JRG: Moved parsing of JSON config file from init() to
!!!   icesat_config_mod.
!!! 2016-06-09 JRG: Class methods need to be private at module level.
!!! 2016-06-14 JRG: Check that parse_icesat_config returned a file list.
!!!
!!! Last SVN commit: $Id: icesat_cism_class_mod.f08 105 2016-07-02 10:46:41Z jguerber $
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE icesat_cism_class_mod

  use icesat_cism_rec_mod      ! defines the record structure

  implicit none
  private
  ! Note: The object methods (init, next_rec, next_file, clean, destroy,
  ! etc.) must be private here but public in the type definition.
  public :: icesat_cism_class
  public :: ICM_INIT_OK, ICM_INIT_NOFILE, ICM_INIT_ERROR
  public :: ICM_NEXT_REC_OK, ICM_NEXT_REC_DONE, ICM_NEXT_REC_NOFILE
  public :: ICM_NEXT_FILE_OK, ICM_NEXT_FILE_END
  public :: ICM_RECORD_GOOD, ICM_RECORD_TOO_HIGH, ICM_RECORD_DEM_DIFFERENCE

  !! Constants
  integer, parameter :: lenstr = 256    ! Length of most strings, eg. filenames

  !! Status constants returned by the methods.  Begin with ICM_
  !! (for Icesat CisM).  0 means successful, as usual.
  integer, parameter :: ICM_INIT_OK = 0
  integer, parameter :: ICM_INIT_NOFILE = 1
  integer, parameter :: ICM_INIT_ERROR = 2

  integer, parameter :: ICM_NEXT_REC_OK = 0
  integer, parameter :: ICM_NEXT_REC_DONE = 1
  integer, parameter :: ICM_NEXT_REC_NOFILE = 2
  integer, parameter :: ICM_NEXT_REC_ERROR = 3

  integer, parameter :: ICM_NEXT_FILE_OK = 0
  integer, parameter :: ICM_NEXT_FILE_END = 1

  !! Record quality codes from clean()
  integer, parameter :: ICM_RECORD_GOOD = 0       ! record passed
  integer, parameter :: ICM_RECORD_TOO_HIGH = 1   ! elevation is too high
  integer, parameter :: ICM_RECORD_DEM_DIFFERENCE = 2  ! elev is too far from DEM

  !!
  !! Base type for the class
  !!
  type :: icesat_cism_class

     ! Note: gfortran 4.6.3 doesn't yet support deferred-len character
     ! components, ie. making these (len=:).
     character(len=lenstr) :: configFile    ! name of the config file
     character(len=lenstr) :: top     ! top-level key in the JSON config file
     character(len=3) :: campaign        ! Icesat campaign ("L1A", etc.)
     character(len=lenstr) :: dir        ! dir for the campaign files
     character(len=lenstr), allocatable :: fileList(:)  ! list of campaign files
     integer :: lun         ! Unit for reading the files
     integer :: cFileIdx=0  ! index of the currently-open file in list

     contains
       ! These procedures are type-bound to this type.  They comprise the
       ! methods of the class.  All are public.
       procedure :: init
       procedure :: next_rec
       procedure :: next_file
       procedure :: clean
       procedure :: destroy
       ! final     :: destroy     ! final not implemented in gfortran 4.6

  end type icesat_cism_class


CONTAINS

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! init: Initializes the object, gets the list of data files from the
  !! config file, and calls next_file to open the first data file.
  !! Type-bound to type(icesat_cism_class).  If already initialized,
  !! destroys self first and reinitializes.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE init( self, configFile, campaign, status, verbosity )

    use icesat_config_mod       ! parse the icesat dataset config file

    !! The class variable.  Don't include as an actual argument.
    class( icesat_cism_class ), intent(inout) :: self

    !! Other dummy args
    character(len=*), intent(in) :: configFile  ! JSON w/ icesat_cism_config
    character(len=*), intent(in) :: campaign    ! Icesat campaign
    integer, intent(out) :: status
    integer, intent(in), optional :: verbosity   ! >0: print extra output

    !! Local variables
    integer :: nfStatus   ! status from next_file
    integer :: verb       ! local copy of verbosity

!!!---------------------------------------------------------------

    print *
    print *, 'In icesat_cism_class_mod%init'

    ! Set verbosity locally rather than for object, so can use different
    ! values when calling other methods.
    verb = 0
    if (present(verbosity))  verb = verbosity

    ! Reset
    call self%destroy

    self%configFile = trim(adjustl(configFile))
    self%top        = 'icesat_cism_config'
    self%campaign   = trim(adjustl(campaign))

    !! Parse the config file.  parse_icesat_config from icesat_config_mod.
    !! No status because fson usually just halts if there's a problem.
    call parse_icesat_config( configfile=self%configFile, top=self%top, &
         campaign=self%campaign, verbosity=verb, &
         dir=self%dir, fileList=self%fileList )
    if (.not. allocated(self%fileList)) then
       print *, 'icesat_cism_class_mod%init: No file list found.'
       status = ICM_INIT_NOFILE
       return
    end if
    print *

    !! Call next_file() to get a LUN and open the first file.
    call self%next_file( status=nfStatus, verbosity=verb )
    if (nfStatus .ne. ICM_NEXT_FILE_OK) then
       print *, 'icesat_cism_class_mod%init: next_file status =', nfStatus
       status = ICM_INIT_NOFILE
       return
    end if

    status = ICM_INIT_OK
    return

  END SUBROUTINE init

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! next_rec: returns the next record, of type(icesat_cism_rec) from the
  !! current file.  If the end of the current file is reached, returns the
  !! first record from the next file.
  !!
  !! NOTE: STATUS only says whether a record was obtained, and if not why.
  !! Caller should call the CLEAN() method on the record and check its
  !! return code to see if the record should be used in calculations.
  !! (There may be situations where the record is desired even if it didn't
  !! pass, eg. for diagnostics on the file.)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  FUNCTION next_rec( self, status, verbosity ) result( retrec )

    use iso_fortran_env, only: iostat_end

    !! The class variable.  Don't include as an actual argument.
    class( icesat_cism_class ), intent(inout) :: self

    !! Dummy args & return value
    type( icesat_cism_rec ) :: retrec
    integer, intent(out) :: status
    integer, intent(in), optional :: verbosity

    !! Local variables
    integer :: readstat    ! i/o status of the read statement
    integer :: nfStatus    ! status returned by next_file
    character(len=lenstr) :: readmsg   ! error message from the read
    integer :: i
    integer :: verb     ! local copy of verbosity

!!!---------------------------------------------------------------

    status = 0

    verb = 0
    if (present(verbosity))  verb = verbosity

    !! Try several times to read a record.  If reach end of current file,
    !! try to open a new file.  Finite loop so that it won't run away (too
    !! far) if something's wrong with next_file.
    READ_LOOP: do i = 1, 10

       if ((verb >= 1) .and. (i /= 1)) then
          print '(a,i0)', 'icesat_cism_class_mod%next_rec: Loop iteration = ', i
       end if

       read (self%lun, iostat=readstat, iomsg=readmsg ) retrec

       READ_STATUS: select case (readstat)
       case (0)
          !! SUCCESS!
          status = ICM_NEXT_REC_OK
          return

       case (IOSTAT_END)

          !! End of this file, so ask for the next one
          if (verb >= 1)  print *, "icesat_cism_class_mod%next_rec: ", &
               "EOF reached, trying next file"

          call self%next_file( status=nfStatus, verbosity=verbosity )

          NEXT_FILE_STATUS: select case (nfStatus)
          case (ICM_NEXT_FILE_OK)
             ! OK: Got a file to read, so cycle the loop to try the read again
             cycle read_loop
          case (ICM_NEXT_FILE_END)
             ! That was the last file, so we're all done!
             if (verb >= 1) print *, 'icesat_cism_class_mod%next_rec: ', &
                  'No more records, no more files'
             status = ICM_NEXT_REC_DONE
             return
          case default
             ! Some other status from next_file
             print '(a,a,i0)', 'icesat_cism_class_mod%next_rec: ', &
                  'Error obtaining next file, status = ', nfStatus
             status = ICM_NEXT_REC_NOFILE
             return
          end select NEXT_FILE_STATUS

       case default
          !! Some other read error
          print '(a,a,i0)', 'icesat_cism_class_mod%next_rec: ', &
               'Error reading record, iostat= ', readstat
          print *, readmsg
          status = ICM_NEXT_REC_ERROR
          return

       end select READ_STATUS

       !! Probably can't get here.  Should this be an error?
       print '(a,a,i0)', 'icesat_cism_class_mod%next_rec: ', &
            'Warning, Unhandled read status, iostat= ', readstat
       print *, readmsg

    end do READ_LOOP

    !!  Should only get here if next_file can't open i_max files in succession
    print *, 'icesat_cism_class_mod%next_rec: Too many read attempts'
    status = ICM_NEXT_REC_ERROR
    return

  END FUNCTION next_rec

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!! subroutine next_file: Keeps track of which file is currently open.
  !!! Open the next one when requested.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine next_file( self, status, verbosity )

    !! The class variable.  Don't include as an actual argument.
    class( icesat_cism_class ), intent(inout) :: self

    !! Dummy args
    integer, intent(out) :: status
    integer, intent(in), optional :: verbosity

    !! Local variables
    integer :: iFile, iLun
    character(len=2*lenstr) :: full   ! full name of the file
    integer :: openstat               ! status of the open
    character(len=lenstr) :: openmsg  ! system msg from an unsuccessful open
    logical :: isOpen                 ! Testing luns, is it open?
    integer :: verb                   ! local verbosity flag

!!!---------------------------------------------------------------

    verb = 0
    if (present(verbosity))  verb = verbosity

    if (self%cFileIdx .eq. 0) then

       !! First time: select a (hopefully) otherwise unused LUN.  (Would
       !! have preferred to use open(newunit=), but that doesn't work if
       !! there's an error opening the file, like not existing.)  What if
       !! we get all the way through and don't find one??  (Possible soln
       !! on fortran90.org.)  Also: _NOLUN error TBD.
       !!
       !! [JRG 2016-05-25: Just why DID I do it this way instead of simply
       !! using newunit=self%lun, and closing the old lun before opening
       !! the new file?  I really don't remember, and can't find any notes
       !! that elucidate my thinking. Perhaps I was thinking someone might
       !! potentially want to have more than one instance simultaneously?
       !! Would that even be a problem?]

       do iLun = 200, 999
          if (verb >= 1) print '(a,i0)', &
               'icesat_cism_class_mod%next_file: Trying LUN ', iLun
          inquire( unit=iLun, opened=isOpen )
          if (.not. isOpen) then
             self%lun = iLun
             exit
          end if
       end do
       if (verb >= 1) print '(a,i0)', &
            'icesat_cism_class_mod%next_file: Got LUN=', self%lun

    end if

    !! Starting at next file (which may be the first), loop up through the
    !! array of file names until we find one we can open.
    do iFile = self%cFileIdx + 1, size(self%fileList)

       ! full filename
       full = trim( self%dir ) // trim( self%fileList(iFile) )

       !! Try to open the next file.  Note we don't explicitly close the
       !! old one, so we can be sure to keep the same LUN.  These files
       !! don't have the record headers that unformatted sequentials do, so
       !! use stream access.  (Direct access requires record numbers to
       !! read... seem to recall used to be able to read them as if
       !! sequential??)  Stream also needs no record length.

       open ( self%lun, file=full, form='unformatted', &
             status='old', access='stream', action='read', &
             iostat=openstat, iomsg=openmsg )

       if ( openstat == 0) then
          print '(a,i0,a)', 'icesat_cism_class_mod%next_file: Reading file ', &
               iFile, ':'
          print *, trim(full)
          self%cFileIdx = iFile
          status = ICM_NEXT_FILE_OK
          return
       else
          !! If couldn't open it, try the next one
          print '(a,i0,a)', &
               'icesat_cism_class_mod%next_file: Error opening file ', &
               iFile, ':'
          print *, trim(full)
          print *, trim(openmsg)
          cycle
       end if
    end do

    !! If we get here, we've exhausted the list of files
    if (verb >= 1) print '(a,i0)', &
         'icesat_cism_class_mod%next_file: File list exhausted. ifile=', ifile
    status = ICM_NEXT_FILE_END

    return
  end subroutine next_file

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!! clean: Apply quality filters to the record, returns code
  !!! ICM_RECORD_GOOD (== 0) or a reason the record should be rejected (/=
  !!! 0).  The filters are based on those in Jack's cism_getdata.pro and
  !!! read_glas_data.pro.
  !!!
  !!! Currently implemented filters are:
  !!!   - Elevation over Greenland < 3300 m
  !!!   - Elevation over Antarctica < 4500 m
  !!!   - Elevation is within 200 m of the GIMP DEM
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  FUNCTION clean( self, rec, verbosity )

    use iso_fortran_env, only: REAL32

    integer :: clean

    !! The class variable.  Don't include as an actual argument.
    class( icesat_cism_class ), intent(inout) :: self

    !! Dummy args
    type( icesat_cism_rec ), intent(in) :: rec
    integer, intent(in), optional :: verbosity

    !! Local variables
    integer :: verb                   ! local verbosity flag

    !! Constants.  These values are taken from Jack's cism_getdata.pro and
    !! read_glas_data.pro.  They should probably be read from configuration.
    real(kind=REAL32), parameter :: GRN_MAX_ELEV = 3300.0
    real(kind=REAL32), parameter :: ANT_MAX_ELEV = 4500.0
    real(kind=REAL32), parameter :: DEM_MAX_DIFF = 200.0

!!!---------------------------------------------------------------

    verb = 0
    if (present(verbosity))  verb = verbosity

    !!! Elevation is whacko
    if ( (rec%lat_degN > 0.0d0) .and. (rec%elev_m >= GRN_MAX_ELEV) ) then
       clean = ICM_RECORD_TOO_HIGH
       if (verb >= 1) then
          print *, 'icesat_cism_class_mod%clean: ', &
               'Elevation exceeds Greenland max'
          print *, rec%elev_m, GRN_MAX_ELEV, ' at ', rec%lat_degN, rec%lon_degE
       end if
       return
    end if

    if ( (rec%lat_degN <= 0.0d0) .and. (rec%elev_m >= ANT_MAX_ELEV) ) then
       clean = ICM_RECORD_TOO_HIGH
       if (verb >= 1) then
          print *, 'icesat_cism_class_mod%clean: ', &
               'Elevation exceeds Antarctica max'
          print *, rec%elev_m, ANT_MAX_ELEV, ' at ', rec%lat_degN, rec%lon_degE
       end if
       return
    end if

    !!! Elevation is too far from the GIMP DEM
    if ( abs(rec%elev_m - rec%gimp_dem) > DEM_MAX_DIFF ) then
       clean = ICM_RECORD_DEM_DIFFERENCE
       if (verb >= 1) then
          print *, 'icesat_cism_class_mod%clean: ', &
               'Elevation is too far from GIMP DEM'
          print *, rec%elev_m, rec%gimp_dem, DEM_MAX_DIFF, &
               ' at ', rec%lat_degN, rec%lon_degE
       end if
       return
    end if

    clean = ICM_RECORD_GOOD
    return

  END FUNCTION clean


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!! destroy - Close file, free allocated arrays.
  !!!
  !!! This should be a finalizer, so the file is closed and arrays
  !!! deallocated if the object goes out of scope, but gfortran doesn't
  !!! support them until 4.9.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE destroy ( self )

    ! self is the class variable, don't include as an actual argument
    class( icesat_cism_class ), intent(inout) :: self


    close ( self%lun )
    self%cFileIdx = 0
    self%lun = -1

    if ( allocated( self%fileList ) )  deallocate ( self%fileList )

    return

  END SUBROUTINE destroy

END MODULE icesat_cism_class_mod
