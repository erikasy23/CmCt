!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! NAME: envisat_config_mod.f08
!!!
!!! PURPOSE: Module to read the envisat datasets' JSON configuration files,
!!! which are organized by campaign.  Uses "fson" to parse the file.
!!!
!!! HISTORY: Author: Jeff Guerber, Sigma Space/GSFC 615, 2016-05-20
!!!    Split out of envisat_cism_class_mod.f08 subroutine init.
!!! 2016-06-14 JRG: Bug fix: fson_get(...,fileListJArrp) used hardcoded
!!!    (and sometimes incorrect) top-level key.  Added some error checking.
!!!    Initialize variables.  Fixed up messages.
!!!
!!! Last SVN commit: $Id: icesat_config_mod.f08 105 2016-07-02 10:46:41Z jguerber $
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE envisat_config_mod


  implicit none
  private
  public parse_envisat_config

  !! Constants
  integer, parameter :: lenstr = 256  ! Length of most strings, eg. filenames

CONTAINS

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! parse_envisat_config: Parse the envisat dataset config files.  Returns
  !! the data directory name and list of filenames.
  !!
  !! Dummy arguments:
  !!   configFile (char in): Name of the JSON configuration file.
  !!   top (char in): Top-level JSON key to look for in the config file.
  !!   campaign (char in): ICESat campaign, expected to be caps.
  !!   verbosity (int in opt): >0 print extra output
  !!   dir (char out): Directory name where data files for this campaign
  !!      are, read from the config file.
  !!   fileList (alloc char array out): List of data file names for this
  !!      campaign, read from the config file.
  !!
  !! Does not return a status because all the work is done by fson, which
  !! if it has an error USUALLY just prints a message and stops.  Callers
  !! probably should check that fileList is allocated, though.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE parse_envisat_config( configfile, top, campaign, verbosity, &
       dir, filelist )

    ! JSON parser:
    use fson, only: fson_value, fson_parse, fson_get, fson_destroy, fson_print
    use fson_value_m, only: fson_value_count, fson_value_get

    !!  Dummy args
    character(len=*), intent(in) :: configFile   ! JSON config file
    character(len=*), intent(in) :: top       ! top level JSON key in configfile
    character(len=*), intent(in) :: campaign     ! envisat campaign
    integer, intent(in), optional :: verbosity   ! >0: print extra output

    character(len=*), intent(out) :: dir       ! directory
    character(len=*), allocatable, intent(out) :: filelist(:)

    !! Local variables
    type( fson_value ), pointer :: configp       ! fson's configFile parse
    type( fson_value ), pointer :: fileListJArrp !JSON array ptr to the file list
    type( fson_value ), pointer :: felemjp       ! JSON ptr into fileListJArrp
    integer :: i
    integer :: nFiles   ! number of files for this camapaign
    character(len=lenstr) :: file
    integer :: verb       ! local copy of verbosity

!!!---------------------------------------------------------------

    print *
    print *, 'In envisat_config_mod%parse_envisat_config'

    verb = 0
    if (present(verbosity))  verb = verbosity

    dir = ''
    if (allocated(filelist))  deallocate (filelist)
    nFiles        = 0
    configp       => null()
    fileListJArrp => null()
    felemjp       => null()

    print *, "   Reading ConfigFile: " // trim(configFile)
    print *, "   for top-level JSON key: " // trim(top)

    ! TBD: VERIFY CAMPAIGN VAIDITY BASED ON VALUES IN FILE.  Bad value
    ! currently makes fson exit with "Unable to resolve path" trying to
    ! read "directory".  Maybe get campaign's value, then check if it's
    ! associated?  fson seems to do this internally.
    print *, "   looking for Campaign: " // trim(campaign)

    configp => fson_parse( configFile )
    if (.not. associated(configp)) then
       print *, 'envisat_config_mod%parse_envisat_config: Could not parse ConfigFile!'
       print *, 'Stopping.'
       stop 1
    end if

    call fson_get( configp, trim(top) // ".campaigns." // &
         trim(campaign) // ".directory", dir )
    if (verb >= 1) print *, "parse_envisat_config: ',  &
         'Dir = " // trim(dir)

    !! Unfortunately fson doesn't have a method to retrieve character arrays
    !! (see interface fson_path_get in fson_path_m.f90; fson_path_get is
    !! renamed to fson_get in use in fson.f90).  So we'll have to do it an
    !! item at a time.  This follows the example in fson's readme.md file.

    !! The fileList for this campaign, as an fson_value array pointer:
    call fson_get( configp, trim(top) // ".campaigns." // &
         trim(campaign) // ".fileList", fileListJArrp )
    if (.not. associated(fileListJArrp)) then
       print *, 'envisat_config_mod%parse_envisat_config: Could not get file list!'
       print *, 'Stopping.'
       stop 1
    end if
    if (verb >= 2) then
       print *,"parse_envisat_config: fileListJArrp="
       call fson_print( fileListJArrp )
    end if

    !! Create array
    nFiles = fson_value_count( fileListJArrp )
    if (nFiles == 0) then
       ! I'm not sure this can even happen; fson may just error and halt.
       print *, 'parse_envisat_config: No files found!'
       return
    end if

    allocate ( fileList(nFiles) )
    if (verb >= 1) print *,"parse_envisat_config: nFiles=", nFiles

    !!  Extract each of the elements of the fileList array
    do i = 1, nFiles
       ! i-th element of the fileList array, as an fson_value pointer:
       felemjp => fson_value_get( fileListJArrp, i )

       ! extract character string from this element:
       call fson_get( this=felemjp, value=file )
       fileList(i) = file
       if (verb >= 1) then
          print '(2a,i0,2a)', "parse_envisat_config: ", &
               "File ", i, " = ", trim(fileList(i))
       end if

    end do

    call fson_destroy( configp )

    return
  END SUBROUTINE parse_envisat_config


END MODULE envisat_config_mod
