!!! Not quite working yet - esp: ldir, lsuffix don't always seem to get the
!!! right lengths
!!!
!!! NAME: temp_file_mod.f08
!!!
!!! PURPOSE: Return the name for a temporary file.  The name will be a
!!! combination of the pid, user-supplied name, and system time since the
!!! program began.
!!!
!!! NOTE: Uses the GNU extension getpid().
!!!
!!! AUTHOR: Jeff Guerber, Sigma Space/GSFC 615, 2015-06-26

module temp_file_mod

  use iso_fortran_env
  implicit none

contains

  function get_temp_name( name, dir, suffix ) result(res)

    character(len=:), allocatable :: res

    !! Dummy argments
    character(len=*), intent(in) :: name
    character(len=*), optional, intent(in) :: dir, suffix

    !! Local variables
    integer (kind=INT64), save :: ftime = -1   ! System time on first call
    integer (kind=INT64) :: stime, dtime       ! Current time, difference
    integer  :: pid                            ! process id
    character(len=512) :: tmpname     ! string used to construct the file name
    integer :: n                      ! num non-blanks in tmpname
    character(len=:), allocatable :: ldir, lsuffix  ! local dir, suffix

    print *, 'In get_temp_name: name>',name,'<  dir>',dir,'<  suffix>',suffix,'<'
    ! Default values for the optional dummy args
    if ( present(dir) ) then
       n = len_trim(dir)
       if ( dir(n:n) /= '/' ) then
          ! Add a / if there isn't one already
          ldir = trim(adjustl(dir)) // '/'
       else
          ldir = trim(adjustl(dir))
       end if
    else
       ldir = ''
    end if
    print *,'ldir='//ldir//'='
    if ( present(suffix) ) then
       lsuffix = trim(adjustl(suffix))
    else
       lsuffix = ''
    end if
    print *,'lsuffix='//lsuffix//'='

    ! Process ID.  Note this is a GNU extension.
    pid = getpid()

    ! Reading of the system clock the first time this routine is called.
    if (ftime < 0) then
       call system_clock( count=ftime )
    end if

    call system_clock( count=stime )
    dtime = abs( stime - ftime )     ! abs just in case it rolls over

    !! Now assemble the result
    tmpname = ''
    write (tmpname, '(a,i0,3a,i0,a)') &
         trim(ldir), pid, '_', trim(adjustl(name)), '_', &
         dtime, trim(lsuffix)
    n = len_trim(tmpname)
    res = tmpname(:n)
    print *, 'trim(tmpname),len = ', '>'//trim(tmpname)//'<', len(trim(tmpname))
    print *, 'res,len = ','>'//res//'<', len(res)

    return
  end function get_temp_name

end module temp_file_mod
