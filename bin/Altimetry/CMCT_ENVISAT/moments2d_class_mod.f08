!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! NAME: moments2d_class_mod.f08
!!!
!!! PURPOSE: Object that calculates incremental sample means and variances
!!! on a 2D grid.  (Can be used for scalars and vectors by specifying
!!! dimensions of (1,1) and (n,1), but results will always be (:,:).)  Each
!!! value is added separately, and cells on the grid may have different
!!! numbers of values.  The current count, mean, variance, and standard
!!! deviation arrays can be retrieved at any time.  Uses the numerically
!!! stable algorithm by B.P.Welford (Technometrics 4(3), 1962), as in
!!! eg. Knuth's The Art of Computer Programming v.2 or
!!! http://en.wikipedia.org/Algorithms_for_calculating_variance (section
!!! "Online algorithm").
!!!
!!! USAGE:
!!!   type(moments2d_class) :: moms
!!!      ...
!!!   call moms%new(5,10)    ! 5x10 arrays
!!!   do          ! loop over data
!!!     ...
!!!     call moms%add_data(x,i,j)
!!!   end do
!!!   call moms%get(n=n, mean=mean, var=var)
!!!
!!! All calculations are done in Double Precision (IEEE Real64).
!!!
!!! NOTE: gfortran 6.4.3 doesn't yet have the ieee_arithmetic module, so
!!! undefined values of the real variables (eg. mean where n<1, variance
!!! where n<2) are flagged with a user-supplied value, which defaults to
!!! the netCDF default _FillValue, 9.96920996838686905e+36 = 15.*(2**119).
!!! (Note that this value converts to real32 without loss of precision.)
!!! May change this to NaN in the future.  Either way, be sure to check
!!! that the value of n(i,j) >= 1 or 2 as appropriate when using the
!!! returned arrays (perhaps with a WHERE statement or construct).
!!!
!!! It is possible to extend this to higher-order moments for skew,
!!! kurtosis, etc.; see the Wikipedia article above, and eg. Pebay (2008)
!!! Sandia Technical Report SAND2008-6212.  If doing this it might be
!!! advantageous to switch from Welford's formula to West's/Terriberry's
!!! (mom2_n = mom2 + ((n-1)/n)(x-mean)^2) which is equivalent, and save
!!! mom1=sum(x-mean) instead of mean.
!!!
!!! TERMINOLOGY: In the comments I've used "moment" to refer to the powers
!!! of (x-mean), but I think that's actually wrong and moments should refer
!!! to mean, variance, etc. themselves.  Needs to be rewritten for clarity.
!!!
!!! AUTHOR:  Jeff Guerber, SigmaSpace/MASA GSFC 615, 2015-07-28
!!! 2015-10-29 JRG: Retrieve invalid flag through get().
!!! 2016-05-11 JRG: Added subr destroy, and call it from new to
!!!   reinitialize an existing object.  Let user specify invflag.  Default
!!!   invflag is now netCDF default fill.  Elements of the class type are
!!!   now private.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE moments2d_class_mod

  use, intrinsic :: iso_fortran_env, only: INT32, REAL64

  implicit none
  private
  public moments2d_class

  !!
  !! BASE TYPE FOR THE CLASS
  !!
  type :: moments2d_class
     private    ! don't want user to change these; retrieve via get()

     ! ubi, ubj: Upper bounds of the arrays in the i and j directions
     integer :: ubi, ubj

     ! invflag: Value used to flag unset elemets.
     real(kind=REAL64) :: invflag

     ! n = number of points in each cell
     ! mean = 1st moment, sum(x)/n
     ! mom2 = 2nd moment, var=mom2/(n-1). Some references call it S_k or M_2.
     integer(kind=INT32), dimension(:,:), allocatable :: n
     real(kind=REAL64), dimension(:,:), allocatable :: mean
     real(kind=REAL64), dimension(:,:), allocatable :: mom2

   contains
     ! Type-bound procedures, which comprise the class methods
     procedure :: new
     procedure :: add_point
     procedure :: get
     procedure :: destroy

  end type moments2d_class

  !! Default flag value for invalid points, if invflag is not given.  This
  !! is the netCDF default fill value, NF90_FILL_DOUBLE = 15.*(2**119) =
  !! 9.96920996838686905e+36 (which has the same binary mantissa in REAL32
  !! so converts without losing precision).  Define it this way so we don't
  !! have to require netCDF for one value.
  real(kind=REAL64), parameter :: DEF_INVFLAG = scale(15.0_REAL64, 119)

CONTAINS

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! new: Initialize a new object or reinitialize an existing one.  Specify
  !! the sizes of and allocate the working arrays.
  !!
  !! Dummy input arguments:
  !!   ubi, ubj (integers): Upper bounds in the i and j directions. (Lower
  !!             bounds are not yet supported.)
  !!   invflag (optional real64): Value to use to flag unset or invalid
  !!             elements (reals only, not integers).  Defaults to
  !!             DEF_INVFLAG.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE new( self, ubi, ubj, invflag )
    ! self is the class variable, don't include as an actual argument
    class( moments2d_class ), intent(inout) :: self
    integer, intent(in) :: ubi, ubj
    real(kind=REAL64), optional, intent(in) :: invflag

    !! If self has already been initialized, deallocate the arrays so we
    !! can reallocate them.
    call self%destroy

    if ( present(invflag) ) then
       self%invflag = invflag
    else
       self%invflag = DEF_INVFLAG
    end if

    self%ubi = ubi
    self%ubj = ubj

    allocate ( self%n(ubi,ubj), self%mean(ubi,ubj), self%mom2(ubi,ubj) )
    self%n = 0
    self%mean = self%invflag
    self%mom2 = self%invflag

    return
  END SUBROUTINE new

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! add_point: Add a data point, updating the arrays.  Implements
  !! Welford's formula.
  !!
  !! Dummy arguments:
  !!   x:  New value to add at point (i,j).  Real64.
  !!   i,j: Indices of the point being added.  Default integer.
  !!   status: 0=success.  Def. integer, optional.
  !!   msg: If status /= 0, an error message explaing why.
  !!        Def. character, optional.
  !!
  !! NOTE: If the point falls outside the array and neither status nor msg
  !! are given, the point is ignored silently.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE add_point( self, x, i, j, status, msg )
    ! Dummy arguments
    ! self is the class variable, don't include as an actual argument
    class( moments2d_class ), intent(inout) :: self
    real(kind=REAL64), intent(in) :: x
    integer, intent(in) :: i, j
    integer, intent(out), optional :: status
    character(len=*), intent(out), optional :: msg

    ! Internal variables
    integer :: n      ! new value for self%n(i,j).
    real(kind=REAL64) :: mean      ! current self%mean(i,j)
    real(kind=REAL64) :: mean_n    ! new (nth) mean(i,j)
    real(kind=REAL64) :: mom2_n    ! new (nth) mom2(i,j)
    integer :: ios   ! iostat of the internal write

    if ( .not. allocated( self%n ) ) then
       if (present(status))  status = 2
       if (present(msg)) then
          msg = "!"         ! just in case string is too short
          ! iostat only so won't terminate program if msg is too short
          write (msg,'(a,6(i0,a))', iostat=ios) &
               "moments2d_class%add_point: Object is not initialized! "
       end if
       return
    end if

    if ( (i < 1) .or. (i > self%ubi) .or. (j < 1) .or. (j > self%ubj) ) then
       if (present(status))  status = 1
       if (present(msg)) then
          msg = "!"         ! just in case string is too short
          ! iostat only so won't terminate program if msg is too short
          write (msg,'(a,6(i0,a))', iostat=ios) &
               "moments2d add_point: Index out of range, (", &
               i, ",", j, ") not in (", 1,"-", self%ubi, ",", 1,"-", self%ubj, ")"
       end if
       return
    end if

    if (self%n(i,j) == 0) then
       ! 1st value added to this cell
       n = 1
       mean_n = x
       mom2_n = 0.0_REAL64
    else
       ! Welford's formula
       mean = self%mean(i,j)
       n = self%n(i,j) + 1
       mean_n = mean + ( (x - mean) / real(n,REAL64) )
       mom2_n = self%mom2(i,j) + (x - mean) * (x - mean_n)
    end if

    ! Update the arrays
    self%n(i,j) = n
    self%mean(i,j) = mean_n
    self%mom2(i,j) = mom2_n

    if (present(status))  status = 0
    if (present(msg))  msg = ""
    return
  END SUBROUTINE add_point

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! get: Retrieve the counts, means, variances, and standard deviations.
  !!
  !! Dummy arguments: All are optional allocatable (:,:) arrays except
  !! invflag.  It's OK if the actual arguments are already allocated, they
  !! will automatically be reallocated.
  !!   n (integer): Array of the number of data points in each cell.
  !!   mean (real64): Means.
  !!   var (real64):  Variances.
  !!   stddev (real64):  Standard Deviations (sqrt(var)).
  !!   invflag (real64): Value used to flag invalid values in mean, var, stddev.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE get (self, n, mean, var, stddev, invflag)
    ! Dummy arguments
    ! self is the class variable, don't include as an actual argument
    class( moments2d_class ), intent(inout) :: self
    integer, allocatable, optional, intent(out) :: n(:,:)
    real(kind=REAL64), allocatable, optional, intent(out) ::  &
         mean(:,:), var(:,:), stddev(:,:)
    real(kind=REAL64), optional, intent(out) :: invflag

    ! If not initialized, print error and return without changing anything.
    if ( .not. allocated( self%n ) ) then
       print *, 'moments2d_class%get: Object is not initialized!'
       return
    end if

    if (present(n)) then
       n = self%n
    end if

    if (present(mean)) then
       mean = self%mean
    end if

    if (present(var)) then
       ! variance
       var = self%mom2
       where (self%n >= 2)
          var = var / (self%n - 1)
       elsewhere
          var = self%invflag
       end where
    end if

    if (present(stddev)) then
       ! standard deviation = sqrt( var )
       stddev = self%mom2
       where (self%n >= 2)
          stddev = sqrt( stddev / (self%n - 1) )
       elsewhere
          stddev = self%invflag
       end where
    end if

    if (present(invflag)) then
       invflag = self%invflag
    end if

    return
  END SUBROUTINE get

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! destroy: Deallocate the arrays and do any other cleanup, if needed.
  !! Usually won't have to call this because the memory is freed
  !! automatically when the object goes out of scope, but it's here if you
  !! want to release it early.  Called by subr new to reinitialize an
  !! existing object.
  !!
  !! This could be a finalizer, but that's not necessary because arrays are
  !! deallocated automatically.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE destroy( self )

    ! self is the class variable, don't include as an actual argument
    class( moments2d_class ), intent(inout) :: self

    if ( allocated( self%n) )     deallocate ( self%n )
    if ( allocated( self%mean) )  deallocate ( self%mean )
    if ( allocated( self%mom2) )  deallocate ( self%mom2 )

    return
  END SUBROUTINE destroy

END MODULE moments2d_class_mod
