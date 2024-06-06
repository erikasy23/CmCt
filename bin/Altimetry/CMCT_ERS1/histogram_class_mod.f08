!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! NAME: histogram_class_mod.f09
!!!
!!! PURPOSE: Object for accumulating a histogram of a series of values.
!!! Calculates the histogram, probability density function (pdf), and
!!! cumulative (actually an empirical) distribution function (cdf).
!!! Currently data values must be accumulated one at a time, not as an
!!! array.  Note that the bins are closed below and open above: Eg, if bins
!!! are (1.0,2.0) and (2.0,3.0), then 2.0 falls in the 2nd bin and 3.0 is
!!! outside the range.  Also, the histogram, bins, pdf, and cdf arrays
!!! begin indexing at 0; this makes bin assignments easier.  The bins array
!!! gives the bounds of each bin: histogram(i) has bounds
!!! bins(i):bins(i+1).  Thus bins has one more element than the other
!!! arrays.
!!!
!!! USAGE: Initialize by calling histogram%new(...) with appropriate bounds
!!! and binsize.  Then for each value x to be inserted into the histogram,
!!! call histogram%accumulate(x).  To retrieve the histogram, call
!!! histogram%get(...).  New values may be added after a call to get().
!!!
!!! CAUTION: Currently NO checking for invalid values is done; they must be
!!! filtered before the value is accumulated.  This is especially true if
!!! the pdf or cdf are desired, since their calculations depend on the
!!! counts below and above the histogram bounds.
!!!
!!! AUTHOR: Jeff Guerber, SigmaSpace/NASA GSFC 615, 2015-12-29
!!! 2016-01-04 JRG: Added cdf.  Pdf calculation includes nbelow+nabove.
!!! 2016-05-11 JRG: Added subroutine destroy, and call it from new to
!!!   reinitialize existing object.  In accumulate and get, check that
!!!   object has been allocated.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE histogram_class_mod

  use, intrinsic :: iso_fortran_env, only: INT32, REAL32, REAL64

  implicit none
  private
  public histogram_class

  !!
  !! BASE TYPE FOR THE CLASS
  !!
  type :: histogram_class

     private    ! don't want user to change these; retrieve via get()
     real(kind=REAL32) :: min, max, binsize
     integer(kind=INT32) :: nbins    ! number of bins

     ! bins and histogram index from 0 to make bin assignment easier.
     ! Boundaries of bin histogram(i) are bins(i):bins(i+1), so bins() has
     ! 1 more element than histogram().
     real(kind=REAL32), allocatable   :: bins(:)      ! bin boundaries, 0:nbins
     integer(kind=INT32), allocatable :: histogram(:)  ! counts, 0:nbins-1
     ! total counts, counts below, counts above hist limits
     integer(kind=INT32) :: n, nbelow, nabove

   contains
     ! Type-bound procedures, which comprise the class methods
     procedure :: new
     procedure :: accumulate
     procedure :: get
     procedure :: destroy

  end type histogram_class

CONTAINS

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! new: Intitialize a new object or reinitialize an already existing one.
  !! Allocate arrays and calculate bins.
  !!
  !! Dummy arguments:
  !!   min, max (in):  Requested range of the histogram
  !!   binsize (in):  Size of each histogram bin
  !!   adjmin, adjmax (opt out): Adjusted range of the histogram, after
  !!      calculating the bins based on binsize.  Adjmax in particular may
  !!      differ from Max.
  !!   nbins (opt out): Number of bins.
  !!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE new( self, min, max, binsize, adjmin, adjmax, nbins )

    ! self is the class variable, don't include as an actual argument
    class( histogram_class ), intent(inout) :: self
    real(kind=REAL32), intent(in)  :: min, max, binsize
    real(kind=REAL32), intent(out), optional :: adjmin, adjmax
    integer(kind=INT32), intent(out), optional  :: nbins

    integer(kind=INT32) :: i

    !! If self has already been initialized, deallocate the arrays so we
    !! can reallocate them.
    call self%destroy

    !! Number of bins.  Uses ceiling() instead of floor() so the top bin isn't
    !! truncated. (Consider floor vs ceiling for eg. min=1.,max=5.,binsize=3.)
    self%nbins = ceiling( (max - min) / binsize )

    self%binsize = binsize
    self%min = min
    self%max = self%min + self%nbins * self%binsize    ! adjusted
    self%n = 0
    self%nbelow = 0
    self%nabove = 0

    !! Allocate and initialize arrays
    allocate ( self%histogram(0:(self%nbins-1)) )
    self%histogram = 0

    allocate ( self%bins(0:self%nbins) )
    self%bins = [ (i * self%binsize + self%min, i=0,self%nbins) ]

    !! If requested, return optional args
    if (present(adjmin)) adjmin = self%min
    if (present(adjmax)) adjmax = self%max
    if (present(nbins))  nbins  = self%nbins

    return
  END SUBROUTINE new


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! accumulate: Put a value into the histogram.
  !!
  !! Dummy arguments:
  !!    value (in): Value to add.
  !!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE accumulate( self, value )

    ! self is the class variable, don't include as an actual argument
    class( histogram_class ), intent(inout) :: self
    real(kind=REAL32), intent(in)  :: value

    integer(kind=INT32) :: ibin

    if ( .not. allocated(self%bins) ) then
       print *, 'histogram_class_mod%accumulate: Object is not initialized!'
       return
    end if

    ibin = floor( (value - self%min) / self%binsize )
    if (ibin < 0) then
       self%nbelow = self%nbelow + 1
    else if (ibin >= self%nbins) then
       self%nabove = self%nabove + 1
    else
       self%histogram(ibin) = self%histogram(ibin) + 1
       self%n = self%n + 1
    end if

    return
  END SUBROUTINE accumulate


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! get: Retrieve the histogram and related items.
  !!
  !! Dummy arguments:
  !!   Most of these are just the object variables.  Histogram and bins
  !!      must be allocatable (:).  In addition:
  !!   pdf (alloc,opt): Probability density function,
  !!      histogram/(n*binsize). Allocatable and same shape as histogram.
  !!   cdf (alloc,opt): Cumulative distribution function.  Allocatable and
  !!      same shape as histogram.  cdf(i) = (sum(histogram(j<=i))+nbelow)/n
  !!
  !! OK for actual argments to already be allocated, they will
  !! automatically be reallocated.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE get( self, min, max, binsize, nbins, bins, histogram, &
       n, nbelow, nabove, pdf, cdf )

    ! self is the class variable, don't include as an actual argument
    class( histogram_class ), intent(inout) :: self
    real(kind=REAL32), intent(out), optional :: min, max, binsize
    real(kind=REAL32), allocatable, intent(out), optional :: bins(:)
    integer(kind=INT32), allocatable, intent(out), optional :: histogram(:)
    integer(kind=INT32), intent(out), optional :: nbins, n, nbelow, nabove
    real(kind=REAL32), allocatable, intent(out), optional :: pdf(:), cdf(:)

    ! local variables
    real(kind=REAL64) :: total
    integer :: lb, ub, i

    if ( .not. allocated(self%bins) ) then
       print *, 'histogram_class_mod%get: Object is not initialized!'
       return
    end if

    lb = lbound(self%histogram,1)
    ub = ubound(self%histogram,1)

    if (present(min)) min = self%min
    if (present(max)) max = self%max
    if (present(binsize)) binsize = self%binsize
    if (present(nbins)) nbins = self%nbins

    if (present(n))      n      = self%n
    if (present(nbelow)) nbelow = self%nbelow
    if (present(nabove)) nabove = self%nabove

    ! Histogram and bins: Shouldn't need to deallocate if they already
    ! exist, automatic reallocation applies.
    if (present(histogram)) histogram = self%histogram
    if (present(bins)) bins = self%bins

    if (present(pdf)) then
       if ( allocated(pdf) )  deallocate ( pdf )
       allocate( pdf(lb:ub) )
       ! Note that the additions are done in integers, mults and divs in dp,
       ! only final result is sp, to avoid any precision loss.
       total = real( ( self%n + self%nabove + self%nbelow ), REAL64)
       pdf = real(self%histogram,REAL64) / (total * real(self%binsize,REAL64))
    end if

    if (present(cdf)) then
       if ( allocated(cdf) )  deallocate ( cdf )
       allocate( cdf(lb:ub) )
       ! Note that the additions are done in integers, mults and divs in dp,
       ! only final result is sp, to avoid any precision loss.
       total = real( (self%n + self%nabove + self%nbelow), REAL64)
       forall (i = lb:ub)
          cdf(i) = real( ( sum(self%histogram(lb:i)) + self%nbelow ), REAL64)
       end forall
       cdf = cdf / total
    end if

    return
  END SUBROUTINE get


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! destroy: Deallocate the arrays and do any other cleanup, if needed.
  !! Usually won't have to call this because the memory is freed
  !! automatically when the object goes out of scope, but it's here if you
  !! want to release it early.  Called by subr new to reinitialize an
  !! existing object.
  !!
  !! This could be a finalizer, but that's not necessary because arrays are
  !! deallocated automatically.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE destroy( self )

    ! self is the class variable, don't include as an actual argument
    class( histogram_class ), intent(inout) :: self

    if ( allocated( self%bins) )       deallocate ( self%bins )
    if ( allocated( self%histogram) )  deallocate ( self%histogram )

    return
  END SUBROUTINE destroy


END MODULE histogram_class_mod
