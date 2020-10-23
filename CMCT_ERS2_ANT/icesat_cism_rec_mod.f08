! NAME: icesat_cism_rec_mod.f03
!
! PURPOSE: Defines the record structure for the subsetted, filtered GLA12
! files that Jack created for the CISM project.  These are the files that
! are in, eg., atlas:/data1/jack/I1_project/GLAS_Data/.  Based on Jack
! Saba's cism_record__define.pro.
!
! The intention here is to provide a record type that can be used to read
! these files, so it is defined as a SEQUENCE type and uses ISO_FORTRAN_ENV
! KIND constants (which makes it f08, but this should be supported by
! gfortran since version 4.5).  Field names are taken from Jack's
! cism_record__define.pro.
!
! AUTHOR: Jeff Guerber, GSFC 615/Sigma Space, 2015-06-01.
!
! HISTORY:

module icesat_cism_rec_mod

  use iso_fortran_env, only: INT32, REAL32, REAL64

  implicit none

  private
  public :: icesat_cism_rec

  type icesat_cism_rec
     sequence
     real( kind=REAL64 )   :: lat_degN
     real( kind=REAL64 )   :: lon_degE
     real( kind=REAL32 )   :: elev_m
     real( kind=REAL32 )   :: refl_uncorr
     real( kind=REAL32 )   :: wf_fit_sigma
     real( kind=REAL32 )   :: glas_dem
     real( kind=REAL32 )   :: gimp_dem
     integer( kind=INT32 ) :: gain
  end type icesat_cism_rec

end module icesat_cism_rec_mod
