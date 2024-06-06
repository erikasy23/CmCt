!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! NAME: icesat_cmct_rec_mod.f08
!!!
!!! PURPOSE: Defines the record structure for the new (March 2016)
!!! subsetted and filtered GLAS elevation files that Jack Saba created for
!!! the CMCT project.  These are the files that are in, eg., (for now)
!!! gs6141-atlas:/data1/jack/CMCT/GLAS_Data/R634_Greenland/, with names
!!! like GLA12_634_2103_001_1135_0_01_0001.CMCT.final.dat.  Based on Jack's
!!! cmct_record__define.pro.
!!!
!!! The intention here is to provide a record type that can be used to read
!!! these files, so it is defined as a SEQUENCE type and uses KIND
!!! constants from ISO_FORTRAN_ENV.  Field names and comments are mostly
!!! taken from Jack's cmct_record__define.pro.
!!!
!!! These definitions are also compatible with the variables in the NetCDF
!!! versions of the CMCT GLAS elevation files (April 2016).
!!!
!!! AUTHOR: Jeff Guerber, GSFC 615/Sigma Space, 2016-03-24.
!!!
!!! HISTORY:
!!! 2016-05-19 JRG: Minor update based on Jack's latest cmct_record__define.pro
!!!    and attributes in the NetCDF versions of the data files.
!!!
!!! Last SVN commit: $Id: icesat_cmct_rec_mod.f08 79 2016-05-25 07:41:05Z jguerber $
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module icesat_cmct_rec_mod

  use iso_fortran_env, only: INT32, REAL32, REAL64

  implicit none

  private
  public :: icesat_cmct_rec

  type icesat_cmct_rec
     sequence

     ! Shot time in seconds since J2000.0 from GLAS record
     real( kind=REAL64 )   :: UTCTime_j2000s

     ! Lat from GLAS record
     real( kind=REAL64 )   :: Lat_degN

     ! Lon from GLAS record
     real( kind=REAL64 )   :: Lon_degE

     ! Elevation relative to Topex/Poseidon ellipsoid from GLAS record
     real( kind=REAL32 )   :: TopexElev_m

     ! Elevation relative to WGS84 ellipsoid from GLAS record
     real( kind=REAL32 )   :: WGS84Elev_m

     ! Elevation relative to EGM2008 mean-tide geoid
     real( kind=REAL32 )   :: EGM08MeanTideElev_m

     ! Elevation relative to EGM2008 tide-free geoid
     real( kind=REAL32 )   :: EGM08TideFreeElev_m

     ! From DiMarzio study of uncertainties
     real( kind=REAL32 )   :: Elev_StdDev_m

     ! From GLAS slope grid
     real( kind=REAL32 )   :: Slope_deg

     ! Uncorrected reflectivity from GLAS record
     real( kind=REAL32 )   :: Refl_Uncorr

     ! Std dev of fit to wf from GLAS record
     real( kind=REAL32 )   :: WF_Fit_Sigma

     ! Gain from GLAS record
     integer( kind=INT32 ) :: Gain

     ! DEM elev from GIMP 90-m data (Grn) or Bamber 1 km DEM (Ant)
     integer( kind=INT32 ) :: DEM_m

     ! Drainage subsystem (1.0-8.2) from 2013 1 km grid, based on Mario
     ! Giovinetto's new DS boundaries and GIMP 90-m ice mask.  For Grn
     ! there are points on isolated ice caps and islands that have DS = 0.
     real( kind=REAL32 )   :: DrainageSystem

     ! Ice type flag.
     !   0 = isolated ice cap (Greenland only)
     !   1 = continental ice
     !   2 = ice shelf (Antarctica only)
     integer( kind=INT32 ) :: IceType

  end type icesat_cmct_rec

end module icesat_cmct_rec_mod
