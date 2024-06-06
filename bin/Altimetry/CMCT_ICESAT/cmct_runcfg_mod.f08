!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! NAME: cmct_runcfg.f08
!!!
!!! PURPOSE: Read and parse the CMCT run configuration JSON file.  Stores
!!! everything in a derived type.
!!!
!!! AUTHOR: Jeff Guerber (JRG), SigmaSpace/GSFC 615, Mar. 2016
!!! HISTORY:
!!! 2016-04-13 JRG: Update for latest json file, esp. adding
!!!   model.model_time_index and run.upload_dir.  Read all items, not just
!!!   those the program expects to use, and include them in the return
!!!   structure.  Simplify by making better use subobject specifiers so
!!!   need fewer fson_value pointers.  Simplified naming of several
!!!   variables.
!!! 2016-04-28 JRG: Add nComps to type::cmct_run_config
!!! 2016-05-11 JRG: parse_runcfg: If cfg%comps already allocated, deallocate.
!!! 2016-05-13 JRG: parse_runcfg: Forgot to call fson_destroy(configp)
!!! 2016-07-02 JRG: Made user_run_title longer.
!!!
!!! Last svn commit: $Id: cmct_runcfg_mod.f08 105 2016-07-02 10:46:41Z jguerber $
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module cmct_runcfg_mod

  use, intrinsic :: iso_fortran_env, only: INT32, REAL32, REAL64
  use :: fson
  use :: fson_value_m

  implicit none
  private
  public cmct_run_config, parse_runcfg

  !!
  !! Type definitions.  These are intended to follow the
  !! cmct_run_config.json file.  (Although we don't necessarily need
  !! everything in there, eg. run.upload_dir.)  cmct_run_config is the main
  !! type.  (Too bad these can't nest!  Or can they?)
  !!
  !! There's nothing sacred about the string lengths here, they just seem
  !! reasonable.

  ! "run" section
  type :: cmct_runcfg_run
     character(len=50)  :: runid
     character(len=25)  :: date
     character(len=150) :: upload_dir
     character(len=100) :: user_run_title
     character(len=25)  :: loginname
     character(len=50)  :: actualusername
     character(len=50)  :: email
  end type cmct_runcfg_run

  ! "comparisons.model" sections
  type :: cmct_runcfg_comparisons_model
     character(len=150) :: modelname  ! model's name
     character(len=150) :: filename   ! model file's name
     character(len=15)  :: format     ! model file fmt, "cism-text" or "netcdf"
     real(kind=REAL64)  :: spacing_km ! spacing of grid cells, in km
     character(len=15)  :: region     ! Region, antarctica or greenland
     character(len=15)  :: variable   ! name of variable to use from model file
     character(len=1000) :: user_comments ! user's comment string
     integer(kind=INT32) :: model_time_index  ! ignored for cism-text
  end type cmct_runcfg_comparisons_model

  ! "comparisons.observations" sections
  type :: cmct_runcfg_comparisons_obs
     character(len=25) :: mission    ! mission name
     character(len=4)  :: campaign   ! for Icesat, campaign
     character(len=50) :: dataset    ! name of the dataset
  end type cmct_runcfg_comparisons_obs

  ! "comparisons" sections
  type :: cmct_runcfg_comparisons
     type(cmct_runcfg_comparisons_model) :: model
     type(cmct_runcfg_comparisons_obs)   :: obs
  end type cmct_runcfg_comparisons

  ! Full file
  type :: cmct_run_config
     type(cmct_runcfg_run) :: run
     integer(kind=INT32)   :: nComps      ! num elems in comps, not a field
     type(cmct_runcfg_comparisons), allocatable :: comps(:)
  end type cmct_run_config

contains

  !!
  !! Function parse_runcfg: Read and parse the run configuration file.
  !! Returns a structure of type cmct_run_config.
  !!

  function parse_runcfg( runCfgFileName ) result( cfg )

    !! Return value and dummy argument
    type(cmct_run_config)        :: cfg               ! return value
    character(len=*), intent(in) :: runCfgFileName    ! input file name

    ! fson pointers to components of the Run Config file
    type( fson_value ), pointer :: &
         configp=>null(), &    ! top level tree
         runInfop=>null(), &   ! "run" section
         compArrp=>null(), &   ! "comparisons" array
         compp=>null()         ! current comparison, element of compArrp

    integer :: i, nCompares

    ! Parse the file
    configp => fson_parse( runCfgFileName )

    ! fson ptr to the RunInfo section
    call fson_get( configp, 'cmct_run_config.run', runInfop )
    call fson_get( runInfop, 'runid', cfg%run%runid )
    call fson_get( runInfop, 'date', cfg%run%date )
    call fson_get( runInfop, 'upload_dir', cfg%run%upload_dir )
    call fson_get( runInfop, 'user_run_title', cfg%run%user_run_title )
    call fson_get( runInfop, 'loginname', cfg%run%loginname )
    call fson_get( runInfop, 'actualusername', cfg%run%actualusername )
    call fson_get( runInfop, 'email', cfg%run%email )

    ! fson ptr to the array of comparisons in the run config file
    call fson_get( configp, 'cmct_run_config.comparisons', compArrp )

    !! Loop over comparisons in this file
    cfg%nComps = fson_value_count( compArrp )
    if ( allocated( cfg%comps ) )  deallocate ( cfg%comps )
    allocate ( cfg%comps(cfg%nComps) )
    do i = 1, cfg%nComps
       compp => fson_value_get( compArrp, i )

       !! Get the gridded model
       call fson_get( compp, 'model.modelname', &
            cfg%comps(i)%model%modelname )
       call fson_get( compp, 'model.filename',&
            cfg%comps(i)%model%filename )
       call fson_get( compp, 'model.format', &
            cfg%comps(i)%model%format )
       call fson_get( compp, 'model.spacing_km', &
            cfg%comps(i)%model%spacing_km )
       call fson_get( compp, 'model.region', &
            cfg%comps(i)%model%region )
       call fson_get( compp, 'model.variable', &
            cfg%comps(i)%model%variable )
       call fson_get( compp, 'model.user_comments', &
            cfg%comps(i)%model%user_comments )
       call fson_get( compp, 'model.model_time_index', &
            cfg%comps(i)%model%model_time_index )

       !! Observation
       !!
       !! Not every mission uses every field.  Probably should initialize
       !! them to "DUMMY", to flag in missions that don't.

       call fson_get( compp, 'observations.mission', cfg%comps(i)%obs%mission )

       select case ( trim(cfg%comps(i)%obs%mission) )
       case ( 'icesat-glas' )
          call fson_get( compp, 'observations.campaign', &
               cfg%comps(i)%obs%campaign )
          call fson_get( compp, 'observations.dataset', &
               cfg%comps(i)%obs%dataset )

       case default
          print *, 'parse_runcfg: Unknown mission=' // &
               trim(cfg%comps(i)%obs%mission) // '. Stopping.'
          stop 1
       end select

    end do

    call fson_destroy( configp )

    return
  end function parse_runcfg

end module cmct_runcfg_mod
