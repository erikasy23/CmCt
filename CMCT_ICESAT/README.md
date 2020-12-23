* CMCT notes

** Requirements 

Basic dependencies:
- gcc-gfortran
- ksh
- jq
- make
- git
- unzip
- wget

CentOS 7 RPMs:
- libaec and libaec-devel (1.0.4-1)
- hdf5 and hdf5-devel (1.8.12-11)
- netcdf and netcdf-devel (4.3.3.1)

from Unidata's github repo:
- netcdf-fortran (4.3.3), provides libnetcdff.so and netcdf.mod

fson from https://github.com/josephalevin/fson
- Alternative json-fortran needs gcc/gfortran 4.9

For further information and history refer to contents of docs subdirectory.

** .mod files

*** TODO Is this the right way to build .mod files?
- gfortran doesn't create new .mod file if interfaces haven't changed.
- Maybe should have their own rules, perhaps using:
-- rm foo.mod                         # force new .mod file
-- gfortran -fsyntax-only foo.f08     # creates .mod file
- App. D of Modern Fortran Explained is about building with modules.
-- Maybe not always recreating the .mod file is the right thing?


** DONE Launch script.
- periodically polls web site
- retrieves new model submissions
- creates working directory
- writes cmct_run_config_*.json
- writes cmct_main_control.json
- spawns cmct_main
- Done (at least for running on ggsghpcc): cmct_main_launch.ksh,
  cmct_main_launch_config.ksh

** http://astroplotlib.stsci.edu/index.html
Potentially useful plot library, some in Python.


** LONGTERM Better way to organize the program?
Perhaps a better way to organize the program would be more like building
blocks.  Have a bunch of modules that do various functions, then for each
comparison type/dataset/etc. a module that calls the appropriate blocks.
This might shorten the individual modules, too.


** TODO Add capability to specify time ranges.
- Especially needed for radar datasets.
- Sophie would like this even for ICEsat, so can combine campaigns.
- Probably need multiple ranges, so can specify eg. all fall data.


========================================================================


** cmct_main.f08
Main program.  Contains the loop over observation records.

gfortran -c -g -O0 cmct_main.f08 -I ./include -I /usr/include
- fson.mod is in ./include
- netcdf.mod in in /usr/include

gfortran -o cmct_main cmct_main.o cmct_runcfg_mod.o  \
  cism_polarstereo_mod.o datasets_mod.o  \
  histogram_class_mod.o icesat_cism_class_mod.o icesat_cism_rec_mod.o   \
  icesat_cmctnc_class_mod.o icesat_cmct_rec_mod.o icesat_config_mod.o   \
  kinds_mod.o model_cism_txt_mod.o model_netcdf_mod.o moments2d_class_mod.o \
  polar_stereographic_mod.o -L ./lib -lfson -lnetcdff

*** DONE Overall statistics
Just make another instance of moments2d, with a 1x1 "grid".

*** DONE Grid-data output file

*** DONE Observation-record output file
**** TODO Write this as a netcdf file with record variables.

*** DONE Clean obs data

*** DONE Output file names.
- Settle on convention.
- If convention includes date, don't put date in filename separately.
- Should dates be ISO, or Unix 1970.0 seconds?
- (user base).CMCT_(comparison#).(file type) ?
- Include run id?  Could make them too long.
- Convention decided on: CMCT_(runid)_(comparison).(filetype).(ext)
**** DONE May change to ..._(runid).(comparison)....

*** TODO Encapsulate grids
Eg. give model objects standard method eg. "find_cell" to find which
cell a lat/lon is in.

*** TODO More comments, esp. in header.

*** DONE Get reclist headers from datasets and model modules
- Most information now from the runcfg, modelInfo, obsInfo structures.

*** DONE meansdgrids.nc should have units attributes

*** TODO meansdgrids.nc should have identifying global attributes

*** TODO meansdgrids.nc should follow CF conventions
- Some potentially useful information on this is at the SeaRISE wiki:
-- http://websrv.cs.umt.edu/isis/index.php/Output_Format
-- http://websrv.cs.umt.edu/isis/index.php/CF_standard_names_for_Glaciology_and_Ice-Sheet_Modeling

*** TODO In meansdgrids.nc, should latitude and longitude have values at all grid points?
- Currently, these have invalid value when outside Greenland too.
- In number array, put 0 where inside model but no points, invalid where
  outside model.  Invalids give good idea of model's extent.

*** DONE meansdgrids.nc should have global attrs for time index, time

**** MAYBE meansdgrids.nc should have time axis
- Instead of time in global attribute
- Would only be 1 element

*** TODO Make meansdgrids.nc arrays conformant with model
- Currently they are 1 element shorter in both directions because are cells
  not grid points.
- Jack suggests adding an extra row and column filled with invalids.

*** DONE Make histograms
**** bilinear_compare.pro (actually cism_results.pro) hardcodes histogram range
to [-700,1000]m (overrides histrange argument) by 25m; also plots central
90% at 1m.
**** IDL: bins = (0,floor((max-min)/binsize)), bin(x) = floor((x-min)/binsize).
For this case, (1000. - -700.)/25. = 68.0, so there are 69 bins.
However, cuts off at max, so last bin may not be as expected: consider
IDL> h=histogram(findgen(10),min=0.,max=6.,binsize=2.,loc=l)
Expect each bin to contain 2?
IDL> print,l,h
      0.00000      2.00000      4.00000      6.00000
           2           2           2           1
**** Better one here: http://badgrads.berkeley.edu/doku.php?id=idl:histograms
**** Probably will accumulate into 1m bins, then sum into 25m.
**** Should this be done in post-processing instead?

**** Implemented as histogram_class_mod.f08

*** DONE Include footprint distance to nearest node.
**** bilinear_compare.pro does this, and scatter-plots dElev vs. D.
Actually done in cism_getdata.pro.  Rounds the footprint x,y to get closest
node, takes distance to there.  Have already filtered out records
not surrounded by 4 good model points.  This is only thing the rounded x,y
is used for.  Jack writes it to 4 decimal places.  Jack does not check if
exceeds expected maximum.

*** DONE Check that iTime is in index range of model_T.
- If someone enters eg. 0, program segfaults.

*** cmct_main_control.json
Current version of main control file.  File name read from command line
argument of cmct_main.  Intended to be created by launch script.  Names
working_dir, cmct_run_config file.

**** cmct_main_control_p03.json
Prototype #3 for main control file.  Version control as
cmct_main_control.json svn Rev r13.

**** cmct_main_control_p04.json
Prototype #4 for main control file.  Adds verbosity, makes comments a
single array.  Version control as cmct_main_control.json svn Rev r14.

**** TODO Add run date.
Print it in log.

========================================================================

** cmct_runcfg_mod.f08
Read and parse the CMCT run configuration JSON file.  Stores everything in
a nested derived type.

*** function parse_runcfg
Function that reads and parses the file.  Returns a type(cmct_run_config)
structure.

*** cmct_run_config.json
Template for the run configuration file.  Name is specified in the main
control file.  Intended to be created by the web page, obtained by launch
script.  Contains info on run (id, email, title), model file (name,
filename, region, variable, etc), requested comparison observation
(mission, dataset, etc).

**** cmct_run_config_proto04.json
Prototype #4 for run configuration file.  Version control as
cmct_run_config.json svn Rev r10.
**** cmct_control_proto*.json: previous versions, changed name.

**** DONE Does comparisons.observations need separate mission, mission-options?
- or can they be combined?
- Have combined them

**** MAYBE Maybe just drop mission, organize by dataset?
- But then wouldn't have just mission for attributes, titles, etc.
- This could be in the datasets.json file.

**** TODO May need a user file title for each comparison

**** TODO user_comments would be better called "description"

**** TODO Have separate submit, run dates.
Run date should probably be in cmct_main_control.json.  Date in this file
should be submission date.

========================================================================

** cism_polarstereo_mod.f03

Initialize polar_stereographic_mod for CISM grid, at given spacing.

*** TODO Convert to standard ISMIP grid?: spacing, scale, units.  Xmin/Ymin?
Scale is "1 unit of x or y is how many units-of-a?".

*** TODO Generalize?: CISM grid at given spacing, ISMIP grid at std spacing.

*** TODO Add subr cism_ps_init(region,spacing,status)
To select parameters constant based on region.

*** type ps_params_def
Derived type for PS grid parameters
**** TODO Add type="...", eg. type="cism_std_grn"

*** subr cism_ps_init_std_grn
Initialize for CISM grid for Greenland.

*** func cism_ps_get_params
Return structure with current parameters.

*** LONGTERM Add Antarctica

========================================================================


** model_cism_txt_mod.f08

Model files in CISM Text format, such as by Steve Price.

Used to be cism_txt_file_mod.f08, renamed because I kept confusing the cism
modules for models and for observations.

gfortran -c -g -O0 model_cism_txt_mod.f08 -I /usr/include

*** subr read_model_cism_txt

Read CISM text format model file, returns data on grid.
Calculates grid based on min lat/lon, CISM grid parameters, supplied spacing.

Used to be read_cism_txt_file.

**** DONE Error handling
**** TODO Check that calculated grid points are close to X,Y in file

*** (invalid value, was CISM_INVALID_R64, no longer public constant)
Value returned for invalid grid points.
**** DONE Change to more standard value such as default netCDF invalid.
Now using NF90_FILL_DOUBLE from module netcdf.
**** CANCELED Revert netCDF invalid?
- Why bring in netCDF just for this?  Since now return value, can be
  anything.  Maybe supply value desired when init'd?
- NCO User's Guide (sec 4.2 footnote 1) strongly discourages using NaN as
  netcdf _FillValue, not sure why.  Availability of IEEE support?
  (gfortran 4.6.3 supports IEEE for I/O and (nan==nan) == .false..)
- Sophie prefers netCDF default invalid [CMCT Weekly 2016-03-13]
- Now calculating value of NF90_DEFAULT_FILL (15*(2**119)) so doesn't need
  the netCDF library any more [2016-05-11]

**** DONE More general name? Or return as a dummy arg of read_cism_txt_file?
So caller doesn't need to know name "CISM_INVALID_R64".
- Changed to intent(out) dummy arg "invflagval".

**** TODO Should x, y be calculated, instead of set from file lats & lons?
Could avoid possibility of invalid values, which seem to be prohibited in
coordinate variables in the netCDF file.

========================================================================


** model_netcdf_mod.f08

Model files in netCDF format.

*** read_model_netcdf

Read netCDF format model file.


========================================================================


** datasets_mod.f08

Available observation datasets.

*** DONE Expand into generic interface to all datasets.
- Could have subr that returns (lat,lon,value) from open dataset.
- Maybe make it a base class for all dataset classes?  Then could use a
  dynamic polymorphic variable.
- Overload a "datasets%get()" that would call appropriate open/get/etc.,
  eg. icesat_cism%next_rec.
- As of r92, also defines class datasets_class with methods open and next_obs.

*** DONE Obey verbosity

*** subr datasets_init
Reads datasets.json and stores the info.

*** subr datasets_info
Returns information about one of the datasets (currently just the config
file name and the dataset description).

*** type datasets_class
Base type for the class.  Access the datasets in a generic way.

*** subr open
Type-bound method.  Opens a dataset in this object instance.

*** subr next_obs
Type-bound method.  Returns the next (lat, lon, value) from the
currently-open dataset.

**** TODO May need to add checks for the various datasets, because they may contain multiple items.

*** datasets.json
- The data file used by datasets_mod.f08.
- Contains list of known datasets, and name of config file and description
  for each.
- Currently includes:
  icesat-cism-elev-grn
  icesat-cmctnc-tpelev-grn
  icesat-cmctnc-wgselev-grn
  icesat-cmctnc-egmmtelev-grn
  icesat-cmctnc-egmtfelev-grn

**** LONGTERM Add more datasets.


========================================================================

** histogram_class_mod.f08

Class for accumulating histograms.  Also does probability density function
(pdf) and empirical cumulative distribution function (cdf).  Cdf should be
useful for estimting percentiles.

*** subr new
Type-bound method. Initializes, calculates bins, allocates working arrays.

*** subr accumulate
Type-bound method. Adds a value into the histogram.

*** subr get
Type-bound method. Retrieve the histogram, bins, pdf, cdf, and other items.

*** subr destroy
Type-bound method. Deallocate arrays.


========================================================================

** icesat_cism_class_mod.f08

Class for reading the ICESat GLA12-subset files Jack created for the CISM
comparison project.

*** TODO error handling.

*** icesat_cism_class
Base type for the class.

*** init
Type-bound method.  Initializes the object, gets the list of data files,
and opens the first data file.

*** next_rec
Type-bound method.  Returns the next record from the observation data set.
Purposely does not call self%clean() itself (why was that?), let the caller
do that.

**** DONE Clean data
- Return a good/reject flag along with the data.
- Maybe make this its own method, called by next_rec.
- But what cleaning needs to be done?  read_glas_data.pro halts if ANY
  abs(elev_m - gimp_dem) > 200; cism_getdata.pro (which calls both
  read_glas_data and read_cism) eliminates points where elev_m > (3300m for
  Grn, 4500m for Ant), point x/y is outside min/max cism grid, or point is
  not surrounded by 4 good grid points.
- IMPLEMENTED as method clean(), see below.

*** next_file
Type-bound method. Closes the current file in the observation data set and
opens the next one.

*** clean
Type-bound method. Checks if the given record meets various criteria,
returns a code indicating if it's good or why it failed.
- Elevation over Greenland < 3300 m
- Elevation over Antarctica < 4500 m
- Elevation is within 200 m of GIMP DEM

*** subr destroy
Type-bound method. Close file, deallocate arrays, free fson pointers.

*** ICM_*
Status constants.

*** icesat_cism_config.json
Configuration file.  Lists all the specific files for each campaign.


========================================================================

** icesat_cism_rec_mod.f08
Defines file record for Jack's GLA12 subset files for CISM.


========================================================================

** icesat_cmct_rec_mod.f08
Defines file records for Jack's new derived ICESat/GLAS files for CMCT.
This is the record format of the flat-binary files, but is also returned by
the netCDF reader icesat_cmctnc_class_mod.f08.

*** icesat_cmct_config.json
- Configuration file for the derived GLAS files for CMCT.  Lists all the
  specific files for each campaign.
- These are the flat-binary versions of these files.
- Greenland only.
- Not currently in use: Superceded by the netCDF version
  icesat_cmctnc_grn_config.json.


========================================================================

** icesat_cmctnc_class_mod.f09

Class for reading the netCDF-format ICESat GLA12-derived files that Jack
created for the CMCT project.  Uses the netcdf nf90 API.

*** icesat_cmctnc_class
Base type for the class.

*** init
Type-bound method.  Initializes the object, gets the list of data files
from the json configuration file, and opens the first data file.

*** next_rec
Type-bound method.  Reads the next value from each of the netCDF variables
and returns them as an icesat_cmct_rec.

*** next_file
Type-bound method.  Opens the next file, stores all the netcdf varids, and
gets the invalid values (_FillValue) for each variable.

**** TODO Find correct chunk cache parameters to read the file efficiently.
- Currently it takes about 20x longer to read these netCDF-4/HDF5 files
  than reading the same campaign from the icesat_cism dataset.
- Alternative: ~/CMCT-datasetsFF has a version of this module that reads
  the entire netCDF file into memory.  Seems to work well (and fast!) for
  the Greenland files, but the Antarctica files are ~ 10x larger.

*** clean
Type-bound method. Checks if the given record meets various criteria,
returns a code indicating if it's good or why it failed.  Need to specify
which value to check, currently only the elevations are supported.
- Value does not have its invalid value
- Elevation over Greenland < 3300 m
- Elevation over Antarctica < 4500 m
- Elevation is within 200 m of GIMP DEM

*** destroy
Type-bound method.  Close the open file, free the allocated arrays, and do
any other cleanup.

*** check_ncdf
Private subroutine for checking the status of the netcdf library calls.

*** icesat_cmctnc_grn_config.json
- Configuration file for the derived GLAS files for CMCT.  Lists all the
  specific files for each campaign.
- These are the netCDF versions of these files.
- Greenland only.

========================================================================

** icesat_config_mod.f08

Module to read the configuration files for the ICESat datasets, which
are organized by campaign.

*** subr parse_icesat_config
Parses the JSON configuration file (eg. icesat_cism_config.json) whose name
is passed in as a dummy argument.  Calls fson.
**** TODO Verify campaign from values in config file.

========================================================================


** kinds_mod.f90
- Symbolic kind names, from GLAS.
- Used by polar_stereographic_mod.f90, cism_polarstereo_mod.f03 (for
  compatibility with polar_stereographic_mod, although that's a bit
  misleading because variables have to be passed back to other procedures).

========================================================================


** moments2d_class_mod.f08
Calculates incremental sample means and variances on a 2D grid.

*** type moments2d_class
*** subr new
*** subr add_point
*** subr get
*** TODO Comments not quite consistent with code. Terminology a bit off too.
*** subr destroy
Type-bound method. Deallocate arrays.

*** CANCELED Add overall statistics option
Use new instance with 1x1 "grid" instead.

*** DONE Use different invalid value
- Maybe -huge(1.0) instead of -huge(1.0d0) so can convert output to real.
- Or, take value to use as input.
- NaN would be best once we have gfortran that supports it. (Actually 4.6.3
  does for I/O and (nan==nan) == .false..)
-- NCO User's Guide (sec 4.2 footnote 1) strongly discourages using NaN as
   netcdf _FillValue, not sure why.
-- Sophie generally prefers netCDF default invalid [CMCT weekly 2016-03-11]
- 2016-05-11: Now allows user to specify invalid, or uses equivalent of
  netCDF NF90_FILL_DOUBLE by default.

========================================================================


** polar_stereographic_mod.f90
Jack's polar stereographic conversions module.

========================================================================


** temp_file_mod.f08
Returns name of a temporary file.  Not used currently.

*** TODO Fix bug.


========================================================================

** MISC NOTES

*** gfortran versions:
**** Atlas: 4.6.3
**** Icesat11: 4.8.4
**** ggsghpcc: 4.9.2
**** gs6141icesat2-dev1: 4.8.2
**** IEEE Arithmetic: I/O yes; ieee modules starting with 5.0
- https://gcc.gnu.org/wiki/Fortran2003Status also says "Input and output of
  IEEE exceptional values: Yes" without giving version.  But without the
  module how would you tell??
- Tested with ~/testnan.f08: I/O works (NaN, Inf, -Inf), as does (Nan ==
  Nan) == .false.  Can generate with (zero/zero), (one/zero), where
  zero=0.0, one=1.0; but compiler complains if divide by constant 0.0
- f2003 allows boz constants as arguments to real(), dble() [MR&C 16.9]
**** json-fortran: needs 4.9
