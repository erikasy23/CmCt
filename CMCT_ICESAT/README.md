#CMCT 

##Build Requirements 

December 2020, JMS

Following software installed in experimental Docker container for CmCt ICESAT:
 
Basic dependencies:
- CentOS 7
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

##More information

For further information, source code notes, todos, history, and caveats, refer to contents of docs subdirectory.
The notes there were compiled in 2016 by CmCt author jguerber.

