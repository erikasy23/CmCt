# CMCT 

To begin with, let's experiment with building and running CMCT ICESAT.

## Build Requirements 


Following software installed in experimental Docker container for CmCt ICESAT:
 
Basic dependencies:
- CentOS 7
- gcc-gfortran 4.8.5 20150623 (Red Hat 4.8.5-44)
- ksh
- jq
- make
- git
- unzip
- wget

CentOS 7 RPMs:
- libaec and libaec-devel (1.0.4-1) -- provides libsz.so.2
- hdf5 and hdf5-devel (1.8.12-11)
- netcdf and netcdf-devel (4.3.3.1)

from Unidata's github repo, https://github.com/Unidata/netcdf-fortran:
- netcdf-fortran (4.3.3) -- provides libnetcdff.so and netcdf.mod

fson from https://github.com/josephalevin/fson
- Alternative json-fortran needs gcc/gfortran 4.9

## TODO

Various steps are desirable/needed to get CmCt containers running at CCR:

- Improve Makefile, perform make install
- Improve Dockerfile, create user, execute build locally
- Convert to Singularity container to run on cluster
- Simplify/convert existing "webby" run scripts for CCR environment

Likewise, useful for the code:

- Standardize input data sources
- Create tests

Note that esimon is separately testing CmCt on the cluster using a combination of modules
(netcdf-gfortran) and user space installs (fson).

## More information

For further information, source code notes, todos, history, and caveats, refer to contents of docs subdirectory.
The notes there were compiled in 2016 by CmCt author jguerber.

December 2020, JMS
