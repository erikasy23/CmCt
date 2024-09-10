The Cryosphere model Comparison tool (CmCt) is a data analysis tool to facilitate ice sheet model comparison, validation, and analysis. Currently, the CmCt is used to compare ice sheet models provided by the user with remotely sensed laser (ICESat (Ice, Cloud, and land Elevation Satellite)), radar (ERS-1, ERS-2 (European Remote-Sensing Satellite), and Envisat (Environmental Satellite)) altimetry and GRACE (Gravity Recovery and Climate Experiment) gravimetry data. The CmCt works for both Greenland and Antarctica.

To install the CmCt:
1. `git clone` this repository
1. `conda create --name CmCt`
1. `conda activate CmCt`
1. `conda install -c conda-forge netcdf-fortran`
1. `mkdir externalpackages`
1. `cd externalpackages`
1. `git clone git@github.com:josephalevin/fson.git`
1. `cd fson`
1. `make`
1. `export LD_LIBRARY_PATH=${HOME}/CmCt/externalpackages/fson/dist:${HOME}/.conda_envs/CmCt/lib` ... NOTE that you may want to add this to your ~./bashrc
1. `mkdir externalpackages/jq`
1. `cd externalpackages/jq`
1. `wget https://github.com/jqlang/jq/releases/download/jq-1.7.1/jq-linux-amd64`
1. For each CMCT code directory (i.e., CMCT_ENVISAT), run `make` inside of that directory

To download CmCt datasets:
1. Install Globus Connect Personal: https://docs.globus.org/globus-connect-personal/install/linux/
1. Start Globus Connect Personal (without GUI): `./globusconnectpersonal -start &`
1. Browse to the web-based Globus File Manager: https://app.globus.org/file-manager
1. Create a datasets directory
1. Transfer data from the GHub-CmCt-Data endpoint to the local endpoint (datasets directory)



To configure and run CmCt 
1. Create the directory for RUNS folder
2. Update RUNDIRS for RUNS folder in cmct_launch_config.ksh file
3. Update the tar variable path in cmct_launch_config.ksh file
   