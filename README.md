# Description
The Cryosphere model Comparison tool (CmCt) compares ice sheet models against remote sensing observations. Currently, the tool supports a comparison of (1) modeled mass change against reconciled mass change from the Ice sheet Mass Balance Inter-comparison Exercise (IMBIE), (2) modeled mass change against observations from satellite gravimetry, and (3) modeled surface elevation against observations from satellite radar and laser altimetry. In the future, the CmCt tool will also be expanded to perform comparisons of other variables such as surface elevation <ins>change</ins> using laser altimetry observations and surface velocity using various observations.

To use the CmCt, you will need the code in this repository and the associated datasets (described below). The mass change comparisons use Jupyter notebooks to configure and run the comparisons, with supporting functions in Python scripts. The surface elevation comparison uses Jupyter notebooks to configure and run the comparisons, with supporting functions in Fortran code that needs to be compiled. **If you do not need the surface elevation comparison, you can skip the Fortran compilation steps below.**

## Getting started (quickly)

### Run the CmCt locally
1. Clone this repository: `git clone git@github.com:NASA-Cryospheric-Sciences-Laboratory/CmCt.git`
2. Create a virtual environment (NOTE that this requires conda): `conda create --name CmCt --file requirements.txt`
3. Activate the environment: `conda activate CmCt`
4. Run the `test/create_lithk_netcdfs.ipynb` notebook to create test model data for both ice sheets.
5. Run the `notebooks/IMBIE/imbie_comparison.ipynb` notebook, which will compare the test input data against IMBIE observations. This can be run for both ice sheets by changing the options in the notebook.

### Launch the CmCt directly on the CryoCloud JupyterHub
If you have an account on [CryoCloud](https://cryointhecloud.com/), you can launch CmCt directly on their JupyterHub using this [link](https://hub.cryointhecloud.com/hub/user-redirect/git-pull?repo=https%3A%2F%2Fgithub.com%2FNASA-Cryospheric-Sciences-Laboratory%2FCmCt&urlpath=lab%2Ftree%2FCmCt%2Fnotebooks%2FIMBIE%2Fimbie_comparison.ipynb&branch=main).


## Going further
TBD

### Ice sheet model data upload requirements
The CmCt requires ice sheet model files to follow the ISMIP6 conventions for model output variables. Currently, the models MUST be defined on a rectangular X-Y grid in the ISMIP6 standard projected polar-stereographic space. Details on this can be found on the [ISMIP6 Wiki](https://theghub.org/groups/ismip6/wiki).


## Contributing

Thanks for your interest in contributing! There are many ways to contribute to this project. Get started [here](CONTRIBUTING.md).

---

# Details on running the surface elevation comparison
The surface elevation comparison calculates differences between modeled surface elevations and observations from satellite laser (ICESat (Ice, Cloud, and land Elevation Satellite)) and radar (ERS-1, ERS-2 (European Remote-Sensing Satellite), and Envisat (Environmental Satellite)) altimetry.

To install the surface elevation comparison code:
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

To download surface elevation datasets:
1. Install Globus Connect Personal: https://docs.globus.org/globus-connect-personal/install/linux/
1. Start Globus Connect Personal (without GUI): `./globusconnectpersonal -start &`
1. Browse to the web-based Globus File Manager: https://app.globus.org/file-manager
1. Create a datasets directory
1. Transfer data from the GHub-CmCt-Data endpoint to the local endpoint (datasets directory)


To configure and run the surface elevation comparison:
1. Create the directory for RUNS folder
1. Update RUNDIRS for RUNS folder in cmct_launch_config.ksh file
1. Update the tar variable path in cmct_launch_config.ksh file
   
