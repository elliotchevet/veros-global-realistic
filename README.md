# veros_global_realistic

## Overview

**veros_global_realistic** is a plugin for Veros that provides a realistic global ocean configuration.

It includes:

- A global realistic model setup  
- Wrappers for hourly ERA5 atmospheric forcings  
- Ocean state initialization from GLORYS12 monthly mean
- Restoring from GLORYS12 climatology
- Pre/Post-processing tools  

The plugin is designed to run realistic global ocean simulations through a simple workflow for data preprocessing and model configuration. By default, it provides a one-degree global setup with 40 vertical levels. The grid and land/water masks used by this configuration are stored in `data/masks.nc`.

This plugin builds on [Veros](https://github.com/team-ocean/veros), the versatile ocean simulator originally developed by Dion Häfner and collaborators. If you use this setup in scientific work, please cite the Veros model description paper: [Häfner et al. (2018)](https://doi.org/10.5194/gmd-11-3299-2018).

## Quick Usage

### 1. Installation

Clone the repository and install the plugin in editable mode:

```bash
git clone https://github.com/elliotchevet/veros_global_realistic.git
cd veros_global_realistic
pip install -e .
```
To use pre/post-processing tools, extra libraries are required. Install plugin with the following:

```bash
pip install -e ".[tools]"
conda install -c conda-forge esmpy
```

### 2. Prepare Forcing, Initialization and Restoring Data

#### (a) ERA5 Forcings

Prepare forcings from ERA5 reanalysis:

- Download ERA5 forcings from the CDS API using `tools/get_era5.py` for the desired range of years, default is 1993 only. See the [CDS API documentation](https://cds.climate.copernicus.eu/how-to-api) for authentication and setup Api key details. One GRIB file per year is saved in `data/ERA5`.

- Install [ecCodes](https://confluence.ecmwf.int/spaces/UDOC/pages/212436594/How%2Bto%2Binstall%2BecCodes%2Bwith%2BPython%2Bbindings%2Bin%2Bconda%2B-%2BecCodes%2BFAQ?utm_source=chatgpt.com) to get grib_to_netcdf converter and convert ERA5 grib files using command:

```bash
grib_to_netcdf -o file_name.nc file_name.grib
```

- Run `tools/preprocessing.py` with the following options:

```bash
python .tools/preprocessing.py --dataset ERA5 
```
Netcdf processed files are generated in `veros/veros_assets/global_1deg_realistic/` directory. They include interpolated hourly ERA5 fields on the Veros simulation grid for the forcing variables required. 
#### (b) GLORYS12 Initialization

Prepare initial conditions from GLORYS12 reanalysis:

- Download the GLORYS12 monthly mean product from the Copernicus Marine Data Store using `tools/get_glorys12.py` for the chosen year and month. The default is January 1993. See the [Copernicus Marine Toolbox API](https://help.marine.copernicus.eu/en/articles/7949409-copernicus-marine-toolbox-introduction#h_9172b5c79a) for authentication and API key setup details.

```bash
python ./tools/get_glorys12.py --mode IC --yyyymm 199301
```

- Run `tools/preprocessing.py` with the GLORYS12 initial-condition option:

```bash
python ./tools/preprocessing.py --dataset GLORYS12_IC
```

NetCDF processed files are generated in the `veros/veros_assets/global_1deg_realistic/` directory. They include filtered and coarse-grained fields on the Veros grid, with masking taken from `data/masks.nc`.

Horizontal filtering uses [GCM-Filters](https://gcm-filters.readthedocs.io/), a diffusion-based spatial filtering package for gridded ocean, weather, and climate data described by [Loose et al. (2022)](https://doi.org/10.21105/joss.03947).

#### (c) GLORYS12 Restoring

Prepare monthly surface restoring data from the GLORYS12 climatology:

- Download the GLORYS12 monthly climatology product from the Copernicus Marine Data Store using `tools/get_glorys12.py`. The restoring download uses the surface level only and writes the file under `data/GLORYS12/Restoring`.

```bash
python ./tools/get_glorys12.py --mode Restoring
```

- Run `tools/preprocessing.py` with the GLORYS12 restoring option:

```bash
python ./tools/preprocessing.py --dataset GLORYS12_Restoring
```

NetCDF processed files are generated in the `veros/veros_assets/global_1deg_realistic/` directory. They include filtered and coarse-grained monthly restoring fields on the Veros grid, with masking taken from `data/masks.nc`.


### 3. Run the Setup

To run the setup:

- Use the copy-setup command of Veros in the directory of your choice:

```bash
veros copy-setup global_1deg_realistic
```

- Finally:

```bash
veros run global_1deg_realistic.py
```

## License

This project is licensed under the MIT License — see the `LICENSE` file for details.
