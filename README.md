# veros_global_realistic

## Overview

**veros_global_realistic** is a plugin for Veros that provides a realistic global ocean configuration.

It includes:

- A global realistic model setup  
- Wrappers for hourly ERA5 atmospheric forcings  
- GLORYS-based ocean initialization  
- Pre/Post-processing tools  

The plugin is designed to run realistic global ocean simulations through a simple workflow for data preprocessing and model configuration. By default, it provides a one-degree global setup with 40 vertical levels. The grid and land/water masks used by this configuration are stored in `data/masks.nc`.

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
```

### 2. Prepare Forcing and Initialization Data

#### (a) ERA5 Forcings

Prepare forcings from ERA5 reanalysis:

- Download ERA5 forcings from the CDS API using `tools/get_era5.py` for the desired range of years, default is 1993 only. See the [CDS API documentation](https://cds.climate.copernicus.eu/how-to-api) for authentication and setup details. One GRIB file per year is saved in `data/ERA5`.

- Install grib_to_netcdf converter and convert ERA5 grib files using command:

```bash
grib_to_netcdf -o file_name.nc file_name.grib
```

- Run `tools/preprocessing.py` with the following options:

```bash
python .tools/preprocessing.py --dataset ERA5 
```
Netcdf processed files are generated in `veros/veros_assets/global_1deg_realistic/` directory. They include interpolated hourly ERA5 fields on the Veros simulation grid for the forcing variables required. 

#### (b) GLORYS Initialization

Prepare initial conditions from GLORYS reanalysis:

- Download GLORYS12 monthly mean product from the Copernicus Marine Data Store using `tools/get_glorys12.py` for the chosen year and month, default is January 1993. See the [Copernicus Marine Toolbox API](https://help.marine.copernicus.eu/en/articles/8283072-copernicus-marine-toolbox-api-subset?pk_vid=ae1ac9aa63e0a5e91786376820793115) for authentication and setup details.

- Run `tools/preprocessing.py` with the following options:

```bash
python ./tools/preprocessing.py --dataset GLORYS12
```

NetCDF processed files are generated in the `veros/veros_assets/global_1deg_realistic/` directory. They include filtered and coarse-grained fields on the Veros grid, with masking taken from `data/masks.nc`.

Horizontal filtering uses [GCM-Filters](https://gcm-filters.readthedocs.io/), a diffusion-based spatial filtering package for gridded ocean, weather, and climate data described by [Loose et al. (2022)](https://doi.org/10.21105/joss.03947).

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
