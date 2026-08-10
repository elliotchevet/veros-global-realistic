# veros_global_realistic

## Overview

**veros_global_realistic** is a plugin for Veros that provides a realistic global ocean configuration.

It includes:

- A global realistic model setup  
- Wrappers for hourly ERA5 atmospheric forcings  
- GLORYS-based ocean initialization  
- Pre/Post-processing tools  

The plugin is designed to run realistic global ocean simulations by providing a simple workflow for data preprocessing and model configuration.

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

- Download ERA5 forcings using get_era5.py for a certain range of years.
A grib file by year.

- Install grib_to_netcdf converter and convert ERA5 grib files using command:

```bash
grib_to_netcdf -o file_name.nc file_name.grib
```

- Run preprocessing.py with the following options:

```bash
python preprocessing.py --single_level=True
```
Netcdf processed files are generated in veros/veros_assets/global_1deg_realistic/ directory. They include interpolated hourly ERA5 fields on the Veros simulation grid for the forcing variables required. 

#### (b) GLORYS Initialization

Prepare initial conditions from GLORYS reanalysis:



### 3. Run the Setup

After installation, the setup becomes available to Veros:

## License

This project is licensed under the MIT License — see the `LICENSE` file for details.
