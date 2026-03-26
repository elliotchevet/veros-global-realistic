# veros_global_realistic

## Overview

**veros_global_realistic** is a plugin for Veros that provides a realistic global ocean configuration.

It includes:

- 🌍 A global realistic model setup  
- 🌬️ Hourly ERA5 atmospheric forcings  
- 🌊 GLORYS-based ocean initialization  
- 🧰 Preprocessing tools to prepare forcing and initialization datasets  

The plugin is designed to simplify realistic ocean simulations by providing a consistent workflow for data preparation and model configuration.

---

## Quick Usage

### 1. Installation

Clone the repository and install the plugin in editable mode:

```bash
git clone https://github.com/yourusername/veros_global_realistic.git
cd veros_global_realistic
pip install -e .
```

Editable installation (`-e`) is recommended so that modifications to the setup or preprocessing tools are immediately available.

---

### 2. Prepare Forcing and Initialization Data

The plugin provides preprocessing utilities to:

- filter raw datasets  
- interpolate data onto the Veros grid  
- generate files directly usable by the setup  

Typical workflow:

#### (a) ERA5 Forcings

Download ERA5 data (hourly resolution) using your preferred method (e.g. CDS API), then run:

```bash
python -m veros_global_realistic.preprocessing.era5 \
    --input /path/to/raw_era5 \
    --grid /path/to/veros_grid.nc \
    --output /path/to/processed_forcing.nc
```

This step:

- filters required atmospheric variables  
- interpolates fields to the Veros grid  
- produces forcing files compatible with the setup  

---

#### (b) GLORYS Initialization

Prepare initial conditions from GLORYS reanalysis:

```bash
python -m veros_global_realistic.preprocessing.glorys \
    --input /path/to/glorys_data \
    --grid /path/to/veros_grid.nc \
    --output /path/to/initial_conditions.nc
```

This preprocessing:

- extracts temperature, salinity, and velocity fields  
- interpolates them onto the Veros grid  
- creates initialization files ready for simulation  

---

### 3. Run the Setup

After installation, the setup becomes available to Veros:

```bash
veros run global_realistic
```

(or run through your standard Veros workflow)

---

## Project Structure

```text
veros_global_realistic/
├── veros_global_realistic/
│   ├── setup/
│   ├── preprocessing/
│   └── ...
├── README.md
├── LICENSE
└── setup.py / pyproject.toml
```

---

## Requirements

- Python ≥ 3.8  
- Veros ≥ 1.4  
- ERA5 datasets  
- GLORYS ocean reanalysis data  

---

## License

This project is licensed under the MIT License — see the `LICENSE` file for details.
