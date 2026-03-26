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

Prepare forcings from ERA5 reanalysis:

This preprocessing:

- filters required atmospheric variables  
- interpolates fields to the Veros grid  
- produces forcing files compatible with the setup  

---

#### (b) GLORYS Initialization

Prepare initial conditions from GLORYS reanalysis:

This preprocessing:

- extracts temperature, salinity, and velocity fields  
- interpolates them onto the Veros grid  
- creates initialization files ready for simulation  

---

### 3. Run the Setup

After installation, the setup becomes available to Veros:

---

## License

This project is licensed under the MIT License — see the `LICENSE` file for details.
