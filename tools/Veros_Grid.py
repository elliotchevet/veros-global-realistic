import numpy as np
import netCDF4 as nc
import xarray as xr
import os
from netCDF4 import Dataset
import warnings

warnings.filterwarnings("ignore")

def Generate_grid_file(veros_in, output_path):
    ds_in = xr.open_dataset(veros_in, chunks={"time":1})
    new_ds = xr.Dataset(coords={
                            "xt": ds_in["xt"],
                            "xu": ds_in["xu"],
                            "yt": ds_in["yt"],
                            "yu": ds_in["yu"],
                            "zt": ds_in["zt"]})
    new_ds.to_netcdf(output_path,mode="w")
    print('Target grid file successfully created')
    return new_ds


veros_in = "/Odyssey/private/e25cheve/simu_veros/runs/global_1deg_glorys/global_1deg.averages.nc"
output_path = "/Odyssey/private/e25cheve/data/interp_grid_60_levels.nc"
Generate_grid_file(veros_in,output_path)

