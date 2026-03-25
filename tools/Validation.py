import matplotlib.pyplot as plt
import numpy as np
import netCDF4 as nc
from cartopy import config
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import os
from pathlib import Path
import sys
from matplotlib import ticker,_cm
from datetime import datetime, timedelta
import cftime
import matplotlib.dates as mdates
import matplotlib.colors as mcolors
from dateutil import tz
import gsw
from matplotlib.ticker import LogFormatter, LogLocator
import matplotlib.ticker as ticker
import scienceplots
import cmasher as cmr
plt.style.use('science')



class Validation:
    def __init__(self,file_path,start_date):

        self.ds = nc.Dataset(file_path)
        self.xt = self.get_coord(["xt", "lon", "longitude"])
        self.xt, self.sort_idx = self.cyclic(self.xt)
        self.yt = self.get_coord(["yt", "lat", "latitude"])
        self.zt = self.get_coord(["zt","time"])
        

        time_var = self.get_coord(["Time","time"])
        unit_var = self.get_units(["Time","time"])


        if unit_var == "days":   
            units = f"days since {start_date:%Y-%m-%d}"
            self.dates = cftime.num2date(
                time_var[:],
                units=units,
                calendar="360_day"
            )
        else:
            self.dates = cftime.num2date(
                time_var[:],
                units=unit_var,
                calendar="365_day"
            )

        
    def get_coord(self, possible_names):
        for name in possible_names:
            if name in self.ds.variables:
                return self.ds.variables[name][:]
        raise KeyError(f"None of these coordinates found: {possible_names}")

    def get_units(self, possible_names):
        for name in possible_names:
            if name in self.ds.variables:
                return self.ds.variables[name].units
        raise KeyError(f"None of these coordinates found: {possible_names}")

    def get_time_indices(self, n_years=10, use_energy=False):
        dates = self.dates
        last_date = dates[-1]
        if last_date.month < 12 or last_date.day < 30:
            last_year = last_date.year - 1
        else:
            last_year = last_date.year
    
        first_year =  last_year - n_years + 1
    
        idx = [i for i, d in enumerate(dates) if first_year <= d.year <= last_year]
    
        return np.array(idx)

    def cyclic(self,x):
        x = ((x + 180) % 360) -180
        sort_idx = np.argsort(x)
        x = x[sort_idx]
        return x,sort_idx


    def find_month(self,month):
        return np.array([i for i, d in enumerate(self.dates) if d.month == month])

    def compute_mld(self,temp,salt,zt,drho=0.03, zref=10.0):
        nz, ny, nx = temp.shape
        kref = np.argmin(np.abs(zt - zref))
        rho = gsw.rho(salt, temp, 0)
        rho_ref = rho[kref, :, :]
        mld = np.full((ny, nx), np.nan)

        for k in range(kref, - 1, -1):
            mask = (rho[k, :, :] - rho_ref) >= drho
            newly_found = mask & np.isnan(mld)
            mld[newly_found] = zt[k]
    
        return -mld

    def mld(self, n_years=10,compute_mld=True):
        idx = self.get_time_indices(n_years=n_years)
        if len(idx)==0:
            idx = [i for i in range(12)]
        zt = self.zt
        lat = self.yt
        idx_winter = [i for i in idx if self.dates[i].month in [1, 2, 3]]
        idx_summer = [i for i in idx if self.dates[i].month in [7, 8, 9]]
        mld_winter_all = []
        mld_summer_all = []
        if compute_mld: 
            if 'temp' in self.ds.variables:
                temp_var = 'temp'
            elif 'thetao' in self.ds.variables:
                temp_var = 'thetao'
            else:
                raise KeyError("No temperature variable found ('temp' or 'thetao')")
            
            if 'salt' in self.ds.variables:
                salt_var = 'salt'
            elif 'so' in self.ds.variables:
                salt_var = 'so'
            else:
                raise KeyError("No salinity variable found ('salt' or 'so')")
    

            for it in idx_winter:
                mld_winter_all.append(
                    self.compute_mld(
                        self.ds.variables[temp_var][it, :, :, :][:, :, self.sort_idx],
                        self.ds.variables[salt_var][it, :, :, :][:, :, self.sort_idx],
                        zt
                    )
                )
    

            for it in idx_summer:
                mld_summer_all.append(
                    self.compute_mld(
                        self.ds.variables[temp_var][it, :, :, :][:, :, self.sort_idx],
                        self.ds.variables[salt_var][it, :, :, :][:, :, self.sort_idx],
                        zt
                    )
                )
    
        else:
            mld_var = self.ds.variables['mlotst'][:]
            for it in idx_winter:
                mld_winter_all.append(mld_var[it,:,:][:,self.sort_idx])
            for it in idx_summer:
                mld_summer_all.append(mld_var[it,:,:][:,self.sort_idx])
        mld_march = np.mean(mld_winter_all, axis=0)
        mld_sept = np.mean(mld_summer_all, axis=0)
    
        north = lat > 0
        north2d = north[:, None]
    
        mld_winter = np.where(north2d, mld_march, mld_sept)
        mld_summer = np.where(north2d, mld_sept, mld_march)
    
        return mld_winter, mld_summer

def plot_mld_comparison(obj1, obj2, n_years=10, vmax_w=500, vmax_s=100):
    mld_w1, mld_s1 = obj1.mld(n_years=n_years)
    mld_w2, mld_s2 = obj2.mld(n_years=n_years,compute_mld=False)
    diff_w = mld_w2 - mld_w1
    diff_s = mld_s2 - mld_s1

    xt = obj1.xt
    yt = obj1.yt

    fig, axs = plt.subplots(
        2, 3,
        figsize=(22, 10),
        subplot_kw={'projection': ccrs.PlateCarree()}
    )

    winter_fields = [mld_w1, mld_w2, diff_w]
    winter_titles = [
        "Winter MLD - File 1",
        "Winter MLD - File 2",
        "Winter MLD Difference (2 - 1)"
    ]

    for j in range(3):
        ax = axs[0, j]

        if j < 2:
            levels = np.linspace(0, vmax_w, 200)
            norm = mcolors.Normalize(vmin=0, vmax=vmax_w)
            cmap = "viridis"
        else:
            vmax_diff = np.nanmax(np.abs(diff_w))
            levels = np.linspace(-vmax_diff, vmax_diff, 200)
            norm = mcolors.TwoSlopeNorm(vmin=-vmax_diff, vcenter=0, vmax=vmax_diff)
            cmap = "RdBu_r"

        cf = ax.contourf(
            xt, yt, winter_fields[j],
            levels=levels,
            cmap=cmap,
            norm=norm,
            transform=ccrs.PlateCarree(),
            extend="both"
        )

        ax.coastlines()
        ax.add_feature(cfeature.GSHHSFeature('low', levels=[1]))
        ax.set_title(winter_titles[j], fontsize=14)

        cbar = fig.colorbar(cf, ax=ax, orientation='horizontal',
                            fraction=0.05, pad=0.05)
        cbar.set_label("Mixed Layer Depth (m)")

    summer_fields = [mld_s1, mld_s2, diff_s]
    summer_titles = [
        "Summer MLD - File 1",
        "Summer MLD - File 2",
        "Summer MLD Difference (2 - 1)"
    ]

    for j in range(3):
        ax = axs[1, j]

        if j < 2:
            levels = np.linspace(0, vmax_s, 200)
            norm = mcolors.Normalize(vmin=0, vmax=vmax_s)
            cmap = "viridis"
        else:
            vmax_diff = np.nanmax(np.abs(diff_s))
            levels = np.linspace(-vmax_diff, vmax_diff, 200)
            norm = mcolors.TwoSlopeNorm(vmin=-vmax_diff, vcenter=0, vmax=vmax_diff)
            cmap = "RdBu_r"

        cf = ax.contourf(
            xt, yt, summer_fields[j],
            levels=levels,
            cmap=cmap,
            norm=norm,
            transform=ccrs.PlateCarree(),
            extend="both"
        )

        ax.coastlines()
        ax.add_feature(cfeature.GSHHSFeature('low', levels=[1]))
        ax.set_title(summer_titles[j], fontsize=14)

        cbar = fig.colorbar(cf, ax=ax, orientation='horizontal',
                            fraction=0.05, pad=0.05)
        cbar.set_label("Mixed Layer Depth (m)")

    plt.tight_layout()
    plt.savefig('/Odyssey/private/e25cheve/transfert/Validation_mld.pdf')

start_date = datetime(1986,1,1) 
veros_file = "/Odyssey/private/e25cheve/simu_veros/runs/global_1deg_glorys/global_1deg.averages.nc"
glorys_file = "/Odyssey/private/e25cheve/veros/veros_assets/global_1deg/Glorys_restoring.nc"

V1 = Validation(veros_file,start_date)
V2 = Validation(glorys_file,start_date)
plot_mld_comparison(V1,V2,n_years=1)
