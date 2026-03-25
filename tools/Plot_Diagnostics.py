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
import cmocean as cmo
import imageio
plt.style.use('science')

class Diagnostics:
    def __init__(self,run_dir,start_date):

        self.run_dir = Path(run_dir)
        self.folder_path = self.run_dir / "plots"
        self.folder_path.mkdir(exist_ok=True)

        self.ds = nc.Dataset(self.run_dir / "global_1deg.averages.nc")
        self.ds_energy = nc.Dataset(self.run_dir / "global_1deg.energy.nc")
        self.ds_snap = nc.Dataset(self.run_dir / "global_1deg.snapshot.nc")

        self.xt = self.ds.variables["xt"][:]
        self.xt, self.sort_idx = self.cyclic(self.xt)
        self.xu = self.ds.variables['xu'][:]
        self.xu, self.sort_idx_u = self.cyclic(self.xu)
        self.yt = self.ds.variables["yt"][:]
        self.yu = self.ds.variables['yu'][:]
        self.zt = self.ds.variables["zt"][:]
        self.zw = self.ds.variables["zw"][:]

        self.start_date = start_date
        units = f"days since {start_date:%Y-%m-%d}"

        time_var = self.ds.variables["Time"]
        self.dates = cftime.num2date(
            time_var[:],
            units=units,
            calendar="360_day"
        )
        
        time_var_nrj = self.ds_energy["Time"] 
        self.nrj_dates = cftime.num2date(
            time_var_nrj[:],
            units=units,
            calendar="360_day"
        )        

        time_var_snap = self.ds_snap["Time"] 
        self.snap_dates = cftime.num2date(
            time_var_snap[:],
            units=units,
            calendar="360_day"
        )        
        self.epsilon = 1e-6

    def get_time_indices(self, n_years=10, snapshot=False):
        dates = self.snap_dates if snapshot else self.dates
        last_date = dates[-1]
        if last_date.year == self.start_date.year:
            return np.arange(len(dates))
        else:
            if last_date.month < 12 or last_date.day < 30:
                last_year = last_date.year - 1
            else:
                last_year = last_date.year
        
            first_year = last_year - n_years + 1
        
            idx = [i for i, d in enumerate(dates) if first_year <= d.year <= last_year]
        
            return np.array(idx)

    def cyclic(self,x):
        x = ((x + 180) % 360) -180
        sort_idx = np.argsort(x)
        x = x[sort_idx]
        return x,sort_idx

    def plot_map(self, lon, lat, field, title, filename,
                 cmap="coolwarm", levels=150,
                 u=None, v=None):

        vmin = np.nanmin(field)
        vmax = np.nanmax(field)+self.epsilon

        plt.figure(figsize=(15, 8))
        ax = plt.axes(projection=ccrs.PlateCarree())

        cf = plt.contourf(
            lon, lat, field,
            np.linspace(vmin, vmax, levels),
            cmap=cmap,
            transform=ccrs.PlateCarree(),
            extend="both"
        )
        cf.set_edgecolor('face')

        plt.colorbar(cf, fraction=0.046, pad=0.04)

        if u is not None and v is not None:
            ax.streamplot(
                lon, lat, u, v,
                color="w",
                linewidth=0.3,
                density=1.5,
                transform=ccrs.PlateCarree()
            )

        ax.add_feature(cfeature.GSHHSFeature("low", levels=[1, 2, 6]))
        ax.coastlines()

        plt.title(title, fontsize=18)
        plt.savefig(self.folder_path / filename,
                    dpi=200, bbox_inches="tight")
        plt.close()

    def plot_map_gif(self, lon, lat, field, title, filename,
                     cmap="coolwarm", levels=150,
                     u=None, v=None, fps=5):
    
        nt = field.shape[0]
    
        vmin = np.nanmin(field)
        vmax = np.nanmax(field) + self.epsilon
    
        frames = []
    
        for t in range(nt):
    
            fig = plt.figure(figsize=(15, 8))
            ax = plt.axes(projection=ccrs.PlateCarree())
    
            cf = plt.contourf(
                lon, lat, field[t],
                np.linspace(vmin, vmax, levels),
                cmap=cmap,
                transform=ccrs.PlateCarree(),
                extend="both"
            )
            cf.set_edgecolor("face")
    
            plt.colorbar(cf, fraction=0.046, pad=0.04)
    
            if u is not None and v is not None:
                ax.streamplot(
                    lon, lat, u[t], v[t],
                    color="w",
                    linewidth=0.3,
                    density=1.5,
                    transform=ccrs.PlateCarree()
                )
    
            ax.add_feature(cfeature.GSHHSFeature("low", levels=[1, 2, 6]))
            ax.coastlines()
    
            plt.title(f"{title} (t={t})", fontsize=18)
    
            frame_file = self.folder_path / f"_frame_{t:04d}.png"
            plt.savefig(frame_file, dpi=200, bbox_inches="tight")
            plt.close()
    
            frames.append(imageio.v2.imread(frame_file))
    
        gif_path = self.folder_path / filename
        imageio.mimsave(gif_path, frames, fps=fps)
    
        for t in range(nt):
            (self.folder_path / f"_frame_{t:04d}.png").unlink()
    
        print(f"GIF saved to {gif_path}")
    

    def vertical_plot_map(self, lat, depth, field, title, filename,
                          cmap="coolwarm", levels=150):
    
        vmin = np.nanmin(field)
        vmax = np.nanmax(field) + self.epsilon
    
        plt.figure(figsize=(10, 6))
    
        cf = plt.contourf(
                lat, -depth, field,
            np.linspace(vmin, vmax, levels),
            cmap=cmap,
            extend="both"
        )
    
        cf.set_edgecolor("face")
    
        plt.colorbar(cf, fraction=0.046, pad=0.04)
    
        plt.gca().invert_yaxis()  # depth increasing downward
    
        plt.xlabel("Latitude")
        plt.ylabel("Depth (m)")
        plt.title(title, fontsize=18)
    
        plt.savefig(
            self.folder_path / filename,
            dpi=200,
            bbox_inches="tight"
        )
    
        plt.close()

    def vertical_plot_map_gif(self, lat, depth, field, title, filename,
                          cmap="coolwarm", levels=150, fps=5):

        nt = field.shape[0]

        vmin = np.nanmin(field)
        vmax = np.nanmax(field) /2

        frames = []

        for t in range(nt):

            plt.figure(figsize=(10, 6))

            cf = plt.contourf(
                lat, -depth, field[t],
                np.linspace(vmin, vmax, levels),
                cmap=cmap,
                extend="both"
            )

            cf.set_edgecolor("face")

            plt.colorbar(cf, fraction=0.046, pad=0.04)

            plt.gca().invert_yaxis()

            plt.xlabel("Longitude")
            plt.ylabel("Depth (m)")
            plt.title(f"{title} (t={t})", fontsize=18)

            frame_file = self.folder_path / f"_frame_{t:04d}.png"
            plt.savefig(frame_file, dpi=200, bbox_inches="tight")
            plt.close()

            frames.append(imageio.v2.imread(frame_file))

        gif_path = self.folder_path / filename
        imageio.mimsave(gif_path, frames, fps=fps)

        for t in range(nt):
            (self.folder_path / f"_frame_{t:04d}.png").unlink()

        print(f"GIF saved to {gif_path}")

    def temp(self, n_years=10, vert=False, snapshot=False, Gif=False, start_gif=None, end_gif=None):
        cmap = cmo.cm.thermal
        if snapshot:
            idx = self.get_time_indices(n_years,snapshot=True)
            year = self.snap_dates[idx[-1]].year
            ds = self.ds_snap
        else:
            idx = self.get_time_indices(n_years)
            year = self.dates[idx[-1]].year
            ds = self.ds

        if Gif:
            idx_gif = idx
            if start_gif is not None or end_gif is not None:
                idx_gif = idx[start_gif:end_gif]
    
            if vert:
                temp = np.mean(
                    ds.variables["temp"][idx_gif, :, :, :],
                    axis=3
                )  
    
                self.vertical_plot_map_gif(
                    self.yt, self.zt, temp,
                    f"Temperature evolution",
                    f"Temp_{year}_vert.gif",
                    cmap=cmap
                )
    
            else:
                temp = ds.variables["temp"][idx_gif, -1, :, :]
                temp = temp[:, :, self.sort_idx]
    
                self.plot_map_gif(
                    self.xt, self.yt, temp,
                    f"SST evolution",
                    f"SST_{year}.gif",
                    cmap=cmap
                )
    
        else:
    
            if vert:
                temp = np.mean(
                    ds.variables["temp"][idx, :, :, :],
                    axis=(0, 3)
                )
    
                self.vertical_plot_map(
                    self.yt, self.zt, temp,
                    f"Temperature (last {n_years} years)",
                    f"Temp_{year}_vert.pdf",
                    cmap=cmap
                )
    
            else:
                temp = np.mean(
                    ds.variables["temp"][idx, -1, :, :],
                    axis=0
                )
    
                temp = temp[:, self.sort_idx]
    
                self.plot_map(
                    self.xt, self.yt, temp,
                    f"SST (last {n_years} years)",
                    f"SST_{year}.pdf",
                    cmap=cmap
                )
    
        if vert:
            print("Vertical temperature done")
        else:
            print("SST done")


    def salt(self, n_years=10, vert=False, snapshot=False, Gif=False, start_gif=None, end_gif=None):
        cmap = cmo.cm.haline
        idx = self.get_time_indices(n_years,snapshot=False)
        if snapshot:
            year = self.snap_dates[idx[-1]].year
            ds = self.ds_snap
        else:
            year = self.dates[idx[-1]].year
            ds = self.ds

        if Gif:
            idx_gif = idx
            if start_gif is not None or end_gif is not None:
                idx_gif = idx[start_gif:end_gif]
    
            if vert:
                salt = np.mean(
                    ds.variables["salt"][idx_gif, :, :, :],
                    axis=3
                )  
    
                self.vertical_plot_map_gif(
                    self.yt, self.zt, salt,
                    f"Salinity evolution",
                    f"Salt_{year}_vert.gif",
                    cmap=cmap
                )
    
            else:
                salt = ds.variables["salt"][idx_gif, -1, :, :]
                salt = salt[:, :, self.sort_idx] 
    
                self.plot_map_gif(
                    self.xt, self.yt, salt,
                    f"SSS evolution",
                    f"SSS_{year}.gif",
                    cmap=cmap
                )
    
        else:
    
            if vert:
                salt = np.mean(
                    ds.variables["salt"][idx, :, :, :],
                    axis=(0, 3)
                )
    
                self.vertical_plot_map(
                    self.yt, self.zt, salt,
                    f"Salinity (last {n_years} years)",
                    f"Salt_{year}_vert.pdf",
                    cmap=cmap
                )
    
            else:
                salt = np.mean(
                    ds.variables["salt"][idx, -1, :, :],
                    axis=0
                )
    
                salt = salt[:, self.sort_idx]
    
                self.plot_map(
                    self.xt, self.yt, salt,
                    f"SSS (last {n_years} years)",
                    f"SSS_{year}.pdf",
                    cmap=cmap
                )
    
        if vert:
            print("Vertical salinity done")
        else:
            print("SSS done")

    def ssh(self, n_years=10,snapshot=False, Gif=False, start_gif=None, end_gif=None):
        idx = self.get_time_indices(n_years,snapshot=snapshot)
        if snapshot:
            ds = self.ds_snap
            year = self.snap_dates[idx[-1]].year
        else:
            ds = self.ds
            year = self.dates[idx[-1]].year
        if Gif:
            if start_gif is not None or end_gif is not None:
                idx_gif = idx[start_gif:end_gif]
            ssh = ds.variables["ssh"][idx_gif, :, :]
            ssh = ssh[:,:,self.sort_idx]
            ssh = ssh - np.nanmean(ssh)
            self.plot_map_gif(
                self.xt, self.yt, ssh,
                f"SSH evolution)",
                f"SSH_{year}.gif"
                )

        else:
            ssh = np.mean(
                ds.variables["ssh"][idx, :, :],
                axis=0
            )
            ssh = ssh[:, self.sort_idx]
            ssh = ssh - np.nanmean(ssh)
            self.plot_map(
                self.xt, self.yt, ssh,
                f"SSH (last {n_years} years)",
                f"SSH_{year}.pdf"
            )
        print("SSH done")

    def heat_flux(self, n_years=10):

        idx = self.get_time_indices(n_years)

        qnet = np.mean(
            self.ds.variables["qnet"][idx,0, :, :],
            axis=0
        )

        qsol = np.mean(
            self.ds.variables["qsol"][idx,0, :, :],
            axis=0
        )

        Q = (qnet + qsol)[:, self.sort_idx]
        year = self.dates[idx[-1]].year

        self.plot_map(
            self.xt, self.yt, Q,
            f"Heat flux (last {n_years} years)",
            f"HeatFlux_{year}.pdf"
        )

        print("Heat flux done")

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

    def seasonal_mld(self, n_years=10):
    
        idx = self.get_time_indices(n_years=n_years)
    
        zt = self.zt
        lat = self.yt
    
        idx_winter = [i for i in idx if self.dates[i].month in [1, 2, 3]]
        idx_summer = [i for i in idx if self.dates[i].month in [7, 8, 9]]
    
        mld_winter_all = []
        for it in idx_winter:
            mld_winter_all.append(
                self.compute_mld(
                    self.ds.variables['temp'][it, :, :, :][:, :, self.sort_idx],
                    self.ds.variables['salt'][it, :, :, :][:, :, self.sort_idx],
                    zt
                )
            )
    
        mld_summer_all = []
        for it in idx_summer:
            mld_summer_all.append(
                self.compute_mld(
                    self.ds.variables['temp'][it, :, :, :][:, :, self.sort_idx],
                    self.ds.variables['salt'][it, :, :, :][:, :, self.sort_idx],
                    zt
                )
            )
    
        mld_march = np.mean(mld_winter_all, axis=0)
        mld_sept = np.mean(mld_summer_all, axis=0)
    
        north = lat > 0
        north2d = north[:, None]
    
        # Hemisphere logic
        mld_winter = np.where(north2d, mld_march, mld_sept)
        mld_summer = np.where(north2d, mld_sept, mld_march)
    
        return mld_winter, mld_summer

    def mld(self, n_years=10):
    
        mld_winter, mld_summer = self.seasonal_mld(n_years)
    
        maxvals = [500,100] #[np.nanmax(mld_winter), np.nanmax(mld_summer)]
        scales = ['linear','linear']
        fig, axs = plt.subplots(
            1, 2,
            figsize=(18, 6),
            subplot_kw={'projection': ccrs.PlateCarree()}
        )
    
        titles = [
            f'Winter MLD (last {n_years} years)',
            f'Summer MLD (last {n_years} years)'
        ]
    
        fields = [mld_winter, mld_summer]
    
        for ax, field, title, vmax, scale in zip(axs, fields, titles, maxvals, scales):
            if scale == 'log':
                norm = mcolors.LogNorm(vmin=1, vmax=vmax)
                levels = np.logspace(np.log10(1), np.log10(vmax), num=200)

            else:
                levels = np.arange(0, vmax, vmax/200)
                norm = mcolors.Normalize(vmin=0, vmax=vmax)

            cf = ax.contourf(
                self.xt, self.yt, field,
                levels=levels,
                cmap='viridis',
                norm=norm,
                extend='max',
                transform=ccrs.PlateCarree()
            )

            cf.set_edgecolor('face')
            ax.add_feature(cfeature.GSHHSFeature('low', levels=[1, 2, 6]))
            ax.coastlines()
            ax.set_title(title, fontsize=16)
        
            cbar = fig.colorbar(
                cf,
                ax=ax,
                orientation='horizontal',
                fraction=0.06,
                pad=0.08
            )
            if scale =='log':
                cbar.locator = LogLocator(base=10.0,subs=[1,2,5], numticks=10)
                cbar.formatter = LogFormatter(base=10, labelOnlyBase=False) 
                cbar.update_ticks()
            cbar.set_label('Mixed Layer Depth (m)', fontsize=12)
        
        plot_name = f"MLD_winter_summer_{self.dates[-1].year}.pdf"
        plt.savefig(self.folder_path / plot_name,
                    dpi=200, bbox_inches='tight')
        plt.close()
    
        print('MLD plot done')

    def velocity(self, n_years=10, vert=False,snapshot=False, Gif=False, start_gif=None, end_gif=None):
        idx = self.get_time_indices(n_years,snapshot)
        cmap = cmo.cm.speed
        if snapshot:
            year = self.snap_dates[idx[-1]].year
            ds = self.ds_snap
        else:
            year = self.dates[idx[-1]].year
            ds = self.ds

        if Gif:
            if start_gif is not None or end_gif is not None:
                idx_gif = idx[start_gif:end_gif]
            if vert:
                u = ds.variables["u"][idx_gif, :, 80, :]
                v = ds.variables["v"][idx_gif, :, 80, :]
                speed = np.sqrt(u**2 + v**2)
                self.vertical_plot_map_gif(
                        self.xt, self.zt, speed,
                        f"Velocity evolution)",
                        f"Velocity_{year}_vert_lat80.gif",
                        cmap=cmap
                    )

            else:
                u = ds.variables["u"][idx_gif, -1, :, :]
                v = ds.variables["v"][idx_gif, -1, :, :]
                u = u[:,:,self.sort_idx]
                v = v[:,:,self.sort_idx]
                speed = np.sqrt(u**2 + v**2)
                self.plot_map_gif(
                        self.xt, self.yt, speed,
                        f"Surface velocity evolution)",
                        f"Velocity_{year}.gif",
                        cmap=cmap
                    )
        else:    
            if vert:
                u = np.mean(
                    ds.variables["u"][idx, :, :, :],
                    axis=(0,3)
                )
                v = np.mean(
                    ds.variables["v"][idx, :, :, :],
                    axis=(0,3)
                )
                speed = np.sqrt(u**2 + v**2)
                self.vertical_plot_map(
                    self.yt, self.zt, speed,
                    f"Velocity (last {n_years} years)",
                    f"Velocity_{year}_vert.pdf",
                    cmap=cmap
                )
            else:
                u = np.mean(
                    ds.variables["u"][idx, -1, :, :],
                    axis=0
                )
                v = np.mean(
                    ds.variables["v"][idx, -1, :, :],
                    axis=0
                )
                u = u[:, self.sort_idx]
                v = v[:, self.sort_idx]
                speed = np.sqrt(u**2 + v**2)
                self.plot_map(
                    self.xt, self.yt, speed,
                    f"Surface velocity (last {n_years} years)",
                    f"Velocity_{year}.pdf",
                    cmap=cmap
                )
        print("Velocity done")

    def energy(self):
        eke_m = self.ds_energy['k_m'][:]
        x = [datetime(d.year, d.month, d.day, d.hour) for d in self.nrj_dates[:]]
        fig, ax = plt.subplots(figsize=(9, 4))
         
        ax.plot(x, eke_m,color='black', lw=0.5)
        ax.set_xlabel("Time")
        ax.set_ylabel(r"$k_m$")
        ax.set_title(r"$k_m$ convergence")
        
        locator = mdates.AutoDateLocator(minticks=5, maxticks=8)
        formatter = mdates.ConciseDateFormatter(locator)
        
        ax.xaxis.set_major_locator(locator)
        ax.xaxis.set_major_formatter(formatter)
        ax.xaxis.set_minor_locator(ticker.AutoMinorLocator(6))
        
        ax.grid(True,alpha=0.3)
        fig.autofmt_xdate()
        plt.tight_layout()
        plt.savefig(self.folder_path/'Pre_Energy_plot_ssh.pdf')
        plt.close()
        print("Energy done")
        return x, eke_m
    

    def run_all(self, n_years_default=10):
        self.temp(n_years_default,vert=True)
        self.salt(n_years_default,vert=True) 
        #self.ssh(n_years_default)
        #self.velocity(n_years=1,vert=True)
        #self.velocity(n_years_default,vert=False)
        #self.energy()
        #self.heat_flux(n_years_default)
        #self.mld(n_years=1)

start = datetime(1986,1,1)
run_dir_bat = "/Odyssey/private/e25cheve/simu_veros/runs/test_bathy/"
run_dir_smooth = "/Odyssey/private/e25cheve/simu_veros/runs/test_bathy/output_bathy_both/"
run_dir_TS = "/Odyssey/private/e25cheve/simu_veros/runs/TS_glorys/"
run_dir_TSUV = "/Odyssey/private/e25cheve/simu_veros/runs/TSUV_glorys/"
run_dir_TSUVSSH = "/Odyssey/private/e25cheve/simu_veros/runs/TSUVSSH_glorys/"
D_bat = Diagnostics(run_dir_bat,start)
#D_TS = Diagnostics(run_dir_TS,start)
D_bat.run_all()
#x_bat, nrj_bat = D_bat.energy()
#x, nrj_TS = D_TS.energy()
#D_TSUV = Diagnostics(run_dir_TSUV,start)
#x, nrj_TSUV = D_TSUV.energy()
#D_TSUVSSH = Diagnostics(run_dir_TSUVSSH,start)
#x, nrj_TSUVSSH = D_TSUVSSH.energy()
#D_smooth = Diagnostics(run_dir_smooth,start)
#x_smooth, nrj_smooth = D_smooth.energy()
#
#fig, ax = plt.subplots(figsize=(9, 4))
#ax.plot(x_bat, nrj_TS,color='black', lw=0.5,label="TS")
#ax.plot(x, nrj_TSUV,color='red', lw=0.5,label="TSUV")
#ax.plot(x, nrj_TSUVSSH,color='green', lw=0.5,label="TSUVSSH")
#ax.plot(x_bat, nrj_bat,color='blue', lw=0.5,label="Smooth + friction")
#ax.plot(x_bat, nrj_smooth,color='purple', lw=0.5,label="Smooth bathy")
#
#ax.set_xlabel("Time")
#ax.set_ylabel(r"$k_m$")
#ax.set_title(r"$k_m$ convergence")
#
#locator = mdates.AutoDateLocator(minticks=5, maxticks=8)
#formatter = mdates.ConciseDateFormatter(locator)
#
#ax.xaxis.set_major_locator(locator)
#ax.xaxis.set_major_formatter(formatter)
#ax.xaxis.set_minor_locator(ticker.AutoMinorLocator(6))
#
#ax.grid(True,alpha=0.3)
#fig.autofmt_xdate()
#plt.tight_layout()
#plt.legend()
#plt.savefig('/Odyssey/private/e25cheve/transfert/'+'Multi_Energy_plot_ssh.pdf')
#plt.close()



