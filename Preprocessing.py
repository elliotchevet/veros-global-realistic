import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
import xarray as xr
import gcm_filters
import os
from netCDF4 import Dataset
import xesmf as xe
from datetime import datetime, timedelta
from dateutil import tz
import warnings
from dask.diagnostics import ProgressBar
from dask.distributed import Client
from glob import glob

warnings.filterwarnings("ignore")

class Preprocessor:
    def __init__(self, src_path, target_path,single_depth=False,verbose=False,ncores=16,filtering=True,weight_path='/Odyssey/private/e25cheve/data/Weights',weight_name='glorys12'):
        self.filtering = filtering
        self.weight_name = weight_name
        self.weight_path = weight_path
        self.ds = xr.open_dataset(src_path, chunks={"time": 1})
        rename_dict = {}
        if "longitude" in self.ds.coords:
            rename_dict["longitude"] = "lon"
        if "latitude" in self.ds.coords:
            rename_dict["latitude"] = "lat"
        if "depth" in self.ds.coords:
            rename_dict["depth"] = "zt"
        if rename_dict:
            self.ds = self.ds.rename(rename_dict)
        if self.ds.lon.min() < 0:
            if verbose: print("Adjusting longitude to 0–360…")
            self.ds = self.ds.assign_coords(lon=(self.ds.lon % 360))
            self.ds = self.ds.sortby("lon")
        if "zt" in self.ds.coords:
            if self.ds.zt.max() > 0:
                if verbose: print("Adjusting depth to negative values")
                self.ds = self.ds.assign_coords(zt=(-self.ds.zt))
                self.ds = self.ds.sortby("zt")
            nlat = self.ds.dims["lat"]
            nlon = self.ds.dims["lon"]
        if self.ds.dims["lat"] > 2000:
            self.ds = self.ds.isel(
                lat=slice(0, -100),
            )

        target_grid = xr.open_dataset(target_path, chunks={"time": 1})
        if "xt" in target_grid:
            if verbose: print("adjusting xt to lon")
            self.target_grid = target_grid.rename({"xt": "lon", "yt": "lat"})
        if "xu" in target_grid:
            if verbose: print("adjusting xu to lon")
            self.target_grid_u = target_grid.rename({"xu": "lon", "yt": "lat"})
        if "yu" in target_grid:
            if verbose: print("adjusting yu to lat")
            self.target_grid_v = target_grid.rename({"xt": "lon", "yu": "lat"})
        if "zt" in self.ds.coords:
            if single_depth:
                self.ds = self.ds.squeeze("zt",drop=True)
        self.z_target = self.target_grid.zt

        dims = self.ds.dims
        self.chunker = {"lon": -1, "lat": -1}
        nt = dims.get("time", None)
        nz = dims.get("zt", None)
        if nt is not None and nt > 1:
            if nt >= ncores:
                self.chunker["time"] = int(np.ceil(nt / ncores))
        
            else:
                self.chunker["time"] = 1
        
                if nz is not None and nz > 1:
                    depth_chunks = max(1, ncores // nt)
                    self.chunker["zt"] = int(np.ceil(nz / depth_chunks))
        elif nz is not None and nz > 1:
            self.chunker["zt"] = int(np.ceil(nz / ncores))

    def run(self, var):
        out = []
        da = self.ds[var]
        vp = VerticalProcessor(self.ds, z_name="zt")
        da = vp.interpolate(da, self.z_target)
        mask_ocean = np.isfinite(da.isel(time=0))
        grid = LatLonGrid(
            self.ds.lat.values,
            self.ds.lon.values,
            mask=mask_ocean
        )
        if self.filtering:
            da_filtered = []
            for k in range(len(self.z_target)):
                hf = HorizontalFilter(grid,k,self.ds.time)
                da_interp_go = da.isel(zt=k)
                var_i, var_filtered = hf.apply_filter(da_interp_go)
                #if k==len(self.z_target)-1:
                #    print("plot")
                #    hf.plot_filtered_field(var_i,var_filtered,var)
                da_filtered.append(var_filtered)
            filtered = xr.concat(da_filtered, dim="zt")
            filtered = filtered.assign_coords(zt=self.z_target)
        else:
            filtered = da
        regrid = HorizontalRegridder(filtered, self.target_grid,weights_path=self.weight_path/'weights_{}.nc'.format(self.weight_name))
        final_var = regrid(filtered)
        final_var = final_var.fillna(0)
        out.append(final_var.to_dataset(name=var))
        return xr.merge(out)

    def run_2D(self, var):
        out = []
        da = self.ds[var]
        mask_ocean = np.isfinite(da.isel(time=0))
        grid = LatLonGrid(
            self.ds.lat.values,
            self.ds.lon.values,
            mask=mask_ocean
        )
        if self.filtering:
            hf = HorizontalFilter(grid,-1,self.ds.time)
            var_i, var_filtered = hf.apply_filter(da)
        else:
            var_filtered = da
        regrid = HorizontalRegridder(var_filtered, self.target_grid,weights_path=self.weight_path/'weights_{}.nc'.format(self.weight_name))
        final_var = regrid(var_filtered)
        final_var = final_var.fillna(0)
        out.append(final_var.to_dataset(name=var))
        return xr.merge(out)
        
    def run_vector(self, var):
        da_u = self.ds[var[0]]
        da_v = self.ds[var[1]]
        vp = VerticalProcessor(self.ds, z_name="zt")
        da_u = vp.interpolate(da_u, self.z_target)
        da_v = vp.interpolate(da_v, self.z_target)
        mask_ocean = np.isfinite(da_u.isel(time=0))
        grid = LatLonGrid(
            self.ds.lat.values,
            self.ds.lon.values,
            mask=mask_ocean
        )
        if self.filtering:
            da_u_filtered = []
            da_v_filtered = []
            for k in range(len(self.z_target)):
                hf = HorizontalFilter(grid,k,self.ds.time)
                da_u_go = da_u.isel(zt=k)
                da_v_go = da_v.isel(zt=k)
                u, u_filtered = hf.apply_filter(da_u_go)
                v, v_filtered = hf.apply_filter(da_v_go)
                #if k==len(self.z_target)-1:
                #    print("plot")
                #    hf.plot_filtered_field(var_i,var_filtered,var)
                da_u_filtered.append(u_filtered)
                da_v_filtered.append(v_filtered)
            
            u_filtered = xr.concat(da_u_filtered, dim="zt")
            v_filtered = xr.concat(da_v_filtered, dim="zt")
            u_filtered = u_filtered.assign_coords(zt=self.z_target)
            v_filtered = v_filtered.assign_coords(zt=self.z_target)
        else:
            u_filtered = da_u
            v_filtered = da_v
        regrid_u = HorizontalRegridder(u_filtered, self.target_grid_u,weights_path=self.weight_path/'weights_u_{}.nc'.format(self.weight_name))
        regrid_v = HorizontalRegridder(v_filtered, self.target_grid_v,weights_path=self.weight_path/'weights_v_{}.nc'.format(self.weight_name))
        final_u = regrid_u(u_filtered)
        final_u = final_u.rename({"lon":"lon_u"})
        final_u = final_u.transpose("time", "zt", "lat", "lon_u")
        final_u = final_u.fillna(0)
        final_v = regrid_v(v_filtered)
        final_v = final_v.rename({"lat":"lat_u"})
        final_v = final_v.transpose("time", "zt", "lat_u", "lon")
        final_v = final_v.fillna(0)
        return xr.Dataset({var[0]:final_u, var[1]:final_v})

    def run_vector_2D(self, var):

        da_u = self.ds[var[0]]
        da_v = self.ds[var[1]]

        mask_ocean = np.isfinite(da_u.isel(time=0))

        grid = LatLonGrid(
            self.ds.lat.values,
            self.ds.lon.values,
            mask=mask_ocean
        )

        if self.filtering:
            hf = HorizontalFilter(grid, -1, self.ds.time)

            u, u_filtered = hf.apply_filter(da_u)
            v, v_filtered = hf.apply_filter(da_v)
        else:
            u_filtered = da_u
            v_filtered = da_v

        regrid_u = HorizontalRegridder(
            u_filtered,
            self.target_grid_u,
            weights_path=self.weight_path/'weights_u_{}.nc'.format(self.weight_name)
        )

        regrid_v = HorizontalRegridder(
            v_filtered,
            self.target_grid_v,
            weights_path=self.weight_path/'weights_v_{}.nc'.format(self.weight_name)
        )

        final_u = regrid_u(u_filtered)
        final_u = final_u.rename({"lon": "lon_u"})
        final_u = final_u.transpose("time", "lat", "lon_u")
        final_u = final_u.fillna(0)

        final_v = regrid_v(v_filtered)
        final_v = final_v.rename({"lat": "lat_u"})
        final_v = final_v.transpose("time", "lat_u", "lon")
        final_v = final_v.fillna(0)

        return xr.Dataset({
            var[0]: final_u,
            var[1]: final_v
        })


    def write(self, variables, output_dir, filename="output.nc"):
        os.makedirs(output_dir, exist_ok=True)
        output_path = os.path.join(output_dir, filename)
         
        for var in variables: 
            if not os.path.exists(output_path):
                if isinstance(var,list):
                    data_var = self.ds[var[0]]
                    ndim = data_var.ndim
                    if ndim == 4:
                        ds_new = self.run_vector(var)
                        ds_u = ds_new[var[0]].chunk({"time":self.chunker["time"],"zt":self.chunker["zt"],"lat":-1,"lon_u":-1,})
                        ds_v = ds_new[var[1]].chunk({"time":self.chunker["time"],"zt":self.chunker["zt"],"lat_u":-1,"lon":-1,})
                    elif ndim ==3:
                        ds_new = self.run_vector_2D(var)
                        ds_u = ds_new[var[0]].chunk({"time":self.chunker["time"],"lat":-1,"lon_u":-1,})
                        ds_v = ds_new[var[1]].chunk({"time":self.chunker["time"],"lat_u":-1,"lon":-1,})

                    with ProgressBar():
                        ds_u.to_netcdf(output_path, mode="w",compute=True)
                    with ProgressBar():
                        ds_v.to_netcdf(output_path, mode="a",compute=True)
                    print(f"Created file and wrote variables '{var[0]}' & '{var[1]}'")
                else:
                    data_var = self.ds[var]
                    ndim = data_var.ndim
                    if ndim == 4:
                        ds_new = self.run(var)
                        ds_new = ds_new.transpose("time", "zt", "lat", "lon")
                        ds_new = ds_new.chunk({"time": self.chunker["time"], "zt": self.chunker["zt"], "lat": -1, "lon": -1})
                    elif ndim == 3:
                        ds_new = self.run_2D(var)
                        ds_new = ds_new.transpose("time", "lat", "lon")
                        ds_new = ds_new.chunk({"time": self.chunker["time"], "lat": -1, "lon": -1})
                    with ProgressBar():
                        ds_new.to_netcdf(output_path, mode="w",compute=True)
                    print(f"Created file and wrote variable '{var}'")
                continue
    
            with xr.open_dataset(output_path) as ds_out:
                if isinstance(var,list):
                    if var[0] and var[1] in ds_out.data_vars:
                        print(f"Variables '{var[0]}' and '{var[1]}' already exist — skipping")
                        continue
                else:
                    if var in ds_out.data_vars:
                        print(f"Variable '{var}' already exists — skipping")
                        continue

            if isinstance(var,list):
                data_var = self.ds[var[0]]
                ndim = data_var.ndim
                if ndim == 4:
                    ds_new = self.run_vector(var)
                    ds_u = ds_new[var[0]].chunk({"time":self.chunker["time"],"zt":self.chunker["zt"],"lat":-1,"lon_u":-1,})
                    ds_v = ds_new[var[1]].chunk({"time":self.chunker["time"],"zt":self.chunker["zt"],"lat_u":-1,"lon":-1,})
                elif ndim ==3:
                    ds_new = self.run_vector_2D(var)
                    ds_u = ds_new[var[0]].chunk({"time":self.chunker["time"],"lat":-1,"lon_u":-1,})
                    ds_v = ds_new[var[1]].chunk({"time":self.chunker["time"],"lat_u":-1,"lon":-1,})

                with ProgressBar():
                    ds_u.to_netcdf(output_path, mode="a",compute=True)
                with ProgressBar():
                    ds_v.to_netcdf(output_path, mode="a",compute=True)
                print(f"Appended variables '{var[0]}' & '{var[1]}'")

            else:
                data_var = self.ds[var]
                ndim = data_var.ndim
                if ndim == 4:
                    ds_new = self.run(var)
                    ds_new = ds_new.transpose("time", "zt", "lat", "lon")
                    ds_new = ds_new.chunk({"time": self.chunker["time"], "zt": self.chunker["zt"], "lat": -1, "lon": -1})
                elif ndim == 3:
                    ds_new = self.run_2D(var)
                    ds_new = ds_new.transpose("time", "lat", "lon")
                    ds_new = ds_new.chunk({"time": self.chunker["time"], "lat": -1, "lon": -1})
                with ProgressBar():
                    ds_new.to_netcdf(output_path, mode="a",compute=True)
                print(f"Appended variable '{var}' to existing file")
    
        return output_path


class LatLonGrid:
    R = 6371000.0

    def __init__(self, lat, lon, mask=None):
        self.lat = lat
        self.lon = lon
        self.mask = mask

        (
            self.area,
            self.dxw,
            self.dyw,
            self.dxs,
            self.dys
        ) = self._compute_metrics()


    def _dx(self, lat, dlon):
        return self.R * np.cos(np.deg2rad(lat)) * np.deg2rad(dlon)

    def _dy(self, dlat):
        return self.R * np.deg2rad(dlat)

    def _compute_metrics(self):
        lat, lon = self.lat, self.lon
        Ny, Nx = len(lat), len(lon)

        dlat = np.diff(lat)
        dlon = np.diff(lon)

        dxw = np.zeros((Ny, Nx))
        dyw = np.zeros((Ny, Nx))
        dxs = np.zeros((Ny, Nx))
        dys = np.zeros((Ny, Nx))

        # dx at west face
        for j in range(1, Nx):
            dxw[:, j] = self._dx(lat, dlon[j - 1])
        dxw[:, 0] = dxw[:, 1]

        # dy at west face
        for i in range(1, Ny):
            dyw[i, :] = self._dy(dlat[i - 1])
        dyw[0, :] = dyw[1, :]

        # dx at south face
        for i in range(1, Ny):
            dxs[i, :] = self._dx(lat[i - 1], dlon.mean())
        dxs[0, :] = dxs[1, :]

        # dy at south face
        for j in range(1, Nx):
            dys[:, j] = self._dy(dlat.mean())
        dys[:, 0] = dys[:, 1]

        area = dxw * dyw

        return (
            np.abs(area),
            np.abs(dxw),
            np.abs(dyw),
            np.abs(dxs),
            np.abs(dys)
        )


    def as_xarray(self,z):
        coords = {"lat": self.lat, "lon": self.lon}
        if z >= 0:
            out_dict = dict(
                dxw=xr.DataArray(self.dxw, coords=coords, dims=("lat", "lon")),
                dyw=xr.DataArray(self.dyw, coords=coords, dims=("lat", "lon")),
                dxs=xr.DataArray(self.dxs, coords=coords, dims=("lat", "lon")),
                dys=xr.DataArray(self.dys, coords=coords, dims=("lat", "lon")),
                area=xr.DataArray(self.area, coords=coords, dims=("lat", "lon")),
                wet_mask=xr.DataArray(self.mask.isel(zt=z), coords=coords, dims=("lat", "lon"))
                if self.mask is not None else None
            )
        else:
            out_dict = dict(
                dxw=xr.DataArray(self.dxw, coords=coords, dims=("lat", "lon")),
                dyw=xr.DataArray(self.dyw, coords=coords, dims=("lat", "lon")),
                dxs=xr.DataArray(self.dxs, coords=coords, dims=("lat", "lon")),
                dys=xr.DataArray(self.dys, coords=coords, dims=("lat", "lon")),
                area=xr.DataArray(self.area, coords=coords, dims=("lat", "lon")),
                wet_mask=xr.DataArray(self.mask, coords=coords, dims=("lat", "lon"))
                if self.mask is not None else None
            )
        return out_dict

class VerticalProcessor:
    def __init__(self, ds, z_name="zt"):
        self.ds = ds
        self.z_name = z_name

    def vertical_filter(self, var, window=3):
        # very light smoothing (optional)
        return var.rolling({self.z_name: window}, center=True).mean()

    def interpolate(self, var, z_target):
        return var.interp({self.z_name: z_target})

class HorizontalFilter:
    def __init__(self, grid,z,time, filter_scale_km=400):
        self.grid = grid
        self.time = time
        self.lat = self.grid.lat
        self.lon = self.grid.lon
        gz = self.grid.as_xarray(z)
        self.filter_scale = filter_scale_km * 1000

        dxw = self._hc(gz["dxw"])
        dyw = self._hc(gz["dyw"])
        dxs = self._hc(gz["dxs"])
        dys = self._hc(gz["dys"])
        mask = self._hc(gz["wet_mask"])
        area = self._hc(gz["area"])
        self.dx_min = min(
            dxw.min(),
            dyw.min(),
            dxs.min(),
            dys.min()
        ).values

        kappa_w = self._hc(xr.ones_like(dxw))
        kappa_s = self._hc(xr.ones_like(dyw))

        
        self.filter = gcm_filters.Filter(
            filter_scale=self.filter_scale,
            dx_min=self.dx_min,
            filter_shape=gcm_filters.FilterShape.GAUSSIAN,
            grid_type=gcm_filters.GridType.IRREGULAR_WITH_LAND,
            grid_vars={
                "dxw": dxw,
                "dyw": dyw,
                "dxs": dxs,
                "dys": dys,
                "wet_mask":mask,
                "area":area,
                "kappa_w": kappa_w,
                "kappa_s": kappa_s,
            }
        )

    def _hc(self,var):
        return var.chunk({"lat":-1,"lon":-1})

    def apply_filter(self,da):
        var = xr.DataArray(da,
                           coords={"time": self.time , "lat": self.lat, "lon": self.lon},
                           dims=("time", "lat", "lon"))
        var = var.chunk({"time": 1,"lat":-1, "lon":-1})
        var_filtered = self.filter.apply(var, dims=["lat", "lon"])
        return var, var_filtered

    def apply_vector_filter(self,u,v):
        u_var = xr.DataArray(u,
                           coords={"time": self.time , "lat": self.lat, "lon": self.lon},
                           dims=("time", "lat", "lon"))
        v_var = xr.DataArray(v,
                           coords={"time": self.time , "lat": self.lat, "lon": self.lon},
                           dims=("time", "lat", "lon"))
        u_var = u_var.chunk({"time": 1,"lat":-1, "lon":-1})
        v_var = v_var.chunk({"time": 1,"lat":-1, "lon":-1})
        (u_var_filtered, v_var_filtered) = self.filter.apply_to_vector(u_var,v_var, dims=["lat", "lon"])
        return u_var_filtered, v_var_filtered

    def plot_filtered_field(self, var_to_filter, var_filtered, var_name, time_index=0,
                            output_dir="/Odyssey/private/e25cheve/transfert/"):

        vmin = float(var_to_filter.isel(time=time_index).min())
        vmax = float(var_to_filter.isel(time=time_index).max())
        units = var_filtered.units
        fig, axs = plt.subplots(1, 2, figsize=(25, 8))

        var_to_filter.isel(time=time_index).plot(
            ax=axs[0],
            vmin=vmin, vmax=vmax,
            cmap="plasma",
            cbar_kwargs={'label': units}
        )
        axs[0].set(title=f'unfiltered field')

        var_filtered.isel(time=time_index).plot(
            ax=axs[1],
            vmin=vmin, vmax=vmax,
            cmap="plasma",
            cbar_kwargs={'label': units}
        )
        axs[1].set(title=f'filtered field',
                   ylabel='')

        #(var_to_filter - var_filtered).isel(time=time_index).plot(
        #    ax=axs[2],
        #    cmap="RdBu_r",
        #    cbar_kwargs={'label': units}
        #)
        #axs[2].set(title=f'Residual ',
        #           ylabel='')

        plt.tight_layout()
        plt.savefig(output_dir+'filtered_{}.png'.format(var_name),format='png',dpi=300)
        plt.close()


class HorizontalRegridder:
    def __init__(self, da_src, ds_tgt,weights_path=None, method="bilinear",verbose=True):
        self.lon_tgt_original = ds_tgt["lon"].values
        grid_src = da_src.to_dataset()[["lat", "lon"]].copy()
        grid_tgt = ds_tgt[["lat", "lon"]].copy()

        grid_tgt["lon"] = grid_tgt["lon"] % 360

        grid_tgt = grid_tgt.sortby("lon")

        self.grid_src = grid_src
        self.grid_tgt = grid_tgt

        if os.path.exists(weights_path):
            if verbose:
                print("Reusing existing weights:", weights_path)
            self.regridder = xe.Regridder(
                self.grid_src,
                self.grid_tgt,
                method,
                periodic=True,
                reuse_weights=True,
                filename=weights_path,
            )
        else:
            if verbose:
                print("Computing weights (one-time cost)…")
            self.regridder = xe.Regridder(
                self.grid_src,
                self.grid_tgt,
                method,
                periodic=True,
                reuse_weights=False,
                filename=weights_path,
            )

    def __call__(self, ds):
        out = self.regridder(ds)
        out = out.assign_coords(lon=self.lon_tgt_original)
        out = out.sortby("lon")
        shift = np.argmin(np.abs(self.grid_tgt.lon.values - self.lon_tgt_original[0]))
        out = out.roll(lon=-shift, roll_coords=False)
        return out

def main():
    client = Client(n_workers=16, threads_per_worker=1)

    variables = [['ewss','nsss'],'tp','e','ssr','str','sshf','slhf']
    #glorys_pattern = "/Odyssey/private/e25cheve/data/Glorys_Climatology/mercatorglorys12v1_gl12_mean_1993_*"
    pattern = '/Odyssey/private/e25cheve/data/ERA5_1deg_Flux/*'

    grid_file = "/Odyssey/private/e25cheve/data/interp_grid_60_levels.nc"

    files = sorted(glob(pattern))

    for f in files:
        print(f"Processing {f}")

        Glo = Preprocessor(
            f,
            grid_file,
            single_depth=True,
            verbose=True,
            filtering=False,
            weight_name='era5'
        )

        output_name = os.path.basename(f).replace(".nc", "_processed.nc")

        Glo.write(
            variables,
            "/Odyssey/private/e25cheve/data/ERA5_1deg_Flux_processed/",
            output_name
        )


if __name__ == "__main__":
    main()

#Lat = Glo.ds.lat
#Lon = Glo.ds.lon
#dxw = grid.dxw
#dyw = grid.dyw
#dxs = grid.dxs
#dys = grid.dys
#area = grid.area
#
#
#fig, axs = plt.subplots(3, 2, figsize=(18, 12))
#cmap = "viridis"
#
#im0 = axs[0, 0].pcolormesh(Lon, Lat, dxw, cmap=cmap)
#axs[0, 0].set_title("dxw (m)")
#fig.colorbar(im0, ax=axs[0, 0])
#
#im1 = axs[0, 1].pcolormesh(Lon, Lat, dyw, cmap=cmap)
#axs[0, 1].set_title("dyw (m)")
#fig.colorbar(im1, ax=axs[0, 1])
#
#im2 = axs[1, 0].pcolormesh(Lon, Lat, dxs, cmap=cmap)
#axs[1, 0].set_title("dxs (m)")
#fig.colorbar(im2, ax=axs[1, 0])
#
#im3 = axs[1, 1].pcolormesh(Lon, Lat, dys, cmap=cmap)
#axs[1, 1].set_title("dys (m)")
#fig.colorbar(im3, ax=axs[1, 1])
#
#im4 = axs[2, 0].pcolormesh(Lon, Lat, area, cmap=cmap)
#axs[2, 0].set_title("Cell Area (m²)")
#fig.colorbar(im4, ax=axs[2, 0])
#
#axs[2, 1].axis("off")
#
#for ax in axs.flat:
#    ax.set_xlabel("Longitude")
#    ax.set_ylabel("Latitude")
#
#plt.tight_layout()
#plt.savefig("/Odyssey/private/e25cheve/transfert/grid_assets_glorys")
#plt.close()
#da_interp = Glo.run(variables)
#
#
#
## --- Extract variables ---
#theta_orig = Glo.ds["thetao"]
#theta_interp = da_interp["thetao"]  # already a DataArray
#
## --- Select one time (important for clarity) ---
#theta_orig = theta_orig.isel(time=0)
#theta_interp = theta_interp.isel(time=0)
#
## --- Compute zonal mean (average over longitude) ---
#theta_orig_zm = theta_orig.mean(dim="lon", skipna=True)
#theta_interp_zm = theta_interp.mean(dim="lon", skipna=True)
#
## --- Coordinates ---
#lat_orig = theta_orig_zm["lat"]
#lat_interp = theta_interp_zm["lat"]
#
#z_orig = theta_orig_zm["depth"] if "depth" in theta_orig_zm.coords else theta_orig_zm["zt"]
#z_interp = theta_interp_zm["zt"]
#
#vmin = np.nanpercentile(theta_orig_zm, 1)
#vmax = np.nanpercentile(theta_orig_zm, 99)
#levels = np.linspace(vmin, vmax, 40)
#
#fig, axes = plt.subplots(
#    nrows=1, ncols=2,
#    figsize=(16, 6),
#    sharey=True,
#    constrained_layout=True
#)
#
## ---- BEFORE vertical interpolation ----
#cf1 = axes[0].contourf(
#    lat_orig, z_orig,
#    theta_orig_zm,
#    levels=levels,
#    cmap="RdBu_r",
#    extend="both"
#)
#axes[0].set_title("Before vertical interpolation")
#axes[0].set_xlabel("Latitude")
#axes[0].set_ylabel("Depth (m)")
#axes[0].invert_yaxis()
#
## ---- AFTER vertical interpolation ----
#cf2 = axes[1].contourf(
#    lat_interp, z_interp,
#    theta_interp_zm,
#    levels=levels,
#    cmap="RdBu_r",
#    extend="both"
#)
#axes[1].set_title("After vertical interpolation")
#axes[1].set_xlabel("Latitude")
#axes[1].invert_yaxis()
#
## ---- Colorbar ----
#cbar = fig.colorbar(cf2, ax=axes, orientation="vertical", shrink=0.9)
#cbar.set_label("Potential temperature (°C)")
#
#plt.savefig("/Odyssey/private/e25cheve/transfert/vert_interp.pdf")
#import cartopy.crs as ccrs
#import cartopy.feature as cfeature
#import scienceplots
#plt.style.use('science')
#
#
#time_index=0
#output_dir = "/Odyssey/private/e25cheve/transfert/"
#ds_glorys = xr.open_dataset(glorys_path, chunks={"time": 1})
#vo_glo = ds_glorys['vo']
#
#res_path ="/Odyssey/private/e25cheve/data/Glorys_CI_test_1993_01.nc"
#ds_res = xr.open_dataset(res_path, chunks={"time": 1})
#ds_res = ds_res.assign_coords(
#        lon=(((ds_res.lon + 180) % 360) - 180)
#).sortby("lon")
#ds_res = ds_res.assign_coords(
#    zt = -ds_res.zt
#).sortby("zt")
#
#vo_res = ds_res['vo']



#thetao_glo_zm = (
#    thetao_glo
#    .isel(time=time_index)
#    .mean(dim="longitude", skipna=True)
#)
#
#uo_res_zm = (
#    uo_res
#    .isel(time=time_index)
#    .mean(dim="lon", skipna=True)
#)
#vmin = float(min(thetao_glo_zm.min(), uo_res_zm.min()))
#vmax = float(max(thetao_glo_zm.max(), uo_res_zm.max()))
#units = uo_glo.units
#
#fig, axs = plt.subplots(1, 2, figsize=(25, 8), sharey=True)
#
#thetao_glo_zm.plot(
#    ax=axs[0],
#    x="latitude",
#    y="depth",
#    vmin=vmin,
#    vmax=vmax,
#    cmap="coolwarm",
#    cbar_kwargs={"label": units}
#)
#axs[0].set(
#    title="Original field – zonal mean",
#    xlabel="Latitude"
#)
#
#thetao_res_zm.plot(
#    ax=axs[1],
#    x="lat",
#    y="zt",
#    vmin=vmin,
#    vmax=vmax,
#    cmap="coolwarm",
#    cbar_kwargs={"label": units}
#)
#axs[1].set(
#    title="Processed field – zonal mean",
#    xlabel="Latitude",
#    ylabel=""
#)
#axs[0].invert_yaxis()












#vmin_glo = float(vo_glo.isel(time=time_index,depth=5).min())
#vmax_glo = float(vo_glo.isel(time=time_index,depth=5).max())
#vmin_res = float(vo_res.isel(time=time_index,zt=0).min())
#vmax_res = float(vo_res.isel(time=time_index,zt=0).max())
#vmin = min(vmin_glo,vmin_res)
#vmax = max(vmax_glo,vmax_res)
#units = vo_glo.units
#levels = np.linspace(vmin, vmax, 100)
#fig, axs = plt.subplots(
#    1, 2,
#    figsize=(25, 7),
#    subplot_kw={'projection': ccrs.PlateCarree()}
#)
#
#cf1 = axs[0].contourf(
#    vo_glo.longitude, vo_glo.latitude,
#    vo_glo.isel(time=time_index, depth=5),
#    levels=levels,
#    cmap='RdBu',
#    extend='both',
#    transform=ccrs.PlateCarree()
#)
#
#axs[0].add_feature(cfeature.GSHHSFeature('low', levels=[1, 2, 6]))
#axs[0].set_title('Original field at −6 m',fontsize=20)
#
#cf2 = axs[1].contourf(
#    vo_res.lon, vo_res.lat_u,
#    vo_res.isel(time=time_index, zt=0),
#    levels=levels,
#    cmap='RdBu',
#    extend='both',
#    transform=ccrs.PlateCarree()
#)
#
#axs[1].add_feature(cfeature.GSHHSFeature('low', levels=[1, 2, 6]))
#axs[1].set_title('Processed field at −6 m',fontsize=20)
#
#cbar = fig.colorbar(
#    cf2,
#    ax=axs,
#    fraction=0.046,
#    pad=0.04
#)
#cbar.set_label(units,size=15)
#cbar.ax.tick_params(labelsize=10)

#(var_to_filter - var_filtered).isel(time=time_index).plot(
#    ax=axs[2],
#    cmap="RdBu_r",
#    cbar_kwargs={'label': units}
#)
#axs[2].set(title=f'Residual ',
#           ylabel='')

