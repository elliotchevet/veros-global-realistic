import argparse
import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
import gcm_filters
import os
import xesmf as xe
import warnings
from dask.diagnostics import ProgressBar
from dask.distributed import Client
from pathlib import Path
from veros_global_realistic.paths import find_veros_assets_dir

warnings.filterwarnings("ignore")

PLUGIN_DIR = Path(__file__).resolve().parents[1]
DATA_DIR = Path(os.environ.get("VEROS_GLOBAL_REALISTIC_DATA_DIR", PLUGIN_DIR / "data")).expanduser()
ERA5_DIR = Path(os.environ.get("ERA5_DIR", DATA_DIR / "ERA5")).expanduser()
GLORYS12_DIR = Path(os.environ.get("GLORYS12_DIR", DATA_DIR / "GLORYS12")).expanduser()
GLORYS12_IC_DIR = GLORYS12_DIR / "Initial_States"
GLORYS12_RESTORING_DIR = GLORYS12_DIR / "Restoring"
WEIGHTS_DIR = Path(os.environ.get("WEIGHTS_DIR", DATA_DIR / "weights")).expanduser()
VEROS_ASSETS_DIR = find_veros_assets_dir()
ASSETS_DIR = VEROS_ASSETS_DIR / "global_1deg_realistic"

GLORYS12_OPTIONS = {
    "filtering": True,
    "weight_name": "glorys12",
}

DATASET_CONFIG = {
    "ERA5": {
        "directory": ERA5_DIR,
        "variables": ["str", ["ewss", "nsss"], "tp", "e", "ssr", "sshf", "slhf"],
        "single_depth": True,
        "filtering": False,
        "weight_name": "era5",
    },
    "GLORYS12_IC": {
        "single_depth": False,
        "directory": GLORYS12_IC_DIR,
        "variables": ["mlotst", "zos", ["uo", "vo"], "thetao", "so"],
        **GLORYS12_OPTIONS,
    },
    "GLORYS12_Restoring": {
        "single_depth": True,
        "directory": GLORYS12_RESTORING_DIR,
        "variables": ["mlotst", "zos", "thetao", "so"],
        **GLORYS12_OPTIONS,
    },
}

class Preprocessor:
    def __init__(
        self,
        src_path,
        target_path,
        single_depth=False,
        use_mask=True,
        verbose=False,
        ncores=16,
        filtering=True,
        filter_scale_km=400,
        weight_path=WEIGHTS_DIR,
        weight_name="glorys12",
    ):
        self.filtering = filtering
        self.filter_scale_km = filter_scale_km
        self.verbose = verbose
        self.weight_name = weight_name
        self.weight_path = Path(weight_path)
        self.weight_path.mkdir(parents=True, exist_ok=True)
        self.use_mask = use_mask
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
        if self.ds.dims["lat"] > 2000:
            self.ds = self.ds.isel(
                lat=slice(0, -100),
            )

        target_grid = xr.open_dataset(target_path, chunks={"time": 1})
        if "xt" in target_grid:
            if verbose: print("adjusting xt to lon")
            self.target_grid = target_grid.rename({"xt": "lon", "yt": "lat", "maskT": "mask"})
        if "xu" in target_grid:
            if verbose: print("adjusting xu to lon")
            self.target_grid_u = target_grid.rename({"xu": "lon", "yt": "lat", "maskU": "mask"})
        if "yu" in target_grid:
            if verbose: print("adjusting yu to lat")
            self.target_grid_v = target_grid.rename({"xt": "lon", "yu": "lat", "maskV": "mask"})
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
            self.chunker["time"] = -1

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
        da_processed = []
        for k in range(len(self.z_target)):
            mask_ocean_k = mask_ocean.isel(zt=k).astype("int8")
            if self.use_mask == False:
                mask_ocean_k = None
            hf = HorizontalFilter(grid, k, self.ds.time, filter_scale_km=self.filter_scale_km)
            da_interp_go = da.isel(zt=k)
            var_i, var_filtered = hf.apply_filter(da_interp_go)
            regrid = HorizontalRegridder(
                    var_filtered, mask_ocean_k, self.target_grid.isel(zt=k),
                    weights_path=self.weight_path/'weights_{}_{}.nc'.format(self.weight_name,k),
                    verbose=self.verbose,
                    )
            var_processed = regrid(var_filtered)
            da_processed.append(var_processed)
        processed = xr.concat(da_processed, dim="zt")
        processed = processed.assign_coords(zt=self.z_target)
        final_var = processed.fillna(0)
        out.append(final_var.to_dataset(name=var))
        return xr.merge(out)

    def run_2D(self, var):
        out = []
        da = self.ds[var]
        mask_ocean = np.isfinite(da.isel(time=0)).astype("int8")
        grid = LatLonGrid(
            self.ds.lat.values,
            self.ds.lon.values,
            mask=mask_ocean
        )
        if self.use_mask == False:
            mask_ocean = None
        if self.filtering:
            hf = HorizontalFilter(grid, -1, self.ds.time, filter_scale_km=self.filter_scale_km)
            var_i, var_filtered = hf.apply_filter(da)
        else:
            var_filtered = da
        regrid = HorizontalRegridder(
            var_filtered,
            mask_ocean,
            self.target_grid.isel(zt=-1),
            weights_path=self.weight_path/'weights_{}.nc'.format(self.weight_name),
            verbose=self.verbose,
        )
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
        da_u_processed = []
        da_v_processed = []
        for k in range(len(self.z_target)):
            mask_ocean_k = mask_ocean.isel(zt=k).astype("int8")
            if self.use_mask == False:
                mask_ocean_k = None
            hf = HorizontalFilter(grid, k, self.ds.time, filter_scale_km=self.filter_scale_km)
            da_u_go = da_u.isel(zt=k)
            da_v_go = da_v.isel(zt=k)
            u, u_filtered = hf.apply_filter(da_u_go)
            v, v_filtered = hf.apply_filter(da_v_go)
            regrid_u = HorizontalRegridder(
                    u_filtered, mask_ocean_k, self.target_grid_u.isel(zt=k),
                    weights_path=self.weight_path/'weights_u_{}_{}.nc'.format(self.weight_name,k),
                    verbose=self.verbose,
                    )
            regrid_v = HorizontalRegridder(
                    v_filtered, mask_ocean_k, self.target_grid_v.isel(zt=k),
                    weights_path=self.weight_path/'weights_v_{}_{}.nc'.format(self.weight_name,k),
                    verbose=self.verbose,
                    )
            var_u_processed = regrid_u(u_filtered)
            var_v_processed = regrid_v(v_filtered)
            da_u_processed.append(var_u_processed)
            da_v_processed.append(var_v_processed)

        u_processed = xr.concat(da_u_processed, dim="zt")
        v_processed = xr.concat(da_v_processed, dim="zt")
        u_processed = u_processed.assign_coords(zt=self.z_target)
        v_processed = v_processed.assign_coords(zt=self.z_target)
        final_u = u_processed.rename({"lon":"lon_u"})
        final_u = final_u.transpose("time", "zt", "lat", "lon_u")
        final_u = final_u.fillna(0)
        final_v = v_processed.rename({"lat":"lat_u"})
        final_v = final_v.transpose("time", "zt", "lat_u", "lon")
        final_v = final_v.fillna(0)
        return xr.Dataset({var[0]:final_u, var[1]:final_v})

    def run_vector_2D(self, var):

        da_u = self.ds[var[0]]
        da_v = self.ds[var[1]]

        mask_ocean = np.isfinite(da_u.isel(time=0)).astype("int8")

        grid = LatLonGrid(
            self.ds.lat.values,
            self.ds.lon.values,
            mask=mask_ocean
        )
        if self.use_mask == False:
            mask_ocean = None

        if self.filtering:
            hf = HorizontalFilter(grid, -1, self.ds.time, filter_scale_km=self.filter_scale_km)

            u, u_filtered = hf.apply_filter(da_u)
            v, v_filtered = hf.apply_filter(da_v)
        else:
            u_filtered = da_u
            v_filtered = da_v

        regrid_u = HorizontalRegridder(
            u_filtered,
            mask_ocean,
            self.target_grid_u.isel(zt=-1),
            weights_path=self.weight_path/'weights_u_{}.nc'.format(self.weight_name),
            verbose=self.verbose,
        )

        regrid_v = HorizontalRegridder(
            v_filtered,
            mask_ocean,
            self.target_grid_v.isel(zt=-1),
            weights_path=self.weight_path/'weights_v_{}.nc'.format(self.weight_name),
            verbose=self.verbose,
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

    def plot_filtered_field(self, var_to_filter, var_filtered, var_name, time_index=0, output_dir=ASSETS_DIR):
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

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

        plt.tight_layout()
        plt.savefig(output_dir / "filtered_{}.png".format(var_name),format='png',dpi=300)
        plt.close()


class HorizontalRegridder:
    def __init__(self, da_src,mask_src, ds_tgt,weights_path=None, method="bilinear",verbose=True):
        self.lon_tgt_original = ds_tgt["lon"].values
        grid_src = da_src.to_dataset()[["lat", "lon"]].copy()
        grid_tgt = ds_tgt[["lat", "lon"]].copy()
        if mask_src is not None:
            grid_src["mask"] = mask_src
            grid_tgt["mask"] = ds_tgt["mask"]
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
                extrap_method='inverse_dist'
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
                extrap_method='inverse_dist'
            )

    def __call__(self, ds):
        out = self.regridder(ds)
        out = out.assign_coords(lon=self.lon_tgt_original)
        out = out.sortby("lon")
        shift = np.argmin(np.abs(self.grid_tgt.lon.values - self.lon_tgt_original[0]))
        out = out.roll(lon=-shift, roll_coords=False)
        return out

def parse_args():
    parser = argparse.ArgumentParser(description="Preprocess ERA5 or GLORYS12 data for the global realistic setup.")
    parser.add_argument("--dataset", choices=DATASET_CONFIG, required=True, help="Dataset to preprocess.")
    parser.add_argument("--data-dir", type=Path, help="Directory containing input files. Defaults to the dataset directory.")
    parser.add_argument("--output-dir", type=Path, default=ASSETS_DIR, help="Directory for processed output files.")
    parser.add_argument("--grid-file", type=Path, default=DATA_DIR / "masks.nc", help="Target Veros grid/mask file.")
    parser.add_argument("--weights-dir", type=Path, default=WEIGHTS_DIR, help="Directory for xESMF weight files.")
    parser.add_argument("--weight-name", help="Weight filename stem. Defaults to a dataset-specific name.")
    parser.add_argument("--filter-scale-km", type=float, default=400, help="Horizontal filter scale in km.")
    parser.add_argument("--filtering", action=argparse.BooleanOptionalAction, default=None, help="Enable or disable filtering.")
    parser.add_argument("--single-depth", action=argparse.BooleanOptionalAction, default=None, help="Treat source zt as singleton.")
    parser.add_argument("--ncores", type=int, default=16, help="Number of Dask workers and chunking target.")
    parser.add_argument("--verbose", action="store_true", help="Print detailed preprocessing progress.")
    return parser.parse_args()


def input_files(data_dir):
    files = sorted(path for path in data_dir.iterdir() if path.is_file() and not path.name.startswith("."))
    if not files:
        raise FileNotFoundError(f"No input files found in {data_dir}")
    return files


def output_filename(input_path):
    if input_path.name.endswith("_processed.nc"):
        return input_path.name
    if input_path.name.endswith(".nc"):
        return input_path.name.replace(".nc", "_processed.nc")
    return f"{input_path.stem}_processed.nc"


def main():
    args = parse_args()
    config = DATASET_CONFIG[args.dataset]

    data_dir = (args.data_dir or config["directory"]).expanduser()
    output_dir = args.output_dir.expanduser()
    weights_dir = args.weights_dir.expanduser()
    output_dir.mkdir(parents=True, exist_ok=True)
    weights_dir.mkdir(parents=True, exist_ok=True)

    grid_file = args.grid_file.expanduser()
    filtering = config["filtering"] if args.filtering is None else args.filtering
    single_depth = config["single_depth"] if args.single_depth is None else args.single_depth
    weight_name = args.weight_name or config["weight_name"]

    if args.verbose:
        print(f"Dataset: {args.dataset}")
        print(f"Input directory: {data_dir}")
        print(f"Output directory: {output_dir}")
        print(f"Weights directory: {weights_dir}")
        print(f"Grid file: {grid_file}")
        print(f"Filtering: {filtering}")
        print(f"Filter scale: {args.filter_scale_km} km")
        print(f"Single depth: {single_depth}")
        print(f"Weight name: {weight_name}")

    client = Client(n_workers=args.ncores, threads_per_worker=1)
    if args.verbose:
        print(client)

    for src_file in input_files(data_dir):
        print(f"Processing {src_file}")
        preprocessor = Preprocessor(
            src_file,
            grid_file,
            single_depth=single_depth,
            use_mask=True,
            filtering=filtering,
            filter_scale_km=args.filter_scale_km,
            verbose=args.verbose,
            ncores=args.ncores,
            weight_path=weights_dir,
            weight_name=weight_name,
        )
        preprocessor.write(config["variables"], output_dir, output_filename(src_file))

if __name__ == "__main__":
    main()

