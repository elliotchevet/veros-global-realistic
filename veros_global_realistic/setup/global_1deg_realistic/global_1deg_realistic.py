__VEROS_VERSION__ = '1.6.0'

if __name__ == "__main__":
    raise RuntimeError(
        "Veros setups cannot be executed directly. "
        f"Try `veros run {__file__}` instead."
    )

# -- end of auto-generated header, original file below --

import os
import h5netcdf
import xarray as xr

from veros.core.utilities import enforce_boundaries
from veros import VerosSetup, tools, time, veros_routine, veros_kernel, KernelOutput
from veros.variables import Variable, allocate
from veros.core.operators import numpy as npx, update, update_multiply, at

BASE_PATH = os.path.dirname(os.path.realpath(__file__))
DATA_FILES = tools.get_assets("global_1deg_realistic", os.path.join(BASE_PATH, "assets.json"))


class GlobalOneDegreeSetup(VerosSetup):
    """Global 1 degree model with era5 forcing and glorys12 initialization.

    `Adapted from pyOM2 <https://wiki.zmaw.de/ifm/TO/pyOM2/1x1%20global%20model>`.

    """
    min_depth = 10.0
    max_depth = 5400.0

    @veros_routine
    def set_parameter(self, state):
        """
        set main parameters
        """
        settings = state.settings

        settings.identifier = "global_1deg_realistic"
        settings.description = "Start from Glorys"
        settings.start_date = "1986-01-01"
        settings.restart_input_filename = "./global_1deg_17280.restart.h5"

        
        settings.nx = 360
        settings.ny = 160
        settings.nz = 60
        settings.dt_mom = 3600
        settings.dt_tracer = 3600
        settings.dt_forcing = 3
        settings.runlen = (30 * 360 - 90) * 24 * 3600

        settings.x_origin = 91.0
        settings.y_origin = -79.0

        settings.coord_degree = True
        settings.enable_cyclic_x = True

        settings.enable_hor_friction = True
        settings.A_h = 5e4
        settings.enable_hor_friction_cos_scaling = True
        settings.hor_friction_cosPower = 1
        settings.enable_tempsalt_sources = True
        settings.enable_implicit_vert_friction = True

        settings.eq_of_state_type = 5

        # isoneutral
        settings.enable_neutral_diffusion = True
        settings.K_iso_0 = 1000.0
        settings.K_iso_steep = 50.0
        settings.iso_dslope = 0.005
        settings.iso_slopec = 0.005
        settings.enable_skew_diffusion = True

        # tke
        settings.enable_tke = True
        settings.c_k = 0.1
        settings.c_eps = 0.7
        settings.alpha_tke = 30.0
        settings.mxl_min = 1e-8
        settings.tke_mxl_choice = 1
        settings.kappaM_min = 2e-4
        settings.kappaH_min = 2e-5
        settings.enable_kappaH_profile = True
        settings.enable_tke_superbee_advection = True

        # eke
        settings.enable_eke = True
        settings.eke_k_max = 1e4
        settings.eke_c_k = 0.4
        settings.eke_c_eps = 0.5
        settings.eke_cross = 2.0
        settings.eke_crhin = 1.0
        settings.eke_lmin = 100.0
        settings.enable_eke_superbee_advection = True
        settings.enable_eke_isopycnal_diffusion = True

        # idemix
        settings.enable_idemix = False
        settings.enable_eke_diss_surfbot = True
        settings.eke_diss_surfbot_frac = 0.2
        settings.enable_idemix_superbee_advection = True
        settings.enable_idemix_hor_diffusion = True
        settings.enable_streamfunction= False

        # custom variables
        state.dimensions["nmonths"] = 12
        state.dimensions["n_window"] = 2
        state.var_meta.update(
            qnet=Variable("qnet", ("xt", "yt", "n_window"), "", "", time_dependent=False),
            qsol=Variable("qsol", ("xt", "yt", "n_window"), "", "", time_dependent=False),
            divpen_shortwave=Variable("divpen_shortwave", ("zt",), "", "", time_dependent=False),
            taux=Variable("taux", ("xt", "yt", "n_window"), "", "", time_dependent=False),
            tauy=Variable("tauy", ("xt", "yt", "n_window"), "", "", time_dependent=False),
            evap=Variable("evap", ("xt", "yt", "n_window"), "", "", time_dependent=False),
            tp =Variable("tp", ("xt", "yt", "n_window"), "", "", time_dependent=False),
            t_star =Variable("t_star", ("xt", "yt", "nmonths"), "", "", time_dependent=False),
            s_star =Variable("s_star", ("xt", "yt", "nmonths"), "", "", time_dependent=False),
        )

        # define ERA5 start time index
        self.n0 = self.ERA5_start_index(DATA_FILES["forcing_era5"], settings.start_date)
        print('n0 = ',self.n0)

    def _read_forcing(self, var):
        from veros.core.operators import numpy as npx
        default_var = ["dz","temperature","sst","salinity","bathymetry","swf","dqdt","tidal_energy","wind_energy"]
        if var in default_var:
            with h5netcdf.File(DATA_FILES["forcing_default"], "r") as infile:
                var = infile.variables[var]
                return npx.asarray(var).T
        else:
            raise ValueError("Variable {} not in default_var list".format(var))

    def _read_init(self,var):
        from veros.core.operators import numpy as npx
        with h5netcdf.File(DATA_FILES["ic_glorys12"], "r") as infile:
            var = infile.variables[var]
            var = npx.asarray(var).T
            return npx.where(npx.isnan(var), 0.0, var)

    def _read_restoring(self,var):
        from veros.core.operators import numpy as npx
        with h5netcdf.File(DATA_FILES["restoring_glorys12"], "r") as infile:
            var = infile.variables[var]
            return npx.asarray(var).T

    def ERA5_start_index(self,forcing_file, start_date):
        from datetime import datetime
        ref = datetime(1900, 1, 1)
        d0 = datetime.strptime(start_date, "%Y-%m-%d")
        with h5netcdf.File(forcing_file, "r") as f:
            time0 = int(f.variables["time"][0])
        h0 = int((d0 - ref).total_seconds() // 3600)
        return h0 - time0

    def _read_ERA5(self,var):
        if not hasattr(self, "_era5_ds"):
            self._era5_ds = xr.open_dataset(
                DATA_FILES["forcing_era5"],
                engine="h5netcdf",
                chunks={"time": 1},  
            )
        return self._era5_ds[var]

    def get_ERA5(self, var, n1, n2):
        from datetime import datetime, timedelta
        var_da = self._read_ERA5(var)
        sliced = var_da.isel(time=[self.n0+n1,self.n0+n2])
        if var == 'ewss':
            t0 = datetime(1900,1,1)
            print('time = ',var_da.time.values[self.n0+n1])
        return npx.asarray(sliced).T

    @veros_routine(dist_safe=False, local_variables=["dxt", "dyt", "dzt"])
    def set_grid(self, state):
        vs = state.variables
        settings = state.settings

        vs.dzt = tools.get_vinokur_grid_steps(settings.nz, self.max_depth, self.min_depth, refine_towards="lower")
        vs.dxt = update(vs.dxt, at[...], 1.0)
        vs.dyt = update(vs.dyt, at[...], 1.0)

    @veros_routine
    def set_coriolis(self, state):
        vs = state.variables
        settings = state.settings
        vs.coriolis_t = update(
            vs.coriolis_t, at[...], 2 * settings.omega * npx.sin(vs.yt[npx.newaxis, :] / 180.0 * settings.pi)
        )

    @veros_routine(dist_safe=False, local_variables=["kbot","xt","yt","zt"])
    def set_topography(self, state):
        import numpy as onp

        vs = state.variables
        settings = state.settings
        t_grid = (vs.xt[2:-2], vs.yt[2:-2], vs.zt)
        bathymetry_data = self._read_forcing("bathymetry")
        salt_raw = self._read_init("so")[:, :,:,0]
        xt_init, yt_init, zt_init = (self._read_init(k) for k in ("lon","lat","zt"))
        salt_data = tools.interpolate((xt_init, yt_init, zt_init), salt_raw, t_grid)

        mask_salt = salt_data == 0.0
        vs.kbot = update(vs.kbot, at[2:-2, 2:-2], 1 + npx.sum(mask_salt.astype("int"), axis=2))

        mask_bathy = bathymetry_data == 0
        vs.kbot = update_multiply(vs.kbot, at[2:-2, 2:-2], ~mask_bathy)
        vs.kbot = vs.kbot * (vs.kbot < settings.nz)

        # close some channels
        i, j = onp.indices((settings.nx, settings.ny))

        mask_channel = (i >= 207) & (i < 214) & (j < 5)  # i = 208,214; j = 1,5
        vs.kbot = update_multiply(vs.kbot, at[2:-2, 2:-2], ~mask_channel)

        # Aleutian islands
        mask_channel = (i == 104) & (j == 134)  # i = 105; j = 135
        vs.kbot = update_multiply(vs.kbot, at[2:-2, 2:-2], ~mask_channel)

        # Engl channel
        mask_channel = (i >= 269) & (i < 271) & (j == 130)  # i = 270,271; j = 131
        vs.kbot = update_multiply(vs.kbot, at[2:-2, 2:-2], ~mask_channel)
        vs.kbot = enforce_boundaries(vs.kbot, settings.enable_cyclic_x, local=True)

    @veros_routine(
        dist_safe=False,
        local_variables=[     
            "qnet",
            "qsol",
            "evap",
            "t_star",
            "s_star",
            "tp",
            "divpen_shortwave",
            "taux",
            "tauy",
            "temp",
            "salt",
            "ssh",
            "u",
            "v",
            "forc_iw_bottom",
            "forc_iw_surface",
            "kbot",
            "maskT",
            "maskU",
            "maskV",
            "maskW",
            "xt",
            "xu",
            "yt",
            "yu",
            "zt",
            "zw",
            "dzt",
        ],
    )
    def set_initial_conditions(self, state):
        vs = state.variables
        settings = state.settings

        rpart_shortwave = 0.58
        efold1_shortwave = 0.35
        efold2_shortwave = 23.0

        # grids

        t_grid = (vs.xt[2:-2], vs.yt[2:-2], vs.zt)
        u_grid = (vs.xu[2:-2], vs.yt[2:-2], vs.zt)
        v_grid = (vs.xt[2:-2], vs.yu[2:-2], vs.zt)

        xt_init, xu_init, yt_init, yu_init, zt_init = (self._read_init(k) for k in ("lon","lon_u","lat","lat_u","zt"))

        # initial conditions
        
        temp_raw = self._read_init("thetao")[...,0]
        temp_data = tools.interpolate((xt_init, yt_init, zt_init), temp_raw, t_grid)
        vs.temp = update(vs.temp, at[2:-2, 2:-2, :, 0], temp_data * vs.maskT[2:-2, 2:-2, :])
        vs.temp = update(vs.temp, at[2:-2, 2:-2, :, 1], temp_data * vs.maskT[2:-2, 2:-2, :])
        vs.temp = update(vs.temp, at[2:-2, 2:-2, :, 2], temp_data * vs.maskT[2:-2, 2:-2, :])

        salt_raw = self._read_init("so")[...,0]
        salt_data = tools.interpolate((xt_init, yt_init, zt_init), salt_raw, t_grid)
        vs.salt = update(vs.salt, at[2:-2, 2:-2, :, 0], salt_data * vs.maskT[2:-2, 2:-2, :])
        vs.salt = update(vs.salt, at[2:-2, 2:-2, :, 1], salt_data * vs.maskT[2:-2, 2:-2, :])
        vs.salt = update(vs.salt, at[2:-2, 2:-2, :, 2], salt_data * vs.maskT[2:-2, 2:-2, :])

        #ssh_raw = self._read_init("zos")[...,0]
        #vs.ssh = update(vs.ssh, at[2:-2, 2:-2], ssh_raw * vs.maskT[2:-2, 2:-2, -1])

        #u_raw = self._read_init("uo")[...,0]
        #u_data = tools.interpolate((xu_init, yt_init, zt_init), u_raw, u_grid)
        #vs.u = update(vs.u, at[2:-2, 2:-2, :, 0], u_data * vs.maskU[2:-2, 2:-2, :])
        #vs.u = update(vs.u, at[2:-2, 2:-2, :, 1], u_data * vs.maskU[2:-2, 2:-2, :])
        #vs.u = update(vs.u, at[2:-2, 2:-2, :, 2], u_data * vs.maskU[2:-2, 2:-2, :])

        #v_raw = self._read_init("vo")[...,0]
        #v_data = tools.interpolate((xt_init, yu_init, zt_init), v_raw, v_grid)
        #vs.v = update(vs.v, at[2:-2, 2:-2, :, 0], v_data * vs.maskV[2:-2, 2:-2, :])
        #vs.v = update(vs.v, at[2:-2, 2:-2, :, 1], v_data * vs.maskV[2:-2, 2:-2, :])
        #vs.v = update(vs.v, at[2:-2, 2:-2, :, 2], v_data * vs.maskV[2:-2, 2:-2, :])

        t_star = self._read_restoring("thetao")
        vs.t_star = update(vs.t_star, at[2:-2, 2:-2,:], t_star * vs.maskT[2:-2, 2:-2, -1, npx.newaxis])
        s_star = self._read_restoring("so")
        vs.s_star = update(vs.s_star, at[2:-2, 2:-2,:], s_star * vs.maskT[2:-2, 2:-2, -1, npx.newaxis])
        self.n_flow = -1

        if settings.enable_idemix:
            tidal_energy_data = self._read_forcing("tidal_energy")
            mask = (
                npx.maximum(0, vs.kbot[2:-2, 2:-2] - 1)[:, :, npx.newaxis]
                == npx.arange(settings.nz)[npx.newaxis, npx.newaxis, :]
            )
            tidal_energy_data *= vs.maskW[2:-2, 2:-2, :][mask].reshape(settings.nx, settings.ny) / settings.rho_0
            vs.forc_iw_bottom = update(vs.forc_iw_bottom, at[2:-2, 2:-2], tidal_energy_data)

            wind_energy_data = self._read_forcing("wind_energy")
            wind_energy_data *= vs.maskW[2:-2, 2:-2, -1] / settings.rho_0 * 0.2
            vs.forc_iw_surface = update(vs.forc_iw_surface, at[2:-2, 2:-2], wind_energy_data)

        """
        Initialize penetration profile for solar radiation and store divergence in divpen
        note that pen is set to 0.0 at the surface instead of 1.0 to compensate for the
        shortwave part of the total surface flux
        """
        swarg1 = vs.zw / efold1_shortwave
        swarg2 = vs.zw / efold2_shortwave
        pen = rpart_shortwave * npx.exp(swarg1) + (1.0 - rpart_shortwave) * npx.exp(swarg2)

        vs.divpen_shortwave = allocate(state.dimensions, ("zt",))
        vs.divpen_shortwave = update(vs.divpen_shortwave, at[1:], (pen[1:] - pen[:-1]) / vs.dzt[1:])
        vs.divpen_shortwave = update(vs.divpen_shortwave, at[0], pen[0] / vs.dzt[0])

    @veros_routine
    def set_forcing(self, state):
        vs = state.variables
        settings = state.settings
        year_in_seconds = time.convert_time(1.0, "years", "seconds")
        n_forcing = 365 * 24 / settings.dt_forcing
        (n1, f1), (n2, f2) = tools.get_periodic_interval(vs.time, year_in_seconds, year_in_seconds / n_forcing, n_forcing)
        if self.n_flow != n1:
            self.n_flow = n1
            n1 = n1*settings.dt_forcing
            n2 = n2*settings.dt_forcing
            print("forcing update n°{}".format(n1))
            taux_data = self.get_ERA5("ewss",n1,n2)/3600
            vs.taux = update(vs.taux, at[2:-2, 2:-2, :], taux_data) 
            tauy_data = self.get_ERA5("nsss",n1,n2)/3600
            vs.tauy = update(vs.tauy, at[2:-2, 2:-2, :], tauy_data) 
            #vs.taux = enforce_boundaries(vs.taux, settings.enable_cyclic_x)
            #vs.tauy = enforce_boundaries(vs.tauy, settings.enable_cyclic_x)
            vs.update(set_forcing_kernel(state))
            str_ = self.get_ERA5("str",n1,n2)
            Qir = str_/3600
            slhf = self.get_ERA5("slhf",n1,n2)
            Qlat = slhf/3600
            sshf = self.get_ERA5("sshf",n1,n2)
            Qh = sshf/3600
            Qnet = Qir + Qlat + Qh
            vs.qnet = update(vs.qnet, at[2:-2, 2:-2, :], Qnet * vs.maskT[2:-2, 2:-2, -1, npx.newaxis])
            ssr = self.get_ERA5("ssr",n1,n2)
            Qsol = ssr/3600
            vs.qsol = update(vs.qsol, at[2:-2, 2:-2, :], Qsol* vs.maskT[2:-2, 2:-2, -1, npx.newaxis])
            
            # Evaporation 
            evap_data = self.get_ERA5("e",n1,n2)/3600
            vs.evap = update(vs.evap, at[2:-2, 2:-2, :], evap_data * vs.maskT[2:-2, 2:-2, -1, npx.newaxis])

            # Total precipitation
            tp_data = self.get_ERA5("tp",n1,n2)/3600
            vs.tp = update(vs.tp, at[2:-2, 2:-2, :], tp_data * vs.maskT[2:-2, 2:-2, -1, npx.newaxis])

        vs.update(set_forcing_kernel(state))

    @veros_routine
    def set_diagnostics(self, state):
        settings = state.settings

        average_vars = [
            "temp",
            "salt",
            "ssh",
            "u",
            "v",
            "w",
            "Nsqr",
            "Hd",
            "K_diss_v",
            "P_diss_v",
            "P_diss_nonlin",
            "P_diss_iso",
            "kappaH",
            "maskT"
        ]
        if settings.enable_skew_diffusion:
            average_vars += ["B1_gm", "B2_gm"]
        if settings.enable_TEM_friction:
            average_vars += ["kappa_gm", "K_diss_gm"]
        if settings.enable_tke:
            average_vars += ["tke", "Prandtlnumber", "mxl", "tke_diss", "forc_tke_surface", "tke_surf_corr"]
        if settings.enable_idemix:
            average_vars += ["E_iw", "forc_iw_surface", "forc_iw_bottom", "iw_diss", "c0", "v0"]
        if settings.enable_eke:
            average_vars += ["eke", "K_gm", "L_rossby", "L_rhines"]

        state.diagnostics["averages"].output_variables = average_vars
        state.diagnostics["cfl_monitor"].output_frequency = 10 * 86400.0
        state.diagnostics["snapshot"].output_frequency = 360 * 86400  
        state.diagnostics["overturning"].output_frequency = 360 * 86400
        state.diagnostics["overturning"].sampling_frequency = 5 * 86400 
        state.diagnostics["energy"].output_frequency = 5*86400
        state.diagnostics["energy"].sampling_frequency = 86400/2
        state.diagnostics["averages"].output_frequency = 360 * 86400 / 12
        state.diagnostics["averages"].sampling_frequency = 86400/2

    @veros_routine
    def after_timestep(self, state):
        pass


@veros_kernel
def set_forcing_kernel(state):
    vs = state.variables
    settings = state.settings

    cp_0 = 3991.86795711963  # J/kg /K
    dQdT = -40 # W/m2 /K (2 months for a 50m mld)
    Ks = -6.37e-6 # piston velocity in m/s (3 months for a 50m mld)

    year_in_seconds = time.convert_time(1.0, "years", "seconds")
    n_forcing = 365 * 24 / settings.dt_forcing
    (n1, f1), (n2, f2) = tools.get_periodic_interval(vs.time, year_in_seconds, year_in_seconds / n_forcing, n_forcing)

    (n1_rest, f1_rest), (n2_rest, f2_rest) = tools.get_periodic_interval(vs.time, year_in_seconds, year_in_seconds / 12.0, 12)
    vs.surface_taux = update(vs.surface_taux, at[:, :], f1 * vs.taux[:, :, 0] + f2 * vs.taux[:, :, 1])
    vs.surface_tauy = update(vs.surface_tauy, at[:, :], f1 * vs.tauy[:, :, 0] + f2 * vs.tauy[:, :, 1])

    if settings.enable_tke:
        vs.forc_tke_surface = update(
            vs.forc_tke_surface,
            at[1:-1, 1:-1],
            npx.sqrt(
                (0.5 * (vs.surface_taux[1:-1, 1:-1] + vs.surface_taux[:-2, 1:-1]) / settings.rho_0) ** 2
                + (0.5 * (vs.surface_tauy[1:-1, 1:-1] + vs.surface_tauy[1:-1, :-2]) / settings.rho_0) ** 2
            )** (3.0 / 2.0),
        )

    qqnet = f1 * vs.qnet[..., 0] + f2 * vs.qnet[..., 1]
    t_star_cur = f1_rest * vs.t_star[..., n1_rest] + f2_rest * vs.t_star[..., n2_rest]
    vs.forc_temp_surface = (qqnet + dQdT * (vs.temp[..., -1,vs.tau] - t_star_cur))  * vs.maskT[..., -1] / cp_0 / settings.rho_0
    evap_cur =  f1 * vs.evap[..., 0]+ f2 * vs.evap[..., 1]
    tp_cur =  f1 * vs.tp[..., 0] + f2 * vs.tp[..., 1]
    s_star_cur = f1_rest * vs.s_star[..., n1_rest] + f2_rest * vs.s_star[..., n2_rest]
    vs.forc_salt_surface = (Ks*(vs.salt[...,-1,vs.tau]-s_star_cur) -(tp_cur - evap_cur) * vs.salt[..., -1,vs.tau]) * vs.maskT[..., -1]

    # apply simple ice mask
    mask1 = vs.temp[:, :, -1, vs.tau] * vs.maskT[:, :, -1] > -1.8
    mask2 = vs.forc_temp_surface > 0
    ice = npx.logical_or(mask1, mask2)
    vs.forc_temp_surface *= ice
    vs.forc_salt_surface *= ice

    # solar radiation
    if settings.enable_tempsalt_sources:
        vs.temp_source = (
            (f1 * vs.qsol[..., 0, None] + f2 * vs.qsol[..., 1, None])
            * vs.divpen_shortwave[None, None, :]
            * ice[..., None]
            * vs.maskT[..., :]
            / cp_0
            / settings.rho_0
        )

    return KernelOutput(
        surface_taux=vs.surface_taux,
        surface_tauy=vs.surface_tauy, 
        temp_source=vs.temp_source,
        forc_tke_surface=vs.forc_tke_surface,
        forc_temp_surface=vs.forc_temp_surface,
        forc_salt_surface=vs.forc_salt_surface,
    )
