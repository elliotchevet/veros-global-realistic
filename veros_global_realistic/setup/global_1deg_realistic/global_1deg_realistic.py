import os
import h5netcdf
import xarray as xr
import veros_global_realistic

from veros import VerosSetup, tools, time, veros_routine, veros_kernel, KernelOutput
from veros.variables import Variable, allocate
from veros.core.operators import numpy as npx, update, at
from veros_global_realistic import forcing
from veros_global_realistic.paths import expanded_asset_manifest


BASE_PATH = os.path.dirname(os.path.realpath(__file__))
ASSET_MANIFEST = expanded_asset_manifest(os.path.join(BASE_PATH, "assets.json"))
DATA_FILES = tools.get_assets("global_1deg_realistic", ASSET_MANIFEST)
ERA5_VARIABLES = ("ewss", "nsss", "str", "slhf", "sshf", "ssr", "e", "tp")


class GlobalOneDegreeSetup(VerosSetup):
    """Global 1 degree model with 40 vertical levels.

    `Adapted from pyOM2 <https://wiki.zmaw.de/ifm/TO/pyOM2/1x1%20global%20model>`.

    """
    min_depth = 10.0
    max_depth = 5400.0
    __veros_plugins__ = (veros_global_realistic,)

    @veros_routine
    def set_parameter(self, state):
        """
        set main parameters
        """
        settings = state.settings

        settings.identifier = "global_1deg_realistic"
        settings.description = "Start from Glorys"
        settings.start_date = "1993-01-01"

        settings.nx = 360
        settings.ny = 160
        settings.nz = 40
        settings.dt_mom = 3600
        settings.dt_tracer = 3600
        settings.dt_forcing = 3
        settings.runlen = (30 * 360) * 24 * 3600

        settings.x_origin = 91.0
        settings.y_origin = -79.0

        settings.coord_degree = True
        settings.enable_cyclic_x = True

        settings.enable_hor_friction = True
        settings.A_h = 5e4
        settings.enable_hor_friction_cos_scaling = True
        settings.hor_friction_cosPower = 1
        settings.enable_quadratic_bottom_friction = False
        settings.r_quad_bot = 1e-3
        settings.enable_ray_friction = False
        settings.r_ray = 1e-5  # weak interior damping
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
        self._era5_year_slice = forcing.era5_year_slice(DATA_FILES["forcing_era5"], settings.start_date)
        self._forcing_year_days = forcing.year_days(self._era5_year_slice)
        state.dimensions["nmonths"] = 12
        if (self._forcing_year_days * 24) % settings.dt_forcing != 0:
            raise ValueError("dt_forcing must divide the hourly forcing year")
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
            forcing_time_shift=Variable(
                "forcing_time_shift", tuple(), "s", "Forcing clock offset from January 1", time_dependent=False
            ),
        )

    def _read_forcing(self, var):
        default_var = (
            "xt",
            "yt",
            "zt",
            "dz",
            "temperature",
            "sst",
            "salinity",
            "bathymetry",
            "swf",
            "dqdt",
            "tidal_energy",
            "wind_energy",
        )
        if var in default_var:
            with h5netcdf.File(DATA_FILES["forcing_default"], "r") as infile:
                var = infile.variables[var]
                return npx.asarray(var).T
        else:
            raise ValueError("Variable {} not in default_var list".format(var))

    def _read_init(self,var):
        with h5netcdf.File(DATA_FILES["ic_glorys12"], "r") as infile:
            var = infile.variables[var]
            var = npx.asarray(var).T
            return npx.where(npx.isnan(var), 0.0, var)

    def _read_restoring(self,var):
        with h5netcdf.File(DATA_FILES["restoring_glorys12"], "r") as infile:
            var = infile.variables[var]
            return npx.asarray(var).T

    def _read_mask(self,var):
        with h5netcdf.File(DATA_FILES["mask"], "r") as infile:
            var = infile.variables[var]
            return npx.asarray(var)


    def _read_ERA5(self,var):
        if not hasattr(self, "_era5_ds"):
            self._era5_ds = xr.open_dataset(
                DATA_FILES["forcing_era5"],
                engine="h5netcdf",
                chunks={"time": 1},
            ).isel(time=self._era5_year_slice)
        return self._era5_ds[var]

    def get_ERA5_window(self, n1, n2):
        window_key = (int(n1), int(n2))
        if getattr(self, "_era5_window_key", None) == window_key:
            return self._era5_window

        print("forcing update n°{}".format(n1))
        self._era5_window_key = window_key
        self._era5_window = {
            var: npx.asarray(self._read_ERA5(var).isel(time=[n1, n2])).T
            for var in ERA5_VARIABLES
        }
        return self._era5_window

    def write_mask_file(self, state, output_path="mask.nc"):
        import numpy as np

        vs = state.variables
        ds = xr.Dataset(
            {
                "maskT": (("xt", "yt", "zt"), np.asarray(vs.maskT[2:-2, 2:-2])),
                "maskU": (("xu", "yt", "zt"), np.asarray(vs.maskU[2:-2, 2:-2])),
                "maskV": (("xt", "yu", "zt"), np.asarray(vs.maskV[2:-2, 2:-2])),
                "maskW": (("xt", "yt", "zw"), np.asarray(vs.maskW[2:-2, 2:-2])),
                "kbot": (("xt", "yt"), np.asarray(vs.kbot[2:-2, 2:-2])),
            },
            coords={
                "xt": np.asarray(vs.xt[2:-2]),
                "xu": np.asarray(vs.xu[2:-2]),
                "yt": np.asarray(vs.yt[2:-2]),
                "yu": np.asarray(vs.yu[2:-2]),
                "zt": np.asarray(vs.zt),
                "zw": np.asarray(vs.zw),
            },
        )
        ds.to_netcdf(output_path)

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

    @veros_routine(dist_safe=False, local_variables=["kbot"])
    def set_topography(self, state):
        vs = state.variables
        kbot_data = self._read_mask("kbot")
        vs.kbot = update(vs.kbot, at[2:-2, 2:-2], kbot_data)

    @veros_routine(
        dist_safe=False,
        local_variables=[
            "t_star",
            "s_star",
            "divpen_shortwave",
            "temp",
            "salt",
            "ssh",
            "u",
            "v",
            "forcing_time_shift",
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

        if settings.write_mask:
            self.write_mask_file(state)

        # initial conditions

        temp_raw = self._read_init("thetao")[:,:,:,0]
        vs.temp = update(vs.temp, at[2:-2, 2:-2, :, 0], temp_raw * vs.maskT[2:-2, 2:-2, :])
        vs.temp = update(vs.temp, at[2:-2, 2:-2, :, 1], temp_raw * vs.maskT[2:-2, 2:-2, :])

        salt_raw = self._read_init("so")[:,:,:,0]
        vs.salt = update(vs.salt, at[2:-2, 2:-2, :, 0], salt_raw * vs.maskT[2:-2, 2:-2, :])
        vs.salt = update(vs.salt, at[2:-2, 2:-2, :, 1], salt_raw * vs.maskT[2:-2, 2:-2, :])

        ssh_raw = self._read_init("zos")[...,0]
        vs.ssh = update(vs.ssh, at[2:-2, 2:-2], ssh_raw * vs.maskT[2:-2, 2:-2, -1])

        u_raw = self._read_init("uo")[...,0]
        vs.u = update(vs.u, at[2:-2, 2:-2, :, 0], u_raw * vs.maskU[2:-2, 2:-2, :])
        vs.u = update(vs.u, at[2:-2, 2:-2, :, 1], u_raw * vs.maskU[2:-2, 2:-2, :])

        v_raw = self._read_init("vo")[...,0]
        vs.v = update(vs.v, at[2:-2, 2:-2, :, 0], v_raw * vs.maskV[2:-2, 2:-2, :])
        vs.v = update(vs.v, at[2:-2, 2:-2, :, 1], v_raw * vs.maskV[2:-2, 2:-2, :])

        t_star = self._read_restoring("thetao")
        vs.t_star = update(vs.t_star, at[2:-2, 2:-2,:], t_star * vs.maskT[2:-2, 2:-2, -1, npx.newaxis])
        s_star = self._read_restoring("so")
        vs.s_star = update(vs.s_star, at[2:-2, 2:-2,:], s_star * vs.maskT[2:-2, 2:-2, -1, npx.newaxis])
        vs.forcing_time_shift = update(
            vs.forcing_time_shift,
            at[...],
            forcing.forcing_time_shift_seconds(
                settings.start_date,
                time.convert_time(1.0, "years", "seconds"),
                self._forcing_year_days,
            ),
        )

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
        n_forcing = self._forcing_year_days * 24 // settings.dt_forcing
        forcing_time = vs.time + vs.forcing_time_shift
        (n1, _), (n2, _) = tools.get_periodic_interval(
            forcing_time,
            year_in_seconds,
            year_in_seconds / n_forcing,
            n_forcing,
        )
        n1 = n1 * settings.dt_forcing
        n2 = n2 * settings.dt_forcing
        era5 = self.get_ERA5_window(n1, n2)
        taux_data = era5["ewss"] / 3600
        vs.taux = update(vs.taux, at[2:-2, 2:-2, :], taux_data)
        tauy_data = era5["nsss"] / 3600
        vs.tauy = update(vs.tauy, at[2:-2, 2:-2, :], tauy_data)
        Qir = era5["str"] / 3600
        Qlat = era5["slhf"] / 3600
        Qh = era5["sshf"] / 3600
        Qnet = Qir + Qlat + Qh
        vs.qnet = update(vs.qnet, at[2:-2, 2:-2, :], Qnet * vs.maskT[2:-2, 2:-2, -1, npx.newaxis])
        Qsol = era5["ssr"] / 3600
        vs.qsol = update(vs.qsol, at[2:-2, 2:-2, :], Qsol* vs.maskT[2:-2, 2:-2, -1, npx.newaxis])

        # Evaporation
        evap_data = era5["e"] / 3600
        vs.evap = update(vs.evap, at[2:-2, 2:-2, :], evap_data * vs.maskT[2:-2, 2:-2, -1, npx.newaxis])

        # Total precipitation
        tp_data = era5["tp"] / 3600
        vs.tp = update(vs.tp, at[2:-2, 2:-2, :], tp_data * vs.maskT[2:-2, 2:-2, -1, npx.newaxis])

        vs.update(set_forcing_kernel(state, n_forcing))

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
            "qsol",
            "qnet",
            "taux",
            "tauy",
            "Nsqr",
            "Hd",
            "K_diss_v",
            "P_diss_v",
            "kbot",
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
        state.diagnostics["cfl_monitor"].output_frequency = 10 * 86400
        state.diagnostics["snapshot"].output_frequency = 360 * 86400 / 12
        state.diagnostics["overturning"].output_frequency = 360 * 86400
        state.diagnostics["overturning"].sampling_frequency = 5 * 86400
        state.diagnostics["energy"].output_frequency = 86400
        state.diagnostics["energy"].sampling_frequency = 3600 * 6
        state.diagnostics["averages"].output_frequency = 360 * 86400 / 12
        state.diagnostics["averages"].sampling_frequency = 86400/2

    @veros_routine
    def after_timestep(self, state):
        pass



@veros_kernel
def set_forcing_kernel(state, n_forcing):
    vs = state.variables
    settings = state.settings

    cp_0 = 3991.86795711963  # J/kg /K
    dQdT = -40 # W/m2 /K (2 months for a 50m mld)
    Ks = -6.37e-6 # piston velocity in m/s (3 months for a 50m mld)

    year_in_seconds = time.convert_time(1.0, "years", "seconds")
    forcing_time = vs.time + vs.forcing_time_shift
    (n1, f1), (n2, f2) = tools.get_periodic_interval(
        forcing_time,
        year_in_seconds,
        year_in_seconds / n_forcing,
        n_forcing,
    )

    (n1_rest, f1_rest), (n2_rest, f2_rest) = tools.get_periodic_interval(
        forcing_time,
        year_in_seconds,
        year_in_seconds / 12.0,
        12,
    )
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

