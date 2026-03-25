import cdsapi

dataset = "reanalysis-era5-single-levels"
for year in range(2009,2025):
    request = {
        "product_type": ["reanalysis"],
        "year": [str(year)],
        "month": [
            "01", "02", "03",
            "04", "05", "06",
            "07", "08", "09",
            "10", "11", "12"
        ],
        "day": [
            "01", "02", "03",
            "04", "05", "06",
            "07", "08", "09",
            "10", "11", "12",
            "13", "14", "15",
            "16", "17", "18",
            "19", "20", "21",
            "22", "23", "24",
            "25", "26", "27",
            "28", "29", "30",
            "31"
        ],
        "time": [
            "00:00", "01:00", "02:00",
            "03:00", "04:00", "05:00",
            "06:00", "07:00", "08:00",
            "09:00", "10:00", "11:00",
            "12:00", "13:00", "14:00",
            "15:00", "16:00", "17:00",
            "18:00", "19:00", "20:00",
            "21:00", "22:00", "23:00"
        ],
        "data_format": "grib",
        "download_format": "unarchived",
        "variable": [
            "Sea surface temperature",
            "total_precipitation",
            "eastward_turbulent_surface_stress",
            "northward_turbulent_surface_stress",
            "evaporation",
            "surface_net_solar_radiation",
            "surface_net_thermal_radiation",
            "surface_sensible_heat_flux",
            "surface_latent_heat_flux"
        ],
        "grid": '1.0/1.0'
    }
    target = "/Odyssey/private/e25cheve/data/ERA5_1deg_{}_Flux.grib".format(year)
    client = cdsapi.Client()
    client.retrieve(dataset, request,target)

