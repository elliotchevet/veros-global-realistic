import copernicusmarine

dataset_id = "cmems_mod_glo_phy_my_0.083deg_P1M-m"

copernicusmarine.get(
    dataset_id=dataset_id,
    filter=f"mercatorglorys12v1_gl12_mean_199301.nc",
    output_directory="../data/GLORYS12"
)
