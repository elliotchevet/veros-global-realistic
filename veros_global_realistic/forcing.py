import re
from datetime import datetime

import numpy as np


def forcing_time_shift_seconds(start_date, model_year_seconds, forcing_year_days):
    date = datetime.strptime(start_date, "%Y-%m-%d")
    year_start = datetime(date.year, 1, 1)
    day_of_year = (date - year_start).days
    if day_of_year >= forcing_year_days:
        raise ValueError(f"start_date {start_date} is outside the {forcing_year_days}-day forcing year")

    return day_of_year / forcing_year_days * model_year_seconds


def year_days(year_slice):
    return (year_slice.stop - year_slice.start) // 24


def datetime_to_file_time(date, units):
    if isinstance(units, bytes):
        units = units.decode()

    match = re.match(r"^(seconds|second|hours|hour|days|day) since (.+)$", units)
    if match is None:
        raise ValueError(f"Unsupported ERA5 time units: {units}")

    unit, ref_date = match.groups()
    ref_date = ref_date.strip().replace("T", " ").rstrip("Z")
    ref_date = ref_date.split("+")[0].split(".")[0].strip()
    for date_format in ("%Y-%m-%d %H:%M:%S", "%Y-%m-%d %H:%M", "%Y-%m-%d"):
        try:
            ref = datetime.strptime(ref_date, date_format)
            break
        except ValueError:
            pass
    else:
        raise ValueError(f"Unsupported ERA5 reference time: {ref_date}")

    seconds = (date - ref).total_seconds()
    if unit.startswith("second"):
        return seconds
    if unit.startswith("hour"):
        return seconds / 3600
    return seconds / 86400


def era5_year_slice(forcing_file, start_date):
    import h5netcdf

    date = datetime.strptime(start_date, "%Y-%m-%d")
    year_start = datetime(date.year, 1, 1)
    next_year_start = datetime(date.year + 1, 1, 1)

    with h5netcdf.File(forcing_file, "r") as f:
        time_var = f.variables["time"]
        time_values = np.asarray(time_var[:])
        units = time_var.attrs.get("units")

    if units is None:
        raise ValueError("ERA5 forcing file time variable must define units")

    year_start_value = datetime_to_file_time(year_start, units)
    next_year_start_value = datetime_to_file_time(next_year_start, units)
    year_start_index = int(np.searchsorted(time_values, year_start_value))
    next_year_start_index = int(np.searchsorted(time_values, next_year_start_value))

    if year_start_index >= time_values.size or not np.isclose(time_values[year_start_index], year_start_value):
        raise ValueError(f"ERA5 forcing file does not contain {year_start:%Y}")

    year_end_index = min(next_year_start_index, time_values.size)
    year_hours = year_end_index - year_start_index
    if year_hours <= 0 or year_hours % 24 != 0:
        raise ValueError("ERA5 forcing year must contain complete hourly days")
    if year_hours // 24 not in (365, 366):
        raise ValueError("ERA5 forcing year must contain 365 or 366 days")

    return slice(year_start_index, year_end_index)

