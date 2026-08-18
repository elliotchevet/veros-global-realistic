import argparse
import os
from calendar import monthrange
from datetime import datetime
from pathlib import Path


BASE_DIR = Path(__file__).resolve().parents[1]
DATA_DIR = Path(os.environ.get("VEROS_GLOBAL_REALISTIC_DATA_DIR", BASE_DIR / "data")).expanduser()
GLORYS12_DIR = Path(os.environ.get("GLORYS12_DIR", DATA_DIR / "GLORYS12")).expanduser()

DEFAULT_IC_DATASET_ID = "cmems_mod_glo_phy_my_0.083deg_P1M-m"
DEFAULT_RESTORING_DATASET_ID = "cmems_mod_glo_phy_my_0.083deg-climatology_P1M-m"
DEFAULT_SURFACE_DEPTH = 0.49402499198913574

IC_VARIABLES = ["so", "thetao", "zos", "mlotst", "uo", "vo"]
RESTORING_VARIABLES = ["so", "thetao", "zos", "mlotst"]


def parse_yyyymm(value):
    try:
        return datetime.strptime(value, "%Y%m")
    except ValueError as exc:
        raise argparse.ArgumentTypeError("expected YYYYMM, for example 200401") from exc


def format_month_start(date):
    return f"{date.year:04d}-{date.month:02d}-01T00:00:00"


def format_month_end(date):
    last_day = monthrange(date.year, date.month)[1]
    return f"{date.year:04d}-{date.month:02d}-{last_day:02d}T23:59:59"


def build_subset_request(args):
    output_root = args.output_root.expanduser()

    if args.mode == "IC":
        month = args.yyyymm
        if month is None:
            raise SystemExit("--yyyymm is required when --mode IC")

        output_directory = output_root / "Initial_States"
        output_filename = f"GLORYS12_1deg_{month:%Y%m}.nc"
        dataset_id = args.ic_dataset_id
        variables = args.ic_variables
        start_datetime = format_month_start(month)
        end_datetime = format_month_end(month)
        minimum_depth = args.minimum_depth
        maximum_depth = args.maximum_depth
    else:
        output_directory = output_root / "Restoring"
        output_filename = "GLORYS12_1deg_Restoring.nc"
        dataset_id = args.restoring_dataset_id
        variables = args.restoring_variables
        start_datetime = f"{args.restoring_year:04d}-01-01T00:00:00"
        end_datetime = f"{args.restoring_year:04d}-12-01T00:00:00"
        minimum_depth = args.surface_depth
        maximum_depth = args.surface_depth

    output_directory.mkdir(parents=True, exist_ok=True)

    request = {
        "dataset_id": dataset_id,
        "variables": variables,
        "minimum_longitude": args.minimum_longitude,
        "maximum_longitude": args.maximum_longitude,
        "minimum_latitude": args.minimum_latitude,
        "maximum_latitude": args.maximum_latitude,
        "start_datetime": start_datetime,
        "end_datetime": end_datetime,
        "output_directory": str(output_directory),
        "output_filename": output_filename,
    }

    if minimum_depth is not None:
        request["minimum_depth"] = minimum_depth
    if maximum_depth is not None:
        request["maximum_depth"] = maximum_depth

    return request


def run_subset(request):
    import copernicusmarine

    return copernicusmarine.subset(**request)


def variable_list(value):
    return [item.strip() for item in value.split(",") if item.strip()]


def parse_args():
    parser = argparse.ArgumentParser(
        description="Download GLORYS12 initial-state or restoring data for the global realistic setup."
    )
    parser.add_argument(
        "--mode",
        choices=["IC", "Restoring"],
        required=True,
        help="Download a monthly initial state (IC) or surface monthly climatology for restoring.",
    )
    parser.add_argument(
        "--yyyymm",
        type=parse_yyyymm,
        help="Initial-state month to download, formatted as YYYYMM. Required with --mode IC.",
    )
    parser.add_argument(
        "--restoring-year",
        type=int,
        default=2004,
        help="Reference climatology year to request with --mode Restoring.",
    )
    parser.add_argument(
        "--output-root",
        type=Path,
        default=GLORYS12_DIR,
        help="Root GLORYS12 directory. Subfolders Initial_States and Restoring are created below it.",
    )
    parser.add_argument("--ic-dataset-id", default=DEFAULT_IC_DATASET_ID)
    parser.add_argument("--restoring-dataset-id", default=DEFAULT_RESTORING_DATASET_ID)
    parser.add_argument("--ic-variables", type=variable_list, default=IC_VARIABLES)
    parser.add_argument("--restoring-variables", type=variable_list, default=RESTORING_VARIABLES)
    parser.add_argument("--minimum-longitude", type=float, default=-180)
    parser.add_argument("--maximum-longitude", type=float, default=179.9166717529297)
    parser.add_argument("--minimum-latitude", type=float, default=-80)
    parser.add_argument("--maximum-latitude", type=float, default=90)
    parser.add_argument(
        "--surface-depth",
        type=float,
        default=DEFAULT_SURFACE_DEPTH,
        help="Depth used for surface-only restoring requests.",
    )
    parser.add_argument(
        "--minimum-depth",
        type=float,
        default=None,
        help="Optional minimum depth for IC requests. Omit to download all depths.",
    )
    parser.add_argument(
        "--maximum-depth",
        type=float,
        default=None,
        help="Optional maximum depth for IC requests. Omit to download all depths.",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    request = build_subset_request(args)
    run_subset(request)


if __name__ == "__main__":
    main()

