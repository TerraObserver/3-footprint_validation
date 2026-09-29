from pathlib import Path

import numpy as np
import pandas as pd
import rasterio
from affine import Affine
from pyproj import CRS, Transformer

import calc_footprint_FFP_climatology as ffp


# ============================================================
# USER CONFIGURATION
# ============================================================

# Keep this script, the Excel input, and calc_footprint_FFP_climatology.py
# in the same folder. This makes the workflow portable to another computer.
BASE_DIR = Path(__file__).resolve().parent
INPUT_FILE = BASE_DIR / "tower_footprint_input.xlsx"
OUTPUT_BASE = BASE_DIR / "footprints_output"

# Edit the sites here. The number of sites is inferred automatically from
# this dictionary, so no separate NUM_SITES variable is needed.
#
# origin_dist : half-width of footprint domain [m]
# latitude/or longitude : EC tower coordinates [decimal degrees, WGS84]
# zm : measurement height above displacement height [m]
# z0 : roughness length [m]
SITE_CONFIG = {
    "OLA": {
        "origin_dist": 700.0,
        "latitude": 36.8013232,
        "longitude": -120.209413,
        "zm": 6.445,
        "z0": 0.45,
    },
    "VAC": {
        "origin_dist": 500.0,
        "latitude": 38.3000364,
        "longitude": -121.9104940,
        "zm": 7.65,
        "z0": 0.35,
    },
    "WWF": {
        "origin_dist": 400.0,
        "latitude": 38.6672946,
        "longitude": -121.8919682,
        "zm": 7.35,
        "z0": 0.35,
    },
}

# Constants used for all sites
H_ATM_BOUND = 1500.0       # atmospheric boundary-layer height [m]
GRID_RESOLUTION = 6.5      # footprint grid resolution [m]

# Daily footprint time window. END_HOUR is INCLUDED.
START_HOUR = 6
END_HOUR = 16

# Minimum number of successfully calculated hourly footprints required
# to create a daily weighted footprint.
MIN_VALID_HOURS = 5

# Retain the cells containing this fraction of each hourly footprint mass.
FOOTPRINT_CUTOFF = 0.90

# False = preserve existing outputs and skip them on reruns.
# IMPORTANT: set True once if you need to regenerate existing daily rasters after
# changing the weighting method, then return it to False.
OVERWRITE = False


# ============================================================
# INTERNAL SETTINGS / REQUIRED INPUTS
# ============================================================

REQUIRED_COLUMNS = {
    "date",
    "SITE_ID",
    "wind_dir",
    "u_star",
    "L",
    "sigma_v",
    "ETo",
}

NUM_SITES = len(SITE_CONFIG)
HOURS = list(range(START_HOUR, END_HOUR + 1))


# ============================================================
# HELPER FUNCTIONS
# ============================================================

def validate_configuration():
    """Validate the editable configuration before processing starts."""
    if not INPUT_FILE.is_file():
        raise FileNotFoundError(f"Input file not found:\n{INPUT_FILE}")

    if not (0 <= START_HOUR <= 23 and 0 <= END_HOUR <= 23):
        raise ValueError("START_HOUR and END_HOUR must be between 0 and 23.")

    if START_HOUR > END_HOUR:
        raise ValueError("START_HOUR must be <= END_HOUR.")

    if GRID_RESOLUTION <= 0:
        raise ValueError("GRID_RESOLUTION must be > 0.")

    if H_ATM_BOUND <= 10:
        raise ValueError("H_ATM_BOUND must be > 10 m for FFP.")

    if not (0 < FOOTPRINT_CUTOFF <= 1):
        raise ValueError("FOOTPRINT_CUTOFF must be in (0, 1].")

    for site, cfg in SITE_CONFIG.items():
        required = {"origin_dist", "latitude", "longitude", "zm", "z0"}
        missing = required.difference(cfg)
        if missing:
            raise ValueError(f"{site}: missing SITE_CONFIG fields: {sorted(missing)}")

        if cfg["origin_dist"] <= 0:
            raise ValueError(f"{site}: origin_dist must be > 0.")
        if not (-90 <= cfg["latitude"] <= 90):
            raise ValueError(f"{site}: invalid latitude.")
        if not (-180 <= cfg["longitude"] <= 180):
            raise ValueError(f"{site}: invalid longitude.")
        if cfg["zm"] <= 0 or cfg["z0"] <= 0:
            raise ValueError(f"{site}: zm and z0 must be > 0.")


def read_input():
    """Read the combined tower footprint input and standardize key columns."""
    df = pd.read_excel(
        INPUT_FILE,
        na_values=["NA", "N/A", "NaN", "nan", "—", "-", "", "null", "None"],
    )

    missing = REQUIRED_COLUMNS.difference(df.columns)
    if missing:
        raise ValueError(
            "Input workbook is missing required columns: "
            + ", ".join(sorted(missing))
        )

    df["date"] = pd.to_datetime(df["date"], errors="coerce")
    df["SITE_ID"] = df["SITE_ID"].astype(str).str.strip().str.upper()

    numeric_cols = ["wind_dir", "u_star", "L", "sigma_v", "ETo"]
    for col in numeric_cols:
        df[col] = pd.to_numeric(df[col], errors="coerce")

    df = df[df["date"].notna()].copy()
    df = df.sort_values(["SITE_ID", "date"]).reset_index(drop=True)

    input_sites = set(df["SITE_ID"].dropna().unique())
    configured_sites = set(SITE_CONFIG)

    missing_sites = configured_sites - input_sites
    extra_sites = input_sites - configured_sites

    if missing_sites:
        print(f"WARNING: configured sites not found in input: {sorted(missing_sites)}")
    if extra_sites:
        print(f"WARNING: input sites not configured and will be ignored: {sorted(extra_sites)}")

    return df


def utm_crs_from_latlon(latitude, longitude):
    """Return local UTM CRS for a WGS84 latitude/longitude."""
    zone = int((longitude + 180) // 6) + 1
    epsg = (32600 if latitude >= 0 else 32700) + zone
    return CRS.from_epsg(epsg)


def tower_xy(latitude, longitude, crs):
    """Project tower WGS84 coordinates to the selected UTM CRS."""
    transformer = Transformer.from_crs("EPSG:4326", crs, always_xy=True)
    return transformer.transform(longitude, latitude)


def snap_station_to_grid(station_x, station_y, dx):
    """
    Preserve the grid-centering logic used in the previous footprint workflow.
    This slightly shifts the computational grid origin so outputs align on a
    consistent half-cell lattice without changing the EC tower metadata.
    """
    half = dx / 2.0
    quarter = dx / 4.0

    rx = station_x % half
    if rx > quarter:
        station_x += half - rx
        if (station_x / half) % 2 == 0:
            station_x -= half
    else:
        station_x -= rx
        if (station_x / half) % 2 == 0:
            station_x += half

    ry = station_y % half
    if ry > quarter:
        station_y += half - ry
        if (station_y / half) % 2 == 0:
            station_y -= half
    else:
        station_y -= ry
        if (station_y / half) % 2 == 0:
            station_y += half

    return station_x, station_y


def mask_fp_cutoff(f_array, cutoff=0.90):
    """Keep cells containing `cutoff` fraction of footprint mass."""
    out = np.asarray(f_array, dtype=float).copy()
    out[~np.isfinite(out)] = 0.0
    out[out < 0] = 0.0

    total = out.sum()
    if total <= 0:
        return np.zeros_like(out)

    flat = out.ravel()
    order = np.argsort(flat)[::-1]
    sorted_vals = flat[order]
    cumulative = np.cumsum(sorted_vals) / total

    idx = np.searchsorted(cumulative, cutoff, side="left")
    idx = min(idx, len(sorted_vals) - 1)
    threshold = sorted_vals[idx]

    out[out < threshold] = 0.0
    return out


def eto_weights(successful_hours):
    """
    Weight each valid hourly footprint in direct proportion to hourly ETo:

        w_i = ETo_i / sum(ETo)

    Non-finite or negative ETo values receive zero weight. If no positive ETo
    remains, fall back to equal weights across the valid hourly footprints.
    Returns weights that sum to 1.
    """
    eto = np.asarray([item["ETo"] for item in successful_hours], dtype=float)

    if len(eto) == 0:
        return np.array([], dtype=float)

    # ETo is a non-negative atmospheric-demand quantity. Invalid or negative
    # values should not contribute to the daily footprint weighting.
    eto[~np.isfinite(eto)] = 0.0
    eto[eto < 0.0] = 0.0

    total = eto.sum()
    if total <= 0:
        return np.repeat(1.0 / len(eto), len(eto))

    return eto / total


def calculate_hourly_footprint(row, cfg, station_x, station_y):
    """Calculate one georeferenced hourly footprint array."""
    temp_ffp = ffp.FFP_climatology(
        domain=[
            -cfg["origin_dist"],
            cfg["origin_dist"],
            -cfg["origin_dist"],
            cfg["origin_dist"],
        ],
        dx=GRID_RESOLUTION,
        dy=GRID_RESOLUTION,
        zm=float(cfg["zm"]),
        h=float(H_ATM_BOUND),
        rs=None,
        z0=float(cfg["z0"]),
        ol=float(row["L"]),
        sigmav=float(row["sigma_v"]),
        ustar=float(row["u_star"]),
        umean=None,  # z0 is known and is therefore used by FFP
        wind_dir=float(row["wind_dir"]),
        crop=0,
        fig=0,
        verbosity=0,
    )

    f_2d = np.asarray(temp_ffp["fclim_2d"], dtype=float)
    x_2d = np.asarray(temp_ffp["x_2d"], dtype=float) + station_x
    y_2d = np.asarray(temp_ffp["y_2d"], dtype=float) + station_y

    # Convert footprint density [m-2] to cell mass, as in the original workflow.
    f_2d = f_2d * (GRID_RESOLUTION ** 2)
    f_2d = mask_fp_cutoff(f_2d, FOOTPRINT_CUTOFF)

    # Normalize each successful hourly footprint to unit mass before ETo weighting.
    total = np.nansum(f_2d)
    if not np.isfinite(total) or total <= 0:
        raise ValueError("hourly footprint has zero/invalid mass")
    f_2d = f_2d / total

    # Preserve the orientation/georeferencing convention from the prior workflow.
    transform = Affine(
        GRID_RESOLUTION,
        0.0,
        float(x_2d[0, 0]),
        0.0,
        GRID_RESOLUTION,
        float(y_2d[0, 0]),
    )

    return f_2d, transform


def write_hourly_stack(path, crs, transform, arrays_by_hour, metadata_by_hour):
    """Write one multiband raster with a fixed band for every configured hour."""
    sample = next(iter(arrays_by_hour.values()))

    with rasterio.open(
        path,
        "w",
        driver="GTiff",
        dtype="float64",
        count=len(HOURS),
        height=sample.shape[0],
        width=sample.shape[1],
        transform=transform,
        crs=crs,
        nodata=0.0,
    ) as dst:
        for band_index, hour in enumerate(HOURS, start=1):
            arr = arrays_by_hour.get(hour)
            if arr is None:
                arr = np.zeros_like(sample)

            dst.write(arr.astype(np.float64), band_index)

            tags = metadata_by_hour.get(
                hour,
                {
                    "hour": f"{hour:02d}:00",
                    "status": "MISSING_OR_INVALID",
                },
            )
            dst.update_tags(band_index, **{k: str(v) for k, v in tags.items()})


def write_daily_raster(path, crs, transform, daily_array):
    """Write the single-band daily ETo-weighted footprint."""
    with rasterio.open(
        path,
        "w",
        driver="GTiff",
        dtype="float64",
        count=1,
        height=daily_array.shape[0],
        width=daily_array.shape[1],
        transform=transform,
        crs=crs,
        nodata=0.0,
    ) as dst:
        dst.write(daily_array.astype(np.float64), 1)
        dst.update_tags(
            1,
            weighting="proportional hourly ETo: ETo_i / sum(ETo)",
            footprint_sum=float(np.nansum(daily_array)),
            start_hour=START_HOUR,
            end_hour=END_HOUR,
        )


# ============================================================
# SITE PROCESSING
# ============================================================

def process_site(all_data, site, cfg):
    """Process every flight day for one site into its own output folder."""
    site_df = all_data[all_data["SITE_ID"] == site].copy()

    if site_df.empty:
        print(f"\n{site}: no rows found; skipping site.")
        return []

    site_output = OUTPUT_BASE / site
    hourly_output = site_output / "hourly"
    daily_output = site_output / "daily_weighted"
    hourly_output.mkdir(parents=True, exist_ok=True)
    daily_output.mkdir(parents=True, exist_ok=True)

    crs = utm_crs_from_latlon(cfg["latitude"], cfg["longitude"])
    station_x, station_y = tower_xy(cfg["latitude"], cfg["longitude"], crs)
    station_x, station_y = snap_station_to_grid(
        station_x, station_y, GRID_RESOLUTION
    )

    print("\n====================================")
    print(f"SITE: {site}")
    print("====================================")
    print(f"Rows: {len(site_df)}")
    print(f"Days: {site_df['date'].dt.date.nunique()}")
    print(f"CRS: {crs.to_string()}")
    print(f"origin_dist: {cfg['origin_dist']} m")
    print(f"hours: {START_HOUR:02d}:00-{END_HOUR:02d}:00 inclusive")

    summaries = []

    # Unique days only. This fixes the previous workflow, which revisited each
    # date once per hourly row.
    dates = sorted(site_df["date"].dt.date.unique())

    for day in dates:
        day_ts = pd.Timestamp(day)
        day_label = day_ts.strftime("%Y_%m_%d")

        hourly_path = hourly_output / f"{site}_hourly_{day_label}.tif"
        daily_path = daily_output / f"{site}_weighted_{day_label}.tif"

        if daily_path.exists() and not OVERWRITE:
            print(f"{site} {day}: daily output already exists; skipping.")
            summaries.append(
                {
                    "SITE_ID": site,
                    "date": day,
                    "status": "SKIPPED_EXISTING",
                    "valid_hours": np.nan,
                    "hourly_raster": str(hourly_path),
                    "daily_raster": str(daily_path),
                }
            )
            continue



        temp_df = site_df[site_df["date"].dt.date == day].copy()
        temp_df = temp_df[
            (temp_df["date"].dt.hour >= START_HOUR)
            & (temp_df["date"].dt.hour <= END_HOUR)
        ].copy()

        # VAC 2023 EC records occur at xx:30 instead of xx:00.
        # Keep the original tower timestamps unchanged and select the appropriate minute.
        if site == "VAC" and day.year == 2023:
            minute_target = 30
        else:
            minute_target = 0

        arrays_by_hour = {}
        metadata_by_hour = {}
        successful_hours = []
        transform = None

        for hour in HOURS:
            timestamp_label = f"{hour:02d}:{minute_target:02d}"

            rows = temp_df[
                (temp_df["date"].dt.hour == hour)
                & (temp_df["date"].dt.minute == minute_target)
            ]

            if rows.empty:
                print(f"  {day} {timestamp_label} -> missing")
                continue

            row = rows.iloc[0]
            required_values = row[["wind_dir", "u_star", "L", "sigma_v", "ETo"]]

            if required_values.isna().any():
                print(f"  {day} {timestamp_label} -> required EC data missing")
                continue

            try:
                arr, hour_transform = calculate_hourly_footprint(
                    row, cfg, station_x, station_y
                )
            except Exception as exc:
                print(f"  {day} {timestamp_label} -> FFP rejected/failed ({exc})")
                continue

            if transform is None:
                transform = hour_transform

            arrays_by_hour[hour] = arr
            metadata_by_hour[hour] = {
                "hour": timestamp_label,
                "status": "VALID",
                "wind_dir": float(row["wind_dir"]),
                "u_star": float(row["u_star"]),
                "L": float(row["L"]),
                "sigma_v": float(row["sigma_v"]),
                "ETo": float(row["ETo"]),
                "footprint_sum": float(np.nansum(arr)),
            }
            successful_hours.append(
                {
                    "hour": hour,
                    "timestamp": timestamp_label,
                    "array": arr,
                    "ETo": float(row["ETo"]),
                }
            )
            print(f"  {day} {timestamp_label} -> OK")

        valid_count = len(successful_hours)

        if valid_count == 0:
            print(f"  {day}: no valid hourly footprints; daily raster not created.")
            summaries.append(
                {
                    "SITE_ID": site,
                    "date": day,
                    "status": "NO_VALID_HOURS",
                    "valid_hours": 0,
                    "hourly_raster": "",
                    "daily_raster": "",
                }
            )
            continue

        # Hourly stack is useful for QA/QC even when the day does not reach the
        # minimum number of valid hours required for a daily footprint.
        if OVERWRITE or not hourly_path.exists():
            write_hourly_stack(
                hourly_path,
                crs,
                transform,
                arrays_by_hour,
                metadata_by_hour,
            )

        if valid_count < MIN_VALID_HOURS:
            print(
                f"  {day}: only {valid_count} valid hours "
                f"(< {MIN_VALID_HOURS}); daily raster not created."
            )
            summaries.append(
                {
                    "SITE_ID": site,
                    "date": day,
                    "status": "INSUFFICIENT_VALID_HOURS",
                    "valid_hours": valid_count,
                    "hourly_raster": str(hourly_path),
                    "daily_raster": "",
                }
            )
            continue

        weights = eto_weights(successful_hours)
        daily_array = np.zeros_like(successful_hours[0]["array"], dtype=float)

        for item, weight in zip(successful_hours, weights):
            daily_array += item["array"] * weight

        daily_sum = np.nansum(daily_array)
        if not np.isfinite(daily_sum) or daily_sum <= 0:
            print(f"  {day}: invalid daily footprint sum; daily raster not created.")
            status = "INVALID_DAILY_SUM"
            daily_file_value = ""
        else:
            daily_array = daily_array / daily_sum
            write_daily_raster(daily_path, crs, transform, daily_array)
            print(
                f"  {day}: daily weighted footprint saved "
                f"({valid_count} valid hours, sum={daily_array.sum():.6f})"
            )
            status = "DAILY_CREATED"
            daily_file_value = str(daily_path)

        summaries.append(
            {
                "SITE_ID": site,
                "date": day,
                "status": status,
                "valid_hours": valid_count,
                "hourly_raster": str(hourly_path),
                "daily_raster": daily_file_value,
            }
        )

    return summaries


# ============================================================
# MAIN
# ============================================================

def main():
    validate_configuration()
    OUTPUT_BASE.mkdir(parents=True, exist_ok=True)

    print("====================================")
    print("FOOTPRINT PROCESSING")
    print("====================================")
    print(f"Input: {INPUT_FILE}")
    print(f"Output: {OUTPUT_BASE}")
    print(f"Configured sites ({NUM_SITES}): {', '.join(SITE_CONFIG.keys())}")

    data = read_input()

    all_summaries = []
    for site, cfg in SITE_CONFIG.items():
        all_summaries.extend(process_site(data, site, cfg))

    summary_df = pd.DataFrame(all_summaries)
    summary_path = OUTPUT_BASE / "footprint_processing_summary.csv"
    summary_df.to_csv(summary_path, index=False)

    print("\n====================================")
    print("PROCESSING COMPLETE")
    print("====================================")
    print(f"Summary: {summary_path}")

    if not summary_df.empty:
        print("\nStatus counts:")
        print(summary_df["status"].value_counts(dropna=False))


if __name__ == "__main__":
    main()
