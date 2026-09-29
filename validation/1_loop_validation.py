from pathlib import Path
import re

import numpy as np
import pandas as pd
import rasterio
from rasterio.warp import reproject, Resampling


# ============================================================
# USER INPUTS
# ============================================================
# ============================================================
# PATHS
# ============================================================

SCRIPT_DIR = Path(__file__).resolve().parent

# 3_footprint_validation/
VALIDATION_FOLDER = SCRIPT_DIR.parent

# Project root: PHD_PAPER1_PROCESSING_v2/
PROJECT_DIR = VALIDATION_FOLDER.parent

# EC tower Excel
TOWER_EXCEL = (
    PROJECT_DIR
    / "1_ectower_extraction"
    / "outputs"
    / "tower_flight_hour.xlsx"
)

# Final output stays inside 3_footprint_validation/validation/
OUTPUT_EXCEL = (
    SCRIPT_DIR
    / "TSEB_EC_hourly_validation.xlsx"
)

# ============================================================
# TSEB BAND DEFINITIONS
# ArcGIS / rasterio physical band numbering starts at 1 here.
# ============================================================

TSEB_BANDS = {
    1: "Rn",
    2: "H",
    3: "LE",
    4: "G",
    5: "CWSI",
}


# ============================================================
# TSEB PHYSICAL QA/QC RANGES
# Conservative screening bounds for footprint validation.
# Values outside these ranges are excluded BEFORE applying
# footprint weights. Units for energy fluxes: W m-2.
#
# These are QA/QC limits, not theoretical physical limits.
# They are intentionally broad to preserve real advection signals.
# ============================================================

TSEB_VALID_RANGES = {
    "Rn": (-200.0, 1200.0),
    "H":  (-700.0, 1000.0),
    "LE": (-300.0, 1500.0),
    "G":  (-300.0, 600.0),
    # CWSI is not range-filtered here; only finite/NoData screening.
    "CWSI": None,
}


# ============================================================
# HELPERS
# ============================================================

def parse_tseb_filename(path):
    """
    Expected:
        OLA_20210729_1500_PT.tif
        VAC_20230710_1230_PM.tif
        WWF_20220710_1300_SW.tif
    """

    pattern = (
        r"^(OLA|VAC|WWF)_"
        r"(\d{8})_"
        r"(\d{4})_"
        r"(PT|PM|SW)\.tif$"
    )

    match = re.match(pattern, path.name, re.IGNORECASE)

    if not match:
        return None

    site, date, time, model = match.groups()

    return {
        "site": site.upper(),
        "date": date,
        "time": time,
        "model": model.upper(),
    }


def get_footprint_band(time_string):
    """
    Band mapping:

    Normal years:
        Band 1 = 06:00
        Band 2 = 07:00
        ...
        Band 11 = 16:00

    VAC 2023:
        Band 1 = 06:30
        Band 2 = 07:30
        ...
        Band 11 = 16:30

    Since both follow the same hourly sequence,
    the band depends on the hour component.
    """

    hour = int(time_string[:2])

    band = hour - 6 + 1

    if band < 1 or band > 11:
        raise ValueError(
            f"Time {time_string} is outside footprint range 06-16."
        )

    return band


def find_footprint(site, date):
    """
    Example:
        footprints_output/OLA/hourly/
        OLA_hourly_2021_07_29.tif
    """

    date_formatted = (
        f"{date[:4]}_{date[4:6]}_{date[6:8]}"
    )

    path = (
        VALIDATION_FOLDER
        / "footprints_output"
        / site
        / "hourly"
        / f"{site}_hourly_{date_formatted}.tif"
    )

    if not path.exists():
        return None

    return path


def read_footprint_on_tseb_grid(
    footprint_path,
    footprint_band,
    tseb_src,
):
    """
    Reproject/resample footprint to exactly match the TSEB raster grid.

    Bilinear interpolation is used because the footprint is a continuous
    weighting surface.
    """

    with rasterio.open(footprint_path) as fp_src:

        if footprint_band > fp_src.count:
            raise ValueError(
                f"Requested footprint band {footprint_band}, "
                f"but raster has only {fp_src.count} bands."
            )

        fp = fp_src.read(footprint_band).astype(np.float64)

        destination = np.zeros(
            (tseb_src.height, tseb_src.width),
            dtype=np.float64,
        )

        reproject(
            source=fp,
            destination=destination,
            src_transform=fp_src.transform,
            src_crs=fp_src.crs,
            dst_transform=tseb_src.transform,
            dst_crs=tseb_src.crs,
            resampling=Resampling.bilinear,
            src_nodata=fp_src.nodata,
            dst_nodata=0.0,
        )

    destination[
        ~np.isfinite(destination)
    ] = 0.0

    destination[
        destination < 0
    ] = 0.0

    return destination


def footprint_weighted_value(
    data,
    footprint,
    nodata=None,
    valid_range=None,
):
    """
    Footprint-weighted spatial mean after QA/QC filtering.

    Processing order:
        1. Keep finite TSEB and footprint pixels.
        2. Remove raster NoData.
        3. Keep footprint > 0.
        4. Apply the variable-specific physical range.
        5. Calculate the weighted product sum:

               weighted_sum = sum(X * FP)

        6. Renormalize by the footprint mass remaining over valid pixels:

               weighted_mean = weighted_sum / sum(FP_valid)

    Renormalization is important after filtering because invalid TSEB pixels
    remove part of the original footprint mass. Without it, the resulting
    flux would be biased downward merely because some pixels were rejected.

    Returns
    -------
    weighted_mean : float
        Final TSEB value used for validation.
    weight_sum : float
        Footprint mass remaining over valid TSEB pixels.
    n_valid : int
        Number of valid pixels used.
    weighted_sum : float
        Raw sum(X * FP) before renormalization.
    n_rejected_range : int
        Number of footprint-covered pixels rejected only by the physical
        range test.
    """

    data = data.astype(np.float64)
    footprint = footprint.astype(np.float64)

    # Base validity: finite data + finite positive footprint.
    base_valid = (
        np.isfinite(data)
        & np.isfinite(footprint)
        & (footprint > 0)
    )

    if nodata is not None and np.isfinite(nodata):
        base_valid &= data != nodata

    valid = base_valid.copy()
    n_rejected_range = 0

    # Variable-specific physical screening.
    if valid_range is not None:
        vmin, vmax = valid_range
        in_range = (data >= vmin) & (data <= vmax)
        n_rejected_range = int(np.sum(base_valid & ~in_range))
        valid &= in_range

    if not np.any(valid):
        return np.nan, np.nan, 0, np.nan, n_rejected_range

    weights = footprint[valid]
    values = data[valid]

    weight_sum = float(np.sum(weights))

    if not np.isfinite(weight_sum) or weight_sum <= 0:
        return np.nan, np.nan, 0, np.nan, n_rejected_range

    weighted_sum = float(np.sum(values * weights))
    weighted_mean = weighted_sum / weight_sum

    return (
        float(weighted_mean),
        weight_sum,
        int(np.sum(valid)),
        weighted_sum,
        n_rejected_range,
    )


# ============================================================
# READ EC TOWER DATA
# ============================================================

tower = pd.read_excel(TOWER_EXCEL)

tower["SITE"] = (
    tower["SITE"]
    .astype(str)
    .str.upper()
    .str.strip()
)

tower["TIMESTAMP_END"] = (
    tower["TIMESTAMP_END"]
    .astype(str)
    .str.replace(".0", "", regex=False)
    .str.zfill(12)
)

tower["KEY"] = (
    tower["SITE"]
    + "_"
    + tower["TIMESTAMP_END"]
)


# ============================================================
# PROCESS TSEB RASTERS
# ============================================================

results = []

for model in ["PT", "PM", "SW"]:

    model_folder = (
        PROJECT_DIR
        / "2_tseb_loop"
        / "outputs"
        / model
    )

    tif_files = sorted(
        model_folder.glob("*.tif")
    )

    print("\n====================================")
    print(f"MODEL: {model}")
    print(f"Files: {len(tif_files)}")
    print("====================================")

    for tseb_path in tif_files:

        info = parse_tseb_filename(tseb_path)

        if info is None:
            print(
                f"Skipping unrecognized file: "
                f"{tseb_path.name}"
            )
            continue

        site = info["site"]
        date = info["date"]
        time = info["time"]

        timestamp = date + time
        key = f"{site}_{timestamp}"

        print(
            f"{site} {date} {time} {model}",
            end=" -> "
        )

        # ----------------------------------------------------
        # Locate footprint
        # ----------------------------------------------------

        footprint_path = find_footprint(
            site,
            date,
        )

        if footprint_path is None:

            print("FOOTPRINT NOT FOUND")

            results.append({
                "SITE": site,
                "TIMESTAMP_END": timestamp,
                "MODEL": model,
                "STATUS": "FOOTPRINT_NOT_FOUND",
            })

            continue

        footprint_band = get_footprint_band(
            time
        )

        # ----------------------------------------------------
        # Open TSEB
        # ----------------------------------------------------

        with rasterio.open(tseb_path) as tseb_src:

            footprint = read_footprint_on_tseb_grid(
                footprint_path,
                footprint_band,
                tseb_src,
            )

            fp_sum_resampled = float(
                np.sum(footprint)
            )

            row = {
                "SITE": site,
                "TIMESTAMP_END": timestamp,
                "DATE": int(date),
                "TIME": int(time),
                "MODEL": model,
                "FOOTPRINT_BAND": footprint_band,
                "FOOTPRINT_SUM_RESAMPLED": fp_sum_resampled,
                "TSEB_FILE": tseb_path.name,
                "FOOTPRINT_FILE": footprint_path.name,
                "STATUS": "OK",
            }

            # ------------------------------------------------
            # Extract each TSEB variable
            # ------------------------------------------------

            for band_number, variable in TSEB_BANDS.items():

                if band_number > tseb_src.count:

                    row[f"TSEB_{variable}"] = np.nan
                    row[f"WEIGHTED_SUM_{variable}"] = np.nan
                    row[f"FP_MASS_{variable}"] = np.nan
                    row[f"NPIX_{variable}"] = 0
                    row[f"NPIX_RANGE_REJECTED_{variable}"] = 0

                    continue

                data = tseb_src.read(
                    band_number
                )

                (
                    value,
                    fp_mass,
                    npix,
                    weighted_sum,
                    n_rejected_range,
                ) = footprint_weighted_value(
                    data,
                    footprint,
                    tseb_src.nodata,
                    TSEB_VALID_RANGES.get(variable),
                )

                row[f"TSEB_{variable}"] = value
                row[f"WEIGHTED_SUM_{variable}"] = weighted_sum
                row[f"FP_MASS_{variable}"] = fp_mass
                row[f"NPIX_{variable}"] = npix
                row[f"NPIX_RANGE_REJECTED_{variable}"] = n_rejected_range

            results.append(row)

        print(
            f"OK | FP band {footprint_band}"
        )


# ============================================================
# RESULTS TABLE
# ============================================================

results_df = pd.DataFrame(results)

results_df["KEY"] = (
    results_df["SITE"]
    + "_"
    + results_df["TIMESTAMP_END"]
)


# ============================================================
# MERGE WITH EC TOWER
# ============================================================

tower_columns = [
    "KEY",
    "ID",
    "TA_1_1_1",
    "VPD_1_1_1",
    "SW_IN",
    "LW_IN",
    "WS",
    "PA",
    "NETRAD",
    "H",
    "LE",
    "G",
    "MATCH_STATUS",
]

tower_subset = tower[
    [
        col
        for col in tower_columns
        if col in tower.columns
    ]
].copy()

final = results_df.merge(
    tower_subset,
    on="KEY",
    how="left",
)


# ============================================================
# RENAME EC VARIABLES
# ============================================================

final = final.rename(
    columns={
        "NETRAD": "EC_Rn",
        "H": "EC_H",
        "LE": "EC_LE",
        "G": "EC_G",
    }
)


# ============================================================
# EC TOWER ENERGY-BALANCE QA/QC
# ============================================================

available_energy = (
    final["EC_Rn"]
    - final["EC_G"]
)

turbulent_flux = (
    final["EC_H"]
    + final["EC_LE"]
)

# Energy Balance Ratio / Closure Ratio
final["EC_EBR"] = (
    turbulent_flux
    / available_energy
)

# Energy balance residual [W m-2]
final["EC_RESIDUAL"] = (
    available_energy
    - turbulent_flux
)

# Bowen ratio
final["EC_BOWEN"] = (
    final["EC_H"]
    / final["EC_LE"]
)


# ============================================================
# DIFFERENCES TSEB - EC
# Useful directly for validation
# ============================================================

for variable in ["Rn", "H", "LE", "G"]:

    tseb_col = f"TSEB_{variable}"
    ec_col = f"EC_{variable}"

    if (
        tseb_col in final.columns
        and ec_col in final.columns
    ):

        final[f"DIFF_{variable}"] = (
            final[tseb_col]
            - final[ec_col]
        )


# ============================================================
# ORDER
# ============================================================

model_order = {
    "PT": 1,
    "PM": 2,
    "SW": 3,
}

final["_MODEL_ORDER"] = (
    final["MODEL"]
    .map(model_order)
)

final = (
    final
    .sort_values(
        [
            "SITE",
            "TIMESTAMP_END",
            "_MODEL_ORDER",
        ]
    )
    .drop(
        columns=[
            "KEY",
            "_MODEL_ORDER",
        ],
        errors="ignore",
    )
    .reset_index(drop=True)
)


# ============================================================
# SAVE EXCEL
# ============================================================

with pd.ExcelWriter(
    OUTPUT_EXCEL,
    engine="openpyxl",
) as writer:

    # Main combined table
    final.to_excel(
        writer,
        sheet_name="Validation",
        index=False,
    )

    # Separate sheets by model
    for model in ["PT", "PM", "SW"]:

        temp = final[
            final["MODEL"] == model
        ]

        temp.to_excel(
            writer,
            sheet_name=model,
            index=False,
        )


# ============================================================
# SUMMARY
# ============================================================

print("\n====================================")
print("PROCESSING COMPLETE")
print("====================================")

print(
    f"Results: {len(final)}"
)

print("\nResults by model:")
print(
    final["MODEL"].value_counts()
)

print("\nStatus:")
print(
    final["STATUS"].value_counts(
        dropna=False
    )
)

print(
    f"\nSaved:\n{OUTPUT_EXCEL}"
)
