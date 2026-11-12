[README_footprint_validation.md](https://github.com/user-attachments/files/34771817/README_footprint_validation.md)
# 🌿 Footprint-weighted validation of UAV TSEB in almond orchards

**From eddy covariance (EC) footprint estimates to tower-representative TSEB fluxes and model evaluation**

![Python](https://img.shields.io/badge/Python-3.10%2B-blue) ![Workflow](https://img.shields.io/badge/workflow-FFP%20%7C%20TSEB%20%7C%20EC-298068) ![Domain](https://img.shields.io/badge/domain-agrometeorology-28724c)

*Micrometeorology · UAV remote sensing · Surface energy balance · Almond orchards*

## Overview

This repository estimates **hourly flux footprints** from EC observations, aligns them with UAV-based **Two-Source Energy Balance (TSEB)** rasters, and evaluates modeled energy fluxes against EC measurements. The example workflow covers three California almond-orchard sites (**OLA, VAC, WWF**) and three TSEB formulations (**PT, PM, SW**).

TSEB GeoTIFFs are produced upstream by [TerraObserver/2-tseb_loop](https://github.com/TerraObserver/2-tseb_loop), which runs models based on [H. Nieto's pyTSEB](https://github.com/hectornieto/pyTSEB). Footprints use the **Kljun et al. (2015)** FFP parameterization, included here as `calc_footprint_FFP_climatology.py`.

## Method

1. **Footprint estimation:** compute one footprint per available hour (06:00–16:00), retain the leading 90% of footprint mass, and normalize each valid hourly footprint. Optional daily composites use normalized hourly ETo weights; **the hourly footprints**, not daily composites, are used in TSEB–EC validation.
2. **Spatial matching:** reproject/resample each hourly footprint onto the corresponding TSEB raster grid, mask invalid TSEB pixels and apply variable-specific QA/QC limits.
3. **Weighted extraction:** calculate a footprint-weighted spatial mean for each variable:

   $$\overline{F}=\frac{\sum_{i\in V}w_iF_i}{\sum_{i\in V}w_i}$$

   Here, $V$ denotes valid pixels and $w_i$ is the footprint weight. The retained weight sum (`FP_MASS_*`) and number of pixels (`NPIX_*`) provide extraction diagnostics.
4. **EC comparison:** match tower records by site and timestamp, calculate energy-balance diagnostics, and visualize TSEB–EC agreement. For the figure, EC sensible and latent heat are corrected using a **Bowen-ratio-preserving** energy-balance closure adjustment; EC net radiation and soil heat flux remain unchanged.

## Repository contents

| File | Purpose |
|---|---|
| `1_footprint_processing.py` | Create hourly footprint GeoTIFF stacks and ETo-weighted daily footprints. |
| `calc_footprint_FFP_climatology.py` | Included Kljun FFP implementation. |
| `tower_footprint_input.xlsx` | Example meteorological footprint inputs. |
| `validation/1_loop_validation.py` | Extract footprint-weighted TSEB fluxes; match with EC data; export Excel. |
| `validation/2_validation_plot.py` | Plot and print PT/PM/SW comparison statistics. |
| `validation/TSEB_EC_hourly_validation.xlsx` | **Included example result** for reproducing the plot without upstream files. |

## Installation and data layout

Python **3.10+** is recommended. From the `3-footprint_validation` folder:

```bash
python -m pip install -r requirements.txt
```

For full end-to-end processing, place the repositories and the EC extraction output under one parent project directory:

```text
PROJECT/
├── 1_ectower_extraction/
│   └── outputs/tower_flight_hour.xlsx   # Required; not included
├── 2-tseb_loop/
│   └── outputs_rst/
│       ├── PT/                       # e.g., OLA_20210729_1500_PT.tif
│       ├── PM/
│       └── SW/
└── 3-footprint_validation/
    ├── tower_footprint_input.xlsx
    ├── 1_footprint_processing.py
    ├── calc_footprint_FFP_climatology.py
    └── validation/
        ├── 1_loop_validation.py
        ├── 2_validation_plot.py
        └── TSEB_EC_hourly_validation.xlsx
```

**Inputs:** `tower_footprint_input.xlsx` needs `date`, `SITE_ID`, `wind_dir`, `u_star`, `L`, `sigma_v`, `ETo`. The EC comparison workbook needs `SITE`, `TIMESTAMP_END`, `NETRAD`, `H`, `LE`, and `G` (other tower variables are retained where available). TSEB primary GeoTIFFs must follow `<SITE>_<YYYYMMDD>_<HHMM>_<PT|PM|SW>.tif`, with bands **1–5 = Rn, H, LE, G, CWSI**. Verify this band mapping for any altered pyTSEB export.

## Run the workflow

```bash
# 1. Create site-specific footprint GeoTIFFs
python 1_footprint_processing.py

# 2. Join hourly footprints to TSEB outputs and EC observations
python validation/1_loop_validation.py

# 3. Produce a validation figure and summary statistics
python validation/2_validation_plot.py
```

**Quick demonstration:** step 3 runs immediately using the included `validation/TSEB_EC_hourly_validation.xlsx`. Steps 1–2 require the corresponding input observations and TSEB outputs; step 2 additionally requires `tower_flight_hour.xlsx`.

## Main outputs and checks

- `footprints_output/<SITE>/hourly/<SITE>_hourly_YYYY_MM_DD.tif`: **11 bands**, corresponding to 06:00–16:00; missing hours contain zero-weight bands.
- `footprints_output/<SITE>/daily_weighted/`: daily ETo-weighted footprint rasters, created when at least five hourly footprints are valid.
- `footprints_output/footprint_processing_summary.csv`: execution status by site and day.
- `validation/TSEB_EC_hourly_validation.xlsx`: merged `Validation`, `PT`, `PM`, and `SW` sheets, including `STATUS`, `FP_MASS_*`, `NPIX_*`, and differences (`TSEB - EC`).
- `validation/TSEB_vs_EC_closure_statistics.png`: pooled Rn/H/LE/G scatterplot and model-level statistics.

**Before interpreting results:** check that timestamps refer to the same observation interval, that footprint and TSEB raster extents overlap, and that `FP_MASS_*` is sufficiently high. The special **VAC 2023** footprint schedule is at **:30**; other example dates use **:00**. The plotted performance statistics **pool four flux components** and should not replace separate Rn/H/LE/G diagnostics for scientific conclusions. Changing spatial QA/QC thresholds changes the values extracted.

## References and attribution

- **Footprint processing:** Original footprint implementation by [Karem Meza Capcha](https://github.com/KaremMezaCapcha), adapted and integrated into a batch-processing loop for this workflow.
- **Flux footprints:** Kljun, N., Calanca, P., Rotach, M. W., & Schmid, H. P. (2015). *A simple two-dimensional parameterisation for Flux Footprint Prediction (FFP).* Geoscientific Model Development, 8, 3695–3713. https://doi.org/10.5194/gmd-8-3695-2015
- **Energy-balance modeling:** [Nieto et al., pyTSEB codebase](https://github.com/hectornieto/pyTSEB) and the project-specific [TSEB loop](https://github.com/TerraObserver/2-tseb_loop).

If redistributing the bundled FFP implementation, retain its original attribution and check the applicable source-code license/terms; the pyTSEB upstream project is GPL-3.0 licensed.
