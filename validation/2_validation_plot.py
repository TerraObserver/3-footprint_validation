from pathlib import Path

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# ============================================================
# PATHS / INPUT
# ============================================================

SCRIPT_DIR = Path(__file__).resolve().parent

excel_path = (
    SCRIPT_DIR
    / "TSEB_EC_hourly_validation.xlsx"
)

output_figure = (
    SCRIPT_DIR
    / "TSEB_vs_EC_closure_statistics.png"
)

df = pd.read_excel(
    excel_path,
    sheet_name="Validation"
)


# ============================================================
# EC TOWER ENERGY BALANCE CLOSURE
# Bowen-ratio preserving correction
# ============================================================

available_energy = df["EC_Rn"] - df["EC_G"]
turbulent_flux = df["EC_H"] + df["EC_LE"]

valid = (
    np.isfinite(available_energy)
    & np.isfinite(turbulent_flux)
    & (np.abs(turbulent_flux) > 1e-6)
)

df["EC_H_CORR"] = np.nan
df["EC_LE_CORR"] = np.nan

df.loc[valid, "EC_H_CORR"] = (
    available_energy[valid]
    * df.loc[valid, "EC_H"]
    / turbulent_flux[valid]
)

df.loc[valid, "EC_LE_CORR"] = (
    available_energy[valid]
    * df.loc[valid, "EC_LE"]
    / turbulent_flux[valid]
)

df["EC_Rn_CORR"] = df["EC_Rn"]
df["EC_G_CORR"] = df["EC_G"]


# ============================================================
# SETTINGS
# ============================================================

components = ["Rn", "H", "LE", "G"]

models = {
    "PT": "o",
    "PM": "s",
    "SW": "^"
}


# ============================================================
# PLOT
# ============================================================

fig, ax = plt.subplots(figsize=(9, 8))

stats_results = []


for model, marker in models.items():

    sub = df[df["MODEL"] == model]

    x_all = []
    y_all = []

    # Combine Rn, H, LE and G
    for variable in components:

        x = sub[f"TSEB_{variable}"]
        y = sub[f"EC_{variable}_CORR"]

        good = (
            np.isfinite(x)
            & np.isfinite(y)
        )

        x_all.extend(x[good])
        y_all.extend(y[good])

    x_all = np.asarray(x_all, dtype=float)
    y_all = np.asarray(y_all, dtype=float)


    # ========================================================
    # ERROR STATISTICS
    # ========================================================

    if len(x_all) < 2:
        print(f"{model}: insufficient valid data, skipped.")
        continue

    error = x_all - y_all

    rmse = np.sqrt(
        np.mean(error ** 2)
    )

    mae = np.mean(
        np.abs(error)
    )

    bias = np.mean(error)


    # ========================================================
    # PEARSON CORRELATION
    # ========================================================

    if np.std(x_all) > 0 and np.std(y_all) > 0:
        pearson_r = np.corrcoef(
            x_all,
            y_all
        )[0, 1]
    else:
        pearson_r = np.nan


    # ========================================================
    # LINEAR REGRESSION R²
    # EC = slope * TSEB + intercept
    # ========================================================

    if np.std(x_all) > 0:
        slope, intercept = np.polyfit(
            x_all,
            y_all,
            1
        )

        y_pred = slope * x_all + intercept

        ss_res = np.sum(
            (y_all - y_pred) ** 2
        )

        ss_tot = np.sum(
            (y_all - np.mean(y_all)) ** 2
        )

        r2 = (
            1 - (ss_res / ss_tot)
            if ss_tot > 0
            else np.nan
        )
    else:
        slope = np.nan
        intercept = np.nan
        r2 = np.nan

    n = len(x_all)


    stats_results.append({
        "MODEL": model,
        "N": n,
        "R2": r2,
        "PEARSON_R": pearson_r,
        "RMSE": rmse,
        "MAE": mae,
        "BIAS": bias,
        "SLOPE": slope,
        "INTERCEPT": intercept
    })


    # ========================================================
    # LABEL
    # ========================================================

    label = (
        f"{model} | "
        f"R²={r2:.2f} | "
        f"Pearson r={pearson_r:.2f} | "
        f"RMSE={rmse:.1f} | "
        f"MAE={mae:.1f} | "
        f"Bias={bias:.1f}"
    )


    ax.scatter(
        x_all,
        y_all,
        marker=marker,
        s=42,
        alpha=0.60,
        label=label
    )


# ============================================================
# 1:1 LINE
# ============================================================

lims = ax.get_xlim() + ax.get_ylim()

minimum = min(lims)
maximum = max(lims)

ax.plot(
    [minimum, maximum],
    [minimum, maximum],
    "--",
    linewidth=1.3,
    label="1:1"
)

ax.set_xlim(minimum, maximum)
ax.set_ylim(minimum, maximum)


# ============================================================
# FORMAT
# ============================================================

ax.set_xlabel(
    "TSEB model (W m$^{-2}$)",
    fontsize=12
)

ax.set_ylabel(
    "EC tower — closure corrected (W m$^{-2}$)",
    fontsize=12
)

ax.set_title(
    "Hourly TSEB vs EC Tower",
    fontsize=14
)

ax.grid(
    True,
    alpha=0.20
)

ax.legend(
    fontsize=8,
    loc="upper left"
)

plt.tight_layout()

plt.savefig(
    output_figure,
    dpi=300,
    bbox_inches="tight"
)

plt.show()


# ============================================================
# PRINT STATISTICS
# ============================================================

stats_df = pd.DataFrame(stats_results)

print("\nMODEL PERFORMANCE")
print("============================================================")

print(
    stats_df.round({
        "R2": 3,
        "PEARSON_R": 3,
        "RMSE": 2,
        "MAE": 2,
        "BIAS": 2,
        "SLOPE": 3,
        "INTERCEPT": 2
    }).to_string(index=False)
)
print(f"\nSaved figure:\n{output_figure}")
