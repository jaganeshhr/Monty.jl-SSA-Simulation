import glob
import os

import pandas as pd

# ---------------------------------------------------------------------------#
try:
    _SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
except NameError:
    _SCRIPT_DIR = os.getcwd()  # fallback when run interactively (e.g. a notebook cell)

RESULTS_DIR = os.path.join(_SCRIPT_DIR, "Results") 
A = 1000.0  # Area of treatment (m^2)
d = 0.15    # Sample core depth (m)
# ---------------------------------------------------------------------------#

MINERALS = [
    "Plagioclase", "Clinopyroxene", "Quartz", "Sericite", "FeTiOxides",
    "Chlorite", "Actinolite", "Olivine", "Orthoclase", "Nepheline",
    "Muscovite", "Apatite", "Calcite",
]

M_ELEMENT = {
    "Plagioclase":   [0.078,  0.0,    0.035,  0.004],
    "Clinopyroxene": [0.053,  0.088,  0.0,    0.0],
    "Quartz":        [0.0,    0.0,    0.0,    0.0],
    "Sericite":      [0.0,    0.0,    0.013,  0.0981],
    "FeTiOxides":    [0.0112, 0.0,    0.0,    0.0],
    "Chlorite":      [0.0,    0.09,   0.0,    0.0],
    "Actinolite":    [0.104,  0.0813, 0.0,    0.0],
    "Olivine":       [0.0,    0.3455, 0.0,    0.0],
    "Orthoclase":    [0.0,    0.0,    0.0,    0.1405],
    "Nepheline":     [0.0,    0.0,    0.1180, 0.0669],
    "Muscovite":     [0.0,    0.0,    0.0,    0.0981],
    "Apatite":       [0.3936, 0.0,    0.0,    0.0],
    "Calcite":       [0.4004, 0.0,    0.0,    0.0],
}

ELEMENTS = ["Ca", "Mg", "Na", "K"]
CATION_FACTORS = {"Ca": 2.196, "Mg": 3.62, "Na": 1.91, "K": 1.12}


def element_concentration(row: pd.Series, elem_idx: int) -> float:
    """Sum(mineral_conc * element_fraction) across all tracked minerals."""
    total = 0.0
    for m in MINERALS:
        if m in row:
            total += row[m] * M_ELEMENT[m][elem_idx]
    return total


def compute_cdr_for_sample(sample_csv_dir: str) -> pd.DataFrame:
    """Load all per-iteration CSVs for one sample and compute cumulative CDR by round."""
    csv_files = glob.glob(os.path.join(sample_csv_dir, "*.csv"))
    if not csv_files:
        raise FileNotFoundError(f"No CSV files found in {sample_csv_dir}")

    df_list = [pd.read_csv(f) for f in csv_files]
    combined = pd.concat(df_list, ignore_index=True)

    value_cols = [c for c in combined.columns if c not in ("round", "control")]
    data_avg = (
        combined.groupby(["round", "control"], as_index=False)[value_cols].mean()
    )

    sorted_rounds = sorted(data_avg["round"].unique())
    baseline_round = sorted_rounds[1] if len(sorted_rounds) > 1 else sorted_rounds[0]

    cdr_rounds = [sorted_rounds[0]]
    cdr_values = [0.0]

    for r_val in sorted_rounds[2:]:
        treat_base = data_avg[(data_avg["control"] == 0) & (data_avg["round"] == baseline_round)]
        treat_curr = data_avg[(data_avg["control"] == 0) & (data_avg["round"] == r_val)]

        if len(treat_base) > 0 and len(treat_curr) > 0:
            base_row = treat_base.iloc[0]
            curr_row = treat_curr.iloc[0]

            loss_sum = 0.0
            for i, elem in enumerate(ELEMENTS):
                base_conc = element_concentration(base_row, i)
                curr_conc = element_concentration(curr_row, i)
                loss = max(0.0, base_conc - curr_conc)
                loss_sum += loss * CATION_FACTORS[elem]

            cdr = (loss_sum * 1e6 / A) * d
            cdr_values.append(cdr)
            cdr_rounds.append(r_val)

    return pd.DataFrame({"round_years": cdr_rounds, "cdr_kg_per_m2": cdr_values})


def main():
    sample_dirs = sorted(glob.glob(os.path.join(RESULTS_DIR, "*_simulation")))
    if not sample_dirs:
        raise FileNotFoundError(
            f"No '*_simulation' folders found under '{RESULTS_DIR}'. "
            "Update RESULTS_DIR to point at your Results folder."
        )

    summary_rows = []
    all_curves = {}

    for sample_dir in sample_dirs:
        sample_name = os.path.basename(sample_dir).replace("_simulation", "")
        csv_dir = os.path.join(sample_dir, "csv")

        try:
            cdr_df = compute_cdr_for_sample(csv_dir)
        except FileNotFoundError as e:
            print(f"Skipping {sample_name}: {e}")
            continue

        all_curves[sample_name] = cdr_df
        total_cdr = cdr_df["cdr_kg_per_m2"].iloc[-1]  # final round = cumulative total
        summary_rows.append({"sample": sample_name, "total_cdr_kg_per_m2": total_cdr})

    summary = pd.DataFrame(summary_rows).sort_values(
        "total_cdr_kg_per_m2", ascending=False
    ).reset_index(drop=True)

    print("\n=== Total CDR drawdown by sample (kg CO2/m^2, final round) ===")
    print(summary.to_string(index=False))

    if not summary.empty:
        best = summary.iloc[0]
        worst = summary.iloc[-1]
        print(f"\nBest performer:  {best['sample']}  ({best['total_cdr_kg_per_m2']:.4f} kg/m^2)")
        print(f"Worst performer: {worst['sample']} ({worst['total_cdr_kg_per_m2']:.4f} kg/m^2)")

    # Save outputs
    summary.to_csv(os.path.join(RESULTS_DIR, "cdr_summary_ranked.csv"), index=False)
    for sample_name, cdr_df in all_curves.items():
        cdr_df.to_csv(
            os.path.join(RESULTS_DIR, f"{sample_name}_cdr_curve.csv"), index=False
        )

    return summary, all_curves


if __name__ == "__main__":
    main()
