import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.colors import LogNorm

RESULTS_DIR = Path(__file__).resolve().parent / "results"
PLOTS_DIR = RESULTS_DIR / "plots"
READV2_CSV = RESULTS_DIR / "readv2_comparison_benchmark.csv"
OUTPUTS_DIR = RESULTS_DIR / "outputs"
FASTPMR_NPZ = OUTPUTS_DIR / "fastpmr" / "fastpmr_results.npz"
READV2_TSV = OUTPUTS_DIR / "readv2" / "Read_Results.tsv"
PMR_TOLERANCE = 1e-6

DEGREE_ORDER = [
    "Identical/Twin",
    "First Degree",
    "Second Degree",
    "Third Degree",
    "Unrelated",
]

READV2_DEGREE_MAP = {
    "IdenticalTwins/SameIndividual": "Identical/Twin",
    "First Degree": "First Degree",
    "Second Degree": "Second Degree",
    "Third Degree": "Third Degree",
    "Unrelated/Consistent with Third Degree": "Unrelated",
    "Unrelated": "Unrelated",
}


def seconds_to_minutes(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df["mean_min"] = df["mean_s"] / 60
    df["stddev_min"] = df["stddev_s"] / 60
    return df


def bytes_to_gb(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df["mean_gb"] = df["mean_bytes"] / (1024**3)
    df["stddev_gb"] = df["stddev_bytes"] / (1024**3)
    return df


def format_duration(seconds: float) -> str:
    if seconds < 60:
        return f"{round(seconds)}s"
    elif seconds < 3600:
        mins = int(seconds // 60)
        secs = int(seconds % 60)
        return f"{mins}min {secs}s" if secs else f"{mins}min"
    else:
        hrs = int(seconds // 3600)
        mins = int((seconds % 3600) // 60)
        return f"{hrs}hr {mins}min"


def format_memory(bytes: float) -> str:
    gb = bytes / (1024**3)
    mb = bytes / (1024**2)
    if gb >= 1:
        return f"{gb:.2f} GB"
    elif mb >= 1:
        return f"{mb:.0f} MB"
    else:
        return f"{bytes:.0f} B"


def save_bar_plot(
    data: pd.DataFrame,
    label_col: str,
    y_col: str,
    err_col: str,
    ylabel: str,
    title: str,
    output_path: Path,
    seconds_col: str | None = None,
    bytes_col: str | None = None,
) -> None:
    assert seconds_col or bytes_col
    fig, ax = plt.subplots(figsize=(5, 5), constrained_layout=True)
    ax.bar(
        data[label_col],
        data[y_col],
        yerr=data[err_col],
        width=0.6,
        capsize=14,
        error_kw={"elinewidth": 2, "capthick": 1},
    )
    y_max = (data[y_col] + data[err_col]).max()
    for x, (y, err) in enumerate(zip(data[y_col], data[err_col], strict=True)):
        if seconds_col:
            annotation = format_duration(data[seconds_col].iloc[x])
        else:
            annotation = format_memory(data[bytes_col].iloc[x])
        ax.text(
            x,
            y + err + y_max * 0.01,
            annotation,
            ha="center",
            va="bottom",
            fontsize=14,
        )
    ax.margins(x=0.2)
    ax.set_title(title, fontsize=16)
    ax.set_ylabel(ylabel, fontsize=16)
    ax.tick_params(axis="both", labelsize=16)
    ax.set_ylim(bottom=0)
    y_max = (data[y_col] + data[err_col]).max()
    ax.set_ylim(bottom=0, top=y_max * 1.08)
    ax.grid(True, axis="y", linewidth=0.8, alpha=0.4)
    sns.despine(ax=ax)
    fig.savefig(output_path, bbox_inches="tight", dpi=600)
    plt.close(fig)


def build_pair_comparison() -> pd.DataFrame:
    # Load fastpmr results
    npz = np.load(FASTPMR_NPZ, allow_pickle=True)
    samples = json.loads(npz["samples.json"])
    degree_labels = json.loads(npz["degree_labels.json"])
    degrees = npz["degrees"]
    mismatch_rates = npz["mismatch_rates"]

    # Map IID -> sample index (READv2 uses IID from the .fam file)
    iid_to_idx = {}
    for idx, name in enumerate(samples):
        iid = name.split(":", 1)[1]
        iid_to_idx[iid] = idx

    # Load READv2 results
    readv2_df = pd.read_csv(READV2_TSV, sep="\t")

    rows = []
    for _, row in readv2_df.iterrows():
        id_a, id_b = row["PairIndividuals"].split(",")
        idx_a = iid_to_idx[id_a]
        idx_b = iid_to_idx[id_b]
        rows.append(
            {
                "fastpmr_degree": degree_labels[degrees[idx_a, idx_b]],
                "READv2_degree": READV2_DEGREE_MAP[row["Rel"]],
                "fastpmr_mismatch_rate": float(mismatch_rates[idx_a, idx_b]),
                "READv2_mismatch_rate": float(row["Nonnormalized_P0"]),
            }
        )
    return pd.DataFrame(rows)


def save_confusion_matrix(comparison: pd.DataFrame, output_path: Path) -> None:
    ct = pd.crosstab(
        comparison["READv2_degree"],
        comparison["fastpmr_degree"],
        dropna=False,
    )
    # Reindex to canonical degree order (fill missing with 0)
    ct = ct.reindex(index=DEGREE_ORDER, columns=DEGREE_ORDER, fill_value=0)

    fig, ax = plt.subplots(figsize=(7, 6), constrained_layout=True)
    sns.heatmap(
        ct,
        annot=True,
        fmt="d",
        cmap="Blues",
        norm=LogNorm(vmin=1, vmax=ct.values.max()),
        linewidths=0.5,
        linecolor="white",
        square=True,
        cbar_kws={"label": "Number of pairs", "shrink": 0.8},
        ax=ax,
    )
    ax.set_xlabel("fastpmr", fontsize=14)
    ax.set_ylabel("READv2", fontsize=14)
    ax.set_title("Degree Classification Confusion Matrix", fontsize=15)
    ax.tick_params(axis="both", labelsize=11)
    # Rotate labels for readability
    ax.set_xticklabels(ax.get_xticklabels(), rotation=30, ha="right")
    ax.set_yticklabels(ax.get_yticklabels(), rotation=0)

    fig.savefig(output_path, bbox_inches="tight", dpi=600)
    plt.close(fig)


def check_pmr_agreement(comparison: pd.DataFrame) -> None:
    diffs = (comparison["fastpmr_mismatch_rate"] - comparison["READv2_mismatch_rate"]).abs()
    max_diff = float(diffs.max())
    if max_diff < PMR_TOLERANCE:
        print(f"All mismatch rates match within {PMR_TOLERANCE} (max diff = {max_diff:.2e}).")
    else:
        print(f"Some mismatch rates differ by more than {PMR_TOLERANCE} (max diff = {max_diff:.2e}).")


def main() -> None:
    PLOTS_DIR.mkdir(parents=True, exist_ok=True)

    pair_comparison = build_pair_comparison()
    save_confusion_matrix(pair_comparison, PLOTS_DIR / "degree_classification_confusion_matrix.pdf")
    check_pmr_agreement(pair_comparison)

    readv2_df = pd.read_csv(READV2_CSV)
    readv2_df = bytes_to_gb(readv2_df)
    readv2_df = seconds_to_minutes(readv2_df)
    save_bar_plot(
        readv2_df,
        "label",
        "mean_min",
        "stddev_min",
        "Mean Runtime (min)",
        "Runtime: fastpmr vs. READv2",
        PLOTS_DIR / "readv2_comparison_benchmark_runtime.pdf",
        seconds_col="mean_s",
    )
    save_bar_plot(
        readv2_df,
        "label",
        "mean_gb",
        "stddev_gb",
        "Peak RSS (GB)",
        "Peak Memory: fastpmr vs. READv2",
        PLOTS_DIR / "readv2_comparison_benchmark_memory.pdf",
        bytes_col="mean_bytes",
    )


if __name__ == "__main__":
    main()
