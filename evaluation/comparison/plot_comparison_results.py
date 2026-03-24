from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

RESULTS_DIR = Path(__file__).resolve().parent / "results"
PLOTS_DIR = RESULTS_DIR / "plots"
READV2_CSV = RESULTS_DIR / "readv2_comparison_benchmark.csv"


def bytes_to_mb(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df["mean_mb"] = df["mean_bytes"] / (1024**2)
    df["stddev_mb"] = df["stddev_bytes"] / (1024**2)
    return df


def save_bar_plot(
    data: pd.DataFrame,
    label_col: str,
    y_col: str,
    err_col: str,
    ylabel: str,
    title: str,
    output_path: Path,
) -> None:
    fig, ax = plt.subplots(figsize=(8, 10), constrained_layout=True)
    ax.bar(
        data[label_col],
        data[y_col],
        yerr=data[err_col],
        width=0.6,
        capsize=14,
        error_kw={"elinewidth": 2.5, "capthick": 1.5},
    )
    ax.margins(x=0.2)
    ax.set_title(title, fontsize=16)
    ax.set_ylabel(ylabel, fontsize=16)
    ax.tick_params(axis="both", labelsize=16)
    ax.set_ylim(bottom=0)
    ax.grid(True, axis="y", linewidth=0.8, alpha=0.4)
    sns.despine(ax=ax)
    fig.savefig(output_path, bbox_inches="tight", dpi=600)
    plt.close(fig)


def main() -> None:
    PLOTS_DIR.mkdir(parents=True, exist_ok=True)

    readv2_df = pd.read_csv(READV2_CSV)
    readv2_df = bytes_to_mb(readv2_df)
    save_bar_plot(
        readv2_df,
        "label",
        "mean_s",
        "stddev_s",
        "Mean Runtime (s)",
        "Runtime: fastpmr vs. READv2",
        PLOTS_DIR / "readv2_comparison_benchmark_runtime.pdf",
    )
    save_bar_plot(
        readv2_df,
        "label",
        "mean_mb",
        "stddev_mb",
        "Peak RSS (MB)",
        "Peak Memory: fastpmr vs. READv2",
        PLOTS_DIR / "readv2_comparison_benchmark_memory.pdf",
    )


if __name__ == "__main__":
    main()
