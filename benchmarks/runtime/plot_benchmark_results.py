from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from benchmark_utils import SCRIPT_DIR
from matplotlib.ticker import MaxNLocator

RESULTS_DIR = SCRIPT_DIR / "results"
PLOTS_DIR = RESULTS_DIR / "plots"

THREADS_CSV = RESULTS_DIR / "thread_count_benchmark.csv"
VARIANTS_CSV = RESULTS_DIR / "variant_count_benchmark.csv"
PAIRS_CSV = RESULTS_DIR / "pair_count_benchmark.csv"


def parse_thread_count(label: str) -> int:
    return int(label.split("threads=", 1)[1])


def parse_variant_count(label: str) -> int:
    spec = label.split("variants=", 1)[1]
    if "-" in spec:
        spec = spec.split("-")[-1]
    return int(spec)


def parse_pair_count(label: str) -> int:
    fields = {}
    for part in label.split("_"):
        if "=" not in part:
            continue
        key, value = part.split("=", 1)
        fields[key] = int(value)
    return fields["samples"]


def save_line_plot(
    data: pd.DataFrame,
    x_col: str,
    y_col: str,
    err_col: str,
    xlabel: str,
    ylabel: str,
    title: str,
    output_path: Path,
    x_max_ticks: int | None = None,
) -> None:
    fig, ax = plt.subplots(figsize=(8, 5), constrained_layout=True)
    ax.errorbar(
        data[x_col],
        data[y_col],
        yerr=data[err_col],
        fmt="-o",
        capsize=9,
        markersize=9,
        elinewidth=2,
    )
    ax.set_title(title, fontsize=16)
    ax.set_xlabel(xlabel, fontsize=16)
    ax.set_ylabel(ylabel, fontsize=16)
    ax.tick_params(axis="both", labelsize=14)
    ax.ticklabel_format(style="plain", axis="x")
    if x_max_ticks is not None:
        ax.xaxis.set_major_locator(MaxNLocator(nbins=x_max_ticks))
    ax.grid(True, linewidth=0.8, alpha=0.4)
    ax.set_xlim(left=0)
    ax.set_ylim(bottom=0)
    sns.despine(ax=ax)
    fig.savefig(output_path, bbox_inches="tight", dpi=600)
    plt.close(fig)


def bytes_to_mb(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df["mean_mb"] = df["mean_bytes"] / (1024**2)
    df["stddev_mb"] = df["stddev_bytes"] / (1024**2)
    return df


def main() -> None:
    PLOTS_DIR.mkdir(parents=True, exist_ok=True)

    threads_df = pd.read_csv(THREADS_CSV)
    threads_df["threads"] = threads_df["label"].apply(parse_thread_count)
    threads_df = threads_df.sort_values("threads")
    threads_df = bytes_to_mb(threads_df)
    save_line_plot(
        threads_df,
        "threads",
        "mean_s",
        "stddev_s",
        "Threads",
        "Mean Runtime (s)",
        "Runtime vs. Thread Count",
        PLOTS_DIR / "thread_count_benchmark_runtime.pdf",
    )
    save_line_plot(
        threads_df,
        "threads",
        "mean_mb",
        "stddev_mb",
        "Threads",
        "Peak RSS (MB)",
        "Peak Memory vs. Thread Count",
        PLOTS_DIR / "thread_count_benchmark_memory.pdf",
    )

    variants_df = pd.read_csv(VARIANTS_CSV)
    variants_df["variants"] = variants_df["label"].apply(parse_variant_count)
    variants_df = variants_df.sort_values("variants")
    variants_df = bytes_to_mb(variants_df)
    save_line_plot(
        variants_df,
        "variants",
        "mean_s",
        "stddev_s",
        "Variant Count",
        "Mean Runtime (s)",
        "Runtime vs. Variant Count",
        PLOTS_DIR / "variant_count_benchmark_runtime.pdf",
    )
    save_line_plot(
        variants_df,
        "variants",
        "mean_mb",
        "stddev_mb",
        "Variant Count",
        "Peak RSS (MB)",
        "Peak Memory vs. Variant Count",
        PLOTS_DIR / "variant_count_benchmark_memory.pdf",
    )

    pairs_df = pd.read_csv(PAIRS_CSV)
    pairs_df["pairs"] = pairs_df["label"].apply(parse_pair_count)
    pairs_df = pairs_df.sort_values("pairs")
    pairs_df = bytes_to_mb(pairs_df)
    save_line_plot(
        pairs_df,
        "pairs",
        "mean_s",
        "stddev_s",
        "Samples",
        "Mean Runtime (s)",
        "Runtime vs. Sample Count",
        PLOTS_DIR / "sample_count_benchmark_runtime.pdf",
    )
    save_line_plot(
        pairs_df,
        "pairs",
        "mean_mb",
        "stddev_mb",
        "Samples",
        "Peak RSS (MB)",
        "Peak Memory vs. Sample Count",
        PLOTS_DIR / "sample_count_benchmark_memory.pdf",
    )


if __name__ == "__main__":
    main()
