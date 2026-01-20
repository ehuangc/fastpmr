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


def parse_thread_count(command: str) -> int:
    return int(command.split("threads=", 1)[1])


def parse_variant_count(command: str) -> int:
    spec = command.split("variants=", 1)[1]
    if "-" in spec:
        spec = spec.split("-")[-1]
    return int(spec)


def parse_pair_count(command: str) -> int:
    fields = {}
    for part in command.split("_"):
        if "=" not in part:
            continue
        key, value = part.split("=", 1)
        fields[key] = int(value)
    return fields["pairs"]


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
    ax.set_xlabel(xlabel, fontsize=16)
    ax.set_ylabel(ylabel, fontsize=16)
    ax.set_title(title, fontsize=16)
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


def main() -> None:
    PLOTS_DIR.mkdir(parents=True, exist_ok=True)

    threads_df = pd.read_csv(THREADS_CSV)
    threads_df["threads"] = threads_df["command"].apply(parse_thread_count)
    threads_df = threads_df.sort_values("threads")
    save_line_plot(
        threads_df,
        "threads",
        "mean",
        "stddev",
        "Threads",
        "Mean Runtime (s)",
        "Runtime vs. Thread Count",
        PLOTS_DIR / "thread_count_benchmark.pdf",
    )

    variants_df = pd.read_csv(VARIANTS_CSV)
    variants_df["variants"] = variants_df["command"].apply(parse_variant_count)
    variants_df = variants_df.sort_values("variants")
    save_line_plot(
        variants_df,
        "variants",
        "mean",
        "stddev",
        "Variant Count",
        "Mean Runtime (s)",
        "Runtime vs. Variant Count",
        PLOTS_DIR / "variant_count_benchmark.pdf",
    )

    pairs_df = pd.read_csv(PAIRS_CSV)
    pairs_df["pairs"] = pairs_df["command"].apply(parse_pair_count)
    pairs_df = pairs_df.sort_values("pairs")
    save_line_plot(
        pairs_df,
        "pairs",
        "mean",
        "stddev",
        "Sample Pairs",
        "Mean Runtime (s)",
        "Runtime vs. Pair Count",
        PLOTS_DIR / "pair_count_benchmark.pdf",
        x_max_ticks=7,
    )


if __name__ == "__main__":
    main()
