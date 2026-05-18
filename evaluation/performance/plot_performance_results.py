import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

from evaluation_utils import PERFORMANCE_DIR, add_panel_label

RESULTS_DIR = PERFORMANCE_DIR / "results"

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


def bytes_to_mb(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df["mean_mb"] = df["mean_bytes"] / (1024**2)
    df["stddev_mb"] = df["stddev_bytes"] / (1024**2)
    return df


def plot_line(
    ax: plt.Axes,
    data: pd.DataFrame,
    x_col: str,
    y_col: str,
    err_col: str,
    xlabel: str,
    ylabel: str,
    title: str,
    x_scale: float = 1.0,
    x_left: float | None = None,
) -> None:
    ax.errorbar(
        data[x_col] / x_scale,
        data[y_col],
        yerr=data[err_col],
        fmt="-o",
        capsize=8,
        markersize=8,
        elinewidth=2,
    )
    ax.set_title(title, fontsize=16)
    ax.set_xlabel(xlabel, fontsize=16)
    ax.set_ylabel(ylabel, fontsize=16)
    ax.tick_params(axis="both", labelsize=14)
    ax.ticklabel_format(style="plain", axis="x")
    ax.grid(True, linewidth=0.8, alpha=0.4)
    if x_left is not None:
        ax.set_xlim(left=x_left)
    y_max = (data[y_col] + data[err_col]).max()
    ax.set_ylim(bottom=0, top=y_max * 1.08)
    sns.despine(ax=ax)


def main() -> None:
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    threads_df = pd.read_csv(THREADS_CSV)
    threads_df["threads"] = threads_df["label"].apply(parse_thread_count)
    threads_df = threads_df.sort_values("threads")
    threads_df = bytes_to_mb(threads_df)

    variants_df = pd.read_csv(VARIANTS_CSV)
    variants_df["variants"] = variants_df["label"].apply(parse_variant_count)
    variants_df = variants_df.sort_values("variants")
    variants_df = bytes_to_mb(variants_df)

    pairs_df = pd.read_csv(PAIRS_CSV)
    pairs_df["pairs"] = pairs_df["label"].apply(parse_pair_count)
    pairs_df = pairs_df.sort_values("pairs")
    pairs_df = bytes_to_mb(pairs_df)

    fig, axes = plt.subplots(3, 2, figsize=(10, 14), constrained_layout=True)
    fig.get_layout_engine().set(hspace=0.04)

    plot_line(
        axes[0, 0],
        threads_df,
        "threads",
        "mean_s",
        "stddev_s",
        "Threads",
        "Mean Runtime (s)",
        "Runtime vs. Thread Count",
        x_left=-20,
    )
    plot_line(
        axes[0, 1],
        threads_df,
        "threads",
        "mean_mb",
        "stddev_mb",
        "Threads",
        "Peak RSS (MB)",
        "Peak Memory vs. Thread Count",
        x_left=-20,
    )
    plot_line(
        axes[1, 0],
        variants_df,
        "variants",
        "mean_s",
        "stddev_s",
        r"Variant Count ($\times 10^5$)",
        "Mean Runtime (s)",
        "Runtime vs. Variant Count",
        x_scale=1e5,
    )
    plot_line(
        axes[1, 1],
        variants_df,
        "variants",
        "mean_mb",
        "stddev_mb",
        r"Variant Count ($\times 10^5$)",
        "Peak RSS (MB)",
        "Peak Memory vs. Variant Count",
        x_scale=1e5,
    )
    plot_line(
        axes[2, 0], pairs_df, "pairs", "mean_s", "stddev_s", "Samples", "Mean Runtime (s)", "Runtime vs. Sample Count"
    )
    plot_line(
        axes[2, 1],
        pairs_df,
        "pairs",
        "mean_mb",
        "stddev_mb",
        "Samples",
        "Peak RSS (MB)",
        "Peak Memory vs. Sample Count",
    )

    for idx, ax in enumerate(axes.flat):
        add_panel_label(ax, chr(ord("A") + idx))

    fig.savefig(RESULTS_DIR / "performance_benchmarks.pdf", bbox_inches="tight", dpi=600)
    plt.close(fig)


if __name__ == "__main__":
    main()
