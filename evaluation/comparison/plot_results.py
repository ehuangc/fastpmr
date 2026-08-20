from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

from evaluation_utils.plot import add_panel_label

RESULTS_DIR = Path(__file__).resolve().parent / "results"
BENCHMARK_CSV = RESULTS_DIR / "comparison_benchmark.csv"
TOOL_COLORS = {"fastpmr": "tab:blue", "READv2": "#CC5500"}
TOOL_MARKERS = {"fastpmr": "o", "READv2": "^"}


def parse_tool(label: str) -> str:
    return label.split("tool=", 1)[1].split("_", 1)[0]


def parse_sample_count(label: str) -> int:
    return int(label.split("samples=", 1)[1])


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


def plot_bar(
    ax: plt.Axes,
    data: pd.DataFrame,
    y_col: str,
    err_col: str,
    ylabel: str,
    title: str,
    seconds_col: str | None = None,
    bytes_col: str | None = None,
) -> None:
    assert seconds_col or bytes_col
    # Error bars are only meaningful when the benchmark ran more than once per config
    errors = data[err_col] if (data[err_col] > 0).any() else None
    ax.bar(
        data["tool"],
        data[y_col],
        yerr=errors,
        color=[TOOL_COLORS[tool] for tool in data["tool"]],
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


def plot_runtime_vs_samples(ax: plt.Axes, data: pd.DataFrame, title: str) -> None:
    for tool, group in data.groupby("tool", sort=False):
        group = group.sort_values("samples")
        # Error bars are only meaningful when the benchmark ran more than once per sample count
        yerr = group["stddev_s"] if (group["stddev_s"] > 0).any() else None
        ax.errorbar(
            group["samples"],
            group["mean_s"],
            yerr=yerr,
            color=TOOL_COLORS[tool],
            fmt=f"-{TOOL_MARKERS[tool]}",
            capsize=8,
            markersize=8,
            elinewidth=2,
            label=tool,
        )
    # fastpmr and READv2 runtimes differ by orders of magnitude, so a linear axis would flatten fastpmr's curve
    ax.set_yscale("log")
    ax.set_title(title, fontsize=16)
    ax.set_xlabel("Samples", fontsize=16)
    ax.set_ylabel("Mean Runtime (s)", fontsize=16)
    ax.tick_params(axis="both", labelsize=14)
    ax.ticklabel_format(style="plain", axis="x")
    ax.grid(True, which="both", linewidth=0.8, alpha=0.4)
    ax.legend(fontsize=14, frameon=False)
    sns.despine(ax=ax)


def main() -> None:
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    benchmark_df = pd.read_csv(BENCHMARK_CSV)
    benchmark_df["tool"] = benchmark_df["label"].apply(parse_tool)
    benchmark_df["samples"] = benchmark_df["label"].apply(parse_sample_count)

    # The bar panels report the full dataset, which is the largest sample count benchmarked
    full_dataset_df = benchmark_df[benchmark_df["samples"] == benchmark_df["samples"].max()]
    full_dataset_df = bytes_to_gb(full_dataset_df)
    full_dataset_df = seconds_to_minutes(full_dataset_df)

    fig = plt.figure(figsize=(10, 10), constrained_layout=True)
    fig.get_layout_engine().set(wspace=0.04, hspace=0.06)
    grid = fig.add_gridspec(2, 2)
    ax_rt = fig.add_subplot(grid[0, 0])
    ax_mem = fig.add_subplot(grid[0, 1])
    # The sample count panel spans the full width of the second row
    ax_samples = fig.add_subplot(grid[1, :])

    plot_bar(
        ax_rt,
        full_dataset_df,
        "mean_min",
        "stddev_min",
        "Mean Runtime (min)",
        "Runtime: fastpmr vs. READv2",
        seconds_col="mean_s",
    )
    add_panel_label(ax_rt, "A")

    plot_bar(
        ax_mem,
        full_dataset_df,
        "mean_gb",
        "stddev_gb",
        "Peak RSS (GB)",
        "Peak Memory: fastpmr vs. READv2",
        bytes_col="mean_bytes",
    )
    add_panel_label(ax_mem, "B")

    plot_runtime_vs_samples(ax_samples, benchmark_df, "Log Runtime vs. Sample Count")
    add_panel_label(ax_samples, "C", x=-0.09)

    fig.savefig(RESULTS_DIR / "comparison_results.pdf", bbox_inches="tight", dpi=600)
    plt.close(fig)


if __name__ == "__main__":
    main()
