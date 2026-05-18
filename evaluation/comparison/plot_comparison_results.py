from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

RESULTS_DIR = Path(__file__).resolve().parent / "results"
PLOTS_DIR = RESULTS_DIR / "plots"
READV2_CSV = RESULTS_DIR / "readv2_comparison_benchmark.csv"


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


def main() -> None:
    PLOTS_DIR.mkdir(parents=True, exist_ok=True)

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
