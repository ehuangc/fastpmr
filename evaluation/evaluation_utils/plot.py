import matplotlib.pyplot as plt


def add_panel_label(ax: plt.Axes, label: str, x: float = -0.18) -> None:
    ax.text(
        x,
        1.08,
        label,
        transform=ax.transAxes,
        fontsize=18,
        fontweight="bold",
        va="top",
        ha="left",
    )
