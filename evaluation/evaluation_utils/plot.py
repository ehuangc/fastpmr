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


# Constrained layout packs the suptitle tight against the panel labels, so widen the padding
# between the figure edges and the panels. Passing an explicit y instead would turn off
# autopositioning, which stops the layout from reserving room for the title at all.
def add_suptitle(fig: plt.Figure, title: str, pad_inches: float = 0.12) -> None:
    fig.get_layout_engine().set(h_pad=pad_inches)
    fig.suptitle(title, fontsize=18)
