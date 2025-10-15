import numpy as np # type: ignore
import matplotlib.pyplot as plt # type: ignore
import pandas as pd
import matplotlib.cm as cm

from mpl_toolkits.mplot3d import Axes3D
from matplotlib.colors import Normalize
from matplotlib.ticker import FuncFormatter
from matplotlib.ticker import ScalarFormatter
from matplotlib.patches import Rectangle
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib import cbook
# import imageio
from scipy.ndimage import gaussian_filter


font = {
        'weight' : 'bold',
        'size'   : 16}
plt.rc('text', usetex=True)
plt.rc('text.latex', preamble=r'\usepackage{amsmath}\usepackage{newtxtext,newtxmath}')
plt.rc('font', **font) 


def plot_phase_diagram(measure, U_vals, V_vals, x_label, y_label, measure_label, fname=None, 
                       x_min=0.0, y_min=0.0, theme="viridis", num_ticks=5, upper_bound=None, highlight=False):

    fig, ax1 = plt.subplots(1, 1, figsize=(8, 6))

    sigma_smoothing = 1.5
    num_contour_levels = 20

    U_min = x_min
    U_max = np.max(U_vals)
    V_min = y_min
    V_max = np.max(V_vals)

    # --- Mask and smoothing ---
    mask_U = (U_vals >= U_min) & (U_vals <= U_max)
    mask_V = (V_vals >= V_min) & (V_vals <= V_max)

    U_vals_new = U_vals[mask_U]
    V_vals_new = V_vals[mask_V]

    masked_measure = measure[np.ix_(mask_V, mask_U)]
    smoothed_measure = gaussian_filter(masked_measure, sigma=sigma_smoothing)

    X, Y = np.meshgrid(U_vals_new, V_vals_new)

    # --- Normalization ---
    vmin = np.nanmin(smoothed_measure)
    vmax = np.nanmax(smoothed_measure) if upper_bound is None else upper_bound
    norm = Normalize(vmin=vmin, vmax=vmax)

    # --- Main image ---
    im0 = ax1.imshow(
        smoothed_measure,
        origin="lower",
        extent=[U_min, U_max, V_min, V_max],
        aspect="auto",
        cmap=theme,
        norm=norm,
    )

    ax1.tick_params(axis='x', labelsize=18)
    ax1.tick_params(axis='y', labelsize=18)

    ax1.set_xlabel(x_label, fontsize=22)
    # ax1.set_ylabel(y_label, fontsize=22)

    # --- COLORBAR ON TOP ---
    pos = ax1.get_position()
    cax = fig.add_axes([pos.x0, -0.06, pos.width, 0.03])

    ticks = np.linspace(vmin, vmax, num_ticks)

    cbar = fig.colorbar(
        im0,
        cax=cax,
        orientation="horizontal",
        ticks=ticks,
        boundaries=np.linspace(vmin, vmax, 256)
    )

    cbar.ax.xaxis.set_ticks_position("bottom")

    # --- Fixed decimal formatting for clean tick labels ---
    def clean_tick(x, _):
        # format to one decimal place and remove trailing .0 if unnecessary
        s = f"{x:.1f}"
        if s.endswith(".0"):
            s = s[:-2]
        return f"${s}$"

    cbar.ax.xaxis.set_major_formatter(FuncFormatter(clean_tick))
    cbar.ax.tick_params(labelsize=14)

    # --- Contours ---
    ax1.contour(X, Y, smoothed_measure, levels=num_contour_levels, colors='black', linewidths=0.5)

    if (x_min == np.min(U_vals)) and (y_min == np.min(V_vals)) and highlight:
        ax1.hlines(y=0.0, xmin=0.0, xmax=np.max(U_vals), color='w', linestyle=':', linewidth=1.5)
        ax1.vlines(x=0.0, ymin=0.0, ymax=np.max(V_vals), color='w', linestyle=':', linewidth=1.5)
    
    if(fname != None):
        fig.savefig(fname, dpi=300, bbox_inches='tight')
    plt.show()

def plot_local_densities(measure, U_vals, V_vals, U_target, V_target, y_label, fname=None, y_min=None, y_max=None):

    i_U = np.argmin(np.abs(U_vals - U_target))
    i_V = np.argmin(np.abs(V_vals - V_target))

    magnetization =  measure[i_U, i_V] 

    if(y_min==None): y_min = np.min(magnetization)
    if(y_max==None): y_max = np.max(magnetization)

    fig, ax = plt.subplots(1, 1, figsize=(6, 4), sharex=True)

    # Connect points with lines
    ax.plot(range(len(magnetization)), magnetization, '-o', markersize=5, linewidth=1.5)

    # Labels and Titles
    ax.set_xlabel(r'Site $j$')
    ax.set_ylabel(y_label, fontsize=22)
    ax.set_ylim(y_min, y_max)

    # Adjust spacing
    fig.subplots_adjust(wspace=0.5)
    if(fname != None):
        plt.savefig(fname, dpi=300, bbox_inches='tight')
    plt.show()

def plot_cut_phase_diagram_dual(
    x_vals, y_vals,
    correlation_1, correlation_2,
    coherence_1, coherence_2,
    target, x_label, y_labels, fname,
    title=None, phase_boundaries=None,
    fixed_axis="y", curve_labels=None,
    x_min=None, x_max=None,
    y_min=None, y_max=None
):
    """
    Plot 1D cuts of four 2D measures (two correlated, two coherent)
    with phase-color bands and masking for chosen x/y limits.
    """

    # ==== Handle limits and masks ====
    if x_min is None: x_min = np.min(x_vals)
    if x_max is None: x_max = np.max(x_vals)
    if y_min is None: y_min = np.min(y_vals)
    if y_max is None: y_max = np.max(y_vals)

    mask_x = (x_vals >= x_min) & (x_vals <= x_max)
    mask_y = (y_vals >= y_min) & (y_vals <= y_max)

    # Apply masks to all data
    x_vals_masked = x_vals[mask_x]
    y_vals_masked = y_vals[mask_y]

    correlation_1_masked = correlation_1[np.ix_(mask_x, mask_y)] if fixed_axis == "x" else correlation_1[np.ix_(mask_y, mask_x)]
    correlation_2_masked = correlation_2[np.ix_(mask_x, mask_y)] if fixed_axis == "x" else correlation_2[np.ix_(mask_y, mask_x)]
    coherence_1_masked   = coherence_1[np.ix_(mask_x, mask_y)] if fixed_axis == "x" else coherence_1[np.ix_(mask_y, mask_x)]
    coherence_2_masked   = coherence_2[np.ix_(mask_x, mask_y)] if fixed_axis == "x" else coherence_2[np.ix_(mask_y, mask_x)]

    # ==== Outer figure ====
    fig_big = plt.figure(figsize=(8, 6))
    inner_pos = [0.12, 0.15, 0.75, 0.6]
    ax_left = fig_big.add_axes(inner_pos)
    ax_right = ax_left.twinx()

    # ==== Helper for slicing ====
    def get_slice(measure, target, fixed_axis):
        if fixed_axis == "y":
            i_fixed = np.argmin(np.abs(y_vals_masked - target))
            return x_vals_masked, measure[:, i_fixed]
        elif fixed_axis == "x":
            i_fixed = np.argmin(np.abs(x_vals_masked - target))
            return y_vals_masked, measure[i_fixed, :]

    # ==== Plot correlation (left y-axis) ====
    x, y1 = get_slice(correlation_1_masked, target, fixed_axis)
    _, y2 = get_slice(correlation_2_masked, target, fixed_axis)
    ax_left.plot(x, y1, '-', color="#bd3053", alpha=1, linewidth=2, label=f"{curve_labels[0]}")
    ax_left.plot(x, y2, '--', color="#bd3053", alpha=0.85, linewidth=1.6, label=f"{curve_labels[1]}")
    ax_left.set_ylim(0.0, 1.0)

    # ==== Plot coherence (right y-axis) ====
    x, y3 = get_slice(coherence_1_masked, target, fixed_axis)
    _, y4 = get_slice(coherence_2_masked, target, fixed_axis)
    ax_right.plot(x, y3, '-', color="#26808e", alpha=1, linewidth=2, label=f"{curve_labels[2]}")
    ax_right.plot(x, y4, '--', color="#26808e", alpha=0.85, linewidth=1.6, label=f"{curve_labels[3]}")

    # ax_right.set_ylim(0.0, max(np.max(coherence_1_masked), np.max(coherence_2_masked)))

    ax_left.set_xlim(x_min, x_max)

    # ==== Phase band ABOVE plot ====
    if phase_boundaries:
        numeric_boundaries = [(val, name) for val, name in phase_boundaries if val is not None]
        numeric_boundaries = sorted(numeric_boundaries, key=lambda x: x[0])

        x_min_, x_max_ = ax_left.get_xlim()
        regions = [x_min_] + [b[0] for b in numeric_boundaries] + [x_max_]
        phase_names = [b[1] for b in phase_boundaries]

        colors = plt.cm.tab20c(np.linspace(0, 1, len(phase_names)))

        fig_big.canvas.draw()
        band_bottom = inner_pos[1] + inner_pos[3]
        band_height = 0.05

        for i, name in enumerate(phase_names):
            x0 = inner_pos[0] + (regions[i] - x_min_) / (x_max_ - x_min_) * inner_pos[2]
            x1 = inner_pos[0] + (regions[i+1] - x_min_) / (x_max_ - x_min_) * inner_pos[2]
            width = x1 - x0
            rect = Rectangle(
                (x0, band_bottom), width, band_height,
                transform=fig_big.transFigure,
                color=colors[i], alpha=0.8, ec="black", lw=1
            )
            fig_big.patches.append(rect)
            fig_big.text(
                x0 + width/2, band_bottom + band_height/2,
                name, ha="center", va="center",
                fontsize=12, weight="bold", color="black"
            )

    # ==== Labels and legend ====
    ax_left.set_xlabel(x_label, fontsize=16)
    ax_left.set_ylabel(y_labels[0], fontsize=13, color='tab:red')
    ax_right.set_ylabel(y_labels[1], fontsize=13, color='tab:blue')

    lines_left, labels_left = ax_left.get_legend_handles_labels()
    lines_right, labels_right = ax_right.get_legend_handles_labels()
    all_lines = lines_left + lines_right
    all_labels = labels_left + labels_right

    fig_big.legend(
        all_lines, all_labels,
        fontsize=14, loc="lower center",
        ncol=4, frameon=False,
        bbox_to_anchor=(0.5, -0.05)
    )

    if title:
        ax_left.set_title(title, fontsize=14)

    fig_big.savefig(fname, dpi=300, bbox_inches="tight")
    plt.show()