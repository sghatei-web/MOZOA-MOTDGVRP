#!/usr/bin/env python3
"""
plot_paper_figures.py
======================
Reads the Excel files produced by run_all_instances_for_paper.m /
run30_mozoa_paper.m and draws every figure needed for the MOZOA paper, plus
several additional figures that are easier to produce well in Python than in
MATLAB and that add visual variety to the results section.

INPUT FILES EXPECTED (in the folder you point --data-dir at; produced by the
MATLAB side, no manual editing needed):
    mozoa_paper_allresults.xlsx          sheets: HV_mean, NDS_mean, CPU_mean
                                         (one row per instance, one column
                                         per algorithm)
    mozoa_paper_results_<instance>.xlsx  per instance, sheets: Raw_HV, Raw_NDS,
                                         Raw_CPU (one row per run) and
                                         PF_<Algorithm> (one row per Pareto
                                         point: Distance/Time/Fuel)

CORE FIGURES (requested):
    1. 3-D Pareto-front scatter plots (distance/time/fuel), several instances
    2. Hypervolume box plots, one panel per instance, grid of all instances
    3. Grouped bar charts: mean HV, mean NDS, mean CPU across all instances

BONUS FIGURES (suggested; richer/more varied than boxplots+bars alone):
    4. Violin+strip plots of HV per instance (shows full distribution shape,
       not just quartiles -- catches bimodal/skewed run distributions that a
       box plot hides)
    5. Win-rate heatmap: for every instance x algorithm pair, which algorithm
       had the best mean HV (green) vs which didn't (graded by distance from
       the winner) -- gives an at-a-glance "who wins where" summary across
       all 13 instances that no single box/bar plot can show
    6. HV-vs-instance-size trend line with a shaded IQR ribbon per algorithm
       -- shows whether an algorithm's relative advantage grows/shrinks as
       problems get larger, which the per-instance panels cannot show at a
       glance
    7. Radar (spider) chart comparing algorithms on normalised HV, NDS, and
       (inverse) CPU simultaneously -- a compact multi-criteria summary
       plot common in algorithm-comparison papers
    8. Pareto-front size (NDS) vs. hypervolume scatter, coloured by
       algorithm, one marker per (instance, algorithm, run) -- exposes the
       quality/quantity trade-off discussed in the paper's NDS-vs-C-metric
       caveat directly as a picture rather than only in two separate tables

USAGE
    python plot_paper_figures.py --data-dir /path/to/excel/files
    python plot_paper_figures.py                      # looks in .
    python plot_paper_figures.py --pf-instances 4      # only 4 PF figures
    python plot_paper_figures.py --box-instances 12    # box-plot grid size
    python plot_paper_figures.py --no-bonus            # core 3 plots only

Requires: pandas, matplotlib, seaborn, openpyxl, numpy
"""

import argparse
import glob
import os
import re
import warnings

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
import pandas as pd
import seaborn as sns

# ----------------------------------------------------------------------------
# Global styling
# ----------------------------------------------------------------------------
# The four literature competitors have fixed, well-known names. The proposed
# method's name has changed during development (ZOA-8op -> MOZOA -> MOZOA), and
# older Excel exports may still use an earlier name. Rather than hardcode one
# name and silently show nothing when it doesn't match, the actual column/
# sheet names found in your Excel files are auto-detected at runtime
# (see detect_algorithms()) and used throughout; --proposed-name lets you
# force a specific label if auto-detection ever picks the wrong column.
KNOWN_COMPETITORS = ["NSGA-2", "MLNSGA-2", "SPEA2", "MOEA/D"]
KNOWN_PROPOSED_ALIASES = ["MOZOA", "MOZOA", "ZOA-8op", "MOZOA"]

BASE_COLORS = {
    "NSGA-2": "#0073bf",
    "MLNSGA-2": "#78ab30",
    "SPEA2": "#7d2e8f",
    "MOEA/D": "#edb120",
}
PROPOSED_COLOR = "#d9541a"  # always drawn in this color, whatever its name is
BASE_MARKERS = {"NSGA-2": "s", "MLNSGA-2": "^", "SPEA2": "D", "MOEA/D": "v"}
PROPOSED_MARKER = "o"

# Populated by detect_algorithms() once the input files are read; every
# plotting function below uses ALGOS/COLORS/MARKERS (not the KNOWN_* constants
# directly), so the whole script adapts to whatever name your Excel files use.
ALGOS: list[str] = []
COLORS: dict[str, str] = {}
MARKERS: dict[str, str] = {}


def detect_algorithms(records, proposed_name_override: str | None = None) -> list[str]:
    """Read the first instance's Raw_HV sheet to find the actual algorithm
    names used in these Excel files, and set up ALGOS/COLORS/MARKERS to match
    -- so the script works whether the proposed method was exported as
    'MOZOA', 'MOZOA', 'ZOA-8op', or anything else, without editing the code."""
    global ALGOS, COLORS, MARKERS

    found_cols: list[str] = []
    for rec in records:
        try:
            df = pd.read_excel(rec["path"], sheet_name="Raw_HV")
        except (ValueError, FileNotFoundError, KeyError):
            continue
        found_cols = [c for c in df.columns if c != "Run"]
        if found_cols:
            break
    if not found_cols:
        raise RuntimeError(
            "Could not find any 'Raw_HV' sheet with algorithm columns in the "
            "given Excel files. Check --data-dir points at the right folder."
        )

    competitors_present = [c for c in KNOWN_COMPETITORS if c in found_cols]

    if proposed_name_override:
        if proposed_name_override not in found_cols:
            raise RuntimeError(
                f"--proposed-name '{proposed_name_override}' was not found as a "
                f"column in Raw_HV. Columns actually present: {found_cols}"
            )
        proposed = proposed_name_override
    else:
        proposed = next((a for a in KNOWN_PROPOSED_ALIASES if a in found_cols), None)
        if proposed is None:
            # fall back to "whatever column isn't one of the 4 known
            # competitors" -- keeps the script working even if the proposed
            # method gets renamed again in the future
            leftover = [c for c in found_cols if c not in competitors_present]
            proposed = leftover[0] if leftover else None

    algos = ([proposed] if proposed else []) + competitors_present
    # de-duplicate while preserving order, in case of odd column sets
    seen = set()
    algos = [a for a in algos if not (a in seen or seen.add(a))]

    colors, markers = {}, {}
    for a in algos:
        if a == proposed:
            colors[a] = PROPOSED_COLOR
            markers[a] = PROPOSED_MARKER
        else:
            colors[a] = BASE_COLORS.get(a, "#777777")
            markers[a] = BASE_MARKERS.get(a, "X")

    ALGOS = algos
    COLORS = colors
    MARKERS = markers

    print(f"Detected algorithm columns: {found_cols}")
    print(f"Using as proposed method: {proposed!r}")
    print(f"Plotting order: {ALGOS}\n")
    return algos


def set_paper_style(font="CMU Serif"):
    """Use a serif font matching the paper if available, else fall back."""
    available = {f.name for f in matplotlib.font_manager.fontManager.ttflist}
    for candidate in [font, "CMU Serif", "Latin Modern Roman", "Times New Roman", "DejaVu Serif"]:
        if candidate in available:
            plt.rcParams["font.family"] = candidate
            break
    else:
        plt.rcParams["font.family"] = "serif"
    plt.rcParams.update({
        "font.size": 10,
        "axes.titlesize": 11,
        "axes.labelsize": 10,
        "legend.fontsize": 9,
        "figure.dpi": 110,
        "savefig.dpi": 220,
        "axes.grid": True,
        "grid.alpha": 0.3,
    })


def sheet_tag(name: str) -> str:
    """Mirror MATLAB's safeSheetTagPaper(): strip forbidden chars, cap length."""
    s = re.sub(r"[\\/\?\*\[\]:]", "", name)
    return s[:27]


def instance_size(tag: str) -> int:
    m = re.search(r"n(\d+)", tag)
    return int(m.group(1)) if m else 0


def save_fig(fig, out_dir: str, base: str):
    fig.savefig(os.path.join(out_dir, base + ".png"), bbox_inches="tight")
    fig.savefig(os.path.join(out_dir, base + ".pdf"), bbox_inches="tight")
    plt.close(fig)


# ----------------------------------------------------------------------------
# Data loading
# ----------------------------------------------------------------------------
def discover_instance_files(data_dir: str):
    pattern = os.path.join(data_dir, "mozoa_paper_results_*.xlsx")
    files = sorted(glob.glob(pattern))
    if not files:
        raise FileNotFoundError(
            f"No mozoa_paper_results_*.xlsx files found in {data_dir!r}. "
            "Run run_all_instances_for_paper.m first, or pass --data-dir."
        )
    records = []
    for f in files:
        tag = re.sub(r"^mozoa_paper_results_|\.xlsx$", "", os.path.basename(f))
        records.append({"tag": tag, "path": f, "n": instance_size(tag)})
    records.sort(key=lambda r: r["n"])
    return records


def load_raw_hv(path: str) -> pd.DataFrame:
    df = pd.read_excel(path, sheet_name="Raw_HV")
    df = df.drop(columns=[c for c in df.columns if c == "Run"], errors="ignore")
    return df


def load_raw_nds(path: str) -> pd.DataFrame:
    return pd.read_excel(path, sheet_name="Raw_NDS").drop(
        columns=["Run"], errors="ignore"
    )


def load_raw_cpu(path: str) -> pd.DataFrame:
    return pd.read_excel(path, sheet_name="Raw_CPU").drop(
        columns=["Run"], errors="ignore"
    )


def load_pf(path: str, algo: str) -> pd.DataFrame | None:
    sn = "PF_" + sheet_tag(algo)
    try:
        return pd.read_excel(path, sheet_name=sn)
    except (ValueError, KeyError):
        return None


def load_summary(data_dir: str):
    path = os.path.join(data_dir, "mozoa_paper_allresults.xlsx")
    if not os.path.exists(path):
        warnings.warn(f"{path} not found; cross-instance bar charts will be skipped.")
        return None
    sheets = {}
    for name in ["HV_mean", "NDS_mean", "CPU_mean"]:
        try:
            sheets[name] = pd.read_excel(path, sheet_name=name)
        except ValueError:
            pass
    return sheets


# ----------------------------------------------------------------------------
# 1) 3-D Pareto-front figures
# ----------------------------------------------------------------------------
def plot_3d_pareto_fronts(records, out_dir, n_instances=6):
    from mpl_toolkits.mplot3d import Axes3D  # noqa: F401  (registers 3-D projection)

    k = min(n_instances, len(records))
    idxs = sorted(set(np.round(np.linspace(0, len(records) - 1, k)).astype(int)))

    for idx in idxs:
        rec = records[idx]
        fig = plt.figure(figsize=(7.5, 6.2))
        ax = fig.add_subplot(111, projection="3d")
        any_plotted = False
        for algo in ALGOS:
            pf = load_pf(rec["path"], algo)
            if pf is None or pf.empty:
                continue
            ax.scatter(
                pf["Distance"], pf["Time"], pf["Fuel"],
                s=55, c=COLORS[algo], marker=MARKERS[algo],
                edgecolors="k", linewidths=0.4, alpha=0.85, label=algo,
            )
            any_plotted = True
        if not any_plotted:
            plt.close(fig)
            continue
        ax.set_xlabel("Distance")
        ax.set_ylabel("Time")
        ax.set_zlabel("Fuel")
        ax.set_title(f"Pareto fronts on {rec['tag']}")
        ax.view_init(elev=22, azim=-50)
        ax.legend(loc="center left", bbox_to_anchor=(1.05, 0.5))
        save_fig(fig, out_dir, f"PF3D_{rec['tag']}")
    print(f"[1/8] 3-D Pareto-front figures: {len(idxs)} instance(s) plotted.")


# ----------------------------------------------------------------------------
# 2) HV box-plot grid across instances
# ----------------------------------------------------------------------------
def plot_hv_boxplot_grid(records, out_dir, n_box=12, seed=None):
    rng = np.random.default_rng(seed)
    k = min(n_box, len(records))
    sel = rng.choice(len(records), size=k, replace=False)
    sel = sorted(sel, key=lambda i: records[i]["n"])

    ncols = 4 if k > 9 else 3
    nrows = int(np.ceil(k / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(4.2 * ncols, 3.4 * nrows))
    axes = np.atleast_1d(axes).flatten()

    for ax, idx in zip(axes, sel):
        rec = records[idx]
        hv = load_raw_hv(rec["path"])
        cols = [a for a in ALGOS if a in hv.columns]
        data = [hv[a].dropna().values for a in cols]
        bp = ax.boxplot(data, patch_artist=True, tick_labels=cols, widths=0.6)
        for patch, a in zip(bp["boxes"], cols):
            patch.set_facecolor(COLORS[a])
            patch.set_alpha(0.55)
        ax.set_title(rec["tag"], fontsize=9)
        ax.tick_params(axis="x", labelrotation=35, labelsize=7)
        ax.tick_params(axis="y", labelsize=7)
    for ax in axes[len(sel):]:
        ax.axis("off")

    fig.suptitle("Hypervolume distribution across instances", fontsize=13)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    save_fig(fig, out_dir, "box_HV_grid")
    print(f"[2/8] HV box-plot grid: {k} instance(s).")


# ----------------------------------------------------------------------------
# 3) Grouped bar charts across all instances (HV / NDS / CPU)
# ----------------------------------------------------------------------------
def plot_grouped_bars(summary, out_dir):
    if summary is None:
        print("[3/8] Grouped bar charts: skipped (no allresults.xlsx).")
        return
    specs = [
        ("HV_mean", "Mean hypervolume", "bar_meanHV", False),
        ("NDS_mean", "Mean # non-dominated solutions", "bar_meanNDS", False),
        ("CPU_mean", "Mean CPU time (s)", "bar_meanCPU", True),
    ]
    for sheet, ylab, base, logscale in specs:
        if sheet not in summary:
            continue
        df = summary[sheet]
        inst_col = df.columns[0]
        cols = [a for a in ALGOS if a in df.columns]
        x = np.arange(len(df))
        width = 0.8 / len(cols)
        fig, ax = plt.subplots(figsize=(10, 4.4))
        for i, a in enumerate(cols):
            ax.bar(x + (i - len(cols) / 2 + 0.5) * width, df[a], width,
                   color=COLORS[a], label=a, edgecolor="black", linewidth=0.3)
        ax.set_xticks(x)
        ax.set_xticklabels(df[inst_col], rotation=30, ha="right")
        if logscale:
            ax.set_yscale("log")
        ax.set_ylabel(ylab)
        ax.set_title(ylab)
        ax.legend(loc="best", ncol=len(cols))
        save_fig(fig, out_dir, base)
    print("[3/8] Grouped bar charts: HV, NDS, CPU.")


# ----------------------------------------------------------------------------
# 4) BONUS: violin + strip plots of HV per instance
# ----------------------------------------------------------------------------
def plot_hv_violin_grid(records, out_dir, n_panels=8, seed=None):
    rng = np.random.default_rng(seed)
    k = min(n_panels, len(records))
    sel = rng.choice(len(records), size=k, replace=False)
    sel = sorted(sel, key=lambda i: records[i]["n"])

    ncols = 4 if k > 6 else 3
    nrows = int(np.ceil(k / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(4.2 * ncols, 3.6 * nrows))
    axes = np.atleast_1d(axes).flatten()

    for ax, idx in zip(axes, sel):
        rec = records[idx]
        hv = load_raw_hv(rec["path"])
        cols = [a for a in ALGOS if a in hv.columns]
        long = hv[cols].melt(var_name="Algorithm", value_name="HV")
        sns.violinplot(data=long, x="Algorithm", y="HV", hue="Algorithm", ax=ax,
                        palette=COLORS, inner=None, cut=0, linewidth=0.8, legend=False)
        sns.stripplot(data=long, x="Algorithm", y="HV", ax=ax,
                       color="black", size=2.5, alpha=0.6, jitter=0.15)
        ax.set_title(rec["tag"], fontsize=9)
        ax.set_xlabel("")
        ax.tick_params(axis="x", labelrotation=35, labelsize=7)
        ax.tick_params(axis="y", labelsize=7)
    for ax in axes[len(sel):]:
        ax.axis("off")

    fig.suptitle("Hypervolume distribution shape (violin + individual runs)", fontsize=13)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    save_fig(fig, out_dir, "violin_HV_grid")
    print(f"[4/8] BONUS violin+strip HV grid: {k} instance(s).")


# ----------------------------------------------------------------------------
# 5) BONUS: win-rate heatmap (best algorithm per instance)
# ----------------------------------------------------------------------------
def plot_winrate_heatmap(summary, out_dir):
    if summary is None or "HV_mean" not in summary:
        print("[5/8] BONUS win-rate heatmap: skipped (no allresults.xlsx).")
        return
    df = summary["HV_mean"].copy()
    inst_col = df.columns[0]
    cols = [a for a in ALGOS if a in df.columns]
    mat = df[cols].values
    # normalise each row (instance) to [0,1] so the heatmap is comparable
    # across instances of very different absolute HV scale
    row_max = mat.max(axis=1, keepdims=True)
    row_min = mat.min(axis=1, keepdims=True)
    norm = (mat - row_min) / np.where(row_max - row_min == 0, 1, row_max - row_min)

    fig, ax = plt.subplots(figsize=(1.3 * len(cols) + 2, 0.42 * len(df) + 2))
    sns.heatmap(norm, annot=mat, fmt=".3g", cmap="RdYlGn", cbar_kws={"label": "Relative HV (row-normalised)"},
                xticklabels=cols, yticklabels=df[inst_col], ax=ax, linewidths=0.5, linecolor="white")
    ax.set_title("Mean hypervolume per instance (green = best on that row)")
    save_fig(fig, out_dir, "heatmap_HV_winrate")
    print("[5/8] BONUS win-rate heatmap.")


# ----------------------------------------------------------------------------
# 6) BONUS: HV vs instance size trend line with IQR ribbon
# ----------------------------------------------------------------------------
def plot_hv_vs_size_trend(records, out_dir):
    rows = []
    for rec in records:
        hv = load_raw_hv(rec["path"])
        for a in ALGOS:
            if a not in hv.columns:
                continue
            col = hv[a].dropna()
            if col.empty:
                continue
            rows.append({
                "n": rec["n"], "tag": rec["tag"], "algo": a,
                "q1": col.quantile(0.25), "median": col.median(), "q3": col.quantile(0.75),
                # min-max normalise within instance so different-scale
                # instances can share one axis and show *relative* trend
            })
    df = pd.DataFrame(rows)
    if df.empty:
        print("[6/8] BONUS HV-vs-size trend: skipped (no data).")
        return
    # normalise HV within each instance (divide by that instance's max median)
    df["norm_median"] = df.groupby("tag")["median"].transform(lambda s: s / s.max())
    df["norm_q1"] = df["q1"] / df.groupby("tag")["median"].transform("max")
    df["norm_q3"] = df["q3"] / df.groupby("tag")["median"].transform("max")

    fig, ax = plt.subplots(figsize=(9, 5))
    for a in ALGOS:
        sub = df[df["algo"] == a].sort_values("n")
        if sub.empty:
            continue
        ax.plot(sub["n"], sub["norm_median"], "-o", color=COLORS[a], label=a, markersize=5)
        ax.fill_between(sub["n"], sub["norm_q1"], sub["norm_q3"], color=COLORS[a], alpha=0.15)
    ax.set_xlabel("Instance size (number of customers, n)")
    ax.set_ylabel("Hypervolume relative to the best algorithm on that instance")
    ax.set_title("Relative hypervolume vs. instance size (median, IQR band)")
    ax.legend(loc="best")
    save_fig(fig, out_dir, "trend_HV_vs_size")
    print("[6/8] BONUS HV-vs-instance-size trend line.")


# ----------------------------------------------------------------------------
# 7) BONUS: radar / spider chart of normalised HV, NDS, 1/CPU
# ----------------------------------------------------------------------------
def plot_radar_chart(summary, out_dir):
    if summary is None or not all(k in summary for k in ["HV_mean", "NDS_mean", "CPU_mean"]):
        print("[7/8] BONUS radar chart: skipped (need HV/NDS/CPU sheets).")
        return
    cols = [a for a in ALGOS if a in summary["HV_mean"].columns]

    def agg(sheet, invert=False):
        vals = summary[sheet][cols].mean(axis=0).values.astype(float)
        if invert:
            vals = 1.0 / np.where(vals == 0, np.nan, vals)
        vmin, vmax = np.nanmin(vals), np.nanmax(vals)
        return (vals - vmin) / (vmax - vmin) if vmax > vmin else np.ones_like(vals) * 0.5

    hv_n = agg("HV_mean")
    nds_n = agg("NDS_mean")
    cpu_n = agg("CPU_mean", invert=True)  # invert: faster (lower CPU) -> higher score

    categories = ["Hypervolume", "Pareto-front size", "Speed (1/CPU)"]
    n_cat = len(categories)
    angles = np.linspace(0, 2 * np.pi, n_cat, endpoint=False).tolist()
    angles += angles[:1]

    fig, ax = plt.subplots(figsize=(6.5, 6.5), subplot_kw=dict(polar=True))
    for i, a in enumerate(cols):
        vals = [hv_n[i], nds_n[i], cpu_n[i]]
        vals += vals[:1]
        ax.plot(angles, vals, color=COLORS[a], linewidth=2, label=a, marker="o")
        ax.fill(angles, vals, color=COLORS[a], alpha=0.08)
    ax.set_xticks(angles[:-1])
    ax.set_xticklabels(categories)
    ax.set_yticklabels([])
    ax.set_title("Multi-criteria comparison (normalised across algorithms, averaged over instances)",
                 y=1.08)
    ax.legend(loc="upper right", bbox_to_anchor=(1.3, 1.1))
    save_fig(fig, out_dir, "radar_multicriteria")
    print("[7/8] BONUS radar/spider chart.")


# ----------------------------------------------------------------------------
# 8) BONUS: NDS vs HV scatter (quality/quantity trade-off), all runs pooled
# ----------------------------------------------------------------------------
def plot_nds_vs_hv_scatter(records, out_dir):
    rows = []
    for rec in records:
        hv = load_raw_hv(rec["path"])
        nds = load_raw_nds(rec["path"])
        for a in ALGOS:
            if a not in hv.columns or a not in nds.columns:
                continue
            for h, d in zip(hv[a], nds[a]):
                rows.append({"algo": a, "HV": h, "NDS": d, "instance": rec["tag"]})
    df = pd.DataFrame(rows)
    if df.empty:
        print("[8/8] BONUS NDS-vs-HV scatter: skipped (no data).")
        return

    fig, ax = plt.subplots(figsize=(8, 6))
    for a in ALGOS:
        sub = df[df["algo"] == a]
        if sub.empty:
            continue
        ax.scatter(sub["NDS"], sub["HV"], s=22, alpha=0.5, color=COLORS[a],
                   marker=MARKERS[a], label=a, edgecolors="none")
    ax.set_xlabel("Number of non-dominated solutions (NDS)")
    ax.set_ylabel("Hypervolume (HV)")
    ax.set_title("Quality (HV) vs. quantity (NDS) across all runs and instances")
    ax.legend(loc="best")
    save_fig(fig, out_dir, "scatter_NDS_vs_HV")
    print("[8/8] BONUS NDS-vs-HV trade-off scatter.")


# ----------------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------------
def _run_step(step_name, fn, *fn_args, **fn_kwargs):
    """Run one plotting function; on error, print it and continue instead of
    aborting the whole script (so one broken figure never hides the rest)."""
    try:
        fn(*fn_args, **fn_kwargs)
    except Exception as exc:  # noqa: BLE001 -- intentionally broad: we want to
        # keep going no matter what kind of error a single figure hits
        print(f"  [FAILED] {step_name}: {type(exc).__name__}: {exc}")


def main():
    p = argparse.ArgumentParser(description="Draw paper figures from Excel results.")
    p.add_argument("--data-dir", default=".", help="Folder containing the Excel files.")
    p.add_argument("--out-dir", default="figures_paper_python", help="Output folder for figures.")
    p.add_argument("--font", default="CMU Serif", help="Preferred font (falls back automatically).")
    p.add_argument("--pf-instances", type=int, default=6, help="How many instances get a 3-D PF figure.")
    p.add_argument("--box-instances", type=int, default=12, help="How many instances in the HV box-plot grid.")
    p.add_argument("--violin-instances", type=int, default=8, help="How many instances in the violin grid.")
    p.add_argument("--seed", type=int, default=None, help="Random seed for instance sub-selection (None = random each run).")
    p.add_argument("--no-bonus", action="store_true", help="Skip the 5 bonus figures; core 3 only.")
    p.add_argument("--proposed-name", default=None,
                    help="Force the exact column/sheet name used for the proposed "
                         "method in your Excel files (e.g. 'MOZOA', 'ZOA-8op', "
                         "'MOZOA'), overriding auto-detection.")
    args = p.parse_args()

    set_paper_style(args.font)
    os.makedirs(args.out_dir, exist_ok=True)

    records = discover_instance_files(args.data_dir)
    detect_algorithms(records, proposed_name_override=args.proposed_name)
    summary = load_summary(args.data_dir)

    print(f"Found {len(records)} instance result file(s) in {args.data_dir!r}.")

    # core 3 -- each wrapped so a failure in one doesn't hide the others
    _run_step("3-D Pareto-front figures", plot_3d_pareto_fronts,
               records, args.out_dir, n_instances=args.pf_instances)
    _run_step("HV box-plot grid", plot_hv_boxplot_grid,
               records, args.out_dir, n_box=args.box_instances, seed=args.seed)
    _run_step("Grouped bar charts", plot_grouped_bars, summary, args.out_dir)

    # bonus 5
    if not args.no_bonus:
        _run_step("BONUS violin+strip HV grid", plot_hv_violin_grid,
                   records, args.out_dir, n_panels=args.violin_instances, seed=args.seed)
        _run_step("BONUS win-rate heatmap", plot_winrate_heatmap, summary, args.out_dir)
        _run_step("BONUS HV-vs-size trend", plot_hv_vs_size_trend, records, args.out_dir)
        _run_step("BONUS radar chart", plot_radar_chart, summary, args.out_dir)
        _run_step("BONUS NDS-vs-HV scatter", plot_nds_vs_hv_scatter, records, args.out_dir)

    print(f"\nAll figures saved to {args.out_dir!r}/")


if __name__ == "__main__":
    main()
