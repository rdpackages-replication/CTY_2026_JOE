from __future__ import annotations

import os
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
os.environ.setdefault("MPLCONFIGDIR", str(SCRIPT_DIR / ".matplotlib-cache"))

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Arc
from scipy.spatial.distance import cdist

OUTPUT_DIR = SCRIPT_DIR / "output"
FIGURES_DIR = SCRIPT_DIR / "figures"
REQUIRED_COLUMNS = [
    "running_saber11",
    "running_sisben",
    "eligible_spp",
    "beneficiary_spp",
    "spadies_any",
    "icfes_educm1",
]


def clean_figures_dir() -> None:
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    for old in FIGURES_DIR.glob("*.png"):
        old.unlink()


def load_spp_data(path: Path = SCRIPT_DIR / "spp.csv") -> pd.DataFrame:
    raw = pd.read_csv(path)
    missing = [name for name in REQUIRED_COLUMNS if name not in raw.columns]
    if missing:
        raise ValueError(f"spp.csv is missing required columns: {', '.join(missing)}")
    dat = raw.loc[:, ["running_saber11", "running_sisben", "spadies_any", "eligible_spp"]].copy()
    dat.columns = ["x.1", "x.2", "y", "d"]
    dat = dat.dropna().copy()
    expected_assignment = ((dat["x.1"] >= 0) & (dat["x.2"] >= 0)).astype(int)
    if not np.all(dat["d"].to_numpy() == expected_assignment.to_numpy()):
        raise ValueError("eligible_spp does not match the quadrant assignment rule.")
    dat["y"] = dat["y"].astype(float)
    dat["d"] = dat["d"].astype(float)
    return dat


def make_eval_grid() -> pd.DataFrame:
    rows = []
    for i in range(1, 21):
        rows.append((0.0, 40 - (i - 1) * 40 / 20))
    for i in range(21, 41):
        rows.append(((i - 20 - 1) * 56 / 20, 0.0))
    return pd.DataFrame(rows, columns=["x.1", "x.2"]).iloc[10:31].reset_index(drop=True)


def scale_data(dat: pd.DataFrame, eval_grid: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    out = dat.copy()
    ev = eval_grid.copy()
    scale2 = out["x.1"].std(ddof=1) / out["x.2"].std(ddof=1)
    out["x.2"] *= scale2
    ev["x.2"] *= scale2
    return out, ev


def signed_distances(dat: pd.DataFrame, eval_grid: pd.DataFrame) -> np.ndarray:
    x = dat[["x.1", "x.2"]].to_numpy()
    b = eval_grid[["x.1", "x.2"]].to_numpy()
    sign = (2 * dat["d"].to_numpy() - 1).reshape(-1, 1)
    return cdist(x, b, metric="euclidean") * sign


def save_png(fig: plt.Figure, file_name: str, width: float = 6, height: float = 5) -> None:
    fig.set_size_inches(width, height)
    fig.tight_layout()
    fig.savefig(FIGURES_DIR / file_name, dpi=300)
    plt.close(fig)


def style_axes(ax: plt.Axes) -> None:
    ax.grid(False)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def theta_function(x: np.ndarray) -> np.ndarray:
    x = np.asarray(x, dtype=float)
    out = np.zeros_like(x)
    left = x <= 0.5
    out[left] = 2 / np.pi * x[left]
    right_x = x[~left]
    out[~left] = (right_x + 0.5) / (np.pi - np.arccos(1 / (2 * right_x)))
    return out


def save_fig1_b() -> None:
    point_x = np.array([0.0, 0.2, 0.4, 0.5, 0.6, 0.8, 1.0])
    point_y = theta_function(point_x)
    xs = np.linspace(0, 1, 1000)
    fig, ax = plt.subplots()
    ax.plot(xs, theta_function(xs), color="dimgray", lw=1.2)
    ax.scatter(point_x, point_y, color="blue", s=24, zorder=3)
    for x, y in zip(point_x, point_y):
        ax.plot([x, x], [0, y], color="lightgray", ls="--", lw=0.8)
        ax.plot([0, x], [y, y], color="lightgray", ls="--", lw=0.8)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 0.75)
    ax.set_xticks(point_x, ["0", "$r_1$", "$r_2$", "$r_3$", "$r_4$", "$r_5$", "$r_6$"])
    ax.set_yticks(point_y, ["$\\theta_{1,\\mathbf{b}}(0)$", "$\\theta_{1,\\mathbf{b}}(r_1)$", "$\\theta_{1,\\mathbf{b}}(r_2)$", "$\\theta_{1,\\mathbf{b}}(r_3)$", "$\\theta_{1,\\mathbf{b}}(r_4)$", "$\\theta_{1,\\mathbf{b}}(r_5)$", "$\\theta_{1,\\mathbf{b}}(r_6)$"])
    style_axes(ax)
    save_png(fig, "fig1-b.png")


def save_fig1_a() -> None:
    cx, cy = 0.5, 0.0
    r_vec = np.array([0.2, 0.4, 0.5, 0.6, 0.8, 1.0])
    theta_vec = np.array([7, 5, 4, 3, 2, 1]) * np.pi / 8
    theta_end_vec = np.array([np.pi, np.pi, np.pi, np.pi - np.arccos(0.5 / 0.6), np.pi - np.arccos(0.5 / 0.8), np.pi - np.arccos(0.5 / 1.0)])
    fig, ax = plt.subplots()
    ax.axvline(0, color="gray", alpha=0.5, lw=3)
    ax.axhline(0, color="gray", alpha=0.5, lw=3)
    ax.scatter([cx], [cy], color="black", s=24)
    ax.text(cx, cy - 0.08, "$\\mathbf{b}$", ha="center", va="top", fontsize=12)
    for i, r in enumerate(r_vec):
        theta = np.linspace(0, theta_end_vec[i], 200)
        ax.plot(cx + r * np.cos(theta), cy + r * np.sin(theta), color="blue", lw=1.0)
        ex = cx + r * np.cos(theta_vec[i])
        ey = cy + r * np.sin(theta_vec[i])
        ax.annotate("", xy=(ex, ey), xytext=(cx, cy), arrowprops={"arrowstyle": "->", "lw": 0.9, "color": "black"})
        ax.text(ex + 0.04, ey + 0.03, f"$r_{i + 1}$", fontsize=10)
    ax.set_xlim(-0.15, 1.55)
    ax.set_ylim(-0.1, 1.05)
    ax.set_aspect("equal")
    ax.set_xlabel("$x_1$")
    ax.set_ylabel("$x_2$")
    style_axes(ax)
    save_png(fig, "fig1-a.png")


def save_scatter(data: pd.DataFrame, eval_grid: pd.DataFrame) -> None:
    fig, ax = plt.subplots()
    control = data["d"] == 0
    treated = data["d"] == 1
    ax.scatter(data.loc[control, "x.1"], data.loc[control, "x.2"], s=3, c="indianred", marker="s", alpha=0.18, label="Control")
    ax.scatter(data.loc[treated, "x.1"], data.loc[treated, "x.2"], s=3, c="#104e8b", marker="o", alpha=0.18, label="Treatment")
    ax.plot([0, 0], [0, 55], color="0.35", lw=1.0, label="Boundary")
    ax.plot([0, 80], [0, 0], color="0.35", lw=1.0)
    ax.scatter(eval_grid["x.1"], eval_grid["x.2"], s=8, c="black", zorder=4)
    ax.set_xlim(-80, 100)
    ax.set_ylim(-40, 60)
    ax.set_xlabel("Saber 11")
    ax.set_ylabel("Sisben")
    ax.legend(loc="upper left", bbox_to_anchor=(0.10, 0.98), frameon=True)
    style_axes(ax)
    save_png(fig, "fig2-a.png")


def read_empapp_output(stem: str) -> pd.DataFrame:
    path = OUTPUT_DIR / f"empapp_{stem}.csv"
    if not path.exists():
        raise FileNotFoundError(f"Missing output file: {path}. Run CTY_2026_JOE--empapp.py first.")
    tab = pd.read_csv(path)
    required = ["row", "h0", "h1", "estimate.p", "ci.lower", "ci.upper", "cb.lower", "cb.upper", "estimate.q", "t.value", "p.value"]
    missing = [name for name in required if name not in tab.columns]
    if missing:
        raise ValueError(f"{path} is missing columns: {', '.join(missing)}")
    return tab.loc[~tab["row"].astype(str).isin(["WBATE", "LBATE"])].reset_index(drop=True)


def save_fig4(tab: pd.DataFrame, file_name: str) -> None:
    x = np.arange(1, len(tab) + 1)
    band_lower = tab["cb.lower"].where(np.isfinite(tab["cb.lower"]), tab["ci.lower"])
    band_upper = tab["cb.upper"].where(np.isfinite(tab["cb.upper"]), tab["ci.upper"])
    fig, ax = plt.subplots()
    ax.axvline(11, color="0.9", lw=0.45)
    ax.fill_between(x, band_lower, band_upper, color="#104e8b", alpha=0.14, label="CB")
    ax.vlines(x, tab["ci.lower"], tab["ci.upper"], color="black", linewidth=0.35, label="CI")
    ax.scatter(x, tab["estimate.p"], s=10, color="black", label="BATEC", zorder=3)
    ax.set_xticks([1, 4, 7, 11, 14, 17, 21], [f"$\\mathbf{{b}}_{{{i}}}$" for i in [1, 4, 7, 11, 14, 17, 21]])
    ax.set_xlim(1, len(tab))
    ax.set_ylim(-0.3, 0.8)
    ax.set_xlabel("Cutoffs on the Boundary")
    ax.set_ylabel("Treatment Effect")
    ax.legend(loc="lower right", frameon=True)
    style_axes(ax)
    save_png(fig, file_name)


def binned_means(x: np.ndarray, y: np.ndarray, bins_each_side: int = 80) -> pd.DataFrame:
    rows = []
    for lo, hi in [(np.nanmin(x), 0.0), (0.0, np.nanmax(x))]:
        if not np.isfinite(lo) or not np.isfinite(hi) or np.isclose(lo, hi):
            continue
        edges = np.linspace(lo, hi, bins_each_side + 1)
        idx = np.digitize(x, edges) - 1
        for j in range(bins_each_side):
            mask = idx == j
            if np.any(mask):
                rows.append({"x": float(np.mean(x[mask])), "y": float(np.mean(y[mask]))})
    return pd.DataFrame(rows).sort_values("x")


def poly_curve(x: np.ndarray, y: np.ndarray, side: str) -> tuple[np.ndarray, np.ndarray]:
    mask = x < 0 if side == "left" else x >= 0
    x_side = x[mask]
    y_side = y[mask]
    if len(x_side) < 4:
        return np.array([]), np.array([])
    deg = min(3, len(np.unique(x_side)) - 1)
    coefs = np.polyfit(x_side, y_side, deg)
    grid = np.linspace(float(np.min(x_side)), -0.01 if side == "left" else float(np.max(x_side)), 250)
    if side == "right":
        grid = np.linspace(0.01, float(np.max(x_side)), 250)
    return grid, np.polyval(coefs, grid)


def save_fig2_detail(y: np.ndarray, distance: np.ndarray, h: float, idx: int, suffix: str) -> None:
    x = distance[:, idx - 1]
    bins = binned_means(x, y)
    in_band = np.abs(bins["x"]) <= h
    colors = np.where(~in_band, "0.65", np.where(bins["x"] < 0, "indianred", "#104e8b"))
    fig, ax = plt.subplots()
    ax.scatter(bins["x"], bins["y"], s=7, c=colors, alpha=0.9)
    for side in ["left", "right"]:
        gx, gy = poly_curve(x, y, side)
        if gx.size:
            ax.plot(gx, gy, color="black", lw=1.1)
    ax.axvline(-h, color="indianred", lw=0.8)
    ax.axvline(h, color="#104e8b", lw=0.8)
    ax.axvline(0, color="black", lw=0.8, ls="--")
    ax.annotate("", xy=(-h, 0.90), xytext=(0, 0.90), arrowprops={"arrowstyle": "<->", "lw": 1.0, "color": "black"})
    ax.annotate("", xy=(0, 0.90), xytext=(h, 0.90), arrowprops={"arrowstyle": "<->", "lw": 1.0, "color": "black"})
    ax.text(-h / 2, 0.94, f"$\\hat{{h}}_{{\\mathbf{{b}}_{{{idx}}}}}$", ha="center", fontsize=11)
    ax.text(h / 2, 0.94, f"$\\hat{{h}}_{{\\mathbf{{b}}_{{{idx}}}}}$", ha="center", fontsize=11)
    ax.set_xlim(-100, 100)
    ax.set_ylim(0, 1)
    ax.set_xlabel(f"Distance to Boundary Point $\\mathbf{{b}}_{{{idx}}}$")
    ax.set_ylabel("College Enrollment")
    style_axes(ax)
    save_png(fig, f"fig2-{suffix}.png")


def check_generated_figures(expected: list[str]) -> None:
    missing = [name for name in expected if not (FIGURES_DIR / name).exists()]
    if missing:
        raise FileNotFoundError(f"Missing expected figures: {', '.join(missing)}")
    too_small = [name for name in expected if (FIGURES_DIR / name).stat().st_size < 1000]
    if too_small:
        raise ValueError(f"Generated figure file(s) look empty: {', '.join(too_small)}")


def main() -> None:
    clean_figures_dir()
    data = load_spp_data()
    eval_raw = make_eval_grid()
    save_fig1_a()
    save_fig1_b()
    save_scatter(data, eval_raw)

    data_scaled, eval_scaled = scale_data(data, eval_raw)
    y = data_scaled["y"].to_numpy()
    distance = signed_distances(data_scaled, eval_scaled)

    smooth = read_empapp_output("smooth_itt")
    adaptive = read_empapp_output("adaptive_itt")
    unknown_kink = read_empapp_output("unknown_kink_itt")
    rdrobust = read_empapp_output("rdrobust_itt")
    save_fig4(smooth, "fig4-smooth.png")
    save_fig4(adaptive, "fig4-adaptive.png")
    save_fig4(unknown_kink, "fig4-unknown_kink.png")
    save_fig4(rdrobust, "fig4-rdrobust.png")

    for idx, suffix in [(1, "b"), (11, "c"), (21, "d")]:
        h = float(smooth.iloc[idx - 1]["h0"])
        save_fig2_detail(y, distance, h, idx, suffix)

    expected = [
        "fig1-a.png",
        "fig1-b.png",
        "fig2-a.png",
        "fig2-b.png",
        "fig2-c.png",
        "fig2-d.png",
        "fig4-smooth.png",
        "fig4-adaptive.png",
        "fig4-unknown_kink.png",
        "fig4-rdrobust.png",
    ]
    check_generated_figures(expected)
    print(f"Wrote empirical figure files to {FIGURES_DIR.relative_to(SCRIPT_DIR)}.")


if __name__ == "__main__":
    main()

