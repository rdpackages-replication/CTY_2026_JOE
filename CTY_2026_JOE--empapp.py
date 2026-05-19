from __future__ import annotations

import argparse
import os
from pathlib import Path
from importlib.metadata import PackageNotFoundError, version

import numpy as np
import pandas as pd
from rdrobust import rdbwselect
from scipy.spatial.distance import cdist


from rd2d import rd2d_distance




def rd2d_version() -> str:
    try:
        return version("rd2d")
    except PackageNotFoundError:
        return "unknown"

REQUIRED_COLUMNS = [
    "running_saber11",
    "running_sisben",
    "spadies_any",
    "eligible_spp",
    "beneficiary_spp",
]


def load_spp_data(path: Path) -> pd.DataFrame:
    raw = pd.read_csv(path)
    missing = [name for name in REQUIRED_COLUMNS if name not in raw.columns]
    if missing:
        raise ValueError(f"spp.csv is missing columns: {', '.join(missing)}")

    dat = raw.loc[:, REQUIRED_COLUMNS].copy()
    dat.columns = ["x.1", "x.2", "y", "d", "w"]
    dat = dat.dropna().copy()
    expected_assignment = ((dat["x.1"] >= 0) & (dat["x.2"] >= 0)).astype(int)
    if not np.all(dat["d"].to_numpy() == expected_assignment.to_numpy()):
        raise ValueError("eligible_spp does not match the quadrant assignment rule.")
    dat["y"] = dat["y"].astype(float)
    dat["d"] = dat["d"].astype(float)
    dat["w"] = dat["w"].astype(float)
    return dat


def make_eval_grid() -> pd.DataFrame:
    rows = []
    for i in range(1, 21):
        rows.append((0.0, 40 - (i - 1) * 40 / 20))
    for i in range(21, 41):
        rows.append(((i - 20 - 1) * 56 / 20, 0.0))
    eval_grid = pd.DataFrame(rows, columns=["x.1", "x.2"])
    return eval_grid.iloc[10:31, :].reset_index(drop=True)


def scale_data(dat: pd.DataFrame, eval_grid: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    out = dat.copy()
    ev = eval_grid.copy()
    scale2 = out["x.1"].std(ddof=1) / out["x.2"].std(ddof=1)
    out["x.2"] = out["x.2"] * scale2
    ev["x.2"] = ev["x.2"] * scale2
    return out, ev


def signed_distances(dat: pd.DataFrame, eval_grid: pd.DataFrame) -> np.ndarray:
    x = dat[["x.1", "x.2"]].to_numpy()
    b = eval_grid[["x.1", "x.2"]].to_numpy()
    sign = (2 * dat["d"].to_numpy() - 1).reshape(-1, 1)
    return cdist(x, b, metric="euclidean") * sign


def rdrobust_bandwidths(y: np.ndarray, distance: np.ndarray) -> np.ndarray:
    bws = np.zeros((distance.shape[1], 2))
    for j in range(distance.shape[1]):
        out = rdbwselect(y, distance[:, j], vce="hc1")
        table = out.bws
        bws[j, 0] = float(table.iloc[0, 0])
        bws[j, 1] = float(table.iloc[0, 1])
    return bws


def as_result_table(result, eval_grid: pd.DataFrame, output: str) -> pd.DataFrame:
    point_res = result[output]
    n = len(point_res)
    summ = result.summary(
        output=output,
        cbands=output,
        WBATE=np.ones(n),
        LBATE=True,
    )
    tab = summ.tables[output].copy()
    point_rows = n

    def col_or_na(name: str) -> pd.Series:
        if name in tab:
            return tab[name]
        return pd.Series(np.repeat(np.nan, len(tab)), index=tab.index)

    row = [str(i) for i in range(1, point_rows + 1)]
    row.extend(str(idx) for idx in tab.index[point_rows:])
    return pd.DataFrame(
        {
            "row": row,
            "b1": list(eval_grid["x.1"]) + [np.nan] * (len(tab) - point_rows),
            "b2": list(eval_grid["x.2"]) + [np.nan] * (len(tab) - point_rows),
            "h0": col_or_na("h0"),
            "h1": col_or_na("h1"),
            "N.Co": col_or_na("N.Co"),
            "N.Tr": col_or_na("N.Tr"),
            "estimate.p": col_or_na("estimate.p"),
            "ci.lower": col_or_na("ci.lower"),
            "ci.upper": col_or_na("ci.upper"),
            "cb.lower": col_or_na("cb.lower"),
            "cb.upper": col_or_na("cb.upper"),
            "estimate.q": col_or_na("estimate.q"),
            "t.value": col_or_na("t.value"),
            "p.value": col_or_na("p.value"),
        }
    )


def check_itt_match(sharp, fuzzy, label: str, tol: float = 1e-10) -> None:
    cols = [
        "estimate.p",
        "std.err.p",
        "estimate.q",
        "std.err.q",
        "t.value",
        "p.value",
        "ci.lower",
        "ci.upper",
        "h0",
        "h1",
        "h0.rbc",
        "h1.rbc",
        "N.Co",
        "N.Tr",
    ]
    diffs = {col: np.nanmax(np.abs(sharp.main[col].to_numpy() - fuzzy.itt[col].to_numpy())) for col in cols}
    worst_col = max(diffs, key=diffs.get)
    if diffs[worst_col] > tol:
        raise ValueError(f"fuzzy itt does not match sharp main for {label}: {worst_col}={diffs[worst_col]:.3e}")


def write_csv(frame: pd.DataFrame, path: Path) -> None:
    frame.to_csv(path, index=False, na_rep="")


def run(data_path: Path, output_dir: Path, repp: int) -> int:
    output_dir.mkdir(parents=True, exist_ok=True)
    for old in output_dir.glob("empapp_*.csv"):
        old.unlink()

    dat = load_spp_data(data_path)
    eval_grid = make_eval_grid()
    dat, eval_grid = scale_data(dat, eval_grid)
    y = dat["y"].to_numpy()
    fuzzy = dat["w"].to_numpy()
    distance = signed_distances(dat, eval_grid)
    b = eval_grid.to_numpy()
    cov_tables = ["main", "itt", "fs"]

    sharp_results = {
        "smooth": rd2d_distance(y, distance, b=b, repp=repp, vce="hc1"),
        "unknown_kink": rd2d_distance(y, distance, b=b, kink_unknown=(True, False), repp=repp, vce="hc1"),
        "adaptive": rd2d_distance(y, distance, b=b, kink_position=10, repp=repp, vce="hc1"),
    }
    fuzzy_results = {
        "smooth": rd2d_distance(y, distance, b=b, repp=repp, vce="hc1", fuzzy=fuzzy, bwparam="itt", params_cov=cov_tables),
        "unknown_kink": rd2d_distance(
            y,
            distance,
            b=b,
            kink_unknown=(True, False),
            repp=repp,
            vce="hc1",
            fuzzy=fuzzy,
            bwparam="itt",
            params_cov=cov_tables,
        ),
        "adaptive": rd2d_distance(
            y,
            distance,
            b=b,
            kink_position=10,
            repp=repp,
            vce="hc1",
            fuzzy=fuzzy,
            bwparam="itt",
            params_cov=cov_tables,
        ),
    }

    bws = rdrobust_bandwidths(y, distance)
    sharp_results["rdrobust"] = rd2d_distance(y, distance, h=bws, b=b, repp=repp, vce="hc1")
    fuzzy_results["rdrobust"] = rd2d_distance(
        y,
        distance,
        h=bws,
        b=b,
        repp=repp,
        vce="hc1",
        fuzzy=fuzzy,
        bwparam="itt",
        params_cov=cov_tables,
    )

    for name in fuzzy_results:
        check_itt_match(sharp_results[name], fuzzy_results[name], name)

    outputs = {"fuzzy": "main", "itt": "itt", "fs": "fs"}
    count = 0
    for name, result in fuzzy_results.items():
        for output_name, output in outputs.items():
            write_csv(as_result_table(result, eval_grid, output), output_dir / f"empapp_{name}_{output_name}.csv")
            count += 1
    return count


def main() -> None:
    os.chdir(Path(__file__).resolve().parent)
    parser = argparse.ArgumentParser()
    parser.add_argument("--data", type=Path, default=Path("spp.csv"))
    parser.add_argument("--output", type=Path, default=Path("output"))
    parser.add_argument("--repp", type=int, default=int(os.getenv("RD2D_EMP_REPP", "5000")))
    args = parser.parse_args()
    if args.repp <= 0:
        raise ValueError("--repp must be a positive integer.")

    count = run(args.data, args.output, args.repp)
    print(f"Wrote {count} empirical output file(s) to {args.output}. Sharp-main/ITT checks passed.")


if __name__ == "__main__":
    main()


