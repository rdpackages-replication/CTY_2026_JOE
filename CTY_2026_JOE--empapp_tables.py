from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent
OUTPUT_DIR = SCRIPT_DIR / "output"
TABLES_DIR = SCRIPT_DIR / "tables"
AGGREGATE_LABELS = ["WBATE", "LBATE"]
METHODS = ["smooth", "adaptive", "unknown_kink", "rdrobust"]
ESTIMANDS = ["fuzzy", "itt", "fs"]
ESTIMAND_HEADERS = {
    "fuzzy": "$\\tau(\\bb)$",
    "itt": "$\\tau_Y(\\bb)$",
    "fs": "$\\tau_W(\\bb)$",
}


def read_empapp_output(method: str, estimand: str) -> pd.DataFrame:
    path = OUTPUT_DIR / f"empapp_{method}_{estimand}.csv"
    if not path.exists():
        raise FileNotFoundError(f"Missing output file: {path}. Run CTY_2026_JOE--empapp.py first.")
    tab = pd.read_csv(path)
    required = [
        "row",
        "h0",
        "h1",
        "N.Co",
        "N.Tr",
        "estimate.p",
        "ci.lower",
        "ci.upper",
        "cb.lower",
        "cb.upper",
        "estimate.q",
        "t.value",
        "p.value",
    ]
    missing = [name for name in required if name not in tab.columns]
    if missing:
        raise ValueError(f"{path} is missing columns: {', '.join(missing)}")
    tab["row"] = tab["row"].astype(str)
    return tab


def clean_tables_dir() -> None:
    TABLES_DIR.mkdir(parents=True, exist_ok=True)
    for old in TABLES_DIR.glob("empapp_*.tex"):
        old.unlink()
    for old in TABLES_DIR.glob("tab-empapp_*.tex"):
        old.unlink()


def scalar(value: object) -> float:
    try:
        out = float(value)
    except (TypeError, ValueError):
        return np.nan
    return out


def fmt_num(value: object, digits: int = 3) -> str:
    x = scalar(value)
    if not np.isfinite(x):
        return ""
    return f"{x:.{digits}f}"


def fmt_pvalue(value: object) -> str:
    x = scalar(value)
    if not np.isfinite(x):
        return ""
    if x < 0.001:
        return "0.000"
    return f"{x:.3f}"


def fmt_ci(lower: object, upper: object) -> str:
    lo = scalar(lower)
    hi = scalar(upper)
    if not np.isfinite(lo) or not np.isfinite(hi):
        return ""
    return f"$({fmt_num(lo)},\\, {fmt_num(hi)})$"


def fmt_row_label(row_id: str) -> str:
    if row_id == "WBATE":
        return "$\\mathtt{WBATE}$"
    if row_id == "LBATE":
        return "$\\mathtt{LBATE}$"
    return f"$\\mathbf{{b}}_{{{row_id}}}$"


def interval_limits(tab: pd.DataFrame) -> tuple[pd.Series, pd.Series]:
    cb_lower = pd.to_numeric(tab["cb.lower"], errors="coerce")
    cb_upper = pd.to_numeric(tab["cb.upper"], errors="coerce")
    use_cb = np.isfinite(cb_lower) & np.isfinite(cb_upper)
    lower = tab["ci.lower"].copy()
    upper = tab["ci.upper"].copy()
    lower.loc[use_cb] = tab.loc[use_cb, "cb.lower"]
    upper.loc[use_cb] = tab.loc[use_cb, "cb.upper"]
    return lower, upper


def display_rows(tab: pd.DataFrame) -> pd.DataFrame:
    point_rows = tab.loc[~tab["row"].isin(AGGREGATE_LABELS)].copy()
    aggregate_rows = tab.loc[tab["row"].isin(AGGREGATE_LABELS)].copy()
    aggregate_rows = aggregate_rows.set_index("row").reindex(AGGREGATE_LABELS).reset_index()
    if aggregate_rows["estimate.p"].isna().all():
        raise ValueError("Empirical output is missing WBATE or LBATE rows.")
    out = pd.concat([point_rows, aggregate_rows], ignore_index=True)
    lower, upper = interval_limits(out)
    return pd.DataFrame(
        {
            "label": [fmt_row_label(row) for row in out["row"]],
            "h": pd.to_numeric(out["h0"], errors="coerce"),
            "N.Co": pd.to_numeric(out["N.Co"], errors="coerce"),
            "N.Tr": pd.to_numeric(out["N.Tr"], errors="coerce"),
            "estimate": pd.to_numeric(out["estimate.p"], errors="coerce"),
            "pvalue": pd.to_numeric(out["p.value"], errors="coerce"),
            "lower": pd.to_numeric(lower, errors="coerce"),
            "upper": pd.to_numeric(upper, errors="coerce"),
        }
    )


def latex_body(rows: pd.DataFrame) -> list[str]:
    lines: list[str] = []
    previous_aggregate = False
    for _, row in rows.iterrows():
        is_aggregate = "WBATE" in row["label"] or "LBATE" in row["label"]
        if is_aggregate and not previous_aggregate:
            lines.append("\\midrule")
        previous_aggregate = is_aggregate
        lines.append(
            row["label"]
            + " & "
            + fmt_num(row["h"])
            + " & "
            + fmt_num(row["N.Co"], 0)
            + " & "
            + fmt_num(row["N.Tr"], 0)
            + " & "
            + fmt_num(row["estimate"])
            + " & "
            + fmt_pvalue(row["pvalue"])
            + " & "
            + fmt_ci(row["lower"], row["upper"])
            + " \\\\"
        )
    return lines


def write_tabular(rows: pd.DataFrame, path: Path, estimand: str) -> None:
    lines = [
        "\\begin{tabular}{@{}crrrrrc@{}}",
        "\\toprule\\toprule",
        " & ".join(
            [
                "\\multicolumn{1}{c}{$\\bb\\in\\B$}",
                "\\multicolumn{1}{c}{$h$}",
                "\\multicolumn{1}{c}{$N_{\\mathrm{Co}}$}",
                "\\multicolumn{1}{c}{$N_{\\mathrm{Tr}}$}",
                f"\\multicolumn{{1}}{{c}}{{{ESTIMAND_HEADERS[estimand]}}}",
                "\\multicolumn{1}{c}{p-value}",
                "\\multicolumn{1}{c}{95\\% RBC CI}",
            ]
        )
        + " \\\\",
        "\\midrule",
        *latex_body(rows),
        "\\bottomrule\\bottomrule",
        "\\end{tabular}",
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def check_generated_tables(expected: list[Path]) -> None:
    missing = [str(path) for path in expected if not path.exists()]
    if missing:
        raise FileNotFoundError(f"Missing expected empirical table: {missing[0]}")
    for path in expected:
        text = path.read_text(encoding="utf-8")
        if "WBATE" not in text or "LBATE" not in text:
            raise ValueError(f"Generated table is missing WBATE or LBATE: {path}")
        if "$N_{\\mathrm{Co}}$" not in text or "$N_{\\mathrm{Tr}}$" not in text:
            raise ValueError(f"Generated table is missing N.Co or N.Tr headers: {path}")


def main() -> None:
    clean_tables_dir()
    outputs = {method: {estimand: read_empapp_output(method, estimand) for estimand in ESTIMANDS} for method in METHODS}
    expected: list[Path] = []
    for method in METHODS:
        for estimand in ESTIMANDS:
            path = TABLES_DIR / f"empapp_{method}_{estimand}.tex"
            write_tabular(display_rows(outputs[method][estimand]), path, estimand)
            expected.append(path)
    check_generated_tables(expected)
    print(f"Wrote {len(expected)} empirical table file(s) to {TABLES_DIR.relative_to(SCRIPT_DIR)}.")


if __name__ == "__main__":
    main()
