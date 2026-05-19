# Cattaneo, Titiunik and Yu (2026, JOE)

Replication files for Cattaneo, Titiunik and Yu (2026), "Estimation and
Inference in Boundary Discontinuity Designs: Distance-Based Methods."

This repository contains the empirical application and simulation scripts for
the distance-based boundary discontinuity design analysis. The computations use
the `rd2d` package for R and Python.

## Website

https://rdpackages.github.io/replication

## Data

The empirical application uses the Ser Pilo Paga data from:

- Londono-Velez, J., C. Rodriguez, and F. Sanchez (2020):
  [Upstream and Downstream Impacts of College Merit-Based Financial Aid for
  Low-Income Students: Ser Pilo Paga in Colombia](https://doi.org/10.1257/pol.20180131),
  _American Economic Journal: Economic Policy_ 12(2): 193-227.

The analysis dataset is [`spp.csv`](spp.csv). The scripts expect these columns:

- `running_saber11`
- `running_sisben`
- `eligible_spp`
- `beneficiary_spp`
- `spadies_any`
- `icfes_educm1`

## Repository Layout

- [`CTY_2026_JOE--empapp.R`](CTY_2026_JOE--empapp.R): runs the R empirical
  application and writes CSV outputs to `output/`.
- [`CTY_2026_JOE--empapp.py`](CTY_2026_JOE--empapp.py): runs the Python
  empirical application and writes CSV outputs to `output/`.
- [`CTY_2026_JOE--empapp_tables.py`](CTY_2026_JOE--empapp_tables.py): converts
  Python empirical CSV outputs into LaTeX table fragments in `tables/`.
- [`CTY_2026_JOE--empapp_figures.py`](CTY_2026_JOE--empapp_figures.py):
  converts Python empirical CSV outputs into figures in `figures/`.
- [`CTY_2026_JOE--empapp_tables.R`](CTY_2026_JOE--empapp_tables.R): converts
  empirical CSV outputs into LaTeX table fragments in `tables/`.
- [`CTY_2026_JOE--empapp_figures.R`](CTY_2026_JOE--empapp_figures.R):
  converts empirical CSV outputs into figures in `figures/`.
- [`CTY_2026_JOE--simuls.R`](CTY_2026_JOE--simuls.R): runs the simulation
  study and writes raw simulation outputs to `output/`.
- [`CTY_2026_JOE--simuls_tables.R`](CTY_2026_JOE--simuls_tables.R): converts
  simulation outputs into LaTeX table fragments in `tables/`.
- [`CTY_2026_JOE--simuls_figures.R`](CTY_2026_JOE--simuls_figures.R): converts
  simulation outputs into figures in `figures/`.
- [`spp.csv`](spp.csv): SPP data used by the empirical application and
  simulation design.

Generated artifacts are reproducible from the scripts. Three local folders are
required:

- `figures/`
- `tables/`
- `output/`

These generated folders are intentionally ignored by Git.

## Requirements

The R scripts require `rd2d` version 0.1.0 or newer, with support for
`params.other`, `params.cov`, and `summary(..., cbands = ...)`. Install the
released R package from CRAN before running the R replication scripts:

```r
install.packages("rd2d")
library(rd2d)
```

The Python empirical script requires Python 3.10 or newer and the Python `rd2d`
package with replication dependencies:

```sh
python -m pip install "rd2d[replication]"
```

## Replication

To replicate the empirical application with R:

```sh
Rscript CTY_2026_JOE--empapp.R
Rscript CTY_2026_JOE--empapp_tables.R
Rscript CTY_2026_JOE--empapp_figures.R
```

To replicate the empirical application with Python:

```sh
python CTY_2026_JOE--empapp.py
python CTY_2026_JOE--empapp_tables.py
python CTY_2026_JOE--empapp_figures.py
```

To replicate the simulation study:

```sh
Rscript CTY_2026_JOE--simuls.R
Rscript CTY_2026_JOE--simuls_tables.R
Rscript CTY_2026_JOE--simuls_figures.R
```

The simulation design considers linear and quadratic mean specifications with
homoskedastic and heteroskedastic variance specifications, for sharp and fuzzy designs.

The full simulation is computationally expensive. These environment variables
can be used for smaller local checks:

- `RD2D_M`: number of simulation replications, default `5000`
- `RD2D_N`: sample size per replication, default `20000`
- `RD2D_REPP`: simulation repetitions for critical values, default `2000`
- `RD2D_SEED`: simulation seed, default `20260510`
- `RD2D_WORKERS`: parallel workers, default based on available cores
- `RD2D_EMP_REPP`: empirical critical-value repetitions, default `5000`; used by both the R and Python empirical scripts

Example quick simulation smoke run:

```sh
RD2D_M=1 RD2D_N=3000 RD2D_REPP=49 RD2D_WORKERS=1 Rscript CTY_2026_JOE--simuls.R
```

## Outputs

- `output/empapp_*.csv`: empirical application data and summary tables.
- `output/simuls_*.csv`: simulation metadata, DGP calibration, targets, and run
  summaries.
- `output/simuls_raw_dgp*.csv`: raw simulation files stacked by DGP, design,
  method, and bandwidth rule.
- `tables/empapp_*.tex`: empirical application LaTeX table fragments.
- `tables/simuls_*.tex`: simulation LaTeX table fragments.
- `figures/*.png`: generated empirical application and simulation figures.

By default, outputs are written to `output/`, `tables/`, and `figures/`. The
simulation scripts also support `RD2D_OUTPUT_DIR`, `RD2D_TABLES_DIR`, and
`RD2D_FIGURES_DIR` for local output-directory overrides.

## Reference

- Cattaneo, M. D., R. Titiunik, and R. R. Yu (2026):
  [Estimation and Inference in Boundary Discontinuity Designs: Distance-Based
  Methods](https://rdpackages.github.io/references/Cattaneo-Titiunik-Yu_2026_JOE.pdf).<br>
  _Journal of Econometrics_, revise and resubmit.<br>
  [Supplemental Appendix](https://rdpackages.github.io/references/Cattaneo-Titiunik-Yu_2026_JOE--Supplement.pdf).


## Acknowledgment

This work was supported in part by the National Science Foundation through
grants [SES-2019432](https://www.nsf.gov/awardsearch/showAward?AWD_ID=2019432),
[DMS-2210561](https://www.nsf.gov/awardsearch/showAward?AWD_ID=2210561),
[SES-2241575](https://www.nsf.gov/awardsearch/showAward?AWD_ID=2241575), and
[SES-2342226](https://www.nsf.gov/awardsearch/showAward?AWD_ID=2342226), and by
the National Institute for Food and Agriculture through grant
[2024-67023-42704](https://www.nifa.usda.gov/data).

