# Beyond Income: Subjective Social Class, Psychological Resources, and Quality of Life in South Korea During the COVID-19 Pandemic

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20960562.svg)](https://doi.org/10.5281/zenodo.20960562)

Analysis code to reproduce all results and figures.

## Contents
- `analysis.R` — full pipeline: preprocessing, CFA of the latent happiness
  construct, intraclass correlations, multicollinearity (VIF), the three
  competing structural equation models (income-only / SSC-only / both), the
  full model, and all robustness analyses (Harman's test, single-factor CFA,
  unmeasured latent method construct, objective-SES composite, overlap-excluded
  model, hierarchical R², FIML, survey-weighted SEM, mediation decomposition,
  cluster-robust SEs, outlier sensitivity). Writes `results/r_revision_results.json`.
- `make_figures.py` — regenerates all manuscript figures (black-and-white) into
  `figures/`.
- `data/` — place the KOSSDA data file here (see `data/README.md`).

## Requirements
- **R** (>= 4.5) with `lavaan` (0.6-21), `psych`, `jsonlite`.
- **Python** (>= 3.9) with `numpy`, `pandas`, `matplotlib`.

## How to run
From the repository root:
```bash
Rscript analysis.R          # -> results/r_revision_results.json
python3 make_figures.py     # -> figures/*.png, *.pdf
```

## Reproducibility notes
- `lavaan` uses the robust maximum likelihood (MLR) estimator; estimates are
  deterministic for a given dataset and package version.
- All headline numbers in the manuscript are computed by `analysis.R`; none are
  hard-coded.
- Newer `lavaan` releases may shift fit indices by trivial amounts without
  changing any substantive conclusion.

## Data availability
Microdata: KOSSDA dataset A1-2021-0003 (v2.0),
https://doi.org/10.22687/KOSSDA-A1-2021-0003-V2.0 . Not redistributed here. No
personally identifiable information is used.

## License
MIT (see `LICENSE`).
