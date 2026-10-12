# PSILAUT_ONE_Paper
Custom code used for publication: Distinct response of resting-state brain networks to psilocybin in autism

# PSILAUT: analysis scripts

R scripts accompanying:

> Whelan, T. P. *et al.* Distinct response of resting-state brain networks to psilocybin in autism. *Nature Mental Health* (2026). DOI: [to add]

## What the code does

- **`parseData.R`** defines `parseData()`, which reads one 100 × 100 ROI-to-ROI connectivity matrix (Fisher z-transformed correlations from CONN v21.a, Schaefer 100-parcel atlas). It assigns each ROI to its Yeo 7-network and returns the mean connectivity for every network pair: 7 within-network and 21 between-network values. Networks are relabelled as DMN (default mode), FPN (frontoparietal), DAN (dorsal attention), VAN (ventral attention), LN (limbic), SMN (somatomotor) and VN (visual).
- **`run.R`** applies `parseData()` to every matrix in the `data/` folder, averages the network-level values across all files, and plots them as a heatmap of mean Fisher z-scores.
- **`ONE_linear_regression_script.R`** tests whether a baseline measure, the Autism Quotient (AQ) total score, predicts the change in a connectivity measure from placebo (0 mg) to 5 mg psilocybin (5 mg minus 0 mg). It adjusts for covariates as appropriate. All variables are standardised. It reports coefficients with 95% confidence intervals and the partial correlation for the baseline predictor.

## Repository contents

```
.
├── README.md
├── LICENSE
├── parseData.R
├── run.R
├── ONE_linear_regression_script.R
└── schaefer_100_labels.csv   # ROI labels, one per row, no header (e.g. LH_Vis_1)
```

## System requirements

- R 2022.07.0.
- R packages: tidyverse, reshape2, broom and RColorBrewer.

## Installation

```r
install.packages(c("tidyverse", "reshape2", "broom", "RColorBrewer"))
```

## Input data

Participant data are not included (see Data availability). The scripts expect:

1. **Connectivity matrices** (`run.R` and `parseData.R`): one CSV file per scan in a folder named `data/`. Each file holds a 100 × 100 matrix of Fisher z-transformed correlations, with no header row, in the same ROI order as `schaefer_100_labels.csv`. Each ROI label must have its network name as the second "_"-separated element (e.g. `LH_Vis_1`).
2. **Long-format table** (`ONE_linear_regression_script.R`): a CSV file with one row per participant per dose, containing:
   - `ID`: participant identifier
   - `Dose`: 0, 2 or 5
   - the connectivity outcome (e.g. DMN–FPN connectivity)
   - the baseline predictor (`AQ`), present on the 0 mg row
   - Covariates such as `Mean_FD` and `Visit` for each session

## Data availability

The data are not publicly available. They are available on reasonable request, subject to approval by the study sponsor and industry partners, participants' consent and ethical approval, and a data-sharing agreement. See the Data Availability statement in the paper.

## Contact

Tobias P. Whelan, King's College London, tobias.whelan@kcl.ac.uk.
