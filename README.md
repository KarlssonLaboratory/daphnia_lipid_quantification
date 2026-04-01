# daphnia_droplets

Statistical analysis pipeline for toxicological experiments using *Daphnia*, examining the effects of chemical and food treatments on lipid droplet accumulation measured by confocal microscopy.

---

## Overview

*Daphnia* were exposed to chemicals (e.g. TBT) across three food regimes:

| Code | Treatment |
|------|-----------|
| NF | No food |
| LF | Low food |
| HF | High food |

After exposure, individuals were imaged by confocal microscopy. Lipid droplets were segmented using deep learning, yielding per-animal summary statistics and individual droplet measurements. This pipeline loads those outputs, fits generalized linear models (GLMs) and linear mixed models (LMMs), and writes diagnostic plots and pairwise comparison tables to disk.

---

## Repository structure

```
daphnia_droplets/
├── bin/
│   ├── analysis_scalable_chems.R     # Main multi-experiment pipeline
│   └── pilot_analysis_only_food.R    # Pilot analysis (food treatment only)
├── data/
│   └── <experiment>/                 # One folder per chemical/experiment
│       ├── measures_export_<N>.csv
│       ├── measures_export_<N>_wells.xlsx
│       └── masks_<N>/
│           └── <well>_droplets.csv   # One file per individual
└── output/
    └── <experiment>/
        ├── graphs/
        └── results/
```

---

## Input data format

All data for an experiment must be placed in `data/<experiment_name>/`. Experiment folders named `pilot` (case-insensitive) are excluded from batch runs.

### Per-well summary table — `measures_export_<N>.csv`

One file per replicate (e.g. `measures_export_1.csv`, `measures_export_2.csv`). Each row is one imaged well. Must contain at least:

| Column | Description |
|---|---|
| `well_name` | Well identifier (e.g. `A1`); first character is the row letter |
| `daphnia_area` | Body area (used to normalise droplet count) |
| `droplet_intensity_total` | Total lipid-droplet fluorescence intensity |
| `droplet_area_total` | Summed area of all droplets |
| `daphnia_size` | Body-size proxy |
| `num_droplets` | Count of detected droplets |

### Treatment map — `measures_export_<N>_wells.xlsx`

One Excel file per replicate, with three columns:

| Column | Description |
|---|---|
| `well` | Capitalised row letter (e.g. `A`, `B`) |
| `food_treatment` | `"no food"`, `"low food"`, or `"high food"` |
| `chem_treatment` | Chemical treatment label (can be a single value if no chemical) |

### Individual droplet data — `masks_<N>/<well>_droplets.csv`

One CSV per individual (well), inside a `masks_<N>/` subfolder for each replicate. Each row is one segmented droplet and must contain an `area` column. Files are optional; if absent the per-droplet LMM is skipped.

---

## Running the analysis

Open R and source the main script:

```r
source("bin/analysis_scalable_chems.R")
```

A dialog box will appear listing all available experiments. Select one, or choose **ALL** to process every experiment sequentially.

For the pilot experiment (food treatment only, no chemical), use:

```r
source("bin/pilot_analysis_only_food.R")
```

### Dependencies

The script will install missing packages automatically. Required packages:

`tidyverse`, `readxl`, `multcomp`, `lme4`, `emmeans`, `DHARMa`, `MASS`, `patchwork`, `car`

---

## Statistical models

### Response variables

| Response | GLM family | Description |
|---|---|---|
| `droplet_intensity_total` | Gamma (log) | Total lipid-droplet fluorescence |
| `droplet_area_total` | Gaussian | Summed droplet area |
| `daphnia_size` | Gaussian | Body-size proxy |
| `num_droplets` | Negative binomial | Count of detected droplets |
| `num_droplets_per_area` | Gamma (log) | Droplet count ÷ body area |

Response variables and their families are defined in the `response` tibble near the top of `analysis_scalable_chems.R` and can be modified there.

### Model structure

Two GLMs are fitted per response variable:

- **Main model** — `response ~ food_treatment * chem_treatment` (or `~ food_treatment` if only one chemical level is present)
- **Replicate-adjusted model** — adds `+ replicate` as a fixed covariate to account for batch effects

Pairwise contrasts are computed with `emmeans` and p-values are adjusted using the Benjamini–Hochberg (BH) method. Compact letter display (CLD) groups are overlaid on boxplots.

For individual droplet area, a **Gamma LMM** (`glmer`) is fitted with `well_name` as a random intercept to account for pseudo-replication within wells.

### Quality control

Wells with five or more zero-valued numeric columns are flagged and removed before modelling (likely imaging failures). A warning listing the affected wells is printed to the console.

---

## Output

Results are written to `output/<experiment>/`.

### `graphs/`

| File | Contents |
|---|---|
| `histograms.pdf` | Distribution of each response variable, overall and by treatment group |
| `<response>_by_treatment.pdf` | Boxplot with jittered points (coloured by replicate) and CLD letters |
| `individual_droplets_by_treatment.pdf` | Violin plot of individual droplet area by treatment |
| `model_diagnostics.pdf` | DHARMa residual diagnostics for each GLM |
| `individual_droplet_area_diagnostics.pdf` | Diagnostics for the droplet-area LMM |

### `results/`

| File | Contents |
|---|---|
| `<response>_anova_treatment.csv` | Type II ANOVA table for the main model |
| `<response>_anova_replicates.csv` | Type II ANOVA table for the replicate-adjusted model |
| `<response>_pairs_treatment.csv` | BH-adjusted pairwise contrasts (main model) |
| `<response>_pairs_replicates.csv` | BH-adjusted pairwise contrasts (replicate-adjusted model) |
| `individual_droplet_area_pairs_treatment.csv` | Pairwise contrasts for droplet area LMM |

---

## Example data layout

The `data/TBT-CL/` folder illustrates the expected structure for a two-replicate experiment:

```
data/TBT-CL/
├── measures_export_1.csv
├── measures_export_1_wells.xlsx
├── measures_export_2.csv
├── measures_export_2_wells.xlsx
├── masks_1/
│   └── *.csv   (one file per individual)
└── masks_2/
    └── *.csv
```