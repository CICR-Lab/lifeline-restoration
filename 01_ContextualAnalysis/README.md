# Module 01 - Contextual Analysis

## 1. Overview

This module implements the contextual analysis of post-earthquake lifeline
restoration. It uses the external canonical dataset placed in Module 00 and supports:

- main-text Figure 2;
- Supplementary Figures 2-5;
- the fitted models, record-set definitions, sensitivity analyses, LOEO-CV
  results and diagnostics underlying those figures.

The module extracts contextual variables, constructs the analysis record sets,
fits the contextual models, evaluates earthquake-held-out predictions and
renders the main and supplementary figures.

## 2. Requirements

- MATLAB with Statistics and Machine Learning Toolbox.
- Optimization Toolbox for beta and restoration-time distribution sensitivities;
  Parallel Computing Toolbox for supplementary bootstrap runs of 100 or more
  replicates.
- The canonical dataset downloaded after author-approved access at
  <https://doi.org/10.5281/zenodo.23113551>, placed at
  `../00_Dataset/standardized_dataset.mat` for analysis workflows.
- The required runtime configuration `contextual_config.json`.
- Sufficient disk space for model fits, bootstrap results and figures.

## 3. Inputs

### Upstream dataset

```text
../00_Dataset/standardized_dataset.mat
```

Before running the analysis, follow the Module 00 instructions to request
access through Zenodo, obtain author approval and download the MAT file to the
path above. The analysis code derives the contextual extraction and analysis
table from it. Figure-only rendering with `DataSource="packaged"` reads the CSV
inputs in `figure_data/main/` and `figure_data/supplementary/`.

### Runtime configuration

```text
contextual_config.json
```

This file records the model settings, random seed and analysis parameters used
by the main and supplementary calculations. It is required for a full run.

### Frozen figure data

`figure_data/main/` and `figure_data/supplementary/` contain the direct CSV
inputs distributed for redrawing Figure 2 and Supplementary Figures 2-5
without refitting the models.

## 4. Directory structure

```text
01_ContextualAnalysis/
├── README.md
├── contextual_config.json
├── run_contextual_analysis.m          # full analysis entry point
├── run_figure_2.m                      # Figure 2 entry point
├── run_supplementary_contextual_analysis.m
│                                      # Supplementary Figures 2-5 entry point
├── Code/                               # internal analysis and rendering functions
├── figure_data/
│   ├── main/                           # frozen Figure 2 inputs
│   └── supplementary/                  # frozen Supplementary Figure inputs
└── outputs/                            # runtime-generated files; empty in release
    ├── intermediate/
    ├── analysis/
    └── supplementary/
```

The public entry points are at the module root. Internal implementation
files are kept in `Code/`; frozen figure inputs are kept outside `outputs/`.

## 5. Running the module

### Full run

From the `01_ContextualAnalysis` directory, run:

```matlab
addpath(pwd);
run_contextual_analysis
```

The full entry point performs the following sequence:

1. extracts the broad contextual input;
2. builds the canonical analysis table and record-set flags;
3. fits the main initial-disruption and restoration-time models;
4. calculates fitted fit statistics and Figure 2 data;
5. calculates LOEO-CV results;
6. runs sensitivity and diagnostic analyses;
7. runs the supplementary sensitivity, prediction and diagnostic analyses.

The default run uses 2,000 bootstrap replicates. A new full run requires that
the target analysis output directory does not already exist.

### Smoke test

For a short execution test, request a small number of bootstrap replicates:

```matlab
run_contextual_analysis([], [], 2)
```

Smoke-test outputs support software validation. The full run with the
configured bootstrap count produces the manuscript results.

The focused CR2 implementation tests can be run independently:

```matlab
results = runtests(fullfile("Code", "tests"));
assertSuccess(results)
```

### Figure-only rendering

The frozen figure inputs can be rendered without rerunning the analysis:

```matlab
run_figure_2("DataSource", "packaged")
run_supplementary_contextual_analysis("DataSource", "packaged")
```

After a local full analysis, the same plotting code can instead read the
generated figure data:

```matlab
run_figure_2("DataSource", "generated")
run_supplementary_contextual_analysis("DataSource", "generated")
```

`packaged` reads only the selected files under `figure_data/`; `generated`
reads only `outputs/analysis/figure_data/` or
`outputs/supplementary/figure_data/`, respectively. Each plotting call uses the
single source selected by `DataSource`. `OutputDir` can be supplied to redirect
the rendered figures without changing the selected data source. `DataRoot` can be supplied
when the selected direct figure-data folder is stored elsewhere.

## 6. Outputs

| Output | Use |
|---|---|
| `outputs/intermediate/` | Extracted contextual data and the canonical analysis table |
| `outputs/analysis/` | Main contextual model results, validation files and generated Figure 2 data |
| `outputs/supplementary/` | Sensitivity, LOEO-CV and diagnostic results |
| `outputs/analysis/figure_data/` | Generated direct inputs for Figure 2 |
| `outputs/supplementary/figure_data/` | Generated direct inputs for Supplementary Figures 2-5 |
| `outputs/figures/main/` | Figure 2 rendered from the selected figure-data source |
| `outputs/figures/supplementary/` | Supplementary Figures 2-5 rendered from the selected figure-data source |
| `figure_data/main/` | Frozen direct inputs for Figure 2 |
| `figure_data/supplementary/` | Frozen direct inputs for Supplementary Figures 2-5 |

The current record-set counts used by the packaged analysis are:

- canonical source: 1,448 records;
- full `D_0` record set: 985 records from 187 earthquakes;
- full `T_{90}` record set: 951 records from 187 earthquakes;
- `D_0`-aligned `T_{90}` record set: 613 records from 120 earthquakes.

## 7. Reproducibility notes

- Earthquake-level bootstrap resampling retains all records associated with a
  sampled earthquake.
- Gaussian restoration-time coefficient and factor intervals and two-sided
  P values use earthquake-cluster-robust CR2 covariance with Satterthwaite
  degrees of freedom, including every fixed-model sensitivity specification.
  Earthquake random-intercept sensitivities use model-based Satterthwaite
  inference. Factor source tables retain the inference method and degrees of
  freedom; bootstrap-count fields contain `NaN` for these coefficient-based
  intervals.
- Fitted R2 statistics use percentile earthquake-cluster bootstrap intervals.
  Each statistic requires the requested number of valid refits (2,000 by
  default); failures are logged and additional samples are drawn independently
  for the affected statistic. The two aligned models share each resample.
  Supplementary checkpoint reuse checks the input, configuration and code
  identity. Main and supplementary primary R2 intervals share the main export.
- The configured bootstrap count and random seed are stored in
  `contextual_config.json`.
- LOEO-CV predictions are generated with the held-out earthquake excluded from
  model fitting and are evaluated from the saved out-of-fold predictions.
  Their intervals resample the completed earthquake-held-out predictions,
  using the same resampled earthquakes for the two aligned models.
  Repeated baseline specifications reuse the primary interval after checking
  record IDs, earthquake IDs, response values, configuration and every
  out-of-fold prediction. Other sensitivity specifications retain their own
  bootstrap intervals.
- `outputs/` contains regenerable analysis results and rendered figures and is
  empty in the release package. The frozen files under `figure_data/` are
  retained as direct inputs for figure rendering.
