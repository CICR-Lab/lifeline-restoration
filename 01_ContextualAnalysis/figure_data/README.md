# Contextual figure data

This directory contains the frozen, direct CSV inputs for the contextual
figures. It is part of the release package and is kept outside `outputs/`,
which contains regenerable run results.

- `main/` contains the four CSV files used by Figure 2a-d.
- `supplementary/` contains the CSV files used by Supplementary Figures 2-5.

The root-level plotting entry points read these files without refitting models
or rerunning bootstrap analyses. The CSV-to-figure mapping is:

| Internal CSV prefix | Published supplementary figure | Content |
|---|---|---|
| `s1a`, `s1b`, `s1c` | S2 | Fitted R2 sensitivity |
| `s3a`, `s3b` | S3 | Adjusted T90 factor sensitivity |
| `s2a`, `s2b` | S4 | LOEO-CV performance and sensitivity |
| `s4` | S5 | D0 calibration and T90 diagnostics |

Column names and units follow the definitions in the module README and the
manuscript Methods.

`s3a_system_factor_sensitivity.csv` and `s3b_D0_factor_sensitivity.csv`
use CR2/Satterthwaite inference for fixed Gaussian models and model-based
Satterthwaite inference for earthquake random-intercept models. The tables
include P values, contrast-specific degrees of freedom and `ci_method`.
For these coefficient-based intervals, bootstrap-count fields contain `NaN`.
D0 sensitivity factors compare D0 = 0.10 with 0.20 for all representations.
These internal s3 files support Supplementary Figure S3 in the final manuscript.

The `log2_d0` and `lognormal` rows in
`s2b_paired_LOEO_delta_sensitivity.csv` repeat the primary paired model.
Their point estimates and confidence intervals share the primary result in
`s2a_primary_LOEO_performance.csv`, reusing its bootstrap interval.
The analysis checks record and prediction identity before reusing an interval.

After a local analysis, generated direct figure inputs are written under
`../outputs/analysis/figure_data/` for Figure 2 and
`../outputs/supplementary/figure_data/` for Supplementary Figures 2-5. The
plotting entries select between retained and generated inputs with
`DataSource="packaged"` or `DataSource="generated"`. Each plotting call uses the
single source selected by `DataSource`.
