# Emergent-simplicity figure data

This directory contains the retained direct inputs for the public plotting
entries. The plotting entries read these files without requiring a preceding
analysis run.

## Main figure

`main/SourceData_Figure3.xlsx` is read by `run_figure_3.m` for Figure 3a-f and
is also the empirical input used by the hierarchical repair-model package.

## Supplementary figures

The supplementary folder contains the direct inputs for Supplementary Figures
S6-S12:

- `SourceData_FigS6.xlsx` for S6;
- `SourceData_FigS7.xlsx` for S7;
- `SourceData_FigS8_S11.xlsx` for the tabular inputs used by S8-S11;
- `unique3_records.mat` and `restoration_shape_ge5_fit_result.mat` for the
  observed trajectories and fitted curves drawn in S8-S11;
- `badfit_weibull_gompertz_loglogistic_ge6_candidates.xlsx` for the candidate
  rows and labels associated with the exception panels; and
- `SourceData_FigS12.xlsx` for S12.

The S8-S11 source data incorporate the confirmed manual classifications stored
in `../Code/exception_classes.xlsx`. The analysis workflows may regenerate
source data under `outputs/supplementary/source_data/` and direct plotting
inputs under `outputs/figure_data/`. The plotting entries select retained or
generated inputs with `DataSource="packaged"` or `DataSource="generated"`.
Each plotting call uses the single source selected by `DataSource`.
