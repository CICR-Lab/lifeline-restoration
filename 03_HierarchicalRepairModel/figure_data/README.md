# Hierarchical repair-model figure data

The `full/` directory contains the direct plotting inputs for the full
publication run. These files are copied from the corresponding validated
calculation outputs and are separate from restartable checkpoints and other
intermediate results under `outputs/full/`.

- `Analytical/Analytical_results.mat`: Figure 4b-d and Supplementary Figures S13-S16.
- `Simulation/source_data/`: Figures 4e-g and Supplementary Figure S17.
- `Interdependent/Interdependent_results.mat`: Figure 4h and Supplementary Figure S18.
- `SupplementaryTables/Supplementary_Table_S8.xlsx`: finite-network input table.

The plotting functions read the frozen files in this directory for the full
publication mode with `DataSource="packaged"`. Each plotting call uses the single
source selected by `DataSource`. Generated calculation results and checkpoints
are stored under `outputs/<mode>/`; generated direct plotting inputs are collected under
`outputs/figure_data/<mode>/`. For quick-mode validation, first run
`run_hierarchical_repair_analysis('quick')`, then plot with `Mode="quick"` and
`DataSource="generated"` using the inputs in `outputs/figure_data/quick/`.
