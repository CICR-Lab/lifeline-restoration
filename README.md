# Code Package - Emergent simplicity in post-earthquake lifeline recovery

This package contains the MATLAB workflows and retained direct figure inputs
accompanying the manuscript *Emergent simplicity in post-earthquake lifeline
recovery*.

The package is organized into four modules. Module 00 documents access to the
canonical standardized dataset and generates Figure 1; Modules 01 and 02 are independent
analysis branches; Module 03 uses the empirical source data exported by Module
02.

## 1. Package Structure

```text
code/
├── README.md
├── 00_Dataset/
│   ├── README.md
│   ├── run_figure_1.m
│   └── outputs/
├── 01_ContextualAnalysis/
│   ├── README.md
│   ├── contextual_config.json
│   ├── run_contextual_analysis.m
│   ├── run_figure_2.m
│   ├── run_supplementary_contextual_analysis.m
│   ├── Code/
│   ├── figure_data/
│   └── outputs/
├── 02_EmergentSimplicity/
│   ├── README.md
│   ├── run_emergent_simplicity_analysis.m
│   ├── run_figure_3.m
│   ├── run_supplementary_emergent_simplicity_analysis.m
│   ├── run_supplementary_emergent_simplicity.m
│   ├── figure_data/
│   ├── Code/
│   └── outputs/
└── 03_HierarchicalRepairModel/
    ├── README.md
    ├── run_hierarchical_repair_analysis.m
    ├── run_figure_4.m
    ├── run_supplementary_hierarchical_repair_model.m
    ├── figure_data/
    ├── Code/
    └── outputs/
```

The module-level READMEs are the detailed operating instructions. This file
defines the project-wide workflow, input boundaries and manuscript mapping.

## 2. Requirements And External Inputs

- MATLAB is required for all modules.
- Module 00 requires Mapping Toolbox for Figure 1.
- Module 01 requires Statistics and Machine Learning Toolbox and Optimization
  Toolbox. Its supplementary bootstrap runs with 100 or more replicates require
  Parallel Computing Toolbox; these runs do not have a serial fallback.
- Module 02 requires Optimization Toolbox and Statistics and Machine Learning
  Toolbox for its fitting and supplementary analyses.
- Module 03 requires Optimization Toolbox and Statistics and Machine Learning
  Toolbox. Parallel Computing Toolbox is optional for parallel full simulations.
- Figure 1 additionally requires Natural Earth Admin 0 - Countries boundary
  data, version 5.1.1. Download the complete shapefile set (`.shp`,
  `.dbf`, `.shx` and `.prj`) and place it in the Module 00 directory
  as described in its README. Obtain the dataset from:
  <https://www.naturalearthdata.com/downloads/10m-cultural-vectors/10m-admin-0-countries/>

The module READMEs document module-specific software and input requirements.

## 3. Inputs And Data Boundaries

### External shared dataset

The canonical `standardized_dataset.mat` is available through Zenodo:

<https://doi.org/10.5281/zenodo.23113551>

Access is restricted. Request access through the Zenodo record and wait for
approval from the authors. After approval, download `standardized_dataset.mat`
and place it at:

```text
00_Dataset/standardized_dataset.mat
```

Figure 1 and the analyses in Modules 01 and 02 read this dataset. The
`DataSource="packaged"` plotting workflows in Modules 01-03 read the supplied
direct inputs under each module's `figure_data/` directory.

### Module 02 to Module 03 input

The retained empirical source workbook is:

```text
02_EmergentSimplicity/figure_data/main/SourceData_Figure3.xlsx
```

Module 03 uses this workbook to parameterize the hierarchical repair model.

### Figure-data retention rule

`figure_data/` contains the supplied direct inputs for drawing the main and
supplementary figures.

`outputs/` is a runtime directory for regenerable calculations, checkpoints,
manifests, figures and tables. It is intentionally empty in the release
package and is populated only when a workflow is run.

The plotting entry points in Modules 01-03 accept `DataSource="packaged"`
(the default) or `DataSource="generated"`. `packaged` reads the retained direct
inputs under the module's `figure_data/` and can redraw the supplied figures
without rerunning analyses. `generated` reads the direct figure inputs written
under `outputs/` by the corresponding local analysis; run the relevant main or
supplementary analysis first.
Each plotting call uses the single source selected by `DataSource`.
Figure 1 reads `00_Dataset/standardized_dataset.mat` and the Natural Earth
boundary data.

## 4. Path Setup

From MATLAB, set the local package root and run one module at a time:

```matlab
codeRoot = "path/to/code";  % replace with the local package path
cd(fullfile(codeRoot, "00_Dataset"));
```

Each module's public entry points are at its root. Module 01 adds its internal
`Code/` directory automatically. Module 02 and Module 03 manage their own
internal code paths through their public entry points.

## 5. Running The Workflows

### 5.1 Module 00: Dataset And Figure 1

After obtaining access to the Zenodo dataset and placing
`standardized_dataset.mat` in `00_Dataset/`, place the Natural Earth
shapefile set in the module directory and run:

```matlab
cd(fullfile(codeRoot, "00_Dataset"));
run_figure_1();
```

This generates Figure 1 under `00_Dataset/outputs/`.

### 5.2 Module 01: Contextual Analysis

To run the analysis and render figures from its generated results:

```matlab
cd(fullfile(codeRoot, "01_ContextualAnalysis"));
addpath(pwd);
run_contextual_analysis;
run_figure_2("DataSource", "generated");
run_supplementary_contextual_analysis("DataSource", "generated");
```

To redraw the retained figures without running the analysis:

```matlab
run_figure_2();
run_supplementary_contextual_analysis();
```

### 5.3 Module 02: Emergent Simplicity

To regenerate the analysis and render figures from its outputs:

```matlab
cd(fullfile(codeRoot, "02_EmergentSimplicity"));
run_emergent_simplicity_analysis();
run_supplementary_emergent_simplicity_analysis();
run_figure_3("DataSource", "generated");
run_supplementary_emergent_simplicity("DataSource", "generated");
```

To redraw the retained figures without running the analyses:

```matlab
run_figure_3();
run_supplementary_emergent_simplicity("DataSource", "packaged");
```

The supplementary analysis workflow uses the manual exception-class workbook
after the researcher has reviewed the candidate trajectories according to the
rules in the Module 02 README. Its generated source data and analysis results
are written under `outputs/`; the supplied direct figure inputs are stored under
`figure_data/`.

### 5.4 Module 03: Hierarchical Repair Model

Module 03 requires the retained Figure 3 source workbook from Module 02.
To run the full analysis and render figures from its generated results:

```matlab
cd(fullfile(codeRoot, "03_HierarchicalRepairModel"));
run_hierarchical_repair_analysis("full");
run_figure_4("Mode", "full", "DataSource", "generated");
run_supplementary_hierarchical_repair_model("Mode", "full", "DataSource", "generated");
```

To redraw the retained full-mode figures without running the analysis:

```matlab
run_figure_4("Mode", "full", "DataSource", "packaged");
run_supplementary_hierarchical_repair_model("Mode", "full", "DataSource", "packaged");
```

For a software and layout check, use the quick mode consistently:

```matlab
run_hierarchical_repair_analysis("quick");
run_figure_4("Mode", "quick", "DataSource", "generated");
run_supplementary_hierarchical_repair_model("Mode", "quick", "DataSource", "generated");
```

## 6. Workflow Dependencies

```text
Zenodo dataset (restricted access; download after author approval)
  -> 00_Dataset/standardized_dataset.mat
     -> Figure 1
     -> 01_ContextualAnalysis
     -> 02_EmergentSimplicity
        -> SourceData_Figure3.xlsx -> 03_HierarchicalRepairModel
```

Modules 01 and 02 can be run independently once the authorized dataset has
been placed in Module 00. Module 03 depends on Module 02's retained Figure 3
source workbook.

## 7. Manuscript Mapping And Outputs

| Manuscript item | Module and public entry point | Direct or principal output |
|---|---|---|
| Figure 1 | Module 00, `run_figure_1.m` | `00_Dataset/outputs/Figure1.*` |
| Figure 2 | Module 01, `run_figure_2.m` | `01_ContextualAnalysis/outputs/figures/main/` |
| Supplementary Figures 2-5 | Module 01, `run_supplementary_contextual_analysis.m` | `01_ContextualAnalysis/outputs/figures/supplementary/` |
| Figure 3 | Module 02, `run_figure_3.m` | `02_EmergentSimplicity/outputs/figures/main/` |
| Supplementary Figures 6-12 | Module 02, `run_supplementary_emergent_simplicity.m` | `02_EmergentSimplicity/outputs/figures/supplementary/` |
| Supplementary Tables S3-S7 | Module 02, `run_supplementary_emergent_simplicity_analysis.m` | `02_EmergentSimplicity/outputs/supplementary/source_data/` |
| Figure 4b-h | Module 03, `run_figure_4.m` | `03_HierarchicalRepairModel/outputs/figures/main/` |
| Supplementary Figures S13-S18 | Module 03, `run_supplementary_hierarchical_repair_model.m` | `03_HierarchicalRepairModel/outputs/figures/supplementary/` |
| Supplementary Table S8, Inputs sheet | Module 03, `run_hierarchical_repair_analysis` | `outputs/<mode>/SupplementaryTables/` |

The computational scope is the plotted results listed above. Figure 4a and
Supplementary Figure S1 are manually prepared schematics. Supplementary Table
S1 is assembled from the relevant analysis-specific sample counts.

Where a figure has a retained direct plotting input, it is stored under the
relevant module `figure_data/` directory. Generated figures, intermediate
results and table workbooks are regenerable runtime files written under
`outputs/`.

## License

This code is released under the MIT License. See the `LICENSE` file for details.
