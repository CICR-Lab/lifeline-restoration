# Module 00 - Dataset

## 1. Overview

This module documents access to the external canonical dataset used by the
downstream analyses and generates Figure 1, which summarizes the spatial and
record-level coverage of the dataset. It is the first stage of the project
workflow:

```text
00_Dataset
  ├── standardized dataset -> 01_ContextualAnalysis
  └── standardized dataset -> 02_EmergentSimplicity
                              └── Figure 3 source data -> 03_HierarchicalRepairModel
```

The module includes the Figure 1 plotting entry point and instructions for
obtaining its dataset and map inputs.

## 2. Requirements

- MATLAB.
- Authorized access to the canonical dataset described below.
- Mapping Toolbox for `shaperead` and `geoshow` when reproducing Figure 1.
- Natural Earth Admin 0 - Countries, version 5.1.1, for the map boundary data.

## 3. Inputs

### External canonical dataset

The canonical `standardized_dataset.mat` is hosted at:

<https://doi.org/10.5281/zenodo.23113551>

Access is restricted. Request access through the Zenodo record. After approval
from the authors, download the file and place it at:

```text
00_Dataset/standardized_dataset.mat
```

This MAT file is the shared input for Figure 1 and the analyses in Modules 01
and 02.

### External map input

Download Natural Earth **Admin 0 - Countries** version 5.1.1 from:

https://www.naturalearthdata.com/downloads/10m-cultural-vectors/10m-admin-0-countries/

Keep the complete shapefile set together under either
`ne_10m_admin_0_countries/` or `natural_earth_vector/` in this module. The
required files include `.shp`, `.dbf`, `.shx` and `.prj`.

## 4. Directory structure

```text
00_Dataset/
├── README.md
├── standardized_dataset.mat        # download from Zenodo after access approval
├── run_figure_1.m
├── ne_10m_admin_0_countries/         # downloaded Natural Earth map input
└── outputs/                        # generated Figure 1 files
```

The authorized downloaded MAT file and the Natural Earth shapefile are inputs. The
`outputs/` directory contains generated files and can be recreated.

## 5. Running the module

### Prepare the canonical dataset

Request access at the Zenodo DOI above. After approval, download
`standardized_dataset.mat` into this module directory and confirm that it is
present before running the workflow.

### Generate Figure 1

From the `00_Dataset` directory, run:

```matlab
run_figure_1()
```

The entry point loads the standardized dataset, summarizes events, regions and
lifeline-system records, reads the Natural Earth boundary data and exports
Figure 1.

## 6. Outputs

| Output | Use |
|---|---|
| `outputs/Figure1.fig` | Editable MATLAB Figure 1 |
| `outputs/Figure1.png` | Raster Figure 1 |
| `outputs/Figure1.pdf` | Vector Figure 1 |

## 7. Reproducibility notes

- Keep one canonical file named `standardized_dataset.mat` for the
  downstream workflow.
- Use the canonical file obtained through the stated Zenodo record after access
  approval.
- Natural Earth boundary data are external inputs and should be downloaded at
  the documented version and retained with the required companion files.
- Figure 1 is generated from `event_ISO` and the standardized record fields in
the MAT file. If the authorized dataset changes, Figure 1 and the dependent
  analyses should be regenerated.
