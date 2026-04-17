# Littoral habitat structure drives freshwater crayfish populations

This repository contains the full reproducible workflow supporting the manuscript:

**"Littoral habitat structure influences freshwater crayfish populations in five volcanic lakes of Aotearoa New Zealand"**

The repository includes all code used for data processing, statistical analyses, and figure generation, derived datasets, generated figures, and peer review documents. Raw field data are not tracked in this repository and will be archived in an open-access data repository upon acceptance.

---

## Repository structure
```
03_natural_habitat_monitoring/
├── scripts/
│   ├── index.qmd              # Main Quarto document (analysis + figures)
│   ├── Site_selection.R       # Site selection script
│   └── _quarto.yml            # Quarto project settings
├── data/
│   ├── raw/                   # Raw field data (tracked)
│   └── derived/               # Processed datasets used in analysis (tracked)
├── outputs/                   # Generated figures and tables (tracked)
├── manuscript/                # Manuscript versions (tracked)
│   ├── Cover_letter.docx
│   ├── Manuscript_v1_preprint.docx  # Before peer review
│   ├── Manuscript_v2_revised.docx   # After peer review
│   ├── Manuscript_v2_revised_track_changes
│   └── Response_to_reviewers.docx              
└── README.md                        # This file
```

---

## Data availability

### Raw data
The primary raw dataset (`data/raw/Natural_habitat.xlsx`) is tracked in this repository. And it will be archived in an open-access repository upon manuscript acceptance.

- **Data repository:** *To be added upon acceptance*
- **DOI:** *To be added upon acceptance*

The raw Excel workbook contains original field observations across multiple worksheets, including site information, environmental variables, kōura catch data, macrophytes, and macroinvertebrates.

### Derived data
Derived CSV files are tracked in this repository and are sufficient to reproduce all statistical analyses and figures without the raw Excel file:

- `Site_info.csv` — site metadata
- `Monitoring_data.csv` — environmental variables per site
- `Weed_data.csv` — aquatic vegetation cover data
- `Fish_data.csv` — fish catch data including kōura
- `Macroinvertebrates.csv` — macroinvertebrate count data
- `Monitoring_CPUE_data.csv` — fully processed dataset used for modelling
- `habitat_classification.csv` — classified habitat types per site

---

## Code availability

All analyses were conducted in **R** and implemented in a reproducible **Quarto** document (`scripts/index.qmd`), which:

- Loads and filters raw data
- Calculates CPUE and BPUE metrics
- Fits GAM statistical models
- Generates all manuscript figures

Supporting site selection code is in `scripts/Site_selection.R`.

All code relies only on publicly available R packages.

---

## How to reproduce the analysis

### Full pipeline (requires raw Excel file)
1. Clone this repository:
```bash
git clone https://github.com/OlivierRaven/Koura_shoreline_habitats.git
```
2. Obtain `Natural_habitat.xlsx` from the data repository (DOI above) and place it in `data/raw/`
3. Open `scripts/index.qmd` in RStudio or Positron
4. Run the Quarto document from the top

### Analysis and figures only (using derived CSVs)
1. Clone this repository
2. Open `scripts/index.qmd` in RStudio or Positron
3. Navigate to the reproducibility note section in the script, uncomment the `read.csv` lines, and run from there — all analyses and figures can be reproduced without the raw Excel file

---

## Outputs

All generated figures and tables are tracked in `outputs/`:

- `Koura_plots.png` — Fig. 2: kōura presence, CPUE and BPUE across lakes
- `Occupancy_panels.png` — Fig. 3: occupancy GAMM results and model performance
- `CPUE_panel_positives.png` — Fig. 4: CPUE GAMM results and model performance
- `BPUE_panel_positives.png` — Fig. 5: BPUE GAMM results and model performance
- `length_weight_plot.png` — Fig. S1: length-weight relationship
- `Raw_data_plots.png` — Fig. S2: raw data underlying modelled relationships
- `EnvBio_summary_table.csv` — Table S1: environmental and biotic summary
- `table_gam_all_models.csv` — GAM model results for all three models

---

## Peer review transparency

Reviewer comments and author responses are included in the `peer_review/` folder in the interest of open and transparent science. Reviewer identities remain anonymous as per journal policy.

- `cover_letter.pdf` — cover letter submitted to Hydrobiologia
- `reviewer_responses.pdf` — responses to reviewer comments

---

## Dependencies

All R packages required are loaded at the top of `scripts/index.qmd`. The analysis was conducted using:

- **R** (≥ 4.3.0)
- **Quarto** (≥ 1.4)

Key packages include: `mgcv`, `gratia`, `glmmTMB`, `pROC`, `ggplot2`, `patchwork`, `tidyverse`

---

## Citation

*To be added upon acceptance.*

---

## Contact

Olivier Raven  
PhD Candidate, Freshwater Ecology  
University of Waikato, Aotearoa New Zealand  
GitHub: [OlivierRaven](https://github.com/OlivierRaven)