# Complex littoral habitats enhance freshwater crayfish populations

This repository contains the full reproducible workflow supporting the manuscript:

**"Complex littoral habitats enhance freshwater crayfish populations in five volcanic lakes of Aotearoa New Zealand."**

The repository includes all code used for data processing, statistical analyses, and figure generation, as well as derived datasets used in the analyses. Raw field data are stored locally and will be archived in an open-access data repository upon acceptance.

---

## Repository structure

```
03_natural_habitat_monitoring/
├── scripts/
│   ├── index.qmd          # Main Quarto document (analysis + figures)
│   ├── Site_selection.R   # Site selection script
│   └── _quarto.yml        # Quarto project settings
├── data/
│   ├── raw/               # Raw field data (not tracked, stored locally/OSF)
│   └── derived/           # Processed datasets used in analysis
├── outputs/               # Generated figures and tables (not tracked)
├── manuscript/            # Manuscript drafts (not tracked)
└── README.md              # This file
```

---

## Data availability

### Raw data
Raw field data are stored locally and will be archived in an open-access repository upon manuscript acceptance.

- **Data repository:** *To be added upon acceptance*
- **DOI:** *To be added upon acceptance*

The primary raw dataset is:
- `data/raw/Natural_habitat.xlsx`  
  An Excel workbook containing original field observations across multiple worksheets, including site information, environmental variables, and kōura catch data, macrophytes, and macroinvertebrates.

### Derived data
Derived datasets used for statistical modelling and figure generation are created by the analysis pipeline and saved to `data/derived/`:

- `Monitoring_CPUE_data.csv` — processed catch-per-unit-effort data
- `habitat_classification.csv` — classified habitat types per site

These files can be regenerated at any time by running `scripts/index.qmd`.

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

1. Clone this repository:
   ```bash
   git clone https://github.com/OlivierRaven/Koura_shoreline_habitats.git
   ```

2. Obtain the raw data from the data repository (DOI above) and place it in `data/raw/`

3. Open `scripts/index.qmd` in RStudio or Positron

4. Run the Quarto document to reproduce all analyses and figures

---

## Dependencies

All R packages required are loaded at the top of `scripts/index.qmd`. The analysis was conducted using:

- **R** (≥ 4.3.0)
- **Quarto** (≥ 1.4)

---

## Citation

*To be added upon acceptance.*

---

## Contact

Olivier Raven  
PhD Candidate, Freshwater Ecology  
University of Waikato, Aotearoa New Zealand  
GitHub: [OlivierRaven](https://github.com/OlivierRaven)
