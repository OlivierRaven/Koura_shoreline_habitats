# Littoral habitat structure influences freshwater crayfish populations in five volcanic lakes of Aotearoa New Zealand

Olivier V. Raven | olivier.raven@icloud.com

**DOI:** To be added upon acceptance

## Overview

This repository contains the full reproducible workflow supporting the manuscript above. It includes all code used for data processing, statistical analyses, and figure generation, derived datasets, generated figures, and peer review documents. Raw field data are not tracked in this repository and will be archived in an open-access data repository upon acceptance.

**Key findings:**

- Kōura occupancy increased with shoreline habitat complexity, particularly coarser substrate types and greater riparian vegetation cover
- Elevated summer surface water temperatures were negatively associated with kōura occupancy and abundance, indicating reduced habitat suitability during warm periods
- pH was positively correlated with kōura abundance and biomass
- Composition and physical structure of littoral habitats are important drivers of kōura distribution and abundance
- Targeted protection and restoration of coarse substrate shoreline habitats may provide a practical conservation response, since climatic drivers (warmer water, reduced dissolved oxygen) are difficult to manage directly

## Repository structure

    .
    +-- data/
    |   +-- raw/          # Raw field data (Natural_habitat.xlsx; archived on acceptance)
    |   +-- derived/      # Processed CSVs used in analysis
    +-- outputs/          # Generated figures and tables
    +-- manuscript/       # Manuscript versions, cover letter, reviewer responses
    +-- images/           # Site photographs
    +-- references/       # Bibliography (.bib) and citation style files
    +-- scripts/          # Site selection and other supporting scripts
    +-- docs/             # Rendered HTML output (GitHub Pages)
    +-- _quarto.yml       # Quarto project configuration
    +-- deploy.R          # Post-render script that publishes docs/ to GitHub Pages
    +-- analysis.qmd      # Full statistical analysis notebook
    +-- index.qmd         # Manuscript

## Reproducing the analysis

This project uses [Quarto](https://quarto.org/) and R.
All analyses are contained in `analysis.qmd`.

### Requirements

- R >= 4.3
- Quarto >= 1.4
- R packages: `mgcv`, `gratia`, `glmmTMB`, `pROC`, `caret`, `performance`,
  `kableExtra`, `emmeans`, `patchwork`, `tidyverse`, `readxl`

### Install R packages

```r
install.packages(c("mgcv", "gratia", "glmmTMB", "pROC", "caret",
                   "performance", "kableExtra", "emmeans",
                   "patchwork", "tidyverse", "readxl"))
```

### Render the analysis

```bash
quarto render analysis.qmd
```

The rendered HTML analysis is available at:
<https://olivierraven.github.io/Koura_shoreline_habitats/>

## Data availability

### Raw data

The primary raw dataset (`data/raw/Natural_habitat.xlsx`) is tracked in this repository and will be archived in an open-access repository upon manuscript acceptance.

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

### Outputs

All generated figures and tables are tracked in `outputs/`:

- `Koura_plots.png` — Fig. 2: kōura presence, CPUE and BPUE across lakes
- `Occupancy_panels.png` — Fig. 3: occupancy GAMM results and model performance
- `CPUE_panel_positives.png` — Fig. 4: CPUE GAMM results and model performance
- `BPUE_panel_positives.png` — Fig. 5: BPUE GAMM results and model performance
- `length_weight_plot.png` — Fig. S1: length-weight relationship
- `Raw_data_plots.png` — Fig. S2: raw data underlying modelled relationships
- `EnvBio_summary_table.csv` — Table S1: environmental and biotic summary
- `table_gam_all_models.csv` — GAM model results for all three models

## Funding

This research was supported by the Fish Futures programme funded through a Ministry of Business, Innovation and Employment grant (CAWX2101) with additional funding provided by the Bay of Plenty Regional Council under the Toihuarewa Waimāori - Bay of Plenty Regional Council Chair in Lake and Freshwater Science programme.

## Acknowledgements

We thank Te Arawa Lakes Trust and Te Komiti Whakahaere for the opportunity to study kōura in the Rotorua Te Arawa lakes. We are grateful to Soweeta Fort-D'ath and William Anaru (Te Arawa Lakes Trust), and Tihini Grant (Ngāti Pikiao), for facilitating access to the lakes and supporting the field programme. We appreciate Joe Butterworth's assistance with fieldwork, and the support of Andy Bruere and the Bay of Plenty Regional Council in helping to fund the fieldwork. We thank Calum MacNeil and three anonymous reviewers for feedback on an earlier version of this manuscript.

## Ethical approval

This study did not involve experimentation on humans or animals. Kōura and fish were captured using standard fisheries sampling methods. Permission to conduct field sampling and handle aquatic fauna was granted by Te Arawa Lakes Trust and Te Komiti Whakahaere. Formal institutional ethical approval was not required.

## Peer review transparency

Reviewer comments and author responses are included in the `manuscript/` folder in the interest of open and transparent science. Reviewer identities remain anonymous as per journal policy.

- `cover_letter.pdf` — cover letter submitted to Hydrobiologia
- `reviewer_responses.pdf` — responses to reviewer comments

## Citation

*To be added upon acceptance.*

## Licence

Code: MIT License
Data: CC BY 4.0 (<https://creativecommons.org/licenses/by/4.0/>)
