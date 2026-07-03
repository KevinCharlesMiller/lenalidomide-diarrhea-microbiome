## Lenalidomide-Related Diarrhea is Associated with Gut Microbiota Dysbiosis and Disruption of the Bile Acid Pool
## Analysis Code

Code to reproduce the figures and tables in *"Lenalidomide-Related Diarrhea is
Associated with Gut Microbiota Dysbiosis and Disruption of the Bile Acid Pool."*
This repository is **code only** — no data are committed.

## Data

Place the de-identified inputs in a local `data/` folder:

- `metadata.csv`, `bile_acids.csv`, `ffq.csv`, `picrust_ko.csv` — journal **Supplementary Data**.
- `len.rds` — de-identified phyloseq object, **available from the corresponding author on request**. Required to run microbiome analyses as written.
- Genus-level microbiome feature table available in journal **Supplementary Data**.
- Raw 16S rRNA reads available at NCBI SRA, BioProject PRJNA1484682

Samples are coded `MM##-<letter>` (`MM##` = subject); the patient-level grouping
variable is `SubjectID`.

## Run

Requires **R ≥ 4.4**.

```r
install.packages(c("readr","dplyr","tidyr","tibble","forcats","stringr","purrr",
  "ggplot2","patchwork","RColorBrewer","scales","colorspace","gghalves","ggpubr",
  "ggnewscale","ggrepel","vegan","permute","lme4","lmerTest","emmeans","car",
  "geepack","gtsummary","gt","conflicted"))
if (!require("BiocManager")) install.packages("BiocManager")
BiocManager::install("phyloseq")
install.packages("FLORAL")

source("run_all.R")          # panels collected into figures/all_figures.pdf
```

`run_all.R` sources `R/00_setup.R` (data + shared objects) then bundles scripts per each
figure. Individual scripts can be run after sourcing `R/00_setup.R` first.
Note, FLORAL (Supp Fig 2C–E) can take some time to run.

## Scripts

| Script | Produces |
|---|---|
| `R/00_setup.R` | shared setup (sourced by all) |
| `R/figure1.R` | Figure 1 (1A–1E); Supp Fig 1A–1D |
| `R/figure2.R` | Figure 2 (2A–2E); Supp Fig 2A–2E |
| `R/figure3.R` | Figure 3 (3A–3C); Supp Fig 3A–3D, Supp Fig 4C/4D, Supplemental Table 9 |
| `R/figure4.R` | Figure 4 (4A–4G); Supp Fig 4A/4B |
| `R/figure5.R` | Figure 5 (5A/5B); Supp Fig 5A–5D |
| `R/gee_models.R` | Supplemental Tables 5–7 (GEE models) |


## License

MIT — see `LICENSE`.

Segments of code for analysis generated with Claude Opus 4.6
