#####
# run_all.R — reproduce all figures and tables in manuscript order.
#####

source("R/00_setup.R")     # data, genus table, alpha-diversity, ratios, PICRUSt CPM, helpers

pdf("figures/all_figures.pdf", width = 9, height = 7, onefile = TRUE)

tryCatch({
  source("R/figure1.R")    # Figure 1 (1C-1E) + Supp Fig 1C-1D
  source("R/figure2.R")    # Figure 2 (2A-2E) + Supp Fig 2A-2E
  source("R/figure3.R")    # Figure 3 (3A-3C) + Supp Fig 3A-3D + Supp Fig 4C-4D + Supplemental Table 9
  source("R/figure4.R")    # Figure 4 (4B-4G) + Supp Fig 4A/4B
  source("R/figure5.R")    # Figure 5 (5A-5B) + Supp Fig 5A-5D
  source("R/gee_models.R") # Supplemental Tables 5-7
}, finally = {
  grDevices::dev.off()
})

message("\nAll scripts complete. Plots written to figures/all_figures.pdf")
