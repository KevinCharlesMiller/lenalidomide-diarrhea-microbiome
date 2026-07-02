#####
# Supplemental Tables 4-6 (GEE models)
#####

if (!exists("df_taxa_genus")) source("R/00_setup.R")

# ===== SUPPLEMENTAL TABLE 4 =====
# GEE: α-diversity ~ SampleCategory + Abx + autoHCT + Sex + Cohort

dat_alpha <- metadata %>%
  dplyr::filter(
    SampleCategory %in% c("Baseline", "OnLen_control", "OnLen_diarrhea"),
    !is.na(invsimpson), invsimpson > 0,
    !is.na(AbxInLast100Days),
    !is.na(autoHCT),
    !is.na(Sex),
    !is.na(Cohort)
  ) %>%
  dplyr::mutate(
    log10_alpha = log10(invsimpson),
    Sex         = factor(Sex, levels = c("Female", "Male")),
    Cohort         = as.factor(Cohort),
    SubjectID         = as.character(SubjectID)
  ) %>%
  dplyr::arrange(SubjectID) %>%
  dplyr::mutate(cluster_id = as.integer(factor(SubjectID, levels = unique(SubjectID))))

cat("\n=== Supplemental Table 4 ===\n")
cat("N samples:", nrow(dat_alpha), "| N patients:", n_distinct(dat_alpha$SubjectID), "\n")

gee_alpha <- geeglm(
  log10_alpha ~ SampleCategory + AbxInLast100Days + autoHCT + Sex + Cohort,
  id     = cluster_id,
  data   = dat_alpha,
  family = gaussian,
  corstr = "exchangeable"
)

tab4 <- tbl_regression(
  gee_alpha,
  intercept = FALSE,
  label = list(
    SampleCategory   ~ "Sample Group",
    AbxInLast100Days ~ "Antibiotics within the last 100 days",
    autoHCT          ~ "Prior AHCT",
    Sex              ~ "Sex",
    Cohort              ~ "Study cohort"
  )
) %>%
  modify_caption("Supplemental Table 4. Multivariable GEE model of α-diversity (log10 Inverse Simpson)") %>%
  bold_p()

print(tab4)
as_gt(tab4)

# ===== SUPPLEMENTAL TABLE 5 =====
# GEE: log10(CA/DCA) ~ log10(α-diversity) + SampleCategory (+/- Abx + autoHCT + Sex + Cohort)

dat_cadca <- df_ratio %>%
  dplyr::left_join(
    metadata %>% dplyr::select(Sample.ID, invsimpson, AbxInLast100Days, autoHCT, Sex, Cohort),
    by = "Sample.ID"
  ) %>%
  dplyr::filter(
    is.finite(invsimpson), invsimpson > 0,
    is.finite(CA_DCA_log10),
    !is.na(SubjectID),
    !is.na(SampleCategory),
    !is.na(AbxInLast100Days),
    !is.na(autoHCT),
    !is.na(Sex),
    !is.na(Cohort)
  ) %>%
  dplyr::mutate(
    log10_alpha    = log10(invsimpson),
    SampleCategory = factor(SampleCategory,
                            levels = c("Baseline", "OnLen_control", "OnLen_diarrhea")),
    Sex            = factor(Sex, levels = c("Female", "Male")),
    Cohort            = as.factor(Cohort),
    ratio          = CA_DCA_log10,
    SubjectID            = as.character(SubjectID)
  ) %>%
  dplyr::arrange(SubjectID) %>%
  dplyr::mutate(cluster_id = as.integer(factor(SubjectID, levels = unique(SubjectID))))

cat("\n=== Supplemental Table 5 ===\n")
cat("N samples:", nrow(dat_cadca), "| N patients:", n_distinct(dat_cadca$SubjectID), "\n")

gee_cadca_unadj <- geeglm(
  ratio ~ log10_alpha + SampleCategory,
  id     = cluster_id,
  data   = dat_cadca,
  family = gaussian,
  corstr = "exchangeable"
)

tab5_unadj <- tbl_regression(
  gee_cadca_unadj,
  intercept = FALSE,
  label = list(
    log10_alpha    ~ "α-diversity (log10-transformed)",
    SampleCategory ~ "Sample Group"
  )
) %>%
  modify_caption("Unadjusted Model") %>%
  bold_p()

gee_cadca_adj <- geeglm(
  ratio ~ log10_alpha + SampleCategory + AbxInLast100Days + autoHCT + Sex + Cohort,
  id     = cluster_id,
  data   = dat_cadca,
  family = gaussian,
  corstr = "exchangeable"
)

tab5_adj <- tbl_regression(
  gee_cadca_adj,
  intercept = FALSE,
  label = list(
    log10_alpha      ~ "α-diversity (log10-transformed)",
    SampleCategory   ~ "Sample Group",
    AbxInLast100Days ~ "Antibiotics within the last 100 days",
    autoHCT          ~ "Prior AHCT",
    Sex              ~ "Sex",
    Cohort              ~ "Study cohort"
  )
) %>%
  modify_caption("Adjusted Model") %>%
  bold_p()

tab5_merged <- tbl_merge(
  tbls = list(tab5_unadj, tab5_adj),
  tab_spanner = c("Unadjusted Model", "Adjusted Model")
) %>%
  modify_caption("Supplemental Table 5. GEE models of CA/DCA ratio (log10-transformed)")

print(tab5_merged)
as_gt(tab5_merged)

# ===== SUPPLEMENTAL TABLE 6 =====
# GEE: log10(CDCA/LCA) ~ log10(α-diversity) + SampleCategory (+/- Abx + autoHCT + Sex + Cohort)

dat_cdcalca <- df_ratio %>%
  dplyr::left_join(
    metadata %>% dplyr::select(Sample.ID, invsimpson, AbxInLast100Days, autoHCT, Sex, Cohort),
    by = "Sample.ID"
  ) %>%
  dplyr::filter(
    is.finite(invsimpson), invsimpson > 0,
    is.finite(CDCA_LCA_log10),
    !is.na(SubjectID),
    !is.na(SampleCategory),
    !is.na(AbxInLast100Days),
    !is.na(autoHCT),
    !is.na(Sex),
    !is.na(Cohort)
  ) %>%
  dplyr::mutate(
    log10_alpha    = log10(invsimpson),
    SampleCategory = factor(SampleCategory,
                            levels = c("Baseline", "OnLen_control", "OnLen_diarrhea")),
    Sex            = factor(Sex, levels = c("Female", "Male")),
    Cohort            = as.factor(Cohort),
    ratio          = CDCA_LCA_log10,
    SubjectID            = as.character(SubjectID)
  ) %>%
  dplyr::arrange(SubjectID) %>%
  dplyr::mutate(cluster_id = as.integer(factor(SubjectID, levels = unique(SubjectID))))

cat("\n=== Supplemental Table 6 ===\n")
cat("N samples:", nrow(dat_cdcalca), "| N patients:", n_distinct(dat_cdcalca$SubjectID), "\n")

gee_cdcalca_unadj <- geeglm(
  ratio ~ log10_alpha + SampleCategory,
  id     = cluster_id,
  data   = dat_cdcalca,
  family = gaussian,
  corstr = "exchangeable"
)

tab6_unadj <- tbl_regression(
  gee_cdcalca_unadj,
  intercept = FALSE,
  label = list(
    log10_alpha    ~ "α-diversity (log10-transformed)",
    SampleCategory ~ "Sample Group"
  )
) %>%
  modify_caption("Unadjusted Model") %>%
  bold_p()

gee_cdcalca_adj <- geeglm(
  ratio ~ log10_alpha + SampleCategory + AbxInLast100Days + autoHCT + Sex + Cohort,
  id     = cluster_id,
  data   = dat_cdcalca,
  family = gaussian,
  corstr = "exchangeable"
)

tab6_adj <- tbl_regression(
  gee_cdcalca_adj,
  intercept = FALSE,
  label = list(
    log10_alpha      ~ "α-diversity (log10-transformed)",
    SampleCategory   ~ "Sample Group",
    AbxInLast100Days ~ "Antibiotics within the last 100 days",
    autoHCT          ~ "Prior AHCT",
    Sex              ~ "Sex",
    Cohort              ~ "Study cohort"
  )
) %>%
  modify_caption("Adjusted Model") %>%
  bold_p()

tab6_merged <- tbl_merge(
  tbls = list(tab6_unadj, tab6_adj),
  tab_spanner = c("Unadjusted Model", "Adjusted Model")
) %>%
  modify_caption("Supplemental Table 6. GEE models of CDCA/LCA ratio (log10-transformed)")

print(tab6_merged)
as_gt(tab6_merged)

cat("\n=== All GEE models with Cohort adjustment complete ===\n")
