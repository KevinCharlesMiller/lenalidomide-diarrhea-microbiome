#####
# Figure 4, Supplemental Figure 4
#####

if (!exists("df_taxa_genus")) source("R/00_setup.R")

# ===== FIGURE 4A =====
# Schematic — baseline sample stratification

# ===== FIGURE 4B =====
# Baseline α-diversity

meta_use <- metadata %>%
  dplyr::mutate(
    Sample.ID     = as.character(Sample.ID),
    LaterDiarrhea = factor(LaterDiarrhea, levels = c("No","Yes"))
  ) %>%
  dplyr::filter(!is.na(LaterDiarrhea), !is.na(invsimpson)) %>%
  dplyr::distinct(Sample.ID, .keep_all = TRUE)

alpha_df <- meta_use %>%
  dplyr::select(Sample.ID, LaterDiarrhea, InvSimpson = invsimpson) %>%
  dplyr::filter(!is.na(LaterDiarrhea), is.finite(InvSimpson)) %>%
  dplyr::mutate(
    LaterDiarrhea    = factor(LaterDiarrhea, levels = c("No","Yes")),
    InvSimpson_log10 = log10(InvSimpson)
  )

stopifnot(nrow(alpha_df) >= 4)

wilcx  <- stats::wilcox.test(InvSimpson_log10 ~ LaterDiarrhea, data = alpha_df, exact = FALSE)
p_val  <- wilcx$p.value
p_lab  <- sprintf("Wilcoxon p%s", ifelse(p_val < .001, "<0.001", sprintf("= %.3f", p_val)))

grp_cols <- c("No" = "#1B9E77", "Yes" = "#D55E00")

fig4_p_alpha_shift <- ggplot(alpha_df,
                             aes(x = LaterDiarrhea, y = InvSimpson_log10,
                                 fill = LaterDiarrhea, colour = LaterDiarrhea)) +
  geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.20, linewidth = 0.7) +
  geom_jitter(width = 0.12, size = 1.6, alpha = 0.7) +
  scale_fill_manual(values = grp_cols, guide = "none") +
  scale_colour_manual(values = grp_cols, guide = "none") +
  labs(x = NULL, y = "log10(Inverse Simpson)", title = "Alpha-diversity") +
  theme_classic(base_size = 13) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1)) +
  annotate("text", x = Inf, y = Inf, hjust = 1.02, vjust = 1.2,
           label = p_lab, fontface = "italic", size = 4) +
  coord_cartesian(clip = "off")

print(fig4_p_alpha_shift)

# ===== FIGURE 4C =====
# Baseline domination

THRESH <- 0.30
TOP_N  <- 10

X <- df_taxa_genus[rownames(df_taxa_genus) %in% metadata$Sample.ID, , drop = FALSE]

meta_dom <- metadata %>%
  dplyr::filter(Sample.ID %in% rownames(X)) %>%
  dplyr::transmute(
    Sample.ID     = as.character(Sample.ID),
    LaterDiarrhea = factor(LaterDiarrhea, levels = c("No","Yes"))
  ) %>%
  dplyr::filter(!is.na(LaterDiarrhea))

keep_ids_dom <- intersect(rownames(X), meta_dom$Sample.ID)
X        <- X[keep_ids_dom, , drop = FALSE]
meta_dom <- meta_dom %>% dplyr::filter(Sample.ID %in% keep_ids_dom)

row_tot <- rowSums(X, na.rm = TRUE)
X_rel   <- sweep(X, 1, ifelse(row_tot > 0, row_tot, 1), "/")

get_dom <- function(v, cn) {
  if (sum(v, na.rm = TRUE) <= 0) {
    return(c(genus = NA_character_, prop = NA_real_))
  }
  j <- which.max(v)
  c(genus = cn[j], prop = as.numeric(v[j]))
}

dom_mat <- t(apply(X_rel, 1, get_dom, cn = colnames(X_rel))) %>%
  as.data.frame()

dom_df <- dom_mat %>%
  dplyr::mutate(
    Sample.ID = rownames(X_rel),
    prop      = as.numeric(prop),
    genus     = ifelse(is.na(genus) | genus == "", "(unlabeled)", genus)
  ) %>%
  dplyr::left_join(meta_dom, by = "Sample.ID") %>%
  dplyr::filter(!is.na(LaterDiarrhea)) %>%
  dplyr::mutate(
    dominated = is.finite(prop) & prop >= THRESH
  )

top_genera <- dom_df %>%
  dplyr::filter(dominated) %>%
  dplyr::count(genus, sort = TRUE) %>%
  dplyr::slice_head(n = TOP_N) %>%
  dplyr::pull(genus)

comp_df <- dom_df %>%
  dplyr::mutate(
    genus2 = dplyr::if_else(
      dominated,
      dplyr::if_else(genus %in% top_genera, genus, "Other"),
      NA_character_
    )
  ) %>%
  dplyr::group_by(LaterDiarrhea) %>%
  dplyr::mutate(n_group = dplyr::n()) %>%
  dplyr::ungroup() %>%
  dplyr::filter(!is.na(genus2)) %>%
  dplyr::count(LaterDiarrhea, genus2, n_group, name = "n_dom_genus") %>%
  dplyr::mutate(frac = n_dom_genus / n_group) %>%
  tidyr::complete(LaterDiarrhea, genus2, fill = list(n_dom_genus = 0, frac = 0)) %>%
  dplyr::mutate(
    LaterDiarrhea = forcats::fct_relevel(LaterDiarrhea, c("No","Yes")),
    genus2        = forcats::fct_infreq(genus2)
  ) %>%
  dplyr::arrange(LaterDiarrhea, dplyr::desc(frac))

lvl       <- levels(comp_df$genus2)
pal_genus <- colorspace::qualitative_hcl(length(lvl), palette = "Set 2")
names(pal_genus) <- lvl
if ("Other" %in% names(pal_genus)) pal_genus["Other"] <- "#BDBDBD"

fig4_p_dom_comp <- ggplot2::ggplot(
  comp_df,
  ggplot2::aes(x = LaterDiarrhea, y = frac, fill = genus2)
) +
  ggplot2::geom_col(width = 0.7, color = "white", linewidth = 0.2) +
  ggplot2::scale_y_continuous(
    breaks = c(0, 0.25, 0.50, 0.75, 1),
    labels = scales::percent_format(accuracy = 1),
    limits = c(0, 1),
    expand = ggplot2::expansion(mult = c(0.02, 0.10))
  ) +
  ggplot2::scale_fill_manual(values = pal_genus, name = "Dominant genus") +
  ggplot2::labs(
    x     = NULL,
    y     = "Fraction of samples",
    title = sprintf("Domination composition (≥%s%% any genus)", THRESH * 100)
  ) +
  ggplot2::theme_classic(base_size = 13) +
  ggplot2::theme(
    axis.text.x     = ggplot2::element_text(angle = 20, hjust = 1),
    legend.position = "right"
  )

tab_frac <- with(dom_df, table(LaterDiarrhea, dominated))

p_frac <- tryCatch(
  stats::fisher.test(tab_frac)$p.value,
  error = function(e) NA_real_
)

if (!is.finite(p_frac) || is.na(p_frac)) {
  if (sum(tab_frac) > 0 && nrow(tab_frac) == 2 && ncol(tab_frac) == 2) {
    tab_cc <- tab_frac + 0.5
    p_frac <- tryCatch(
      stats::fisher.test(tab_cc)$p.value,
      error = function(e) NA_real_
    )
  }
}

if (!is.finite(p_frac) || is.na(p_frac)) {
  p_frac <- 1
}

lab_frac <- paste0(
  "Fisher (dominated %) p ",
  ifelse(p_frac < .001, "<0.001", sprintf("= %.3f", p_frac))
)

fig4_p_dom_comp_annot <- fig4_p_dom_comp +
  ggplot2::annotate(
    "text",
    x        = Inf,
    y        = Inf,
    hjust    = 1.02,
    vjust    = 1.2,
    label    = lab_frac,
    fontface = "italic",
    size     = 4
  ) +
  ggplot2::coord_cartesian(ylim = c(0, 1.08), clip = "off")

print(fig4_p_dom_comp_annot)

# ===== FIGURE 4D =====
# Baseline β-diversity (Bray–Curtis PCoA + PERMANOVA)

meta_beta <- metadata %>%
  dplyr::mutate(
    Sample.ID     = as.character(Sample.ID),
    LaterDiarrhea = factor(LaterDiarrhea, levels = c("No","Yes"))
  ) %>%
  dplyr::filter(!is.na(LaterDiarrhea)) %>%
  dplyr::distinct(Sample.ID, .keep_all = TRUE)

keep_ids_beta <- intersect(sample_names(phy), meta_beta$Sample.ID)
stopifnot(length(keep_ids_beta) >= 4)

phy_beta <- prune_samples(keep_ids_beta, phy)

bray_beta <- phyloseq::distance(phy_beta, method = "bray")
ord_beta  <- ordinate(phy_beta, method = "PCoA", distance = bray_beta)

meta_perm <- meta_beta %>%
  dplyr::filter(Sample.ID %in% rownames(ord_beta$vectors)) %>%
  dplyr::select(Sample.ID, LaterDiarrhea)
rownames(meta_perm) <- meta_perm$Sample.ID

adon_tab <- as.data.frame(vegan::adonis2(bray_beta ~ LaterDiarrhea, data = meta_perm))

get_perm_label <- function(tab, preferred = "LaterDiarrhea") {
  rn   <- rownames(tab)
  idx  <- which(rn == preferred)
  if (!length(idx)) idx <- which(!(rn %in% c("Residual","Total")))[1]
  fcol <- if ("F" %in% names(tab)) "F" else if ("F.Model" %in% names(tab)) "F.Model" else names(tab)[2]
  Fv   <- as.numeric(tab[idx, fcol])
  Pv   <- as.numeric(tab[idx, "Pr(>F)"])
  sprintf("PERMANOVA  F=%.2f, p %s",
          Fv,
          ifelse(is.finite(Pv) && Pv < .001, "<0.001", sprintf("= %.3g", Pv)))
}
ann_perm <- get_perm_label(adon_tab)
print(ann_perm)

pcv <- ord_beta$values$Relative_eig
x_lab <- if (!is.null(pcv) && length(pcv) >= 1 && is.finite(pcv[1])) sprintf("PCoA1 (%.1f%%)", 100*pcv[1]) else "PCoA1"
y_lab <- if (!is.null(pcv) && length(pcv) >= 2 && is.finite(pcv[2])) sprintf("PCoA2 (%.1f%%)", 100*pcv[2]) else "PCoA2"

cat_cols2 <- c("No" = "#1B9E77", "Yes" = "#D55E00")

df_ordination <- ord_beta$vectors %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Sample.ID") %>%
  dplyr::select(Sample.ID, Axis.1, Axis.2) %>%
  dplyr::left_join(meta_perm, by = "Sample.ID")

fig4_p_beta <- ggplot(df_ordination, aes(Axis.1, Axis.2)) +
  stat_ellipse(aes(group = LaterDiarrhea, fill = LaterDiarrhea),
               geom = "polygon", level = 0.95, alpha = 0.1,
               colour = NA, show.legend = FALSE) +
  stat_ellipse(aes(color = LaterDiarrhea),
               level = 0.95, linewidth = 1, show.legend = TRUE) +
  geom_point(aes(color = LaterDiarrhea), size = 2, alpha = 0.6) +
  scale_color_manual(values = cat_cols2, name = "") +
  scale_fill_manual(values = cat_cols2, guide = "none") +
  labs(x = x_lab, y = y_lab, title = "Beta-diversity (Bray-Curtis)") +
  scale_x_continuous(expand = expansion(mult = 0.12)) +
  scale_y_continuous(expand = expansion(mult = 0.12)) +
  coord_equal(expand = TRUE, clip = "off") +
  theme(plot.margin = margin(6, 10, 6, 6)) +
  theme_classic(base_size = 14) +
  theme(legend.position = "right",
        legend.direction = "vertical",
        panel.grid = element_blank()) +
  annotate("text", x = Inf, y = -Inf,
           hjust = 1.02, vjust = -0.6, size = 3.5,
           label = ann_perm)

print(fig4_p_beta)

# ===== FIGURE 4E =====
# Baseline PICRUSt2 bile acid enzymes

plot_marker_by_group <- function(marker_df, metadata,
                                 marker_label = "BSH",
                                 filename     = NULL) {
  stopifnot(all(c("Sample.ID","log_CPM") %in% names(marker_df)),
            all(c("Sample.ID","LaterDiarrhea") %in% names(metadata)))
  
  dat <- marker_df %>%
    dplyr::select(Sample.ID, log_CPM) %>%
    dplyr::inner_join(
      metadata %>%
        dplyr::transmute(Sample.ID = as.character(Sample.ID),
                         LaterDiarrhea = factor(LaterDiarrhea, levels = c("No","Yes"))),
      by = "Sample.ID"
    ) %>%
    dplyr::filter(!is.na(LaterDiarrhea), is.finite(log_CPM))
  
  p_val <- tryCatch(
    wilcox.test(log_CPM ~ LaterDiarrhea, data = dat)$p.value,
    error = function(e) NA_real_
  )
  p_lab <- paste0(
    "Wilcoxon p",
    ifelse(
      !is.finite(p_val), "NA",
      ifelse(p_val < .001, "<0.001", sprintf("=%.3f", p_val))
    )
  )
  
  p <- ggplot(dat, aes(x = LaterDiarrhea, y = log_CPM,
                       fill = LaterDiarrhea, colour = LaterDiarrhea)) +
    geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.30, linewidth = 0.7) +
    geom_jitter(width = 0.12, size = 1.6, alpha = 0.65) +
    scale_fill_manual(values = cat_cols, guide = "none") +
    scale_colour_manual(values = cat_cols, guide = "none") +
    labs(
      x     = NULL,
      y     = "log10(CPM)",
      title = sprintf("%s abundance by later diarrhea status", marker_label)
    ) +
    theme_classic(base_size = 13) +
    theme(axis.text.x = element_text(angle = 20, hjust = 1)) +
    annotate("text", x = Inf, y = Inf, hjust = 1.02, vjust = 1.2,
             label = p_lab, fontface = "italic", size = 4) +
    coord_cartesian(clip = "off")

  p
}

p_bsh_baseline  <- plot_marker_by_group(
  bsh_df,  metadata,
  marker_label = "BSH",
  filename     = "Fig_BSH_logCPM_vs_LaterDiarrhea.svg"
)
p_hsdh_baseline <- plot_marker_by_group(
  hsdh_df, metadata,
  marker_label = "7α-HSDH",
  filename     = "Fig_7aHSDH_logCPM_vs_LaterDiarrhea.svg"
)
p_bai_baseline  <- plot_marker_by_group(
  bai_df,  metadata,
  marker_label = "bai operon",
  filename     = "Fig_baiOperon_logCPM_vs_LaterDiarrhea.svg"
)

print(p_bsh_baseline)
print(p_hsdh_baseline)
print(p_bai_baseline)

# ===== FIGURE 4F =====
# Baseline taxa relative abundances

phy_base      <- prune_samples(
  metadata$Sample.ID[metadata$SampleCategory == "Baseline"], phy)
phy_genus     <- tax_glom(phy_base, taxrank = "genus", NArm = FALSE)
phy_genus_rel <- transform_sample_counts(phy_genus, function(x) x / sum(x))
rel_mat <- as(otu_table(phy_genus_rel), "matrix")
if (taxa_are_rows(phy_genus_rel)) rel_mat <- t(rel_mat)
tax_df <- as.data.frame(tax_table(phy_genus)) |> tibble::rownames_to_column("otu_id")
colnames(rel_mat) <- tax_df$genus[match(colnames(rel_mat), tax_df$otu_id)]

dat_genus_bl <- data.frame(
  Sample.ID        = rownames(rel_mat),
  Romboutsia       = rel_mat[, "Romboutsia"],
  Faecalibacterium = rel_mat[, "Faecalibacterium"]
) |>
  dplyr::left_join(metadata |> dplyr::select(Sample.ID, LaterDiarrhea), by = "Sample.ID") |>
  dplyr::filter(!is.na(LaterDiarrhea)) |>
  dplyr::mutate(LaterDiarrhea = factor(LaterDiarrhea, levels = c("No", "Yes")))

plot_genus_fig4 <- function(df, y, title) {
  p_val <- tryCatch(
    wilcox.test(reformulate("LaterDiarrhea", y), data = df, exact = FALSE)$p.value,
    error = function(e) NA_real_)
  p <- ggplot(df, aes(x = LaterDiarrhea, y = .data[[y]],
                      fill = LaterDiarrhea, colour = LaterDiarrhea)) +
    geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.30, linewidth = 0.7) +
    geom_jitter(width = 0.12, size = 1.6, alpha = 0.65) +
    scale_fill_manual(values = cat_cols, guide = "none") +
    scale_colour_manual(values = cat_cols, guide = "none") +
    labs(x = NULL, y = "Relative abundance", title = title) +
    theme_classic(base_size = 13) +
    theme(axis.text.x = element_text(angle = 20, hjust = 1),
          plot.title  = element_text(face = "italic")) +
    annotate("text", x = Inf, y = Inf, hjust = 1.02, vjust = 1.2,
             label = wilcox_label(p_val), fontface = "italic", size = 4) +
    coord_cartesian(clip = "off")
  print(p); invisible(p)
}

p_romboutsia <- plot_genus_fig4(dat_genus_bl, "Romboutsia", "Romboutsia")
p_faecalibacterium <- plot_genus_fig4(dat_genus_bl, "Faecalibacterium", "Faecalibacterium")

# ===== FIGURE 4G + SUPPLEMENTAL FIGURE 4A =====
# Baseline bile acids
#   CA / CDCA / DCA / LCA and the two ratios = Figure 4G
#   GCA / TCA                                = Supplemental Figure 4A

acid_cols <- c(
  log10_CA   = "Cholate (CA)",              
  log10_CDCA = "Chenodeoxycholate (CDCA)", 
  log10_DCA  = "Deoxycholate (DCA)",        
  log10_LCA  = "Lithocholate (LCA)",    
  log10_GCA  = "Glycocholate (GCA)",     
  log10_TCA  = "Taurocholate (TCA)"    
)

ratio_cols <- c(
  CA_DCA_log10   = "CA/DCA ratio",    
  CDCA_LCA_log10 = "CDCA/LCA ratio" 
)

wilcox_label <- function(p) {
  paste0("Wilcoxon p",
         ifelse(!is.finite(p), "NA",
                ifelse(p < .001, "<0.001", sprintf("=%.3f", p))))
}

plot_box_jitter <- function(df, y, ylab, title, file_stub) {
  stopifnot(all(c("LaterDiarrhea", y) %in% names(df)))
  ok <- df %>% dplyr::count(LaterDiarrhea) %>% dplyr::filter(n >= 2) %>% nrow() == 2
  p_val <- tryCatch(
    if (ok) wilcox.test(reformulate("LaterDiarrhea", y), data = df, exact = FALSE)$p.value else NA_real_,
    error = function(e) NA_real_
  )
  p_lab <- wilcox_label(p_val)
  
  fig4_p_boxplots <- ggplot(df, aes(x = LaterDiarrhea, y = .data[[y]],
                      fill = LaterDiarrhea, colour = LaterDiarrhea)) +
    geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.30, linewidth = 0.7) +
    geom_jitter(width = 0.12, size = 1.6, alpha = 0.65) +
    scale_fill_manual(values = cat_cols, guide = "none") +
    scale_colour_manual(values = cat_cols, guide = "none") +
    labs(x = NULL, y = ylab, title = title) +
    theme_classic(base_size = 13) +
    theme(axis.text.x = element_text(angle = 20, hjust = 1)) +
    annotate("text", x = Inf, y = Inf, hjust = 1.02, vjust = 1.2,
             label = p_lab, fontface = "italic", size = 4) +
    coord_cartesian(clip = "off")
  
  print(fig4_p_boxplots)
  invisible(fig4_p_boxplots)
}

need_cols <- c(
  "Sample.ID","SampleCat","LaterDiarrhea",
  names(acid_cols), names(ratio_cols)
)
missing <- setdiff(need_cols, names(df_ratio))
if (length(missing)) {
  stop("df_ratio is missing: ", paste(missing, collapse = ", "))
}

dat_base <- df_ratio %>%
  dplyr::transmute(
    Sample.ID,
    SampleCat      = as.character(SampleCat),
    LaterDiarrhea  = factor(LaterDiarrhea, levels = c("No","Yes")),
    dplyr::across(all_of(c(names(acid_cols), names(ratio_cols))), ~ suppressWarnings(as.numeric(.x)))
  ) %>%
  dplyr::filter(SampleCat == "Baseline", !is.na(LaterDiarrhea)) %>%
  dplyr::distinct(Sample.ID, .keep_all = TRUE)

for (nm in names(acid_cols)) {
  pretty <- acid_cols[[nm]]
  df <- dat_base %>% dplyr::filter(is.finite(.data[[nm]]))
  if (nrow(df) >= 3) {
    plot_box_jitter(
      df, y = nm,
      ylab = "log10(AUC)",
      title = pretty,
      file_stub = paste0("BL_", gsub("^log10_", "", nm), "_logAUC_vs_LaterDiarrhea")
    )
  } else {
    message("Skipping ", pretty, ": too few baseline samples with finite values.")
  }
}

for (nm in names(ratio_cols)) {
  pretty <- ratio_cols[[nm]]
  df <- dat_base %>% dplyr::filter(is.finite(.data[[nm]]))
  if (nrow(df) >= 3) {
    plot_box_jitter(
      df, y = nm,
      ylab = paste0(
        "log10(",
        gsub("_", "/", gsub("_log10$", "", nm)),
        ")"
      ),
      title = pretty,
      file_stub = paste0("BL_", gsub("_log10$", "", nm), "_logRatio_vs_LaterDiarrhea")
    )
  } else {
    message("Skipping ", pretty, ": too few baseline samples with finite values.")
  }
}


# ===== SUPPLEMENTAL FIGURE 4B =====
# Baseline BA ratios (CA/DCA, CDCA/LCA) vs α-diversity

dat_bl_cadca <- df_ratio %>%
  dplyr::select(Sample.ID, SampleCat, CA_DCA_log10) %>%
  dplyr::left_join(
    metadata %>%
      dplyr::select(Sample.ID, SampleCategory, invsimpson, LaterDiarrhea),
    by = "Sample.ID"
  ) %>%
  dplyr::mutate(
    Group         = dplyr::coalesce(as.character(SampleCat), as.character(SampleCategory)),
    LaterDiarrhea = factor(LaterDiarrhea, levels = c("No", "Yes"))
  ) %>%
  dplyr::filter(Group == "Baseline",
                is.finite(CA_DCA_log10),
                is.finite(invsimpson),
                !is.na(LaterDiarrhea)) %>%
  dplyr::mutate(x_alpha = log10(invsimpson), y_ratio = CA_DCA_log10)

if (nrow(dat_bl_cadca) < 3) stop("Too few baseline samples for CA/DCA vs alpha correlation.")

ct_cadca    <- stats::cor.test(dat_bl_cadca$x_alpha, dat_bl_cadca$y_ratio, method = "pearson")
p_txt_cadca <- ifelse(ct_cadca$p.value < .001, "<0.001", sprintf("=%.3f", ct_cadca$p.value))
ann_cadca   <- sprintf("R² = %.2f, p%s", unname(ct_cadca$estimate)^2, p_txt_cadca)

fig4_p_ca_dca_alpha_bl <- ggplot(dat_bl_cadca, aes(x = x_alpha, y = y_ratio)) +
  geom_smooth(method = "lm", se = TRUE,
              color = "black", fill = "grey70", alpha = 0.18, linewidth = 1) +
  geom_point(aes(color = LaterDiarrhea), size = 2.2, alpha = 0.85) +
  scale_color_manual(values = cat_cols, name = "Diarrhea later") +
  labs(x = "log10(Inverse Simpson)", y = "log10(ratio)", title = "CA/DCA vs. alpha-diversity") +
  theme_classic(base_size = 13) +
  theme(legend.position = "right", panel.grid = element_blank()) +
  annotate("text", x = Inf, y = Inf, hjust = 1.02, vjust = 1.2,
           size = 3.6, fontface = "italic", label = ann_cadca) +
  coord_cartesian(clip = "off")

print(fig4_p_ca_dca_alpha_bl)

dat_bl_cdcalca <- df_ratio %>%
  dplyr::select(Sample.ID, SampleCat, CDCA_LCA_log10) %>%
  dplyr::left_join(
    metadata %>%
      dplyr::select(Sample.ID, SampleCategory, invsimpson, LaterDiarrhea),
    by = "Sample.ID"
  ) %>%
  dplyr::mutate(
    Group         = dplyr::coalesce(as.character(SampleCat), as.character(SampleCategory)),
    LaterDiarrhea = factor(LaterDiarrhea, levels = c("No", "Yes"))
  ) %>%
  dplyr::filter(Group == "Baseline",
                is.finite(CDCA_LCA_log10),
                is.finite(invsimpson),
                !is.na(LaterDiarrhea)) %>%
  dplyr::mutate(x_alpha = log10(invsimpson), y_ratio = CDCA_LCA_log10)

if (nrow(dat_bl_cdcalca) < 3) stop("Too few baseline samples for CDCA/LCA vs alpha correlation.")

ct_cdcalca    <- stats::cor.test(dat_bl_cdcalca$x_alpha, dat_bl_cdcalca$y_ratio, method = "pearson")
p_txt_cdcalca <- ifelse(ct_cdcalca$p.value < .001, "<0.001", sprintf("=%.3f", ct_cdcalca$p.value))
ann_cdcalca   <- sprintf("R² = %.2f, p%s", unname(ct_cdcalca$estimate)^2, p_txt_cdcalca)

fig4_p_cdca_lca_alpha_bl <- ggplot(dat_bl_cdcalca, aes(x = x_alpha, y = y_ratio)) +
  geom_smooth(method = "lm", se = TRUE,
              color = "black", fill = "grey70", alpha = 0.18, linewidth = 1) +
  geom_point(aes(color = LaterDiarrhea), size = 2.2, alpha = 0.85) +
  scale_color_manual(values = cat_cols, name = "Diarrhea later") +
  labs(x = "log10(Inverse Simpson)", y = "log10(ratio)", title = "CDCA/LCA vs. alpha-diversity") +
  theme_classic(base_size = 13) +
  theme(legend.position = "right", panel.grid = element_blank()) +
  annotate("text", x = Inf, y = Inf, hjust = 1.02, vjust = 1.2,
           size = 3.6, fontface = "italic", label = ann_cdcalca) +
  coord_cartesian(clip = "off")

print(fig4_p_cdca_lca_alpha_bl)
