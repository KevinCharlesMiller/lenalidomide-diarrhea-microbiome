#####
# Figure 5, Supplemental Figure 5
#####

if (!exists("df_taxa_genus")) source("R/00_setup.R")

dens_labels <- c(
  DT_FIBER        = "Total fiber",
  DT_FIBER_INSOL = "Insoluble fiber",
  DT_FIBER_SOL   = "Soluble fiber",
  STARCH         = "Starch",
  DT_SUG_T       = "Total sugars",
  DT_SFAT        = "Saturated fat",
  DT_MFAT        = "Monounsaturated fat",
  DT_PFAT        = "Polyunsaturated fat",
  DT_TFAT        = "Total fat",
  DT_CARB        = "Total carbohydrate",
  DT_PROT        = "Total protein"
)

macro_codes   <- c("DT_TFAT", "DT_CARB", "DT_PROT")
subnutr_codes <- c(
  "DT_FIBER", "DT_FIBER_SOL", "DT_FIBER_INSOL", "STARCH", "DT_SUG_T",
  "DT_SFAT", "DT_MFAT", "DT_PFAT"
)

diet_base <- readr::read_csv(
  "data/ffq.csv",
  show_col_types = FALSE
) %>%
  mutate(
    Sample.ID = as.character(Sample.ID),
    SubjectID       = as.character(SubjectID),
    SubjectID       = str_replace(SubjectID, "^0+(?=\\d+$)", "")
  )

dens_long_all <- diet_base %>%
  select(SubjectID, Sample.ID, LaterDiarrhea, ends_with("_dens")) %>%
  pivot_longer(
    cols      = ends_with("_dens"),
    names_to  = "nut_code",
    values_to = "density"
  ) %>%
  mutate(
    base_code = sub("_dens$", "", nut_code),
    nutrient  = factor(dens_labels[base_code], levels = unname(dens_labels))
  ) %>%
  drop_na(density)

# ===== FIGURE 5A =====
# Macronutrient densities

lab_pq <- function(p, q = NULL, digits = 3) {
  if (is.null(q)) {
    paste0("p", fmt_p(p))
  } else {
    paste0("p", fmt_p(p), ", q", fmt_p(q))
  }
}
fmt_p <- function(p) {
  ifelse(
    is.na(p), "=NA",
    ifelse(p < 0.001, "<0.001", paste0("=", sprintf("%.3f", p)))
  )
}

lab_pq <- function(p, q) {
  paste0("p", fmt_p(p), ", q", fmt_p(q))
}

dens_long_A <- dens_long_all %>%
  filter(base_code %in% macro_codes) %>%
  mutate(LaterDiarrhea = factor(LaterDiarrhea, levels = c("No", "Yes")))

wilcox_tbl_A <- dens_long_A %>%
  group_by(nutrient) %>%
  summarise(
    p_raw = tryCatch(
      stats::wilcox.test(density ~ LaterDiarrhea, exact = FALSE)$p.value,
      error = function(e) NA_real_
    ),
    .groups = "drop"
  ) %>%
  mutate(
    q_bh  = p.adjust(p_raw, method = "BH"),
    p_lab = lab_pq(p_raw, q_bh)
  )

print(wilcox_tbl_A %>% arrange(q_bh))

ann_df_A <- wilcox_tbl_A %>% select(nutrient, p_lab)

p_ffq_macros <- ggplot(
  dens_long_A,
  aes(
    x = LaterDiarrhea, y = density,
    fill = LaterDiarrhea, colour = LaterDiarrhea
  )
) +
  geom_boxplot(
    width = 0.6, outlier.shape = NA,
    alpha = 0.30, linewidth = 0.7
  ) +
  geom_jitter(
    width = 0.12, size = 1.6, alpha = 0.65
  ) +
  scale_fill_manual(values = cat_cols, guide = "none") +
  scale_colour_manual(values = cat_cols, guide = "none") +
  scale_x_discrete(drop = FALSE, limits = c("No", "Yes")) +
  scale_y_continuous(expand = expansion(mult = c(0.04, 0.18))) +
  facet_wrap(~ nutrient, scales = "free_y") +
  geom_text(
    data = ann_df_A,
    mapping = aes(label = p_lab),
    x = Inf, y = Inf,
    hjust = 1.02, vjust = 1.4,
    fontface = "italic", size = 5.0,
    inherit.aes = FALSE
  ) +
  labs(x = NULL, y = "Nutrient density (g/1000 kcal)") +
  coord_cartesian(clip = "off") +
  theme_classic(base_size = 13) +
  theme(
    axis.text.x      = element_text(angle = 20, hjust = 1),
    strip.background = element_blank(),
    strip.text       = element_text(size = 15, face = "bold"),
    plot.title       = element_blank()
  )

print(p_ffq_macros)

# ===== FIGURE 5B =====
# Sub-nutrient densities

dens_long_B <- dens_long_all %>%
  filter(base_code %in% subnutr_codes)

flav_long_B <- diet_base %>%
  dplyr::filter(!is.na(EnergyAdj_Total_flavonoids_milligrams)) %>%
  dplyr::transmute(
    SubjectID, Sample.ID, LaterDiarrhea,
    nut_code  = "FLAV", base_code = "FLAV",
    density   = EnergyAdj_Total_flavonoids_milligrams,
    nutrient  = "Total flavonoids"
  )
dens_long_B <- dplyr::bind_rows(dens_long_B, flav_long_B)

facet_order <- c(
  "Total fiber","Soluble fiber", "Insoluble fiber", "Starch", "Total sugars",
  "Total flavonoids",
  "Saturated fat", "Monounsaturated fat", "Polyunsaturated fat"
)
present   <- unique(as.character(dens_long_B$nutrient))
facet_use <- intersect(facet_order, present)

dens_long_B <- dens_long_B %>%
  mutate(
    nutrient      = factor(as.character(nutrient), levels = facet_use),
    LaterDiarrhea = factor(LaterDiarrhea, levels = c("No", "Yes"))
  )

wilcox_tbl_B <- dens_long_B %>%
  group_by(nutrient) %>%
  summarise(
    p_raw = tryCatch(
      stats::wilcox.test(density ~ LaterDiarrhea, exact = FALSE)$p.value,
      error = function(e) NA_real_
    ),
    .groups = "drop"
  ) %>%
  mutate(
    q_bh = p.adjust(p_raw, method = "BH"),
    lab  = paste0("p", fmt_p(p_raw), ", q", fmt_p(q_bh))
  )

ann_df_B <- wilcox_tbl_B %>% select(nutrient, lab)

p_ffq_subnutr_long <- ggplot(
  dens_long_B,
  aes(
    x = LaterDiarrhea, y = density,
    fill = LaterDiarrhea, colour = LaterDiarrhea
  )
) +
  geom_boxplot(
    width = 0.6, outlier.shape = NA,
    alpha = 0.30, linewidth = 0.7
  ) +
  geom_jitter(
    width = 0.12, size = 1.6, alpha = 0.65
  ) +
  scale_fill_manual(values = cat_cols, guide = "none") +
  scale_colour_manual(values = cat_cols, guide = "none") +
  scale_x_discrete(drop = FALSE, limits = c("No", "Yes")) +
  scale_y_continuous(expand = expansion(mult = c(0.04, 0.18))) +
  facet_wrap(~ nutrient, scales = "free_y") +  
  geom_text(
    data = ann_df_B,
    mapping = aes(label = lab),
    x = Inf, y = Inf,
    hjust = 1.02, vjust = 1.4,
    fontface = "italic", size = 5.0,
    inherit.aes = FALSE
  ) +
  labs(x = NULL, y = "Nutrient density (g/1000 kcal)") +
  coord_cartesian(clip = "off") +
  theme_classic(base_size = 13) +
  theme(
    axis.text.x      = element_text(angle = 20, hjust = 1),
    strip.background = element_blank(),
    strip.text       = element_text(size = 15, face = "bold"),
    plot.title       = element_blank()
  )

print(p_ffq_subnutr_long)

# ===== SUPPLEMENTAL FIGURE 5A/5B =====
# Diet vs baseline α-diversity
fmt_p_simple <- function(p) ifelse(is.na(p), "NA", ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))
macro_alpha_labels <- c(DT_TFAT = "Total fat", DT_CARB = "Total carbohydrate", DT_PROT = "Total protein")

dat_macro <- diet_base %>%
  dplyr::left_join(metadata %>% dplyr::select(Sample.ID, invsimpson), by = "Sample.ID") %>%
  dplyr::filter(is.finite(invsimpson), invsimpson > 0) %>%
  dplyr::mutate(x_alpha = log10(invsimpson),
                LaterDiarrhea = factor(LaterDiarrhea, levels = c("No", "Yes"))) %>%
  dplyr::select(SubjectID, Sample.ID, LaterDiarrhea, x_alpha, ends_with("_dens")) %>%
  tidyr::pivot_longer(cols = ends_with("_dens"), names_to = "nut_code", values_to = "density") %>%
  dplyr::mutate(base_code = sub("_dens$", "", nut_code)) %>%
  dplyr::filter(base_code %in% names(macro_alpha_labels)) %>%
  dplyr::mutate(nutrient = factor(macro_alpha_labels[base_code], levels = unname(macro_alpha_labels))) %>%
  tidyr::drop_na(density)

for (nut in levels(dat_macro$nutrient)) {
  sub <- dat_macro %>% dplyr::filter(as.character(nutrient) == nut)
  if (nrow(sub) < 3) next
  ct <- tryCatch(cor.test(sub$x_alpha, sub$density, method = "pearson"), error = function(e) NULL)
  r  <- if (!is.null(ct)) unname(ct$estimate) else NA_real_
  p  <- if (!is.null(ct)) ct$p.value else NA_real_
  p_gp <- ggplot(sub, aes(x = x_alpha, y = density)) +
    geom_smooth(method = "lm", se = TRUE, color = "black", fill = "grey70", alpha = 0.18, linewidth = 1) +
    geom_point(aes(color = LaterDiarrhea), size = 2.2, alpha = 0.85) +
    scale_color_manual(values = cat_cols, name = "Diarrhea later") +
    labs(x = "log10(Inverse Simpson)", y = "Nutrient density (g/1000 kcal)", title = as.character(nut)) +
    theme_classic(base_size = 13) +
    theme(legend.position = "right", panel.grid = element_blank(),
          axis.line = element_line(color = "black")) +
    annotate("text", x = Inf, y = Inf, hjust = 1.05, vjust = 1.5,
             label = sprintf("R² = %.2f, p=%s", r^2, fmt_p_simple(p)), size = 3.6)
  print(p_gp)
}

fiber_labels <- c(
  DT_FIBER       = "Total fiber",
  DT_FIBER_SOL   = "Soluble fiber",
  DT_FIBER_INSOL = "Insoluble fiber"
)
fiber_codes <- names(fiber_labels)

dat_fiber <- diet_base %>%
  left_join(metadata %>% select(Sample.ID, invsimpson), by = "Sample.ID") %>%
  filter(is.finite(invsimpson), invsimpson > 0) %>%
  mutate(
    x_alpha = log10(invsimpson),
    LaterDiarrhea = factor(LaterDiarrhea, levels = c("No", "Yes"))
  ) %>%
  select(SubjectID, Sample.ID, LaterDiarrhea, x_alpha, ends_with("_dens")) %>%
  pivot_longer(
    cols      = ends_with("_dens"),
    names_to  = "nut_code",
    values_to = "density"
  ) %>%
  mutate(base_code = sub("_dens$", "", nut_code)) %>%
  filter(base_code %in% fiber_codes) %>%
  mutate(
    nutrient = factor(fiber_labels[base_code], levels = unname(fiber_labels))
  ) %>%
  drop_na(density)

nut_list <- levels(dat_fiber$nutrient)
if (is.null(nut_list)) nut_list <- unique(as.character(dat_fiber$nutrient))

fmt_p_simple <- function(p) {
  ifelse(is.na(p), "NA", ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))
}

plot_list <- list()

for (nut in nut_list) {
  sub <- dat_fiber %>% filter(as.character(nutrient) == nut)
  if (nrow(sub) < 3) next

  ct <- tryCatch(cor.test(sub$x_alpha, sub$density, method = "pearson"),
                 error = function(e) NULL)
  r  <- if (!is.null(ct)) unname(ct$estimate) else NA_real_
  p  <- if (!is.null(ct)) ct$p.value else NA_real_

  ann <- sprintf("R\u00b2 = %.2f, p=%s", r^2, fmt_p_simple(p))

  x_min <- min(sub$x_alpha,  na.rm = TRUE); x_max <- max(sub$x_alpha,  na.rm = TRUE)
  y_min <- min(sub$density,  na.rm = TRUE); y_max <- max(sub$density,  na.rm = TRUE)
  xpos  <- x_min + 0.02 * (x_max - x_min)
  ypos  <- y_max - 0.05 * (y_max - y_min)

  p_gp <- ggplot(sub, aes(x = x_alpha, y = density)) +
    geom_smooth(
      method = "lm", se = TRUE,
      color = "black", fill = "grey70", alpha = 0.18, linewidth = 1
    ) +
    geom_point(aes(color = LaterDiarrhea), size = 2.2, alpha = 0.85) +
    scale_color_manual(values = cat_cols, name = "Diarrhea later") +
    labs(
      x = "log10(Inverse Simpson)",
      y = "Nutrient density (g/1000 kcal)",
      title = as.character(nut)
    ) +
    theme_classic(base_size = 13) +
    theme(
      legend.position = "right",
      panel.grid      = element_blank(),
      axis.line       = element_line(color = "black")
    ) +
    annotate("text", x = xpos, y = ypos, hjust = 0, vjust = 1,
             label = ann, size = 3.6)

  print(p_gp)
  plot_list[[as.character(nut)]] <- p_gp
}
  
# ===== SUPPLEMENTAL FIGURE 5C =====
# Diet PCoA overlays

macro_fiber_codes <- c("DT_TFAT", "DT_CARB", "DT_PROT", "DT_FIBER")
macro_fiber_dens_cols <- paste0(macro_fiber_codes, "_dens")

meta_baseline <- metadata %>%
  mutate(
    Sample.ID      = as.character(Sample.ID),
    SubjectID            = as.character(SubjectID),
    SampleCategory = factor(SampleCategory,
                            levels = c("Baseline", "OnLen_control", "OnLen_diarrhea"))
  ) %>%
  filter(SampleCategory == "Baseline") %>%
  distinct(Sample.ID, .keep_all = TRUE)

diet_bl_macro_fiber <- meta_baseline %>%
  left_join(
    diet_base %>%
      select(Sample.ID, SubjectID, LaterDiarrhea, all_of(macro_fiber_dens_cols)),
    by = c("Sample.ID", "SubjectID", "LaterDiarrhea")
  )

ids_keep <- intersect(diet_bl_macro_fiber$Sample.ID, rownames(df_taxa_genus))
X_genus_b <- df_taxa_genus[ids_keep, , drop = FALSE]
stopifnot(nrow(X_genus_b) >= 3L)

diet_bl_macro_fiber <- diet_bl_macro_fiber %>%
  filter(Sample.ID %in% ids_keep)

row_sums <- rowSums(X_genus_b)
row_sums[row_sums == 0] <- 1
X_rel     <- sweep(as.matrix(X_genus_b), 1, row_sums, "/")
bray_base <- vegan::vegdist(X_rel, method = "bray")

diet_df_macro_fiber <- diet_bl_macro_fiber %>%
  select(Sample.ID, SubjectID, LaterDiarrhea, all_of(macro_fiber_dens_cols)) %>%
  column_to_rownames("Sample.ID")
diet_df_macro_fiber <- diet_df_macro_fiber[attr(bray_base, "Labels"), , drop = FALSE]

dens_labels_mac_fiber <- c(
  DT_CARB  = "Total carbohydrate (g/1000 kcal)",
  DT_TFAT  = "Total fat (g/1000 kcal)",
  DT_PROT  = "Total protein (g/1000 kcal)",
  DT_FIBER = "Total fiber (g/1000 kcal)"
)

perm_one <- function(var) {
  v     <- diet_df_macro_fiber[[var]]
  keep  <- is.finite(v)
  n_use <- sum(keep)
  if (n_use < 3L) {
    return(tibble::tibble(
      variable  = var,
      F         = NA_real_,
      R2        = NA_real_,
      p_value   = NA_real_,
      n_samples = n_use
    ))
  }
  Dmat  <- as.matrix(bray_base)
  D_sub <- stats::as.dist(Dmat[keep, keep, drop = FALSE])
  df_sub <- data.frame(x = v[keep])
  a2 <- vegan::adonis2(D_sub ~ x, data = df_sub,
                       permutations = 999, by = "margin")
  tibble::tibble(
    variable  = var,
    F         = unclass(a2$F)[1],
    R2        = unclass(a2$R2)[1],
    p_value   = unclass(a2$`Pr(>F)`)[1],
    n_samples = n_use
  )
}

res_tbl <- purrr::map_dfr(macro_fiber_dens_cols, perm_one) %>%
  mutate(
    q_value = p.adjust(p_value, method = "BH"),
    nice    = dplyr::recode(
      variable,
      !!!setNames(unname(dens_labels_mac_fiber),
                  paste0(names(dens_labels_mac_fiber), "_dens")),
      .default = variable
    )
  ) %>%
  select(variable, nice, n_samples, F, R2, p_value, q_value) %>%
  arrange(q_value)

print(res_tbl, n = Inf)

ids_keep_phy <- intersect(ids_keep, sample_names(phy))
phy_sub      <- prune_samples(ids_keep_phy, phy)

bray_sub <- phyloseq::distance(phy_sub, method = "bray")
ord_sub  <- ordinate(phy_sub, method = "PCoA", distance = bray_sub)

diet_bl_macro_fiber_pcoa <- diet_bl_macro_fiber %>%
  filter(Sample.ID %in% ids_keep_phy)

plot_overlay_with_permanova <- function(
    marker_col, 
    marker_label,
    cat_cols    = cat_cols,
    point_size  = 4.5,
    point_alpha = 0.85,
    zoom_mult   = 0.14,
    nperm       = 999
) {
  df_ord <- ord_sub$vectors %>%
    as.data.frame() %>%
    tibble::rownames_to_column("Sample.ID") %>%
    select(Sample.ID, Axis.1, Axis.2) %>%
    left_join(diet_bl_macro_fiber_pcoa, by = "Sample.ID") %>%
    mutate(marker_val = .data[[marker_col]]) %>%
    filter(is.finite(marker_val), !is.na(LaterDiarrhea)) %>%
    droplevels()
  
  stopifnot(nrow(df_ord) >= 6L)
  
  Dmat  <- as.matrix(bray_sub)
  D_sub <- stats::as.dist(Dmat[df_ord$Sample.ID, df_ord$Sample.ID, drop = FALSE])
  
  ctrl <- if (any(table(df_ord$SubjectID) > 1)) {
    h <- permute::how(nperm = nperm)
    permute::setBlocks(h) <- df_ord$SubjectID
    h
  } else nperm
  
  fit <- vegan::adonis2(
    D_sub ~ marker_val + LaterDiarrhea,
    data        = df_ord,
    permutations = ctrl,
    by          = "margin"
  )
  
  tab   <- as.data.frame(fit)
  
  f_col <- dplyr::case_when(
    "F" %in% names(tab)       ~ "F",
    "F.Model" %in% names(tab) ~ "F.Model",
    TRUE                      ~ names(tab)[2]
  )
  r2_col <- dplyr::case_when(
    "R2" %in% names(tab) ~ "R2",
    TRUE                 ~ names(tab)[grep("R2", names(tab), fixed = TRUE)[1]]
  )
  
  F_mrk  <- tab["marker_val", f_col]
  R2_mrk <- tab["marker_val", r2_col]
  P_mrk  <- tab["marker_val", "Pr(>F)"]
  
  ann <- sprintf(
    "%s PERMANOVA  F=%.2f, p %s",
    marker_label, 
    F_mrk, 
    fmt_p(P_mrk)
  )
  
  pcv   <- ord_sub$values$Relative_eig
  x_lab <- if (is.finite(pcv[1])) sprintf("PCoA 1 (%.1f%%)", 100 * pcv[1]) else "PCoA 1"
  y_lab <- if (is.finite(pcv[2])) sprintf("PCoA 2 (%.1f%%)", 100 * pcv[2]) else "PCoA 2"
  
  p <- ggplot(df_ord, aes(Axis.1, Axis.2)) +
    geom_point(
      aes(color = marker_val, shape = LaterDiarrhea),
      size  = point_size,
      alpha = point_alpha,
      na.rm = TRUE
    ) +
    scale_shape_manual(
      values = c("No" = 16, "Yes" = 17), 
      name = NULL
    ) +
    scale_color_viridis_c(
      option    = "inferno",   
      direction = -1,          
      name      = "Nutrient density\n(g/1000 kcal)",
      na.value  = "grey85"
    ) +
    labs(
      x = x_lab,
      y = y_lab,
      title = paste0("Bray–Curtis PCoA: ", marker_label)
    ) +
    scale_x_continuous(expand = expansion(mult = zoom_mult)) +
    scale_y_continuous(expand = expansion(mult = zoom_mult)) +
    coord_equal(expand = TRUE, clip = "off") +
    theme_classic(base_size = 13) +   
    theme(
      legend.position  = "right",
      panel.grid       = element_blank(),
      plot.margin      = margin(6, 10, 6, 6)
    ) +
    annotate(
      "text",
      x = Inf, y = -Inf,
      hjust = 1.02, vjust = -0.6,
      label = ann, size = 3.6   
    )
  
  list(
    plot = p, 
    permanova = fit, 
    data = df_ord,
    stats = tibble::tibble(
      marker = marker_label,
      n      = nrow(df_ord),
      F      = F_mrk,
      R2     = R2_mrk,
      p      = P_mrk
    )
  )
}

res_carb <- plot_overlay_with_permanova(
  marker_col   = "DT_CARB_dens",
  marker_label = "Total carbohydrate",
  cat_cols     = cat_cols
)
print(res_carb$plot)
print(res_carb$permanova)

res_fat <- plot_overlay_with_permanova(
  marker_col   = "DT_TFAT_dens",
  marker_label = "Total fat",
  cat_cols     = cat_cols
)
print(res_fat$plot)
print(res_fat$permanova)

res_prot <- plot_overlay_with_permanova(
  marker_col   = "DT_PROT_dens",
  marker_label = "Total protein",
  cat_cols     = cat_cols
)
print(res_prot$plot)
print(res_prot$permanova)

summ_tbl <- dplyr::bind_rows(
  res_carb$stats,
  res_fat$stats,
  res_prot$stats,

) %>%
  dplyr::mutate(
    p_fmt = dplyr::if_else(
      !is.finite(p),
      NA_character_,
      dplyr::if_else(p < .001, "<0.001", sprintf("%.3f", p))
    )
  ) %>%
  dplyr::select(marker, n, F, R2, p = p_fmt)

print(summ_tbl)

# ===== SUPPLEMENTAL FIGURE 5D =====
# Diet vs baseline bile-acid ratios (CA/DCA, CDCA/LCA)

dens_labels_ba <- c(
  DT_CARB = "Total carbohydrate (g/1000 kcal)",
  DT_TFAT = "Total fat (g/1000 kcal)",
  DT_PROT = "Total protein (g/1000 kcal)"
)
dens_vars   <- names(dens_labels_ba)
nut_cols    <- paste0(dens_vars, "_dens")
nut_labels  <- unname(dens_labels_ba)

ffq_dens <- diet_base %>%
  mutate(SubjectID = as.character(SubjectID)) %>%
  distinct(SubjectID, .keep_all = TRUE) %>%
  select(SubjectID, all_of(nut_cols))

meta_bl <- metadata %>%
  dplyr::mutate(
    SubjectID            = as.character(SubjectID),
    SampleCategory = factor(SampleCategory,
                            levels = c("Baseline","OnLen_control","OnLen_diarrhea")),
    LaterDiarrhea  = factor(LaterDiarrhea, levels = c("No","Yes"))
  ) %>%
  dplyr::filter(SampleCategory == "Baseline") %>%
  dplyr::arrange(SubjectID) %>%
  dplyr::distinct(SubjectID, .keep_all = TRUE) %>%
  dplyr::select(Sample.ID, SubjectID, LaterDiarrhea)

ratio_bl <- df_ratio %>%
  dplyr::mutate(
    Sample.ID = as.character(Sample.ID),
    SubjectID       = as.character(SubjectID)
  ) %>%
  dplyr::filter(SampleCat %in% c("Baseline","baseline")) %>%
  dplyr::select(
    Sample.ID, SubjectID, SampleCat,
    CA_DCA_log10, CDCA_LCA_log10,
    dplyr::any_of(c("CA","DCA","CDCA","LCA"))
  )

need_ratio <- c("CA_DCA_log10","CDCA_LCA_log10")
stopifnot(all(need_ratio %in% names(ratio_bl)))

dat_ba <- ratio_bl %>%
  dplyr::inner_join(meta_bl,  by = c("Sample.ID","SubjectID")) %>%
  dplyr::inner_join(ffq_dens, by = "SubjectID") %>%
  dplyr::filter(
    dplyr::if_all(dplyr::all_of(need_ratio), is.finite),
    dplyr::if_any(dplyr::all_of(nut_cols),   is.finite)
  )

make_ratio_plot <- function(data, ratio_col, nutrient_col, nutrient_label, ratio_title) {
  df <- data %>%
    dplyr::select(
      LaterDiarrhea,
      x = dplyr::all_of(ratio_col),
      y = dplyr::all_of(nutrient_col)
    ) %>%
    dplyr::filter(is.finite(x), is.finite(y))
  
  if (nrow(df) < 5L) return(NULL)
  
  ct <- suppressWarnings(stats::cor.test(df$x, df$y, method = "pearson"))
  r  <- unname(ct$estimate)
  p  <- ct$p.value
  ann <- sprintf("R\u00b2 = %.2f, p%s", r^2, fmt_p(p))
  
  ggplot(df, aes(x = x, y = y)) +
    geom_smooth(
      method = "lm", se = TRUE,
      color  = "black", fill = "grey70",
      alpha  = 0.18, linewidth = 1
    ) +
    geom_point(aes(color = LaterDiarrhea), size = 2.2, alpha = 0.85) +
    scale_color_manual(values = cat_cols, name = "Diarrhea later") +
    labs(
      x     = "log10(ratio)",
      y     = nutrient_label,
      title = ratio_title
    ) +
    theme_classic(base_size = 13) +
    theme(legend.position = "right", panel.grid = element_blank()) +
    coord_cartesian(clip = "off") +
    annotate(
      "text",
      x = Inf, y = Inf,
      hjust = 1.02, vjust = 1.2,
      label = ann, size = 3.9, fontface = "italic"
    )
}

pairs <- expand.grid(
  ratio = c("CA_DCA_log10","CDCA_LCA_log10"),
  nut_i = seq_along(nut_cols),
  stringsAsFactors = FALSE
)

plots_ba <- vector("list", nrow(pairs))
names(plots_ba) <- apply(
  pairs, 1,
  function(r) paste0(r["ratio"], "_vs_", nut_cols[as.integer(r["nut_i"])])
)

for (i in seq_len(nrow(pairs))) {
  ratio_col  <- pairs$ratio[i]
  nut_col    <- nut_cols[pairs$nut_i[i]]
  nut_lab    <- nut_labels[pairs$nut_i[i]]
  ratio_name <- if (ratio_col == "CA_DCA_log10") "CA/DCA" else "CDCA/LCA"
  ratio_title <- sprintf(
    "%s vs. %s",
    sub(" \\(.*\\)$", "", nut_lab),
    ratio_name
  )
  plots_ba[[i]] <- make_ratio_plot(dat_ba, ratio_col, nut_col, nut_lab, ratio_title)
}

invisible(lapply(plots_ba, print))
