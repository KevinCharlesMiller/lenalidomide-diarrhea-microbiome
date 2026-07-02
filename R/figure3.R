#####
# Figure 3, Supplemental Figure 3, Supplemental Figure 4C-4D
#####

if (!exists("df_taxa_genus")) source("R/00_setup.R")

# ===== FIGURE 3A =====
# Ordered α-diversity and 6 BA heatmap (OnLen samples only)

group_levels <- c("OnLen_control", "OnLen_diarrhea")
group_labels <- c(
  OnLen_control  = "Len (control)",
  OnLen_diarrhea = "Len (diarrhea)"
)

pal_dark2 <- RColorBrewer::brewer.pal(8, "Dark2")
cols_group <- c(
  OnLen_control  = pal_dark2[2],
  OnLen_diarrhea = pal_dark2[3]
)

plot_dat <- df_ratio |>
  dplyr::transmute(
    Sample.ID,
    SampleCat = forcats::fct_relevel(as.character(SampleCat), group_levels),
    CA   = suppressWarnings(as.numeric(CA)),
    CDCA = suppressWarnings(as.numeric(CDCA)),
    DCA  = suppressWarnings(as.numeric(DCA)),
    LCA  = suppressWarnings(as.numeric(LCA)),
    GCA  = suppressWarnings(as.numeric(GCA)),
    TCA  = suppressWarnings(as.numeric(TCA)
    )) |>
  dplyr::filter(SampleCat %in% group_levels) |>
  dplyr::left_join(
    df_alpha,
    by = "Sample.ID"
  ) |>
  dplyr::filter(is.finite(invsimpson)) |>
  dplyr::mutate(
    dplyr::across(
      c(CA, CDCA, DCA, LCA, GCA, TCA),
      ~ log10(pmax(., 0) + LOG10_OFFSET),
      .names = "{.col}_log10"
    ),
    alpha_log10    = log10(invsimpson)
  ) |>
  dplyr::arrange(SampleCat, dplyr::desc(invsimpson)) |>
  dplyr::group_by(SampleCat) |>
  dplyr::mutate(rank_in_group = dplyr::row_number()) |>
  dplyr::ungroup() |>
  dplyr::mutate(x = dplyr::row_number())

gbounds <- plot_dat |>
  dplyr::count(SampleCat, name = "n") |>
  dplyr::mutate(
    right = cumsum(n),
    left  = right - n + 1,
    mid   = (left + right) / 2
  )

fig3_p_group <- ggplot(plot_dat, aes(x = x, y = 1, fill = SampleCat)) +
  geom_tile(height = 1) +
  scale_fill_manual(
    values = cols_group,
    labels = group_labels,
    name   = NULL
  ) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  {
    if (nrow(gbounds) > 1) {
      geom_vline(
        data        = gbounds[-nrow(gbounds), ],
        aes(xintercept = right + 0.5),
        inherit.aes = FALSE,
        color       = "white",
        linewidth   = 1.1
      )
    } else {
      NULL
    }
  } +
  theme_void(base_size = 11) +
  guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
  theme(
    legend.position   = "top",
    legend.box.margin = margin(t = 4, r = 0, b = 10, l = 0),
    legend.margin     = margin(2, 6, 2, 6),
    plot.margin       = margin(8, 10, 2, 10)
  )

fig3_p_alpha <- ggplot(plot_dat, aes(x = x, y = alpha_log10)) +
  geom_point(size = 2) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(breaks = c(0.5, 1.0, 1.5)) +
  labs(y = "log10(Inverse Simpson)", x = NULL) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x         = element_blank(),
    panel.grid.major.x  = element_blank(),
    panel.grid.minor    = element_blank(),
    plot.margin         = margin(0, 10, 0, 10)
  )

ba_order <- c("CA", "CDCA", "DCA", "LCA", "GCA", "TCA")
ba_labs  <- c(
  CA   = "Cholate (CA)",
  CDCA = "Chenodeoxycholate (CDCA)",
  DCA  = "Deoxycholate (DCA)",
  LCA  = "Lithocholate (LCA)",
  GCA  = "Glycocholate (GCA)",
  TCA  = "Taurocholate (TCA)"
)
ba_cols <- paste0(ba_order, "_log10")

hm_log <- plot_dat |>
  dplyr::select(x, dplyr::all_of(ba_cols)) |>
  tidyr::pivot_longer(
    ba_cols,
    names_to  = "Analyte_raw",
    values_to = "log10_value"
  ) |>
  dplyr::mutate(
    code    = sub("_log10$", "", Analyte_raw),
    Analyte = factor(ba_labs[code], levels = rev(ba_labs[ba_order]))
  )

ca_limits <- range(plot_dat$CA_log10, na.rm = TRUE)
if (!all(is.finite(ca_limits))) ca_limits <- c(-3, 3)

fig3_p_heat_log <- ggplot(hm_log, aes(x = x, y = Analyte, fill = log10_value)) +
  geom_tile(width = 1, height = 1, colour = NA) +
  scale_fill_viridis_c(
    option   = "cividis",
    name     = "log10(AUC) (CA scale)",
    limits   = ca_limits,
    oob      = scales::squish,
    na.value = "grey90"
  ) +
  scale_x_continuous(expand = c(0, 0), breaks = NULL) +
  labs(y = NULL, x = NULL) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid   = element_blank(),
    axis.title.y = element_blank(),
    plot.margin  = margin(0, 10, 0, 10)
  )

FIGURE_3A <- fig3_p_group / fig3_p_alpha / fig3_p_heat_log +
  plot_layout(heights = c(0.10, 1.00, 0.70))

print(FIGURE_3A)

# ===== FIGURE 3B + SUPPLEMENTAL FIGURE 3A =====
# BA/PICRUSt vs α-diversity (LME)

rownames(metadata) <- metadata$Sample.ID

analysis_data <- bile_focus |>
  dplyr::left_join(
    metadata |> dplyr::select(Sample.ID, invsimpson, SampleCategory),
    by = "Sample.ID"
  ) |>
  dplyr::mutate(
    SampleCat = dplyr::coalesce(as.character(SampleCat), as.character(SampleCategory)),
    SampleCat = factor(
      SampleCat,
      levels = c("Baseline", "OnLen_control", "OnLen_diarrhea")
    ),
    log_alpha      = log10(invsimpson),
    CA_DCA_log10   = logdiff(CA,   DCA),
    CDCA_LCA_log10 = logdiff(CDCA, LCA)
  ) |>
  dplyr::filter(
    SampleCat %in% c("OnLen_control", "OnLen_diarrhea"),
    is.finite(log_alpha)
  )

bile_acids <- c("CA", "CDCA", "DCA", "LCA", "GCA", "TCA")

results_list <- list()

for (acid in bile_acids) {
  cat("\n=== Running LME for", acid, "===\n")
  
  y_col <- paste0("log10_", acid)
  
  dat <- analysis_data |>
    dplyr::select(SubjectID, Sample.ID, SampleCat, log_alpha, y = all_of(y_col)) |>
    dplyr::filter(is.finite(y))
  
  cat("Sample size:", nrow(dat), "\n")
  cat("Number of patients:", length(unique(dat$SubjectID)), "\n")
  
  mod <- lmer(y ~ log_alpha + (1 | SubjectID), data = dat)
  
  print(summary(mod))
  
  results_list[[acid]] <- list(
    model = mod,
    data = dat,
    acid = acid
  )
}

save_svg <- function(plot, fname) {
  invisible(NULL)
}

for (acid in bile_acids) {
  cat("\n=== Creating plot for", acid, "===\n")
  
  res <- results_list[[acid]]
  mod <- res$model
  dat <- res$data
  
  coef_alpha <- fixef(mod)["log_alpha"]
  summ <- summary(mod)
  coef_table <- coef(summ)
  p_val <- coef_table["log_alpha", "Pr(>|t|)"]
  
  p_text <- ifelse(p_val < .001, "<0.001", sprintf("= %.3f", p_val))
  label_text <- sprintf("B = %.2f, p%s", coef_alpha, p_text)
  
  x_range <- range(dat$log_alpha, na.rm = TRUE)
  y_range <- range(dat$y, na.rm = TRUE)
  pred_df <- data.frame(
    log_alpha = seq(x_range[1], x_range[2], length.out = 100)
  )
  pred_df$y_pred <- predict(mod, newdata = pred_df, re.form = NA)
  
  xpos <- x_range[1] + 0.02 * diff(x_range)
  ypos <- y_range[2] - 0.06 * diff(y_range)
  
  p <- ggplot(dat, aes(x = log_alpha, y = y)) +
    geom_line(
      aes(group = SubjectID),
      color = "grey60",
      alpha = 0.2,
      linewidth = 0.5
    ) +
    geom_line(
      data = pred_df,
      aes(x = log_alpha, y = y_pred),
      color = "black",
      linewidth = 0.8,
      alpha = 0.9
    ) +
    geom_point(
      aes(
        color = SampleCat,
        shape = SampleCat,
        size = SampleCat,
        alpha = SampleCat
      )
    ) +
    scale_color_manual(
      values = c(
        "OnLen_control" = "#D95F02",
        "OnLen_diarrhea" = "#7570B3"
      ),
      labels = c(
        "OnLen_control" = "Len (control)",
        "OnLen_diarrhea" = "Len (diarrhea)"
      ),
      name = NULL
    ) +
    scale_shape_manual(
      values = c(
        "OnLen_control" = 16,
        "OnLen_diarrhea" = 17
      ),
      labels = c(
        "OnLen_control" = "Len (control)",
        "OnLen_diarrhea" = "Len (diarrhea)"
      ),
      name = NULL
    ) +
    scale_size_manual(
      values = c(
        "OnLen_control" = 2.0,
        "OnLen_diarrhea" = 2.6
      ),
      guide = "none"
    ) +
    scale_alpha_manual(
      values = c(
        "OnLen_control" = 0.9,
        "OnLen_diarrhea" = 0.9
      ),
      guide = "none"
    ) +
    labs(
      x = "log10(Inverse Simpson)",
      y = paste0("log10(", acid, ")"),
      title = paste(acid, "vs. Alpha-diversity (LME)")
    ) +
    annotate(
      "text",
      x = xpos,
      y = ypos,
      label = label_text,
      hjust = 0,
      size = 3.6
    ) +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "right",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid = element_blank(),
      axis.line = element_line(color = "black", linewidth = 0.4),
      axis.ticks = element_line(color = "black", linewidth = 0.4)
    )
  
  print(p)
  
  save_svg(p, paste0("Fig3B_", acid, "_vs_alpha_LME.svg"))
}

analysis_data <- bile_focus |>
  dplyr::left_join(
    metadata |> dplyr::select(Sample.ID, invsimpson, SampleCategory),
    by = "Sample.ID"
  ) |>
  dplyr::mutate(
    SampleCat = dplyr::coalesce(as.character(SampleCat), as.character(SampleCategory)),
    SampleCat = factor(
      SampleCat,
      levels = c("Baseline", "OnLen_control", "OnLen_diarrhea")
    ),
    log_alpha      = log10(invsimpson),
    CA_DCA_log10   = logdiff(CA,   DCA),
    CDCA_LCA_log10 = logdiff(CDCA, LCA)
  ) |>
  dplyr::filter(
    SampleCat %in% c("OnLen_control", "OnLen_diarrhea"),
    is.finite(log_alpha)
  )

ratios <- c("CA_DCA", "CDCA_LCA")

for (ratio in ratios) {
  cat("\n=== Running LME for", ratio, "===\n")
  
  y_col <- paste0(ratio, "_log10")
  
  dat <- analysis_data |>
    dplyr::select(SubjectID, Sample.ID, SampleCat, log_alpha, y = all_of(y_col)) |>
    dplyr::filter(is.finite(y))
  
  cat("Sample size:", nrow(dat), "\n")
  cat("Number of patients:", length(unique(dat$SubjectID)), "\n")
  
  mod <- lmer(y ~ log_alpha + (1 | SubjectID), data = dat)
  
  print(summary(mod))
  
  results_list[[ratio]] <- list(
    model = mod,
    data = dat,
    ratio = ratio
  )
}

ratio_titles <- c(
  CA_DCA = "CA/DCA ratio",
  CDCA_LCA = "CDCA/LCA ratio"
)

for (ratio in ratios) {
  cat("\n=== Creating plot for", ratio, "===\n")
  
  res <- results_list[[ratio]]
  mod <- res$model
  dat <- res$data
  
  coef_alpha <- fixef(mod)["log_alpha"]
  summ <- summary(mod)
  coef_table <- coef(summ)
  p_val <- coef_table["log_alpha", "Pr(>|t|)"]
  
  p_text <- ifelse(p_val < .001, "< 0.001", sprintf("= %.3f", p_val))
  label_text <- sprintf("B = %.2f, p %s", coef_alpha, p_text)
  
  x_range <- range(dat$log_alpha, na.rm = TRUE)
  y_range <- range(dat$y, na.rm = TRUE)
  pred_df <- data.frame(
    log_alpha = seq(x_range[1], x_range[2], length.out = 100)
  )
  pred_df$y_pred <- predict(mod, newdata = pred_df, re.form = NA)
  
  xpos <- x_range[1] + 0.02 * diff(x_range)
  ypos <- y_range[2] - 0.06 * diff(y_range)
  
  ratio_label <- gsub("_", "/", ratio)
  
  p <- ggplot(dat, aes(x = log_alpha, y = y)) +
    geom_line(
      aes(group = SubjectID),
      color = "grey60",
      alpha = 0.2,
      linewidth = 0.5 ,
      alpha = 0.9
    ) +
    geom_line(
      data = pred_df,
      aes(x = log_alpha, y = y_pred),
      color = "black",
      linewidth = 0.8
    ) +
    geom_point(
      aes(
        color = SampleCat,
        shape = SampleCat,
        size = SampleCat,
        alpha = SampleCat
      )
    ) +
    scale_color_manual(
      values = c(
        "OnLen_control" = "#D95F02",
        "OnLen_diarrhea" = "#7570B3"
      ),
      labels = c(
        "OnLen_control" = "Len (control)",
        "OnLen_diarrhea" = "Len (diarrhea)"
      ),
      name = NULL
    ) +
    scale_shape_manual(
      values = c(
        "OnLen_control" = 16,
        "OnLen_diarrhea" = 17
      ),
      labels = c(
        "OnLen_control" = "Len (control)",
        "OnLen_diarrhea" = "Len (diarrhea)"
      ),
      name = NULL
    ) +
    scale_size_manual(
      values = c(
        "OnLen_control" = 2.0,
        "OnLen_diarrhea" = 2.6
      ),
      guide = "none"
    ) +
    scale_alpha_manual(
      values = c(
        "OnLen_control" = 0.9,
        "OnLen_diarrhea" = 0.90
      ),
      guide = "none"
    ) +
    labs(
      x = "log10(Inverse Simpson)",
      y = paste0("log10(", ratio_label, ")"),
      title = paste(ratio_titles[ratio], "vs. Alpha-diversity (LME)")
    ) +
    annotate(
      "text",
      x = xpos,
      y = ypos,
      label = label_text,
      hjust = 0,
      size = 3.6
    ) +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "right",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid = element_blank(),
      axis.line = element_line(color = "black", linewidth = 0.4),
      axis.ticks = element_line(color = "black", linewidth = 0.4)
    )
  
  print(p)
  
  save_svg(p, paste0("Fig3B_", ratio, "_vs_alpha_LME.svg"))
}

add_alpha_onlen <- function(d) {
  d |>
    dplyr::left_join(metadata |> dplyr::select(Sample.ID, invsimpson), by = "Sample.ID") |>
    dplyr::mutate(log_alpha = log10(invsimpson)) |>
    dplyr::filter(SampleCat %in% c("OnLen_control", "OnLen_diarrhea"),
                  is.finite(log_alpha))
}

ko_list <- list(
  BSH  = list(data = add_alpha_onlen(bsh_df),  label = "BSH (K01442)"),
  HSDH = list(data = add_alpha_onlen(hsdh_df), label = "7α-HSDH (K00076)"),
  BAI  = list(data = add_alpha_onlen(bai_df),  label = "bai operon")
)

ko_results <- list()

for (ko_name in names(ko_list)) {
  cat("\n=== Running LME for", ko_list[[ko_name]]$label, "===\n")
  
  dat <- ko_list[[ko_name]]$data
  
  cat("Sample size:", nrow(dat), "\n")
  cat("Number of patients:", length(unique(dat$SubjectID)), "\n")
  
  mod <- lmer(log_CPM ~ log_alpha + (1 | SubjectID), data = dat)
  
  print(summary(mod))
  
  ko_results[[ko_name]] <- list(
    model = mod,
    data = dat,
    label = ko_list[[ko_name]]$label
  )
}

for (ko_name in names(ko_results)) {
  cat("\n=== Creating plot for", ko_name, "===\n")
  
  res <- ko_results[[ko_name]]
  mod <- res$model
  dat <- res$data
  ko_label <- res$label
  
  coef_alpha <- fixef(mod)["log_alpha"]
  summ <- summary(mod)
  coef_table <- coef(summ)
  p_val <- coef_table["log_alpha", "Pr(>|t|)"]
  
  p_text <- ifelse(p_val < .001, "< 0.001", sprintf("= %.3f", p_val))
  label_text <- sprintf("B = %.2f, p %s", coef_alpha, p_text)
  
  x_range <- range(dat$log_alpha, na.rm = TRUE)
  y_range <- range(dat$log_CPM, na.rm = TRUE)
  pred_df <- data.frame(
    log_alpha = seq(x_range[1], x_range[2], length.out = 100)
  )
  pred_df$y_pred <- predict(mod, newdata = pred_df, re.form = NA)
  
  xpos <- x_range[1] + 0.02 * diff(x_range)
  ypos <- y_range[2] - 0.06 * diff(y_range)
  
  p <- ggplot(dat, aes(x = log_alpha, y = log_CPM)) +
    geom_line(
      aes(group = SubjectID),
      color = "grey60",
      alpha = 0.2,
      linewidth = 0.5
    ) +
    geom_line(
      data = pred_df,
      aes(x = log_alpha, y = y_pred),
      color = "black",
      linewidth = 0.8
    ) +
    geom_point(
      aes(
        color = SampleCat,
        shape = SampleCat,
        size = SampleCat,
        alpha = SampleCat
      )
    ) +
    scale_color_manual(
      values = c(
        "OnLen_control" = "#D95F02",
        "OnLen_diarrhea" = "#7570B3"
      ),
      labels = c(
        "OnLen_control" = "Len (control)",
        "OnLen_diarrhea" = "Len (diarrhea)"
      ),
      name = NULL
    ) +
    scale_shape_manual(
      values = c(
        "OnLen_control" = 16,
        "OnLen_diarrhea" = 17
      ),
      labels = c(
        "OnLen_control" = "Len (control)",
        "OnLen_diarrhea" = "Len (diarrhea)"
      ),
      name = NULL
    ) +
    scale_size_manual(
      values = c(
        "OnLen_control" = 2.0,
        "OnLen_diarrhea" = 2.6
      ),
      guide = "none"
    ) +
    scale_alpha_manual(
      values = c(
        "OnLen_control" = 0.9,
        "OnLen_diarrhea" = 0.90
      ),
      guide = "none"
    ) +
    labs(
      x = "log10(Inverse Simpson)",
      y = "log10(CPM)",
      title = paste(ko_label, "vs. Alpha-diversity (LME)")
    ) +
    annotate(
      "text",
      x = xpos,
      y = ypos,
      label = label_text,
      hjust = 0,
      size = 3.6
    ) +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "right",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid = element_blank(),
      axis.line = element_line(color = "black", linewidth = 0.4),
      axis.ticks = element_line(color = "black", linewidth = 0.4)
    )
  
  print(p)
  
  save_svg(p, paste0("Fig3C_", ko_name, "_vs_alpha_LME.svg"))
}

# ===== FIGURE 3C =====
# Bray-Curtis PCoA overlays, PICRUSt

bray_all <- phyloseq::distance(phy, method = "bray")
ord_all  <- phyloseq::ordinate(phy, method = "PCoA", distance = bray_all)

plot_onlen_locked_axes_with_marker <- function(
    ord_all,
    bray_all,
    metadata,
    marker_df,
    marker_label  = "BSH",
    colorbar_name = sprintf("%s log10(CPM)", marker_label),
    point_size    = 4.5,
    point_alpha   = 0.85,
    show_ellipse  = FALSE,
    zoom_mult = 0.14,
    nperm         = 999
) {
  onlen_levels <- c("OnLen_control", "OnLen_diarrhea")
  
  meta_onlen <- metadata |>
    dplyr::transmute(
      Sample.ID,
      SampleCategory = factor(
        SampleCategory,
        levels = c("Baseline", onlen_levels)
      ),
      SubjectID = factor(SubjectID)
    ) |>
    dplyr::filter(SampleCategory %in% onlen_levels)
  
  keep_ids <- Reduce(
    intersect,
    list(
      rownames(ord_all$vectors),
      rownames(as.matrix(bray_all)),
      meta_onlen$Sample.ID,
      marker_df$Sample.ID
    )
  )
  if (length(keep_ids) < 6L) {
    stop("Too few overlapping OnLen samples.")
  }
  
  df_ord <- ord_all$vectors |>
    as.data.frame() |>
    tibble::rownames_to_column("Sample.ID") |>
    dplyr::select(Sample.ID, Axis.1, Axis.2) |>
    dplyr::filter(Sample.ID %in% keep_ids) |>
    dplyr::left_join(meta_onlen, by = "Sample.ID") |>
    dplyr::left_join(
      marker_df |>
        dplyr::transmute(Sample.ID, marker_logCPM = log_CPM),
      by = "Sample.ID"
    ) |>
    dplyr::filter(is.finite(marker_logCPM)) |>
    droplevels()
  
  Dmat   <- as.matrix(bray_all)
  keep2  <- df_ord$Sample.ID
  D_onln <- stats::as.dist(Dmat[keep2, keep2, drop = FALSE])
  
  ctrl <- if (any(table(df_ord$SubjectID) > 1)) {
    h <- permute::how(nperm = nperm)
    permute::setBlocks(h) <- df_ord$SubjectID
    h
  } else {
    nperm
  }
  
  fit <- vegan::adonis2(
    D_onln ~ marker_logCPM + SampleCategory,
    data        = df_ord,
    permutations = ctrl,
    by          = "margin"
  )
  
  tab <- as.data.frame(fit)
  
  f_col  <- dplyr::case_when(
    "F"       %in% names(tab) ~ "F",
    "F.Model" %in% names(tab) ~ "F.Model",
    TRUE                      ~ names(tab)[2]
  )
  r2_col <- dplyr::case_when(
    "R2" %in% names(tab) ~ "R2",
    TRUE ~ names(tab)[grep("R2", names(tab), fixed = TRUE)[1]]
  )
  
  F_mrk <- as.numeric(tab["marker_logCPM", f_col])
  R2_mrk<- as.numeric(tab["marker_logCPM", r2_col])
  P_mrk <- as.numeric(tab["marker_logCPM", "Pr(>F)"])
  
  ann <- sprintf(
    "%s PERMANOVA  F=%.2f, p %s",
    marker_label,
    F_mrk,
    ifelse(is.finite(P_mrk) && P_mrk < .001,
           "<0.001",
           sprintf("= %.3f", P_mrk))
  )
  
  pc_var <- ord_all$values$Relative_eig
  x_lab  <- if (!is.null(pc_var) && length(pc_var) >= 1)
    sprintf("PCoA 1 (%.1f%%)", 100 * pc_var[1]) else "PCoA 1"
  y_lab  <- if (!is.null(pc_var) && length(pc_var) >= 2)
    sprintf("PCoA 2 (%.1f%%)", 100 * pc_var[2]) else "PCoA 2"
  
  shp <- c(
    OnLen_control  = 16,
    OnLen_diarrhea = 17
  )
  lab <- c(
    OnLen_control  = "Len (control)",
    OnLen_diarrhea = "Len (diarrhea)"
  )
  
  p <- ggplot2::ggplot(df_ord, ggplot2::aes(Axis.1, Axis.2)) +
    ggplot2::geom_point(
      ggplot2::aes(
        color = marker_logCPM,
        shape = SampleCategory
      ),
      size  = point_size,
      alpha = point_alpha,
      na.rm = TRUE
    ) +
    ggplot2::scale_color_viridis_c(
      option   = "inferno",
      name     = colorbar_name,
      na.value = "grey85"
    ) +
    ggplot2::scale_shape_manual(
      values = shp,
      labels = lab,
      name   = NULL
    ) +
    ggplot2::labs(
      x     = x_lab,
      y     = y_lab,
      title = sprintf(
        "Bray-Curtis PCoA (OnLen) — %s intensity",
        marker_label
      )
    ) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = zoom_mult)) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = zoom_mult)) +
    ggplot2::coord_equal(expand = TRUE, clip = "off") +
    ggplot2::theme_classic(base_size = 13) +
    ggplot2::theme(
      legend.position = "right",
      panel.grid      = ggplot2::element_blank()
    ) +
    ggplot2::annotate(
      "text",
      x     = Inf,
      y     = -Inf,
      hjust = 1.02,
      vjust = -0.6,
      label = ann,
      size  = 3.6
    )
  
  list(
    plot      = p,
    permanova = fit,
    data      = df_ord,
    stats     = tibble::tibble(
      marker = marker_label,
      F      = F_mrk,
      R2     = R2_mrk,
      p      = P_mrk,
      n      = nrow(df_ord)
    )
  )
}

res_bsh  <- plot_onlen_locked_axes_with_marker(
  ord_all, bray_all, metadata, bsh_df,
  marker_label  = "BSH",
  colorbar_name = "log10(CPM)\nBSH"
)
res_hsdh <- plot_onlen_locked_axes_with_marker(
  ord_all, bray_all, metadata, hsdh_df,
  marker_label  = "7α-HSDH",
  colorbar_name = "log10(CPM)\n7α-HSDH"
)
res_bai  <- plot_onlen_locked_axes_with_marker(
  ord_all, bray_all, metadata, bai_df,
  marker_label  = "bai operon",
  colorbar_name = "log10(CPM)\nbai operon"
)

print(res_bsh$plot)
print(res_hsdh$plot)
print(res_bai$plot)

stats_tbl <- dplyr::bind_rows(
  res_bsh$stats,
  res_hsdh$stats,
  res_bai$stats
) |>
  dplyr::mutate(
    p_fmt = dplyr::if_else(
      !is.finite(p),
      NA_character_,
      dplyr::if_else(p < .001, "<0.001", sprintf("%.3f", p))
    )
  ) |>
  dplyr::select(marker, n, F, R2, p = p_fmt)

print(stats_tbl)

# ===== SUPPLEMENTAL FIGURE 3B =====
# Baseline PCoA, AHCT covariate

meta_beta <- metadata %>%
  dplyr::mutate(
    Sample.ID     = as.character(Sample.ID),
    LaterDiarrhea = factor(LaterDiarrhea, levels = c("No", "Yes")),
    autoHCT_fac   = factor(autoHCT, levels = c(0, 1), labels = c("No AHCT", "AHCT"))
  ) %>%
  dplyr::filter(!is.na(LaterDiarrhea), !is.na(autoHCT)) %>%
  dplyr::distinct(Sample.ID, .keep_all = TRUE)

keep_ids_beta <- intersect(sample_names(phy), meta_beta$Sample.ID)
stopifnot(length(keep_ids_beta) >= 4)

phy_beta  <- prune_samples(keep_ids_beta, phy)
bray_beta <- phyloseq::distance(phy_beta, method = "bray")
ord_beta  <- ordinate(phy_beta, method = "PCoA", distance = bray_beta)

meta_perm <- meta_beta %>%
  dplyr::filter(Sample.ID %in% rownames(ord_beta$vectors)) %>%
  dplyr::select(Sample.ID, LaterDiarrhea, autoHCT_fac)
rownames(meta_perm) <- meta_perm$Sample.ID

adon_tab <- as.data.frame(
  vegan::adonis2(
    bray_beta ~ autoHCT_fac + LaterDiarrhea,
    data         = meta_perm,
    permutations = 999,
    by           = "margin"
  )
)

extract_perm_stats <- function(tab, term) {
  fcol <- if ("F" %in% names(tab)) "F" else if ("F.Model" %in% names(tab)) "F.Model" else names(tab)[2]
  r2col <- if ("R2" %in% names(tab)) "R2" else names(tab)[grep("R2", names(tab), fixed = TRUE)[1]]
  Fv  <- as.numeric(tab[term, fcol])
  R2v <- as.numeric(tab[term, r2col])
  Pv  <- as.numeric(tab[term, "Pr(>F)"])
  list(F = Fv, R2 = R2v, p = Pv)
}

stats_ahct <- extract_perm_stats(adon_tab, "autoHCT_fac")
stats_lrd  <- extract_perm_stats(adon_tab, "LaterDiarrhea")

fmt_p <- function(p) ifelse(is.finite(p) && p < .001, "<0.001", sprintf("= %.3g", p))

ann_lrd  <- sprintf("LRD: F=%.2f, R\u00B2=%.3f, p%s", stats_lrd$F, stats_lrd$R2, fmt_p(stats_lrd$p))
ann_ahct <- sprintf("AHCT: F=%.2f, R\u00B2=%.3f, p%s", stats_ahct$F, stats_ahct$R2, fmt_p(stats_ahct$p))

cat("\n--- PERMANOVA (marginal): autoHCT + LaterDiarrhea ---\n")
print(adon_tab)
cat("\n", ann_lrd, "\n", ann_ahct, "\n")

pcv   <- ord_beta$values$Relative_eig
x_lab <- if (!is.null(pcv) && length(pcv) >= 1 && is.finite(pcv[1])) sprintf("PCoA1 (%.1f%%)", 100 * pcv[1]) else "PCoA1"
y_lab <- if (!is.null(pcv) && length(pcv) >= 2 && is.finite(pcv[2])) sprintf("PCoA2 (%.1f%%)", 100 * pcv[2]) else "PCoA2"

df_ordination <- ord_beta$vectors %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Sample.ID") %>%
  dplyr::select(Sample.ID, Axis.1, Axis.2) %>%
  dplyr::left_join(meta_perm, by = "Sample.ID")

cat_cols_lrd <- c("No" = "#1B9E77", "Yes" = "#D55E00")
shp_ahct     <- c("No AHCT" = 16, "AHCT" = 17)  

fig4c_revised <- ggplot(df_ordination, aes(Axis.1, Axis.2)) +
  stat_ellipse(
    aes(group = LaterDiarrhea, fill = LaterDiarrhea),
    geom = "polygon", level = 0.95, alpha = 0.1,
    colour = NA, show.legend = FALSE
  ) +
  stat_ellipse(
    aes(color = LaterDiarrhea),
    level = 0.95, linewidth = 1, show.legend = TRUE
  ) +
  geom_point(
    aes(color = LaterDiarrhea, shape = autoHCT_fac),
    size = 2.4, alpha = 0.7
  ) +
  scale_color_manual(
    values = cat_cols_lrd,
    name   = "Len-related\ndiarrhea",
    labels = c("No" = "No", "Yes" = "Yes")
  ) +
  scale_fill_manual(values = cat_cols_lrd, guide = "none") +
  scale_shape_manual(
    values = shp_ahct,
    name   = "Prior AHCT"
  ) +
  labs(x = x_lab, y = y_lab, title = "Beta-diversity (Bray-Curtis)") +
  scale_x_continuous(expand = expansion(mult = 0.12)) +
  scale_y_continuous(expand = expansion(mult = 0.12)) +
  coord_equal(expand = TRUE, clip = "off") +
  theme_classic(base_size = 14) +
  theme(
    legend.position  = "right",
    legend.direction = "vertical",
    legend.spacing.y = grid::unit(4, "pt"),
    panel.grid       = element_blank(),
    plot.margin      = margin(6, 10, 6, 6)
  ) +
  annotate(
    "text", x = Inf, y = -Inf,
    hjust = 1.02, vjust = -1.8,
    size = 3.3, fontface = "italic",
    label = ann_lrd
  ) +
  annotate(
    "text", x = Inf, y = -Inf,
    hjust = 1.02, vjust = -0.4,
    size = 3.3, fontface = "italic",
    label = ann_ahct
  )

print(fig4c_revised)

# ===== SUPPLEMENTAL FIGURE 3B =====
# Baseline PCoA, cohort covariate

meta_beta_irb <- metadata %>%
  dplyr::mutate(
    Sample.ID     = as.character(Sample.ID),
    LaterDiarrhea = factor(LaterDiarrhea, levels = c("No", "Yes")),
    Cohort_fac       = factor(Cohort)
  ) %>%
  dplyr::filter(!is.na(LaterDiarrhea), !is.na(Cohort)) %>%
  dplyr::distinct(Sample.ID, .keep_all = TRUE)

keep_ids_beta_irb <- intersect(sample_names(phy), meta_beta_irb$Sample.ID)
stopifnot(length(keep_ids_beta_irb) >= 4)

phy_beta_irb  <- prune_samples(keep_ids_beta_irb, phy)
bray_beta_irb <- phyloseq::distance(phy_beta_irb, method = "bray")
ord_beta_irb  <- ordinate(phy_beta_irb, method = "PCoA", distance = bray_beta_irb)

meta_perm_irb <- meta_beta_irb %>%
  dplyr::filter(Sample.ID %in% rownames(ord_beta_irb$vectors)) %>%
  dplyr::select(Sample.ID, LaterDiarrhea, Cohort_fac)
rownames(meta_perm_irb) <- meta_perm_irb$Sample.ID

adon_tab_irb <- as.data.frame(
  vegan::adonis2(
    bray_beta_irb ~ Cohort_fac + LaterDiarrhea,
    data         = meta_perm_irb,
    permutations = 999,
    by           = "margin"
  )
)

extract_perm_stats <- function(tab, term) {
  fcol  <- if ("F" %in% names(tab)) "F" else if ("F.Model" %in% names(tab)) "F.Model" else names(tab)[2]
  r2col <- if ("R2" %in% names(tab)) "R2" else names(tab)[grep("R2", names(tab), fixed = TRUE)[1]]
  Fv  <- as.numeric(tab[term, fcol])
  R2v <- as.numeric(tab[term, r2col])
  Pv  <- as.numeric(tab[term, "Pr(>F)"])
  list(F = Fv, R2 = R2v, p = Pv)
}

stats_irb <- extract_perm_stats(adon_tab_irb, "Cohort_fac")
stats_lrd_irb <- extract_perm_stats(adon_tab_irb, "LaterDiarrhea")

fmt_p <- function(p) ifelse(is.finite(p) && p < .001, "<0.001", sprintf("= %.3g", p))

ann_lrd_irb <- sprintf("LRD: F=%.2f, R\u00B2=%.3f, p%s",
                        stats_lrd_irb$F, stats_lrd_irb$R2, fmt_p(stats_lrd_irb$p))
ann_irb     <- sprintf("Study: F=%.2f, R\u00B2=%.3f, p%s",
                        stats_irb$F, stats_irb$R2, fmt_p(stats_irb$p))

cat("\n--- PERMANOVA (marginal): Cohort + LaterDiarrhea ---\n")
print(adon_tab_irb)
cat("\n", ann_lrd_irb, "\n", ann_irb, "\n")

pcv_irb <- ord_beta_irb$values$Relative_eig
x_lab_irb <- if (!is.null(pcv_irb) && length(pcv_irb) >= 1 && is.finite(pcv_irb[1])) {
  sprintf("PCoA1 (%.1f%%)", 100 * pcv_irb[1])
} else "PCoA1"
y_lab_irb <- if (!is.null(pcv_irb) && length(pcv_irb) >= 2 && is.finite(pcv_irb[2])) {
  sprintf("PCoA2 (%.1f%%)", 100 * pcv_irb[2])
} else "PCoA2"

df_ordination_irb <- ord_beta_irb$vectors %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Sample.ID") %>%
  dplyr::select(Sample.ID, Axis.1, Axis.2) %>%
  dplyr::left_join(meta_perm_irb, by = "Sample.ID")

cat_cols_lrd <- c("No" = "#1B9E77", "Yes" = "#D55E00")

irb_levels <- levels(df_ordination_irb$Cohort_fac)
shp_irb    <- setNames(c(16, 17, 15, 18)[seq_along(irb_levels)], irb_levels)

fig_beta_irb <- ggplot(df_ordination_irb, aes(Axis.1, Axis.2)) +
  stat_ellipse(
    aes(group = LaterDiarrhea, fill = LaterDiarrhea),
    geom = "polygon", level = 0.95, alpha = 0.1,
    colour = NA, show.legend = FALSE
  ) +
  stat_ellipse(
    aes(color = LaterDiarrhea),
    level = 0.95, linewidth = 1, show.legend = TRUE
  ) +
  geom_point(
    aes(color = LaterDiarrhea, shape = Cohort_fac),
    size = 2.4, alpha = 0.7
  ) +
  scale_color_manual(
    values = cat_cols_lrd,
    name   = "Len-related\ndiarrhea",
    labels = c("No" = "No", "Yes" = "Yes")
  ) +
  scale_fill_manual(values = cat_cols_lrd, guide = "none") +
  scale_shape_manual(
    values = shp_irb,
    name   = "Study"
  ) +
  labs(x = x_lab_irb, y = y_lab_irb, title = "Beta-diversity (Bray-Curtis)") +
  scale_x_continuous(expand = expansion(mult = 0.12)) +
  scale_y_continuous(expand = expansion(mult = 0.12)) +
  coord_equal(expand = TRUE, clip = "off") +
  theme_classic(base_size = 14) +
  theme(
    legend.position  = "right",
    legend.direction = "vertical",
    legend.spacing.y = grid::unit(4, "pt"),
    panel.grid       = element_blank(),
    plot.margin      = margin(6, 10, 6, 6)
  ) +
  annotate(
    "text", x = Inf, y = -Inf,
    hjust = 1.02, vjust = -1.8,
    size = 3.3, fontface = "italic",
    label = ann_lrd_irb
  ) +
  annotate(
    "text", x = Inf, y = -Inf,
    hjust = 1.02, vjust = -0.4,
    size = 3.3, fontface = "italic",
    label = ann_irb
  )

print(fig_beta_irb)

# ===== SUPPLEMENTAL FIGURE 3C + SUPPLEMENTAL TABLE 8 =====
# envfit biplot (genus contributors to baseline β-diversity)

cat_cols <- c("No" = "#1B9E77", "Yes" = "#D55E00")

meta_beta <- metadata %>%
  dplyr::filter(SampleCategory == "Baseline") %>%
  dplyr::mutate(
    LaterDiarrhea = factor(LaterDiarrhea, levels = c("No", "Yes"))
  ) %>%
  dplyr::filter(!is.na(LaterDiarrhea)) %>%
  dplyr::distinct(Sample.ID, .keep_all = TRUE)

keep_ids <- intersect(sample_names(phy), meta_beta$Sample.ID)
phy_beta <- prune_samples(keep_ids, phy)

cat("Baseline samples:\n")
print(table(meta_beta$LaterDiarrhea))

bray_beta <- phyloseq::distance(phy_beta, method = "bray")
ord_beta  <- ordinate(phy_beta, method = "PCoA", distance = bray_beta)

meta_perm <- meta_beta %>%
  dplyr::filter(Sample.ID %in% rownames(ord_beta$vectors)) %>%
  dplyr::select(Sample.ID, LaterDiarrhea)
rownames(meta_perm) <- meta_perm$Sample.ID

adon_tab <- as.data.frame(vegan::adonis2(bray_beta ~ LaterDiarrhea, data = meta_perm))
Fv <- as.numeric(adon_tab["LaterDiarrhea", "F"])
Pv <- as.numeric(adon_tab["LaterDiarrhea", "Pr(>F)"])
ann_perm <- sprintf("PERMANOVA F=%.2f, p %s", Fv,
                    ifelse(is.finite(Pv) && Pv < .001, "<0.001", sprintf("= %.3g", Pv)))

pcv   <- ord_beta$values$Relative_eig
x_lab <- sprintf("PCoA1 (%.1f%%)", 100 * pcv[1])
y_lab <- sprintf("PCoA2 (%.1f%%)", 100 * pcv[2])

df_ord <- ord_beta$vectors %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Sample.ID") %>%
  dplyr::select(Sample.ID, Axis.1, Axis.2) %>%
  dplyr::left_join(meta_perm, by = "Sample.ID")

phy_genus     <- tax_glom(phy_beta, taxrank = "genus", NArm = FALSE)
phy_genus_rel <- transform_sample_counts(phy_genus, function(x) x / sum(x))

rel_mat <- as(otu_table(phy_genus_rel), "matrix")
if (taxa_are_rows(phy_genus_rel)) rel_mat <- t(rel_mat)

tax_df <- as.data.frame(tax_table(phy_genus)) |> rownames_to_column("otu_id")
colnames(rel_mat) <- tax_df$genus[match(colnames(rel_mat), tax_df$otu_id)]

bad <- is.na(colnames(rel_mat)) | colnames(rel_mat) %in% c("", "g", "NA", "Unclassified")
rel_mat <- rel_mat[, !bad]
colnames(rel_mat) <- make.unique(colnames(rel_mat))
rel_mat <- rel_mat[, apply(rel_mat, 2, var) > 0]
rel_mat <- rel_mat[rownames(ord_beta$vectors), , drop = FALSE]

prevalence <- colMeans(rel_mat > 0)
rel_mat <- rel_mat[, prevalence >= 0.20]
cat("Genera after prevalence filter (≥20%):", ncol(rel_mat), "\n")

cat("Genera retained for envfit:", ncol(rel_mat), "\n")

pcoa_pts <- ord_beta$vectors[, c("Axis.1", "Axis.2")]
set.seed(123)
ef <- envfit(pcoa_pts, rel_mat, permutations = 9999)

contrib_df <- data.frame(
  genus = names(ef$vectors$r),
  Axis1 = unname(ef$vectors$arrows[, 1]),
  Axis2 = unname(ef$vectors$arrows[, 2]),
  r2    = unname(ef$vectors$r),
  pval  = unname(ef$vectors$pvals)
) |>
  mutate(padj = p.adjust(pval, "BH")) |>
  arrange(desc(r2))

contrib_top <- contrib_df |> filter(padj < 0.05) |> slice_head(n = 15)

cat("\nSignificant genus contributors (BH-adj p < 0.05):\n")
print(contrib_top[, c("genus", "r2", "padj")])

grp <- meta_perm$LaterDiarrhea[match(rownames(rel_mat), meta_perm$Sample.ID)]

grp_means <- as.data.frame(rel_mat) |>
  rownames_to_column("Sample.ID") |>
  pivot_longer(-Sample.ID, names_to = "genus", values_to = "abund") |>
  left_join(data.frame(Sample.ID = rownames(rel_mat), outcome = grp), by = "Sample.ID") |>
  group_by(genus, outcome) |>
  summarise(mean_abund = mean(abund), .groups = "drop") |>
  pivot_wider(names_from = outcome, values_from = mean_abund) |>
  mutate(enriched_in = ifelse(Yes > No, "LRD", "No LRD"))

contrib_top <- contrib_top |> left_join(grp_means, by = "genus")

arrow_df <- contrib_top |>
  mutate(
    Axis1_scaled = Axis1 * sqrt(r2),
    Axis2_scaled = Axis2 * sqrt(r2)
  )

scale_factor <- max(abs(c(df_ord$Axis.1, df_ord$Axis.2))) * 1 /
                max(abs(c(arrow_df$Axis1_scaled, arrow_df$Axis2_scaled)))

arrow_df <- arrow_df |>
  mutate(
    x_end = Axis1_scaled * scale_factor,
    y_end = Axis2_scaled * scale_factor
  )

arrow_cols <- c("No LRD" = "#1B9E77", "LRD" = "#D55E00")

p_biplot <- ggplot(df_ord, aes(Axis.1, Axis.2)) +
  stat_ellipse(aes(group = LaterDiarrhea, fill = LaterDiarrhea),
               geom = "polygon", level = 0.95, alpha = 0.03,
               colour = NA, show.legend = FALSE) +
  stat_ellipse(aes(color = LaterDiarrhea),
               level = 0.95, linewidth = 0.5, show.legend = FALSE) +
  geom_point(aes(color = LaterDiarrhea), size = 2, alpha = 0.6) +
  geom_segment(data = arrow_df,
               aes(x = 0, y = 0, xend = x_end, yend = y_end),
               arrow = arrow(length = unit(0.22, "cm")),
               colour = "black", linewidth = 0.5, inherit.aes = FALSE) +
  geom_text_repel(data = arrow_df,
                  aes(x = x_end, y = y_end, label = genus),
                  size = 3.8, colour = "black", fontface = "italic",
                  inherit.aes = FALSE,
                  max.overlaps = 20,
                  box.padding = 0.4,
                  point.padding = 0.3,
                  force = 2,
                  min.segment.length = 0,
                  segment.size = 0) +
  scale_color_manual(
    values = cat_cols,
    name = "Len-related\ndiarrhea"
  ) +
  scale_fill_manual(values = cat_cols, guide = "none") +
  labs(x = x_lab, y = y_lab,
       title = "Baseline beta-diversity: top genus contributors") +
  scale_x_continuous(expand = expansion(mult = 0.15)) +
  scale_y_continuous(expand = expansion(mult = 0.15)) +
  coord_equal(expand = TRUE, clip = "off") +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "right",
    legend.direction = "vertical",
    panel.grid = element_blank(),
    plot.margin = margin(6, 10, 6, 6)
  ) 

print(p_biplot)

# ===== SUPPLEMENTAL FIGURE 3D =====
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
  Ruminococcoides  = rel_mat[, "Ruminococcoides"],
  Blautia          = rel_mat[, "Blautia"]
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

p_ruminococcoides <- plot_genus_fig4(dat_genus_bl, "Ruminococcoides", "Ruminococcoides")
p_blautia <- plot_genus_fig4(dat_genus_bl, "Blautia", "Blautia")

# ===== SUPPLEMENTAL FIGURE 4C =====
# BA-ratio PCoA overlays, baseline (CA/DCA, CDCA/LCA)
meta_use <- metadata %>%
  dplyr::mutate(
    Sample.ID     = as.character(Sample.ID),
    LaterDiarrhea = factor(LaterDiarrhea, levels = c("No", "Yes")),
    SubjectID           = dplyr::coalesce(as.character(SubjectID), NA_character_)
  ) %>%
  dplyr::filter(!is.na(LaterDiarrhea)) %>%
  dplyr::distinct(Sample.ID, .keep_all = TRUE)

keep_ids <- intersect(phyloseq::sample_names(phy), meta_use$Sample.ID)
stopifnot(length(keep_ids) >= 4)

phy_sub  <- phyloseq::prune_samples(keep_ids, phy)
meta_sub <- meta_use %>% dplyr::filter(Sample.ID %in% keep_ids)

bray_sub <- phyloseq::distance(phy_sub, method = "bray")
ord_sub  <- phyloseq::ordinate(phy_sub, method = "PCoA", distance = bray_sub)

plot_overlay_like_previous <- function(ord_sub,
                                       bray_sub,
                                       metadata,
                                       marker_df,
                                       marker_col    = "value",
                                       marker_label  = "Marker",
                                       colorbar_name = marker_label,
                                       cat_cols      = c("No" = "#1B9E77", "Yes" = "#D55E00"),
                                       point_size    = 4.5,
                                       point_alpha   = 0.85,
                                       zoom_mult     = 0.14,  
                                       nperm         = 999) {
  meta_use <- metadata %>%
    dplyr::transmute(
      Sample.ID     = as.character(Sample.ID),
      LaterDiarrhea = factor(LaterDiarrhea, levels = c("No", "Yes")),
      SubjectID           = dplyr::coalesce(as.character(SubjectID), NA_character_)
    ) %>%
    dplyr::filter(!is.na(LaterDiarrhea)) %>%
    dplyr::distinct(Sample.ID, .keep_all = TRUE)
  
  keep_ids <- Reduce(
    intersect,
    list(
      rownames(ord_sub$vectors),
      labels(bray_sub),
      meta_use$Sample.ID,
      marker_df$Sample.ID
    )
  )
  if (length(keep_ids) < 6L) stop("Too few overlapping samples to plot.")
  
  df_ord <- ord_sub$vectors %>%
    as.data.frame() %>%
    tibble::rownames_to_column("Sample.ID") %>%
    dplyr::select(Sample.ID, Axis.1, Axis.2) %>%
    dplyr::filter(Sample.ID %in% keep_ids) %>%
    dplyr::left_join(meta_use, by = "Sample.ID") %>%
    dplyr::left_join(
      marker_df %>%
        dplyr::transmute(Sample.ID, marker_val = .data[[marker_col]]),
      by = "Sample.ID"
    ) %>%
    dplyr::filter(is.finite(marker_val)) %>%
    droplevels()
  
  Dmat  <- as.matrix(bray_sub)
  D_sub <- stats::as.dist(Dmat[df_ord$Sample.ID, df_ord$Sample.ID, drop = FALSE])
  
  ctrl <- nperm
  if (!all(is.na(df_ord$SubjectID)) && any(table(df_ord$SubjectID) > 1)) {
    h <- permute::how(nperm = nperm)
    permute::setBlocks(h) <- df_ord$SubjectID
    ctrl <- h
  }
  
  fit <- vegan::adonis2(
    D_sub ~ marker_val + LaterDiarrhea,
    data        = df_ord,
    permutations = ctrl,
    by          = "margin"
  )
  
  tab   <- as.data.frame(fit)
  f_col <- if ("F" %in% names(tab)) {
    "F"
  } else if ("F.Model" %in% names(tab)) {
    "F.Model"
  } else {
    names(tab)[2]
  }
  
  r2_col <- dplyr::case_when(
    "R2" %in% names(tab) ~ "R2",
    TRUE                 ~ names(tab)[grep("R2", names(tab), fixed = TRUE)[1]]
  )
  
  F_mrk <- suppressWarnings(as.numeric(tab["marker_val", f_col]))
  R2_mrk <- suppressWarnings(as.numeric(tab["marker_val", r2_col]))
  P_mrk <- suppressWarnings(as.numeric(tab["marker_val", "Pr(>F)"]))
  
  ann <- sprintf(
    "%s PERMANOVA  F=%.2f, p %s",   
    marker_label,
    F_mrk,
    ifelse(is.finite(P_mrk) && P_mrk < .001, "<0.001", sprintf("= %.3f", P_mrk))
  )
  
  pcv   <- ord_sub$values$Relative_eig
  x_lab <- if (!is.null(pcv) && length(pcv) >= 1 && is.finite(pcv[1])) {
    sprintf("PCoA 1 (%.1f%%)", 100 * pcv[1])   
  } else {
    "PCoA 1"
  }
  y_lab <- if (!is.null(pcv) && length(pcv) >= 2 && is.finite(pcv[2])) {
    sprintf("PCoA 2 (%.1f%%)", 100 * pcv[2])
  } else {
    "PCoA 2"
  }
  
  shp <- c("No" = 16, "Yes" = 17)
  
  p <- ggplot2::ggplot(df_ord, ggplot2::aes(Axis.1, Axis.2)) +
    ggplot2::geom_point(
      ggplot2::aes(color = marker_val, shape = LaterDiarrhea),
      size  = point_size,
      alpha = point_alpha,
      na.rm = TRUE
    ) +
    ggplot2::scale_shape_manual(
      values = shp,
      name   = NULL
    ) +
    ggplot2::scale_color_viridis_c(
      option   = "inferno",  
      direction = -1,    
      name     = colorbar_name,
      na.value = "grey85"
    ) +
    ggplot2::labs(
      x     = x_lab,
      y     = y_lab,
      title = sprintf("Bray-Curtis PCoA: %s", marker_label)
    ) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = zoom_mult)) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = zoom_mult)) +
    ggplot2::coord_equal(expand = TRUE, clip = "off") +
    ggplot2::theme_classic(base_size = 13) +  
    ggplot2::theme(
      legend.position   = "right",
      panel.grid        = ggplot2::element_blank(),
      plot.margin       = margin(6, 10, 6, 6)
    ) +
    ggplot2::annotate(
      "text",
      x     = Inf,
      y     = -Inf,
      hjust = 1.02,
      vjust = -0.6,
      label = ann,
      size  = 3.6  
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

res_ca_dca <- plot_overlay_like_previous(
  ord_sub,
  bray_sub,
  metadata,
  marker_df = df_ratio %>%
    dplyr::select(Sample.ID, CA_DCA_log10) %>%
    dplyr::filter(is.finite(CA_DCA_log10)),
  marker_col    = "CA_DCA_log10",
  marker_label  = "CA/DCA ratio",
  colorbar_name = "log10(CA/DCA)",
  cat_cols      = cat_cols
)

print(res_ca_dca$plot)
print(res_ca_dca$permanova)

res_cdca_lca <- plot_overlay_like_previous(
  ord_sub,
  bray_sub,
  metadata,
  marker_df = df_ratio %>%
    dplyr::select(Sample.ID, CDCA_LCA_log10) %>%
    dplyr::filter(is.finite(CDCA_LCA_log10)),
  marker_col    = "CDCA_LCA_log10",
  marker_label  = "CDCA/LCA ratio",
  colorbar_name = "log10(CDCA/LCA)",
  cat_cols      = cat_cols
)

print(res_cdca_lca$plot)
print(res_cdca_lca$permanova)

stats_tbl <- dplyr::bind_rows(
  res_ca_dca$stats,
  res_cdca_lca$stats
) |>
  dplyr::mutate(
    p_fmt = dplyr::if_else(
      !is.finite(p),
      NA_character_,
      dplyr::if_else(p < .001, "<0.001", sprintf("%.3f", p))
    )
  ) |>
  dplyr::select(marker, n, F, R2, p = p_fmt)

print(stats_tbl)

# ===== SUPPLEMENTAL FIGURE 4D =====
# BA-ratio PCoA overlays, on-treatment (CA/DCA, CDCA/LCA)
if (!exists("bray_all") || !exists("ord_all")) {
  bray_all <- phyloseq::distance(phy, method = "bray")
  ord_all  <- phyloseq::ordinate(phy, method = "PCoA", distance = bray_all)
}

find_ratio_col <- function(df, candidates) {
  hit <- intersect(candidates, names(df))
  if (!length(hit)) stop("None of the candidate columns found: ", paste(candidates, collapse = ", "))
  hit[1]
}

plot_onlen_locked_axes_with_ratio <- function(ord_all,
                                              bray_all,
                                              metadata,
                                              df_ratio,
                                              ratio_candidates,
                                              marker_label,
                                              colorbar_name = NULL,
                                              point_size    = 4.5,
                                              point_alpha   = 0.85,
                                              zoom_mult     = 0.14,   
                                              nperm         = 999) {
  onlen_levels <- c("OnLen_control", "OnLen_diarrhea")
  
  meta_onlen <- metadata %>%
    dplyr::transmute(
      Sample.ID     = as.character(Sample.ID),
      SampleCategory = factor(
        SampleCategory,
        levels = c("Baseline", onlen_levels)
      ),
      SubjectID = factor(SubjectID)
    ) %>%
    dplyr::filter(SampleCategory %in% onlen_levels)
  
  ratio_col <- find_ratio_col(df_ratio, ratio_candidates)
  
  keep_ids <- Reduce(
    intersect,
    list(
      rownames(ord_all$vectors),
      rownames(as.matrix(bray_all)),
      meta_onlen$Sample.ID,
      df_ratio$Sample.ID
    )
  )
  if (length(keep_ids) < 6L) stop("Too few overlapping OnLen samples.")
  
  df_ord <- ord_all$vectors %>%
    as.data.frame() %>%
    tibble::rownames_to_column("Sample.ID") %>%
    dplyr::select(Sample.ID, Axis.1, Axis.2) %>%
    dplyr::filter(Sample.ID %in% keep_ids) %>%
    dplyr::left_join(meta_onlen, by = "Sample.ID") %>%
    dplyr::left_join(
      df_ratio %>%
        dplyr::transmute(Sample.ID, marker_val = .data[[ratio_col]]),
      by = "Sample.ID"
    ) %>%
    dplyr::filter(is.finite(marker_val)) %>%
    droplevels()
  
  Dmat  <- as.matrix(bray_all)
  D_sub <- stats::as.dist(Dmat[df_ord$Sample.ID, df_ord$Sample.ID, drop = FALSE])
  
  ctrl <- if (any(table(df_ord$SubjectID) > 1)) {
    h <- permute::how(nperm = nperm)
    permute::setBlocks(h) <- df_ord$SubjectID
    h
  } else {
    nperm
  }
  
  fit <- vegan::adonis2(
    D_sub ~ marker_val + SampleCategory,
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
  
  F_mrk <- as.numeric(tab["marker_val", f_col])
  R2_mrk <- as.numeric(tab["marker_val", r2_col])
  P_mrk <- as.numeric(tab["marker_val", "Pr(>F)"])
  
  ann <- sprintf(
    "%s PERMANOVA  F=%.2f, p %s",   
    marker_label,
    F_mrk,
    ifelse(is.finite(P_mrk) && P_mrk < .001, "<0.001", sprintf("= %.3f", P_mrk))
  )
  
  pc_var <- ord_all$values$Relative_eig
  x_lab <- if (!is.null(pc_var) && length(pc_var) >= 1) {
    sprintf("PCoA 1 (%.1f%%)", 100 * pc_var[1])
  } else {
    "PCoA 1"
  }
  y_lab <- if (!is.null(pc_var) && length(pc_var) >= 2) {
    sprintf("PCoA 2 (%.1f%%)", 100 * pc_var[2])
  } else {
    "PCoA 2"
  }
  
  shp <- c(OnLen_control = 16, OnLen_diarrhea = 17)
  lab <- c(
    OnLen_control  = "Len (control)",
    OnLen_diarrhea = "Len (diarrhea)"
  )
  
  if (is.null(colorbar_name)) {
    colorbar_name <- marker_label
  }
  
  p <- ggplot2::ggplot(df_ord, ggplot2::aes(Axis.1, Axis.2)) +
    ggplot2::geom_point(
      ggplot2::aes(color = marker_val, shape = SampleCategory),
      size  = point_size,
      alpha = point_alpha,
      na.rm = TRUE
    ) +
    ggplot2::scale_shape_manual(
      values = shp,
      labels = lab,
      name   = NULL
    ) +
    ggplot2::scale_color_viridis_c(
      option   = "inferno",   
      direction = -1,         
      name     = colorbar_name,
      na.value = "grey85",
    ) +
    ggplot2::labs(
      x     = x_lab,
      y     = y_lab,
      title = sprintf("Bray–Curtis PCoA (OnLen) — %s", marker_label)
    ) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = zoom_mult)) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = zoom_mult)) +
    ggplot2::coord_equal(expand = TRUE, clip = "off") +
    ggplot2::theme_classic(base_size = 13) +
    ggplot2::theme(
      legend.position  = "right",
      panel.grid       = ggplot2::element_blank(),
      plot.margin      = margin(6, 10, 6, 6)
    ) +
    ggplot2::annotate(
      "text",
      x     = Inf,
      y     = -Inf,
      hjust = 1.02,
      vjust = -0.6,
      label = ann,
      size  = 3.6
    )
  
  list(
    plot  = p,
    permanova = fit,
    stats = tibble::tibble(
      marker = marker_label,
      n      = nrow(df_ord),
      F      = F_mrk,
      R2     = R2_mrk,
      p      = P_mrk
    ),
    data  = df_ord
  )
}

res_cadca <- plot_onlen_locked_axes_with_ratio(
  ord_all,
  bray_all,
  metadata,
  df_ratio,
  ratio_candidates = c("CA_DCA_log10", "log10_CA_DCA", "CACDCA_log10"),
  marker_label     = "CA/DCA",
  colorbar_name    = "log10(CA/DCA)"
)

res_cdcalca <- plot_onlen_locked_axes_with_ratio(
  ord_all,
  bray_all,
  metadata,
  df_ratio,
  ratio_candidates = c("CDCA_LCA_log10", "log10_CDCA_LCA"),
  marker_label     = "CDCA/LCA",
  colorbar_name    = "log10(CDCA/LCA)"
)

if (!is.null(res_cadca)) {
  print(res_cadca$plot)
  cat("\n--- PERMANOVA: CA/DCA ---\n")
  print(res_cadca$permanova)
}

if (!is.null(res_cdcalca)) {
  print(res_cdcalca$plot)
  cat("\n--- PERMANOVA: CDCA/LCA ---\n")
  print(res_cdcalca$permanova)
}

stats_tbl <- dplyr::bind_rows(
  res_cadca$stats,
  res_cdcalca$stats
) |>
  dplyr::mutate(
    p_fmt = dplyr::if_else(
      !is.finite(p),
      NA_character_,
      dplyr::if_else(p < .001, "<0.001", sprintf("%.3f", p))
    )
  ) |>
  dplyr::select(marker, n, F, R2, p = p_fmt)

print(stats_tbl)
