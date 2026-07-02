#####
# Figure 1, Supplemental Figure 1
#####

if (!exists("df_taxa_genus")) source("R/00_setup.R")

# ===== FIGURE 1A =====
# Swimmer plot 

# ===== FIGURE 1B =====
# Schematic of bile acid metabolism

# ===== FIGURE 1C =====
# Total bile acids
single_totalBA_plot <- function(font_base = 14,
                                bracket_top_pad = 0.10,
                                bracket_gap     = 0.06,
                                y_limits = c(7,10)     ,
                                y_breaks = waiver(),
                                y_label  = "log10(AUC)",
                                axis_text_size = 12,
                                axis_title_size = 12) {

  dat <- bile_long %>%
    dplyr::filter(!is.na(SubjectID)) %>%
    dplyr::group_by(SubjectID, Sample.ID,
                    SampleCat = factor(SampleCat, levels = c("Baseline","OnLen_control","OnLen_diarrhea"))) %>%
    dplyr::summarise(AUC = sum(Amount, na.rm = TRUE), .groups = "drop") %>%
    dplyr::mutate(log_AUC = log10(AUC + LOG10_OFFSET))

  m0 <- lme4::lmer(log_AUC ~ SampleCat + (1|SubjectID), data = dat)
  p_global <- car::Anova(m0, type = 3, test.statistic = "F")["SampleCat","Pr(>F)"]

  group_summ <- dat %>%
    dplyr::group_by(SampleCat) %>%
    dplyr::summarise(mean = mean(log_AUC),
                     se   = stats::sd(log_AUC)/sqrt(dplyr::n()),
                     .groups = "drop") %>%
    dplyr::mutate(ymin = mean - 1.96*se, ymax = mean + 1.96*se)

  emm <- emmeans::emmeans(m0, "SampleCat")
  pw  <- emmeans::contrast(emm, method = "pairwise", adjust = "tukey") %>%
    as.data.frame() %>%
    tidyr::separate(contrast, c("group1","group2"), " - ") %>%
    dplyr::mutate(label = cut(p.value, c(-Inf,.001,.01,.05,Inf),
                              c("***","**","*","ns"), right = FALSE))

  if (is.null(y_limits)) {
    yr <- diff(range(dat$log_AUC, na.rm = TRUE))
    anchor <- max(group_summ$ymax, na.rm = TRUE) + bracket_top_pad * yr
    y_pos_dir <- +1
    ylim_up <- NA_real_
  } else {
    yr <- y_limits[2] - y_limits[1]
    anchor <- y_limits[2] - bracket_top_pad * yr
    y_pos_dir <- -1
    ylim_up <- y_limits[2]
  }

  ord <- tibble::tibble(
    group1 = c("Baseline","OnLen_control","Baseline"),
    group2 = c("OnLen_control","OnLen_diarrhea","OnLen_diarrhea")
  ) %>%
    dplyr::left_join(pw, by = c("group1","group2")) %>%
    dplyr::mutate(y.position = anchor + y_pos_dir * bracket_gap * yr * (dplyr::row_number()-1))

  if (is.na(ylim_up)) {
    ylim_up <- max(ord$y.position, na.rm = TRUE) + 0.06 * yr
  }

  ggplot2::ggplot() +
    ggplot2::geom_line(data = dat, ggplot2::aes(SampleCat, log_AUC, group = SubjectID),
                       colour = "grey70", linewidth = .3, alpha = .4) +
    gghalves::geom_half_violin(
      data = dat, ggplot2::aes(SampleCat, log_AUC, fill = SampleCat),
      side = "l", width = .9, alpha = .5, trim = FALSE
    ) +
    ggplot2::geom_jitter(data = dat, ggplot2::aes(SampleCat, log_AUC, colour = SampleCat),
                         width = .1, size = 1.4, alpha = .60) +
    ggplot2::geom_ribbon(data = group_summ, ggplot2::aes(SampleCat, ymin = ymin, ymax = ymax, group = 1),
                         inherit.aes = FALSE, fill = "#377eb8", alpha = .20) +
    ggplot2::geom_line(data = group_summ, ggplot2::aes(SampleCat, mean, group = 1),
                       colour = "#377eb8", linewidth = 1.2) +
    ggplot2::scale_x_discrete(
      labels = c(Baseline = "Baseline",
                 OnLen_control = "Len (control)",
                 OnLen_diarrhea = "Len (diarrhea)")
    ) +
    ggplot2::scale_y_continuous(limits = y_limits, breaks = y_breaks,
                                expand = expansion(mult = c(0.02, if (is.null(y_limits)) 0.08 else 0))) +
    ggplot2::scale_fill_brewer(palette = "Dark2") +
    ggplot2::scale_colour_brewer(palette = "Dark2") +
    ggplot2::labs(x = NULL, y = y_label, title = "Total Bile Acids") +
    ggplot2::theme_classic(base_size = font_base) +
    ggplot2::theme(
      legend.position = "none",
      axis.text.x  = ggplot2::element_text(size = axis_text_size, angle = 20, hjust = 1),
      axis.text.y  = ggplot2::element_text(size = axis_text_size),
      axis.title.y = ggplot2::element_text(size = axis_title_size, margin = ggplot2::margin(r = 8))
    ) +
    ggpubr::stat_pvalue_manual(
      ord, label = "label", hide.ns = FALSE,
      tip.length = .01, bracket.size = .6, step.increase = 0,
      size = font_base * 0.3
    ) +
    ggplot2::annotate(
      "text", x = Inf, y = ylim_up, hjust = 1.02, vjust = 1.2,
      label = sprintf("Overall p=%s",
                      ifelse(is.na(p_global), "NA",
                             ifelse(p_global < .001, "<0.001", sprintf("%.3f", p_global)))),
      fontface = "italic", size = axis_title_size * 0.35
    ) +
    ggplot2::coord_cartesian(ylim = c(NA, ylim_up), clip = "off")
}

p_total <- single_totalBA_plot()
print(p_total)

# ===== FIGURE 1D + SUPPLEMENTAL FIGURE 1C =====
# Bile-acid heatmap (1D); mean composition by class/conjugation (Supp 1C)
p_thresh  <- 0.05
col_order <- c("Len (control) vs. Baseline",
               "Len (diarrhea) vs. Baseline",
               "Len (diarrhea) vs. Len (control)")

dat_ba <- bile_long %>%
  filter(!is.na(SubjectID)) %>%
  mutate(SampleCat = factor(SampleCat, levels = c("Baseline","OnLen_control","OnLen_diarrhea"))) %>%
  transmute(SubjectID, SampleCat, bilenames, log10_AUC = log10(Amount + LOG10_OFFSET)) %>%
  group_by(SubjectID, SampleCat, bilenames) %>%
  summarise(log10_AUC = mean(log10_AUC, na.rm = TRUE), .groups = "drop")

pairs_tbl <- tribble(
  ~group1,         ~group2,            ~comparison,
  "Baseline",      "OnLen_control",    "Len (control) vs. Baseline",
  "Baseline",      "OnLen_diarrhea",   "Len (diarrhea) vs. Baseline",
  "OnLen_control", "OnLen_diarrhea",   "Len (diarrhea) vs. Len (control)"
)

get_wilcox <- function(acid) {
  d <- filter(dat_ba, bilenames == acid)
  if (n_distinct(d$SampleCat) < 2L) return(NULL)

  pmap_dfr(pairs_tbl, function(group1, group2, comparison) {
    x1 <- d$log10_AUC[d$SampleCat == group1]
    x2 <- d$log10_AUC[d$SampleCat == group2]
    if (!length(x1) || !length(x2)) return(NULL)

    tibble(
      bilenames  = acid,
      comparison = comparison,
      logFC10    = mean(x2, na.rm = TRUE) - mean(x1, na.rm = TRUE),
      p.value    = stats::wilcox.test(x2, x1, exact = FALSE)$p.value
    )
  })
}

df_all <- map_dfr(unique(dat_ba$bilenames), get_wilcox) %>%
  group_by(comparison) %>%
  mutate(q.value = p.adjust(p.value, method = "BH")) %>%
  ungroup()

df_wide  <- df_all %>%
  select(bilenames, comparison, logFC10) %>%
  pivot_wider(names_from = comparison, values_from = logFC10) %>%
  select(bilenames, all_of(col_order))

mat_logFC <- df_wide %>% select(-bilenames) %>% as.matrix()
rownames(mat_logFC) <- df_wide$bilenames

df_wide_q <- df_all %>%
  select(bilenames, comparison, q.value) %>%
  pivot_wider(names_from = comparison, values_from = q.value) %>%
  select(bilenames, all_of(col_order))

qval_mat <- df_wide_q %>% select(-bilenames) %>% as.matrix()
rownames(qval_mat) <- df_wide_q$bilenames

row_anno <- bile_long %>%
  distinct(bilenames, Classification, Conjugation) %>%
  mutate(
    Classification = fct_relevel(Classification, "primary/host", "secondary/microbe"),
    Conjugation   = fct_relevel(Conjugation, "unconjugated", "glycine conjugated",
                                "taurine conjugated", "sulfated")
  )

row_order <- row_anno %>%
  arrange(Classification, Conjugation, bilenames) %>%
  pull(bilenames) %>%
  intersect(rownames(mat_logFC))

anno_colors <- list(
  Classification = c("primary/host" = "darkgoldenrod1",
                     "secondary/microbe" = "#8BBF50"),
  Conjugation   = c("unconjugated"       = "#87CEFA",
                    "glycine conjugated" = "#9370DB",
                    "taurine conjugated" = "#FFB6C1",
                    "sulfated"           = "#A9A9A9")
)

rng <- na.omit(as.numeric(mat_logFC))
lim <- max(quantile(abs(rng), 0.95, names = FALSE), 0.5)

hm_long <- as.data.frame(mat_logFC) %>%
  tibble::rownames_to_column("bilenames") %>%
  tidyr::pivot_longer(-bilenames, names_to = "comparison", values_to = "logFC10") %>%
  dplyr::left_join(
    as.data.frame(qval_mat) %>%
      tibble::rownames_to_column("bilenames") %>%
      tidyr::pivot_longer(-bilenames, names_to = "comparison", values_to = "q.value"),
    by = c("bilenames", "comparison")
  ) %>%
  dplyr::filter(bilenames %in% row_order) %>%
  dplyr::mutate(
    bilenames  = factor(bilenames,  levels = rev(row_order)),
    comparison = factor(comparison, levels = col_order),
    star       = ifelse(!is.na(q.value) & q.value < p_thresh, "*", "")
  )

anno_long <- row_anno %>%
  dplyr::filter(bilenames %in% row_order) %>%
  dplyr::mutate(bilenames = factor(bilenames, levels = rev(row_order))) %>%
  tidyr::pivot_longer(c(Classification, Conjugation),
                      names_to = "anno", values_to = "value") %>%
  dplyr::mutate(anno = factor(anno, levels = c("Classification", "Conjugation")))

p_anno <- ggplot(anno_long, aes(anno, bilenames, fill = value)) +
  geom_tile(color = "white") +
  scale_fill_manual(values = c(anno_colors$Classification, anno_colors$Conjugation),
                    name = NULL) +
  scale_x_discrete(position = "top") +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 10) +
  theme(axis.text.x = element_text(angle = 40, hjust = 0),
        panel.grid  = element_blank(),
        legend.position = "right")

p_hmap <- ggplot(hm_long, aes(comparison, bilenames, fill = logFC10)) +
  geom_tile(color = "grey92") +
  geom_text(aes(label = star), vjust = 0.78, size = 5) +
  scale_fill_gradient2(low = "#3B6FB6", mid = "white", high = "#D54E4E",
                       midpoint = 0, limits = c(-lim, lim),
                       oob = scales::squish, name = "log10 FC") +
  scale_x_discrete(position = "top") +
  labs(x = NULL, y = NULL, caption = sprintf("* FDR q < %.2f", p_thresh)) +
  theme_minimal(base_size = 10) +
  theme(axis.text.x = element_text(angle = 40, hjust = 0),
        axis.text.y = element_blank(),
        panel.grid  = element_blank(),
        legend.position = "right")

p_ht <- p_anno + p_hmap + patchwork::plot_layout(widths = c(1, 2.2))
print(p_ht)

pal_class <- c(
  "primary/host"      = "darkgoldenrod1",
  "secondary/microbe" = "#8BBF50"
)

lab_map <- c(
  Baseline       = "Baseline",
  OnLen_control  = "Len (control)",
  OnLen_diarrhea = "Len (diarrhea)"
)

comp_class <- bile_long %>%
  dplyr::group_by(Sample.ID) %>%
  dplyr::mutate(total = sum(Amount, na.rm = TRUE)) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(Sample.ID, SampleCat, Classification) %>%
  dplyr::summarise(
    pct = sum(Amount, na.rm = TRUE) / dplyr::first(total),
    .groups = "drop"
  ) %>%
  dplyr::group_by(SampleCat, Classification) %>%
  dplyr::summarise(mean_pct = mean(pct, na.rm = TRUE), .groups = "drop")

comp_class_lab <- comp_class %>%
  dplyr::mutate(
    SampleCat = factor(SampleCat, levels = c("Baseline","OnLen_control","OnLen_diarrhea")),
    label     = ifelse(mean_pct < 0.03, "", scales::percent(mean_pct, accuracy = 1))
  )

comp_class_lab <- comp_class %>%
  mutate(
    SampleCat = factor(SampleCat, levels = c("Baseline","OnLen_control","OnLen_diarrhea")),
    label = ifelse(mean_pct < 0.03, "", percent(mean_pct, accuracy = 1))
  )

p_pie <- ggplot(comp_class_lab, aes(x = 1, y = mean_pct, fill = Classification)) +
  geom_col(width = 1, color = "white", size = 0.4) +
  coord_polar(theta = "y", clip = "off") +
  facet_wrap(~ SampleCat, nrow = 1, labeller = as_labeller(lab_map)) +
  geom_text(aes(label = label),
            position = position_stack(vjust = 0.5),
            color = "black", size = 4.8) +
  scale_fill_manual(
    values = pal_class,
    breaks = c("primary/host","secondary/microbe"),
    name   = "Classification"
  ) +
  labs(x = NULL, y = NULL, title = "Primary/Host vs. Secondary/Microbe-Derived Bile Acids") +
  theme_void(base_size = 15) +
  theme(
    strip.text      = element_text(face = "bold", size = 14),
    plot.title      = element_text(face = "bold", size = 16, margin = margin(b = 8)),
    legend.position = "bottom",
    legend.direction= "horizontal",
    legend.title    = element_text(size = 12),
    legend.text     = element_text(size = 11),
    legend.box      = "vertical",
    plot.margin     = margin(t = 6, r = 8, b = 6, l = 8)
  ) +
  guides(fill = guide_legend(nrow = 1, byrow = TRUE))

print(p_pie)

conj_colors <- c(
  "unconjugated"       = "#87CEFA",
  "glycine conjugated" = "#9370DB",
  "taurine conjugated" = "#FFB6C1",
  "sulfated"           = "#A9A9A9"
)

lab_map <- c(
  Baseline       = "Baseline",
  OnLen_control  = "Len (control)",
  OnLen_diarrhea = "Len (diarrhea)"
)

comp_conj <- bile_long %>%
  group_by(Sample.ID) %>%
  mutate(total = sum(Amount, na.rm = TRUE)) %>%
  ungroup() %>%
  group_by(Sample.ID, SampleCat, Conjugation) %>%
  summarise(pct = sum(Amount, na.rm = TRUE) / first(total), .groups = "drop") %>%
  group_by(SampleCat, Conjugation) %>%
  summarise(mean_pct = mean(pct, na.rm = TRUE), .groups = "drop")

comp_conj_lab <- comp_conj %>%
  mutate(
    SampleCat   = factor(SampleCat, levels = c("Baseline","OnLen_control","OnLen_diarrhea")),
    Conjugation = factor(Conjugation, levels = names(conj_colors)),
    label       = ifelse(mean_pct < 0.03, "", percent(mean_pct, accuracy = 1))
  )

p_pie_conj <- ggplot(comp_conj_lab, aes(x = 1, y = mean_pct, fill = Conjugation)) +
  geom_col(width = 1, color = "white", size = 0.4) +
  coord_polar(theta = "y", clip = "off") +
  facet_wrap(~ SampleCat, nrow = 1, labeller = as_labeller(lab_map)) +
  geom_text(aes(label = label),
            position = position_stack(vjust = 0.5),
            color = "black", size = 4.8) +
  scale_fill_manual(
    values = conj_colors,
    breaks = names(conj_colors),
    name   = "Conjugation"
  ) +
  labs(x = NULL, y = NULL, title = "Bile acid composition by Conjugation") +
  theme_void(base_size = 15) +
  theme(
    strip.text      = element_text(face = "bold", size = 14),
    plot.title      = element_text(face = "bold", size = 16, margin = margin(b = 8)),
    legend.position = "bottom",
    legend.direction= "horizontal",
    legend.title    = element_text(size = 12),
    legend.text     = element_text(size = 11),
    legend.box      = "vertical",
    plot.margin     = margin(t = 6, r = 8, b = 6, l = 8)
  ) +
  guides(fill = guide_legend(nrow = 1, byrow = TRUE))

print(p_pie_conj)

# ===== FIGURE 1E =====
# Cholate (CA)
single_acid_plot_CA <- function(font_base = 14,
                                bracket_top_pad = 0.10,
                                bracket_gap     = 0.06,
                                y_limits = c(3, 10.5),
                                y_breaks = waiver(),
                                y_label  = "log10(AUC)",
                                axis_text_size = 12,
                                axis_title_size = 12) {

  nm <- c("cholate (CA)", "cholic acid")

  dat <- bile_long %>%
    dplyr::filter(bilenames %in% nm, !is.na(SubjectID)) %>%
    dplyr::group_by(SubjectID, Sample.ID,
                    SampleCat = factor(SampleCat, levels = c("Baseline","OnLen_control","OnLen_diarrhea"))) %>%
    dplyr::summarise(AUC = sum(Amount, na.rm = TRUE), .groups = "drop") %>%
    dplyr::mutate(log_AUC = log10(AUC + LOG10_OFFSET))

  m0 <- lme4::lmer(log_AUC ~ SampleCat + (1|SubjectID), data = dat)
  p_global <- car::Anova(m0, type = 3, test.statistic = "F")["SampleCat","Pr(>F)"]

  group_summ <- dat %>%
    dplyr::group_by(SampleCat) %>%
    dplyr::summarise(mean = mean(log_AUC), se = stats::sd(log_AUC)/sqrt(dplyr::n()), .groups = "drop") %>%
    dplyr::mutate(ymin = mean - 1.96*se, ymax = mean + 1.96*se)

  emm <- emmeans::emmeans(m0, "SampleCat")
  pw  <- emmeans::contrast(emm, method = "pairwise", adjust = "tukey") %>%
    as.data.frame() %>%
    tidyr::separate(contrast, c("group1","group2"), " - ") %>%
    dplyr::mutate(label = cut(p.value, c(-Inf,.001,.01,.05,Inf),
                              c("***","**","*","ns"), right = FALSE))

  if (is.null(y_limits)) {
    yr <- diff(range(dat$log_AUC, na.rm = TRUE))
    anchor <- max(group_summ$ymax, na.rm = TRUE) + bracket_top_pad * yr
    y_pos_dir <- +1
    ylim_up <- NA_real_
  } else {
    yr <- y_limits[2] - y_limits[1]
    anchor <- y_limits[2] - bracket_top_pad * yr
    y_pos_dir <- -1
    ylim_up <- y_limits[2]
  }

  ord <- tibble::tibble(
    group1 = c("Baseline","OnLen_control","Baseline"),
    group2 = c("OnLen_control","OnLen_diarrhea","OnLen_diarrhea")
  ) %>%
    dplyr::left_join(pw, by = c("group1","group2")) %>%
    dplyr::mutate(y.position = anchor + y_pos_dir * bracket_gap * yr * (dplyr::row_number()-1))

  if (is.na(ylim_up)) {
    ylim_up <- max(ord$y.position, na.rm = TRUE) + 0.06 * yr
  }

  ggplot2::ggplot() +
    ggplot2::geom_line(data = dat, ggplot2::aes(SampleCat, log_AUC, group = SubjectID),
                       colour = "grey70", linewidth = .3, alpha = .4) +
    gghalves::geom_half_violin(data = dat, ggplot2::aes(SampleCat, log_AUC, fill = SampleCat),
                               side = "l", width = .9, alpha = .5, trim = FALSE) +
    ggplot2::geom_jitter(data = dat, ggplot2::aes(SampleCat, log_AUC, colour = SampleCat),
                         width = .1, size = 1.4, alpha = .60) +
    ggplot2::geom_ribbon(data = group_summ, ggplot2::aes(SampleCat, ymin = ymin, ymax = ymax, group = 1),
                         inherit.aes = FALSE, fill = "#377eb8", alpha = .20) +
    ggplot2::geom_line(data = group_summ, ggplot2::aes(SampleCat, mean, group = 1),
                       colour = "#377eb8", linewidth = 1.2) +
    ggplot2::scale_x_discrete(
      labels = c(Baseline = "Baseline",
                 OnLen_control = "Len (control)",
                 OnLen_diarrhea = "Len (diarrhea)")
    ) +
    ggplot2::scale_y_continuous(limits = y_limits, breaks = y_breaks,
                                expand = expansion(mult = c(0.02, if (is.null(y_limits)) 0.08 else 0))) +
    ggplot2::scale_fill_brewer(palette = "Dark2") +
    ggplot2::scale_colour_brewer(palette = "Dark2") +
    ggplot2::labs(x = NULL, y = y_label, title = "Cholate (CA)") +
    ggplot2::theme_classic(base_size = font_base) +
    ggplot2::theme(
      legend.position = "none",
      axis.text.x  = ggplot2::element_text(size = axis_text_size, angle = 20, hjust = 1),
      axis.text.y  = ggplot2::element_text(size = axis_text_size),
      axis.title.y = ggplot2::element_text(size = axis_title_size, margin = ggplot2::margin(r = 8))
    ) +
    ggpubr::stat_pvalue_manual(
      ord, label = "label", hide.ns = FALSE,
      tip.length = .01, bracket.size = .6, step.increase = 0,
      size = font_base * 0.3
    ) +
    ggplot2::annotate(
      "text", x = Inf, y = ylim_up, hjust = 1.02, vjust = 1.2,
      label = sprintf("Overall p%s",
                      ifelse(is.na(p_global), "NA",
                             ifelse(p_global < .001, "<0.001", sprintf("%.3f", p_global)))),
      fontface = "italic", size = axis_title_size * 0.35
    ) +
    ggplot2::coord_cartesian(ylim = c(NA, ylim_up), clip = "off")
}

p_CA <- single_acid_plot_CA()

print(p_CA)

# ===== FIGURE 1E =====
# Chenodeoxycholate (CDCA)
single_acid_plot_CDCA <- function(font_base = 14,
                                bracket_top_pad = 0.10,
                                bracket_gap     = 0.06,
                                y_limits = c(3, 10.5),
                                y_breaks = waiver(),
                                y_label  = "log10 (AUC)",
                                axis_text_size = 12,
                                axis_title_size = 12) {

  nm <- c("chenodeoxycholate (CDCA)", "chenodeoxycholic acid")

  dat <- bile_long %>%
    dplyr::filter(bilenames %in% nm, !is.na(SubjectID)) %>%
    dplyr::group_by(SubjectID, Sample.ID,
                    SampleCat = factor(SampleCat, levels = c("Baseline","OnLen_control","OnLen_diarrhea"))) %>%
    dplyr::summarise(AUC = sum(Amount, na.rm = TRUE), .groups = "drop") %>%
    dplyr::mutate(log_AUC = log10(AUC + LOG10_OFFSET))

  m0 <- lme4::lmer(log_AUC ~ SampleCat + (1|SubjectID), data = dat)
  p_global <- car::Anova(m0, type = 3, test.statistic = "F")["SampleCat","Pr(>F)"]

  group_summ <- dat %>%
    dplyr::group_by(SampleCat) %>%
    dplyr::summarise(mean = mean(log_AUC), se = stats::sd(log_AUC)/sqrt(dplyr::n()), .groups = "drop") %>%
    dplyr::mutate(ymin = mean - 1.96*se, ymax = mean + 1.96*se)

  emm <- emmeans::emmeans(m0, "SampleCat")
  pw  <- emmeans::contrast(emm, method = "pairwise", adjust = "tukey") %>%
    as.data.frame() %>%
    tidyr::separate(contrast, c("group1","group2"), " - ") %>%
    dplyr::mutate(label = cut(p.value, c(-Inf,.001,.01,.05,Inf),
                              c("***","**","*","ns"), right = FALSE))

  if (is.null(y_limits)) {
    yr <- diff(range(dat$log_AUC, na.rm = TRUE))
    anchor <- max(group_summ$ymax, na.rm = TRUE) + bracket_top_pad * yr
    y_pos_dir <- +1
    ylim_up <- NA_real_
  } else {
    yr <- y_limits[2] - y_limits[1]
    anchor <- y_limits[2] - bracket_top_pad * yr
    y_pos_dir <- -1
    ylim_up <- y_limits[2]
  }

  ord <- tibble::tibble(
    group1 = c("Baseline","OnLen_control","Baseline"),
    group2 = c("OnLen_control","OnLen_diarrhea","OnLen_diarrhea")
  ) %>%
    dplyr::left_join(pw, by = c("group1","group2")) %>%
    dplyr::mutate(y.position = anchor + y_pos_dir * bracket_gap * yr * (dplyr::row_number()-1))

  if (is.na(ylim_up)) {
    ylim_up <- max(ord$y.position, na.rm = TRUE) + 0.06 * yr
  }

  ggplot2::ggplot() +
    ggplot2::geom_line(data = dat, ggplot2::aes(SampleCat, log_AUC, group = SubjectID),
                       colour = "grey70", linewidth = .3, alpha = .4) +
    gghalves::geom_half_violin(data = dat, ggplot2::aes(SampleCat, log_AUC, fill = SampleCat),
                               side = "l", width = .9, alpha = .5, trim = FALSE) +
    ggplot2::geom_jitter(data = dat, ggplot2::aes(SampleCat, log_AUC, colour = SampleCat),
                         width = .1, size = 1.4, alpha = .60) +
    ggplot2::geom_ribbon(data = group_summ, ggplot2::aes(SampleCat, ymin = ymin, ymax = ymax, group = 1),
                         inherit.aes = FALSE, fill = "#377eb8", alpha = .20) +
    ggplot2::geom_line(data = group_summ, ggplot2::aes(SampleCat, mean, group = 1),
                       colour = "#377eb8", linewidth = 1.2) +
    ggplot2::scale_x_discrete(
      labels = c(Baseline = "Baseline",
                 OnLen_control = "Len (control)",
                 OnLen_diarrhea = "Len (diarrhea)")
    ) +
    ggplot2::scale_y_continuous(limits = y_limits, breaks = y_breaks,
                                expand = expansion(mult = c(0.02, if (is.null(y_limits)) 0.08 else 0))) +
    ggplot2::scale_fill_brewer(palette = "Dark2") +
    ggplot2::scale_colour_brewer(palette = "Dark2") +
    ggplot2::labs(x = NULL, y = y_label, title = "Chenodeoxycholate (CDCA)") +
    ggplot2::theme_classic(base_size = font_base) +
    ggplot2::theme(
      legend.position = "none",
      axis.text.x  = ggplot2::element_text(size = axis_text_size, angle = 20, hjust = 1),
      axis.text.y  = ggplot2::element_text(size = axis_text_size),
      axis.title.y = ggplot2::element_text(size = axis_title_size, margin = ggplot2::margin(r = 8))
    ) +
    ggpubr::stat_pvalue_manual(
      ord, label = "label", hide.ns = FALSE,
      tip.length = .01, bracket.size = .6, step.increase = 0,
      size = font_base * 0.3
    ) +
    ggplot2::annotate(
      "text", x = Inf, y = ylim_up, hjust = 1.02, vjust = 1.2,
      label = sprintf("Overall p=%s",
                      ifelse(is.na(p_global), "NA",
                             ifelse(p_global < .001, "<0.001", sprintf("%.3f", p_global)))),
      fontface = "italic", size = axis_title_size * 0.35
    ) +
    ggplot2::coord_cartesian(ylim = c(NA, ylim_up), clip = "off")
}

p_CDCA <- single_acid_plot_CDCA()
print(p_CDCA)

# ===== FIGURE 1E =====
# Deoxycholate (DCA)
single_acid_plot_DCA <- function(font_base = 14,
                                  bracket_top_pad = 0.10,
                                  bracket_gap     = 0.06,
                                  y_limits = c(3, 10.5),
                                  y_breaks = waiver(),
                                  y_label  = "log10 (AUC)",
                                  axis_text_size = 12,
                                  axis_title_size = 12) {

  nm <- c("deoxycholate (DCA)", "deoxycholic acid")

  dat <- bile_long %>%
    dplyr::filter(bilenames %in% nm, !is.na(SubjectID)) %>%
    dplyr::group_by(SubjectID, Sample.ID,
                    SampleCat = factor(SampleCat, levels = c("Baseline","OnLen_control","OnLen_diarrhea"))) %>%
    dplyr::summarise(AUC = sum(Amount, na.rm = TRUE), .groups = "drop") %>%
    dplyr::mutate(log_AUC = log10(AUC + LOG10_OFFSET))

  m0 <- lme4::lmer(log_AUC ~ SampleCat + (1|SubjectID), data = dat)
  p_global <- car::Anova(m0, type = 3, test.statistic = "F")["SampleCat","Pr(>F)"]

  group_summ <- dat %>%
    dplyr::group_by(SampleCat) %>%
    dplyr::summarise(mean = mean(log_AUC), se = stats::sd(log_AUC)/sqrt(dplyr::n()), .groups = "drop") %>%
    dplyr::mutate(ymin = mean - 1.96*se, ymax = mean + 1.96*se)

  emm <- emmeans::emmeans(m0, "SampleCat")
  pw  <- emmeans::contrast(emm, method = "pairwise", adjust = "tukey") %>%
    as.data.frame() %>%
    tidyr::separate(contrast, c("group1","group2"), " - ") %>%
    dplyr::mutate(label = cut(p.value, c(-Inf,.001,.01,.05,Inf),
                              c("***","**","*","ns"), right = FALSE))

  if (is.null(y_limits)) {
    yr <- diff(range(dat$log_AUC, na.rm = TRUE))
    anchor <- max(group_summ$ymax, na.rm = TRUE) + bracket_top_pad * yr
    y_pos_dir <- +1
    ylim_up <- NA_real_
  } else {
    yr <- y_limits[2] - y_limits[1]
    anchor <- y_limits[2] - bracket_top_pad * yr
    y_pos_dir <- -1
    ylim_up <- y_limits[2]
  }

  ord <- tibble::tibble(
    group1 = c("Baseline","OnLen_control","Baseline"),
    group2 = c("OnLen_control","OnLen_diarrhea","OnLen_diarrhea")
  ) %>%
    dplyr::left_join(pw, by = c("group1","group2")) %>%
    dplyr::mutate(y.position = anchor + y_pos_dir * bracket_gap * yr * (dplyr::row_number()-1))

  if (is.na(ylim_up)) {
    ylim_up <- max(ord$y.position, na.rm = TRUE) + 0.06 * yr
  }

  ggplot2::ggplot() +
    ggplot2::geom_line(data = dat, ggplot2::aes(SampleCat, log_AUC, group = SubjectID),
                       colour = "grey70", linewidth = .3, alpha = .4) +
    gghalves::geom_half_violin(data = dat, ggplot2::aes(SampleCat, log_AUC, fill = SampleCat),
                               side = "l", width = .9, alpha = .5, trim = FALSE) +
    ggplot2::geom_jitter(data = dat, ggplot2::aes(SampleCat, log_AUC, colour = SampleCat),
                         width = .1, size = 1.4, alpha = .60) +
    ggplot2::geom_ribbon(data = group_summ, ggplot2::aes(SampleCat, ymin = ymin, ymax = ymax, group = 1),
                         inherit.aes = FALSE, fill = "#377eb8", alpha = .20) +
    ggplot2::geom_line(data = group_summ, ggplot2::aes(SampleCat, mean, group = 1),
                       colour = "#377eb8", linewidth = 1.2) +
    ggplot2::scale_x_discrete(
      labels = c(Baseline = "Baseline",
                 OnLen_control = "Len (control)",
                 OnLen_diarrhea = "Len (diarrhea)")
    ) +
    ggplot2::scale_y_continuous(limits = y_limits, breaks = y_breaks,
                                expand = expansion(mult = c(0.02, if (is.null(y_limits)) 0.08 else 0))) +
    ggplot2::scale_fill_brewer(palette = "Dark2") +
    ggplot2::scale_colour_brewer(palette = "Dark2") +
    ggplot2::labs(x = NULL, y = y_label, title = "Deoxycholate (DCA)") +
    ggplot2::theme_classic(base_size = font_base) +
    ggplot2::theme(
      legend.position = "none",
      axis.text.x  = ggplot2::element_text(size = axis_text_size, angle = 20, hjust = 1),
      axis.text.y  = ggplot2::element_text(size = axis_text_size),
      axis.title.y = ggplot2::element_text(size = axis_title_size, margin = ggplot2::margin(r = 8))
    ) +
    ggpubr::stat_pvalue_manual(
      ord, label = "label", hide.ns = FALSE,
      tip.length = .01, bracket.size = .6, step.increase = 0,
      size = font_base * 0.3
    ) +
    ggplot2::annotate(
      "text", x = Inf, y = ylim_up, hjust = 1.02, vjust = 1.2,
      label = sprintf("Overall p=%s",
                      ifelse(is.na(p_global), "NA",
                             ifelse(p_global < .001, "<0.001", sprintf("%.3f", p_global)))),
      fontface = "italic", size = axis_title_size * 0.35
    ) +
    ggplot2::coord_cartesian(ylim = c(NA, ylim_up), clip = "off")
}

p_DCA <- single_acid_plot_DCA()
print(p_DCA)

# ===== FIGURE 1E =====
# Lithocholate (LCA)
single_acid_plot_LCA <- function(font_base = 14,
                                 bracket_top_pad = 0.10,
                                 bracket_gap     = 0.06,
                                 y_limits = c(3, 10.5),
                                 y_breaks = waiver(),
                                 y_label  = "log10 (AUC)",
                                 axis_text_size = 12,
                                 axis_title_size = 12) {

  nm <- c("lithocholate (LCA)", "lithocholic acid")

  dat <- bile_long %>%
    dplyr::filter(bilenames %in% nm, !is.na(SubjectID)) %>%
    dplyr::group_by(SubjectID, Sample.ID,
                    SampleCat = factor(SampleCat, levels = c("Baseline","OnLen_control","OnLen_diarrhea"))) %>%
    dplyr::summarise(AUC = sum(Amount, na.rm = TRUE), .groups = "drop") %>%
    dplyr::mutate(log_AUC = log10(AUC + LOG10_OFFSET))

  m0 <- lme4::lmer(log_AUC ~ SampleCat + (1|SubjectID), data = dat)
  p_global <- car::Anova(m0, type = 3, test.statistic = "F")["SampleCat","Pr(>F)"]

  group_summ <- dat %>%
    dplyr::group_by(SampleCat) %>%
    dplyr::summarise(mean = mean(log_AUC), se = stats::sd(log_AUC)/sqrt(dplyr::n()), .groups = "drop") %>%
    dplyr::mutate(ymin = mean - 1.96*se, ymax = mean + 1.96*se)

  emm <- emmeans::emmeans(m0, "SampleCat")
  pw  <- emmeans::contrast(emm, method = "pairwise", adjust = "tukey") %>%
    as.data.frame() %>%
    tidyr::separate(contrast, c("group1","group2"), " - ") %>%
    dplyr::mutate(label = cut(p.value, c(-Inf,.001,.01,.05,Inf),
                              c("***","**","*","ns"), right = FALSE))

  if (is.null(y_limits)) {
    yr <- diff(range(dat$log_AUC, na.rm = TRUE))
    anchor <- max(group_summ$ymax, na.rm = TRUE) + bracket_top_pad * yr
    y_pos_dir <- +1
    ylim_up <- NA_real_
  } else {
    yr <- y_limits[2] - y_limits[1]
    anchor <- y_limits[2] - bracket_top_pad * yr
    y_pos_dir <- -1
    ylim_up <- y_limits[2]
  }

  ord <- tibble::tibble(
    group1 = c("Baseline","OnLen_control","Baseline"),
    group2 = c("OnLen_control","OnLen_diarrhea","OnLen_diarrhea")
  ) %>%
    dplyr::left_join(pw, by = c("group1","group2")) %>%
    dplyr::mutate(y.position = anchor + y_pos_dir * bracket_gap * yr * (dplyr::row_number()-1))

  if (is.na(ylim_up)) {
    ylim_up <- max(ord$y.position, na.rm = TRUE) + 0.06 * yr
  }

  ggplot2::ggplot() +
    ggplot2::geom_line(data = dat, ggplot2::aes(SampleCat, log_AUC, group = SubjectID),
                       colour = "grey70", linewidth = .3, alpha = .4) +
    gghalves::geom_half_violin(data = dat, ggplot2::aes(SampleCat, log_AUC, fill = SampleCat),
                               side = "l", width = .9, alpha = .5, trim = FALSE) +
    ggplot2::geom_jitter(data = dat, ggplot2::aes(SampleCat, log_AUC, colour = SampleCat),
                         width = .1, size = 1.4, alpha = .60) +
    ggplot2::geom_ribbon(data = group_summ, ggplot2::aes(SampleCat, ymin = ymin, ymax = ymax, group = 1),
                         inherit.aes = FALSE, fill = "#377eb8", alpha = .20) +
    ggplot2::geom_line(data = group_summ, ggplot2::aes(SampleCat, mean, group = 1),
                       colour = "#377eb8", linewidth = 1.2) +
    ggplot2::scale_x_discrete(
      labels = c(Baseline = "Baseline",
                 OnLen_control = "Len (control)",
                 OnLen_diarrhea = "Len (diarrhea)")
    ) +
    ggplot2::scale_y_continuous(limits = y_limits, breaks = y_breaks,
                                expand = expansion(mult = c(0.02, if (is.null(y_limits)) 0.08 else 0))) +
    ggplot2::scale_fill_brewer(palette = "Dark2") +
    ggplot2::scale_colour_brewer(palette = "Dark2") +
    ggplot2::labs(x = NULL, y = y_label, title = "Lithocholate (LCA)") +
    ggplot2::theme_classic(base_size = font_base) +
    ggplot2::theme(
      legend.position = "none",
      axis.text.x  = ggplot2::element_text(size = axis_text_size, angle = 20, hjust = 1),
      axis.text.y  = ggplot2::element_text(size = axis_text_size),
      axis.title.y = ggplot2::element_text(size = axis_title_size, margin = ggplot2::margin(r = 8))
    ) +
    ggpubr::stat_pvalue_manual(
      ord, label = "label", hide.ns = FALSE,
      tip.length = .01, bracket.size = .6, step.increase = 0,
      size = font_base * 0.3
    ) +
    ggplot2::annotate(
      "text", x = Inf, y = ylim_up, hjust = 1.02, vjust = 1.2,
      label = sprintf("Overall p=%s",
                      ifelse(is.na(p_global), "NA",
                             ifelse(p_global < .001, "<0.001", sprintf("%.3f", p_global)))),
      fontface = "italic", size = axis_title_size * 0.35
    ) +
    ggplot2::coord_cartesian(ylim = c(NA, ylim_up), clip = "off")
}

p_LCA <- single_acid_plot_LCA()
print(p_LCA)

# ===== FIGURE 1E =====
# CA/DCA ratio
single_ratio_plot_CA_DCA <- function(
font_base = 14,
bracket_top_pad = 0.10,
bracket_gap     = 0.06,
y_limits = c(-5,7),
y_breaks = waiver(),
y_label  = "log10 (CA/DCA)",
axis_text_size = 12,
axis_title_size = 12,
backtransform_ticks = FALSE
) {
  ca_names  <- c("cholate (CA)", "cholic acid")
  dca_names <- c("deoxycholate (DCA)", "deoxycholic acid")

  dat <- bile_long %>%
    dplyr::filter(bilenames %in% c(ca_names, dca_names), !is.na(SubjectID)) %>%
    dplyr::mutate(acid = dplyr::case_when(
      bilenames %in% ca_names  ~ "CA",
      bilenames %in% dca_names ~ "DCA"
    )) %>%
    dplyr::group_by(SubjectID, Sample.ID,
                    SampleCat = factor(SampleCat, levels = c("Baseline","OnLen_control","OnLen_diarrhea")),
                    acid) %>%
    dplyr::summarise(AUC = sum(Amount, na.rm = TRUE), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = acid, values_from = AUC, values_fill = 0) %>%
    dplyr::mutate(
      log_ratio = log10(CA  + LOG10_OFFSET) - log10(DCA + LOG10_OFFSET)
    )

  m0 <- lme4::lmer(log_ratio ~ SampleCat + (1|SubjectID), data = dat)
  p_global <- car::Anova(m0, type = 3, test.statistic = "F")["SampleCat","Pr(>F)"]

  group_summ <- dat %>%
    dplyr::group_by(SampleCat) %>%
    dplyr::summarise(mean = mean(log_ratio),
                     se   = stats::sd(log_ratio)/sqrt(dplyr::n()),
                     .groups = "drop") %>%
    dplyr::mutate(ymin = mean - 1.96*se, ymax = mean + 1.96*se)

  emm <- emmeans::emmeans(m0, "SampleCat")
  pw  <- emmeans::contrast(emm, method = "pairwise", adjust = "tukey") %>%
    as.data.frame() %>%
    tidyr::separate(contrast, c("group1","group2"), " - ") %>%
    dplyr::mutate(label = cut(p.value, c(-Inf,.001,.01,.05,Inf),
                              c("***","**","*","ns"), right = FALSE))

  if (is.null(y_limits)) {
    yr <- diff(range(dat$log_ratio, na.rm = TRUE))
    anchor <- max(group_summ$ymax, na.rm = TRUE) + bracket_top_pad * yr
    y_dir <- +1
    ylim_up <- NA_real_
  } else {
    yr <- y_limits[2] - y_limits[1]
    anchor <- y_limits[2] - bracket_top_pad * yr
    y_dir <- -1
    ylim_up <- y_limits[2]
  }

  ord <- tibble::tibble(
    group1 = c("Baseline","OnLen_control","Baseline"),
    group2 = c("OnLen_control","OnLen_diarrhea","OnLen_diarrhea")
  ) %>%
    dplyr::left_join(pw, by = c("group1","group2")) %>%
    dplyr::mutate(y.position = anchor + y_dir * bracket_gap * yr * (dplyr::row_number()-1))

  if (is.na(ylim_up)) {
    ylim_up <- max(ord$y.position, na.rm = TRUE) + 0.06 * yr
  }

  y_scale <- if (isTRUE(backtransform_ticks)) {
    ggplot2::scale_y_continuous(
      limits = y_limits,
      breaks = if (is.null(y_limits)) y_breaks else seq(ceiling(y_limits[1]), floor(y_limits[2]), by = 1),
      labels = scales::math_format(10^.x),
      expand = ggplot2::expansion(mult = c(0.02, if (is.null(y_limits)) 0.08 else 0))
    )
  } else {
    ggplot2::scale_y_continuous(
      limits = y_limits,
      breaks = y_breaks,
      expand = ggplot2::expansion(mult = c(0.02, if (is.null(y_limits)) 0.08 else 0))
    )
  }

  ggplot2::ggplot() +
    ggplot2::geom_line(data = dat, ggplot2::aes(SampleCat, log_ratio, group = SubjectID),
                       colour = "grey70", linewidth = .3, alpha = .4) +
    gghalves::geom_half_violin(
      data = dat, ggplot2::aes(SampleCat, log_ratio, fill = SampleCat),
      side = "l", width = .9, alpha = .5, trim = FALSE
    ) +
    ggplot2::geom_jitter(data = dat, ggplot2::aes(SampleCat, log_ratio, colour = SampleCat),
                         width = .1, size = 1.4, alpha = .60) +
    ggplot2::geom_ribbon(data = group_summ,
                         ggplot2::aes(SampleCat, ymin = ymin, ymax = ymax, group = 1),
                         inherit.aes = FALSE, fill = "#377eb8", alpha = .20) +
    ggplot2::geom_line(data = group_summ, ggplot2::aes(SampleCat, mean, group = 1),
                       colour = "#377eb8", linewidth = 1.2) +
    ggplot2::scale_x_discrete(
      labels = c(Baseline = "Baseline",
                 OnLen_control = "Len (control)",
                 OnLen_diarrhea = "Len (diarrhea)")
    ) +
    y_scale +
    ggplot2::scale_fill_brewer(palette = "Dark2") +
    ggplot2::scale_colour_brewer(palette = "Dark2") +
    ggplot2::labs(
      x = NULL,
      y = if (isTRUE(backtransform_ticks)) "CA/DCA" else y_label,
      title = "CA/DCA ratio"
    ) +
    ggplot2::theme_classic(base_size = font_base) +
    ggplot2::theme(
      legend.position = "none",
      axis.text.x  = ggplot2::element_text(size = axis_text_size, angle = 20, hjust = 1),
      axis.text.y  = ggplot2::element_text(size = axis_text_size),
      axis.title.y = ggplot2::element_text(size = axis_title_size, margin = ggplot2::margin(r = 8))
    ) +
    ggpubr::stat_pvalue_manual(
      ord, label = "label", hide.ns = FALSE,
      tip.length = .01, bracket.size = .6, step.increase = 0,
      size = font_base * 0.3
    ) +
    ggplot2::annotate(
      "text", x = Inf, y = ylim_up, hjust = 1.02, vjust = 1.2,
      label = sprintf("Overall p=%s",
                      ifelse(is.na(p_global), "NA",
                             ifelse(p_global < .001, "<0.001", sprintf("%.3f", p_global)))),
      fontface = "italic", size = axis_title_size * 0.35
    ) +
    ggplot2::coord_cartesian(ylim = c(NA, ylim_up), clip = "off")
}

p_CA_DCA <- single_ratio_plot_CA_DCA()

print(p_CA_DCA)

# ===== FIGURE 1E =====
# CDCA/LCA ratio
    single_ratio_plot_CDCA_LCA <- function(
    font_base = 14,
    bracket_top_pad = 0.10,
    bracket_gap     = 0.06,
    y_limits = c(-4, 7),
    y_breaks = waiver(),
    y_label  = "log10 (CDCA/LCA)",
    axis_text_size = 12,
    axis_title_size = 12,
    backtransform_ticks = FALSE
        ) {
          cdca_names <- c("chenodeoxycholate (CDCA)", "chenodeoxycholic acid")
          lca_names  <- c("lithocholate (LCA)", "lithocholic acid")

          dat <- bile_long %>%
            dplyr::filter(bilenames %in% c(cdca_names, lca_names), !is.na(SubjectID)) %>%
            dplyr::mutate(acid = dplyr::case_when(
              bilenames %in% cdca_names ~ "CDCA",
              bilenames %in% lca_names  ~ "LCA"
            )) %>%
            dplyr::group_by(SubjectID, Sample.ID,
                            SampleCat = factor(SampleCat, levels = c("Baseline","OnLen_control","OnLen_diarrhea")),
                            acid) %>%
            dplyr::summarise(AUC = sum(Amount, na.rm = TRUE), .groups = "drop") %>%
            tidyr::pivot_wider(names_from = acid, values_from = AUC, values_fill = 0) %>%
            dplyr::mutate(
              log_ratio = log10(CDCA + LOG10_OFFSET) - log10(LCA + LOG10_OFFSET)
            )

          m0 <- lme4::lmer(log_ratio ~ SampleCat + (1|SubjectID), data = dat)
          p_global <- car::Anova(m0, type = 3, test.statistic = "F")["SampleCat","Pr(>F)"]

          group_summ <- dat %>%
            dplyr::group_by(SampleCat) %>%
            dplyr::summarise(mean = mean(log_ratio),
                             se   = stats::sd(log_ratio)/sqrt(dplyr::n()),
                             .groups = "drop") %>%
            dplyr::mutate(ymin = mean - 1.96*se, ymax = mean + 1.96*se)

          emm <- emmeans::emmeans(m0, "SampleCat")
          pw  <- emmeans::contrast(emm, method = "pairwise", adjust = "tukey") %>%
            as.data.frame() %>%
            tidyr::separate(contrast, c("group1","group2"), " - ") %>%
            dplyr::mutate(label = cut(p.value, c(-Inf,.001,.01,.05,Inf),
                                      c("***","**","*","ns"), right = FALSE))

          if (is.null(y_limits)) {
            yr <- diff(range(dat$log_ratio, na.rm = TRUE))
            anchor <- max(group_summ$ymax, na.rm = TRUE) + bracket_top_pad * yr
            y_dir <- +1
            ylim_up <- NA_real_
          } else {
            yr <- y_limits[2] - y_limits[1]
            anchor <- y_limits[2] - bracket_top_pad * yr
            y_dir <- -1
            ylim_up <- y_limits[2]
          }

          ord <- tibble::tibble(
            group1 = c("Baseline","OnLen_control","Baseline"),
            group2 = c("OnLen_control","OnLen_diarrhea","OnLen_diarrhea")
          ) %>%
            dplyr::left_join(pw, by = c("group1","group2")) %>%
            dplyr::mutate(y.position = anchor + y_dir * bracket_gap * yr * (dplyr::row_number()-1))

          if (is.na(ylim_up)) {
            ylim_up <- max(ord$y.position, na.rm = TRUE) + 0.06 * yr
          }

          y_scale <- if (isTRUE(backtransform_ticks)) {
            ggplot2::scale_y_continuous(
              limits = y_limits,
              breaks = if (is.null(y_limits)) y_breaks else seq(ceiling(y_limits[1]), floor(y_limits[2]), by = 1),
              labels = scales::math_format(10^.x),
              expand = ggplot2::expansion(mult = c(0.02, if (is.null(y_limits)) 0.08 else 0))
            )
          } else {
            ggplot2::scale_y_continuous(
              limits = y_limits,
              breaks = y_breaks,
              expand = ggplot2::expansion(mult = c(0.02, if (is.null(y_limits)) 0.08 else 0))
            )
          }

          ggplot2::ggplot() +
            ggplot2::geom_line(data = dat, ggplot2::aes(SampleCat, log_ratio, group = SubjectID),
                               colour = "grey70", linewidth = .3, alpha = .4) +
            gghalves::geom_half_violin(
              data = dat, ggplot2::aes(SampleCat, log_ratio, fill = SampleCat),
              side = "l", width = .9, alpha = .5, trim = FALSE
            ) +
            ggplot2::geom_jitter(data = dat, ggplot2::aes(SampleCat, log_ratio, colour = SampleCat),
                                 width = .1, size = 1.4, alpha = .60) +
            ggplot2::geom_ribbon(data = group_summ,
                                 ggplot2::aes(SampleCat, ymin = ymin, ymax = ymax, group = 1),
                                 inherit.aes = FALSE, fill = "#377eb8", alpha = .20) +
            ggplot2::geom_line(data = group_summ, ggplot2::aes(SampleCat, mean, group = 1),
                               colour = "#377eb8", linewidth = 1.2) +
            ggplot2::scale_x_discrete(
              labels = c(Baseline = "Baseline",
                         OnLen_control = "Len (control)",
                         OnLen_diarrhea = "Len (diarrhea)")
            ) +
            y_scale +
            ggplot2::scale_fill_brewer(palette = "Dark2") +
            ggplot2::scale_colour_brewer(palette = "Dark2") +
            ggplot2::labs(
              x = NULL,
              y = if (isTRUE(backtransform_ticks)) "CDCA/LCA" else y_label,
              title = "CDCA/LCA ratio"
            ) +
            ggplot2::theme_classic(base_size = font_base) +
            ggplot2::theme(
              legend.position = "none",
              axis.text.x  = ggplot2::element_text(size = axis_text_size, angle = 20, hjust = 1),
              axis.text.y  = ggplot2::element_text(size = axis_text_size),
              axis.title.y = ggplot2::element_text(size = axis_title_size, margin = ggplot2::margin(r = 8))
            ) +
            ggpubr::stat_pvalue_manual(
              ord, label = "label", hide.ns = FALSE,
              tip.length = .01, bracket.size = .6, step.increase = 0,
              size = font_base * 0.3
            ) +
            ggplot2::annotate(
              "text", x = Inf, y = ylim_up, hjust = 1.02, vjust = 1.2,
              label = sprintf("Overall p=%s",
                              ifelse(is.na(p_global), "NA",
                                     ifelse(p_global < .001, "<0.001", sprintf("%.3f", p_global)))),
              fontface = "italic", size = axis_title_size * 0.35
            ) +
            ggplot2::coord_cartesian(ylim = c(NA, ylim_up), clip = "off")
        }

p_CDCA_LCA <- single_ratio_plot_CDCA_LCA()

print(p_CDCA_LCA)

# ===== FIGURE 1E =====
# Glycocholate (GCA)
single_acid_plot_GCA <- function(font_base = 14,
                                 bracket_top_pad = 0.10,
                                 bracket_gap     = 0.06,
                                 y_limits = c(3, 10.5),
                                 y_breaks = waiver(),
                                 y_label  = "log10 (AUC)",
                                 axis_text_size = 12,
                                 axis_title_size = 12) {

  nm <- c("glycocholate (GCA)", "glycocholic acid")

  dat <- bile_long %>%
    dplyr::filter(bilenames %in% nm, !is.na(SubjectID)) %>%
    dplyr::group_by(SubjectID, Sample.ID,
                    SampleCat = factor(SampleCat, levels = c("Baseline","OnLen_control","OnLen_diarrhea"))) %>%
    dplyr::summarise(AUC = sum(Amount, na.rm = TRUE), .groups = "drop") %>%
    dplyr::mutate(log_AUC = log10(AUC + LOG10_OFFSET))

  m0 <- lme4::lmer(log_AUC ~ SampleCat + (1|SubjectID), data = dat)
  p_global <- car::Anova(m0, type = 3, test.statistic = "F")["SampleCat","Pr(>F)"]

  group_summ <- dat %>%
    dplyr::group_by(SampleCat) %>%
    dplyr::summarise(mean = mean(log_AUC),
                     se   = stats::sd(log_AUC)/sqrt(dplyr::n()),
                     .groups = "drop") %>%
    dplyr::mutate(ymin = mean - 1.96*se, ymax = mean + 1.96*se)

  emm <- emmeans::emmeans(m0, "SampleCat")
  pw  <- emmeans::contrast(emm, method = "pairwise", adjust = "tukey") %>%
    as.data.frame() %>%
    tidyr::separate(contrast, c("group1","group2"), " - ") %>%
    dplyr::mutate(label = cut(p.value, c(-Inf,.001,.01,.05,Inf),
                              c("***","**","*","ns"), right = FALSE))

  if (is.null(y_limits)) {
    yr <- diff(range(dat$log_AUC, na.rm = TRUE))
    anchor <- max(group_summ$ymax, na.rm = TRUE) + bracket_top_pad * yr
    y_pos_dir <- +1
    ylim_up <- NA_real_
  } else {
    yr <- y_limits[2] - y_limits[1]
    anchor <- y_limits[2] - bracket_top_pad * yr
    y_pos_dir <- -1
    ylim_up <- y_limits[2]
  }

  ord <- tibble::tibble(
    group1 = c("Baseline","OnLen_control","Baseline"),
    group2 = c("OnLen_control","OnLen_diarrhea","OnLen_diarrhea")
  ) %>%
    dplyr::left_join(pw, by = c("group1","group2")) %>%
    dplyr::mutate(y.position = anchor + y_pos_dir * bracket_gap * yr * (dplyr::row_number()-1))

  if (is.na(ylim_up)) {
    ylim_up <- max(ord$y.position, na.rm = TRUE) + 0.06 * yr
  }

  ggplot2::ggplot() +
    ggplot2::geom_line(data = dat, ggplot2::aes(SampleCat, log_AUC, group = SubjectID),
                       colour = "grey70", linewidth = .3, alpha = .4) +
    gghalves::geom_half_violin(
      data = dat, ggplot2::aes(SampleCat, log_AUC, fill = SampleCat),
      side = "l", width = .9, alpha = .5, trim = FALSE
    ) +
    ggplot2::geom_jitter(data = dat, ggplot2::aes(SampleCat, log_AUC, colour = SampleCat),
                         width = .1, size = 1.4, alpha = .60) +
    ggplot2::geom_ribbon(data = group_summ,
                         ggplot2::aes(SampleCat, ymin = ymin, ymax = ymax, group = 1),
                         inherit.aes = FALSE, fill = "#377eb8", alpha = .20) +
    ggplot2::geom_line(data = group_summ, ggplot2::aes(SampleCat, mean, group = 1),
                       colour = "#377eb8", linewidth = 1.2) +
    ggplot2::scale_x_discrete(
      labels = c(Baseline = "Baseline",
                 OnLen_control = "Len (control)",
                 OnLen_diarrhea = "Len (diarrhea)")
    ) +
    ggplot2::scale_y_continuous(limits = y_limits, breaks = y_breaks,
                                expand = ggplot2::expansion(mult = c(0.02, if (is.null(y_limits)) 0.08 else 0))) +
    ggplot2::scale_fill_brewer(palette = "Dark2") +
    ggplot2::scale_colour_brewer(palette = "Dark2") +
    ggplot2::labs(x = NULL, y = y_label, title = "Glycocholate (GCA)") +
    ggplot2::theme_classic(base_size = font_base) +
    ggplot2::theme(
      legend.position = "none",
      axis.text.x  = ggplot2::element_text(size = axis_text_size, angle = 20, hjust = 1),
      axis.text.y  = ggplot2::element_text(size = axis_text_size),
      axis.title.y = ggplot2::element_text(size = axis_title_size, margin = ggplot2::margin(r = 8))
    ) +
    ggpubr::stat_pvalue_manual(
      ord, label = "label", hide.ns = FALSE,
      tip.length = .01, bracket.size = .6, step.increase = 0,
      size = font_base * 0.3
    ) +
    ggplot2::annotate(
      "text", x = Inf, y = ylim_up, hjust = 1.02, vjust = 1.2,
      label = sprintf("Overall p=%s",
                      ifelse(is.na(p_global), "NA",
                             ifelse(p_global < .001, "<0.001", sprintf("%.3f", p_global)))),
      fontface = "italic", size = axis_title_size * 0.35
    ) +
    ggplot2::coord_cartesian(ylim = c(NA, ylim_up), clip = "off")
}

p_GCA <- single_acid_plot_GCA()
print(p_GCA)

# ===== FIGURE 1E =====
# Taurocholate (TCA)
single_acid_plot_TCA <- function(font_base = 14,
                                 bracket_top_pad = 0.10,
                                 bracket_gap     = 0.06,
                                 y_limits = c(2.5, 11.5),
                                 y_breaks = waiver(),
                                 y_label  = "log10 (AUC)",
                                 axis_text_size = 12,
                                 axis_title_size = 12) {

  nm <- c("taurocholate (TCA)", "taurocholic acid")

  dat <- bile_long %>%
    dplyr::filter(bilenames %in% nm, !is.na(SubjectID)) %>%
    dplyr::group_by(SubjectID, Sample.ID,
                    SampleCat = factor(SampleCat, levels = c("Baseline","OnLen_control","OnLen_diarrhea"))) %>%
    dplyr::summarise(AUC = sum(Amount, na.rm = TRUE), .groups = "drop") %>%
    dplyr::mutate(log_AUC = log10(AUC + LOG10_OFFSET))

  m0 <- lme4::lmer(log_AUC ~ SampleCat + (1|SubjectID), data = dat)
  p_global <- car::Anova(m0, type = 3, test.statistic = "F")["SampleCat","Pr(>F)"]

  group_summ <- dat %>%
    dplyr::group_by(SampleCat) %>%
    dplyr::summarise(mean = mean(log_AUC),
                     se   = stats::sd(log_AUC)/sqrt(dplyr::n()),
                     .groups = "drop") %>%
    dplyr::mutate(ymin = mean - 1.96*se, ymax = mean + 1.96*se)

  emm <- emmeans::emmeans(m0, "SampleCat")
  pw  <- emmeans::contrast(emm, method = "pairwise", adjust = "tukey") %>%
    as.data.frame() %>%
    tidyr::separate(contrast, c("group1","group2"), " - ") %>%
    dplyr::mutate(label = cut(p.value, c(-Inf,.001,.01,.05,Inf),
                              c("***","**","*","ns"), right = FALSE))

  if (is.null(y_limits)) {
    yr <- diff(range(dat$log_AUC, na.rm = TRUE))
    anchor <- max(group_summ$ymax, na.rm = TRUE) + bracket_top_pad * yr
    y_pos_dir <- +1
    ylim_up <- NA_real_
  } else {
    yr <- y_limits[2] - y_limits[1]
    anchor <- y_limits[2] - bracket_top_pad * yr
    y_pos_dir <- -1
    ylim_up <- y_limits[2]
  }

  ord <- tibble::tibble(
    group1 = c("Baseline","OnLen_control","Baseline"),
    group2 = c("OnLen_control","OnLen_diarrhea","OnLen_diarrhea")
  ) %>%
    dplyr::left_join(pw, by = c("group1","group2")) %>%
    dplyr::mutate(y.position = anchor + y_pos_dir * bracket_gap * yr * (dplyr::row_number()-1))

  if (is.na(ylim_up)) {
    ylim_up <- max(ord$y.position, na.rm = TRUE) + 0.06 * yr
  }

  ggplot2::ggplot() +
    ggplot2::geom_line(data = dat, ggplot2::aes(SampleCat, log_AUC, group = SubjectID),
                       colour = "grey70", linewidth = .3, alpha = .4) +
    gghalves::geom_half_violin(
      data = dat, ggplot2::aes(SampleCat, log_AUC, fill = SampleCat),
      side = "l", width = .9, alpha = .5, trim = FALSE
    ) +
    ggplot2::geom_jitter(data = dat, ggplot2::aes(SampleCat, log_AUC, colour = SampleCat),
                         width = .1, size = 1.4, alpha = .60) +
    ggplot2::geom_ribbon(data = group_summ,
                         ggplot2::aes(SampleCat, ymin = ymin, ymax = ymax, group = 1),
                         inherit.aes = FALSE, fill = "#377eb8", alpha = .20) +
    ggplot2::geom_line(data = group_summ, ggplot2::aes(SampleCat, mean, group = 1),
                       colour = "#377eb8", linewidth = 1.2) +
    ggplot2::scale_x_discrete(
      labels = c(Baseline = "Baseline",
                 OnLen_control = "Len (control)",
                 OnLen_diarrhea = "Len (diarrhea)")
    ) +
    ggplot2::scale_y_continuous(limits = y_limits, breaks = y_breaks,
                                expand = ggplot2::expansion(mult = c(0.02, if (is.null(y_limits)) 0.08 else 0))) +
    ggplot2::scale_fill_brewer(palette = "Dark2") +
    ggplot2::scale_colour_brewer(palette = "Dark2") +
    ggplot2::labs(x = NULL, y = y_label, title = "Taurocholate (TCA)") +
    ggplot2::theme_classic(base_size = font_base) +
    ggplot2::theme(
      legend.position = "none",
      axis.text.x  = ggplot2::element_text(size = axis_text_size, angle = 20, hjust = 1),
      axis.text.y  = ggplot2::element_text(size = axis_text_size),
      axis.title.y = ggplot2::element_text(size = axis_title_size, margin = ggplot2::margin(r = 8))
    ) +
    ggpubr::stat_pvalue_manual(
      ord, label = "label", hide.ns = FALSE,
      tip.length = .01, bracket.size = .6, step.increase = 0,
      size = font_base * 0.3
    ) +
    ggplot2::annotate(
      "text", x = Inf, y = ylim_up, hjust = 1.02, vjust = 1.2,
      label = sprintf("Overall p=%s",
                      ifelse(is.na(p_global), "NA",
                             ifelse(p_global < .001, "<0.001", sprintf("%.3f", p_global)))),
      fontface = "italic", size = axis_title_size * 0.35
    ) +
    ggplot2::coord_cartesian(ylim = c(NA, ylim_up), clip = "off")
}

p_TCA <- single_acid_plot_TCA()
print(p_TCA)

# ===== SUPPLEMENTAL FIGURE 1A =====
# Swimmer plot with QOL surveys

# ===== SUPPLEMENTAL FIGURE 1B =====
# Venn diagram of sample overlap

# ===== SUPPLEMENTAL FIGURE 1D =====
# Bile-acid dot plot

stopifnot(exists("bile_long"))

row_anno <- bile_long %>%
  distinct(bilenames, Classification, Conjugation) %>%
  mutate(
    Classification = fct_relevel(Classification, "primary/host", "secondary/microbe"),
    Conjugation   = fct_relevel(Conjugation, "unconjugated",
                                "glycine conjugated", "taurine conjugated", "sulfated")
  )

row_order <- row_anno %>%
  arrange(Classification, Conjugation, bilenames) %>%
  pull(bilenames)

ba_group_mean <- bile_long %>%
  filter(!is.na(SampleCat)) %>%
  mutate(SampleCat = factor(SampleCat,
                            levels = c("Baseline","OnLen_control","OnLen_diarrhea"))) %>%
  group_by(SubjectID, SampleCat, bilenames) %>%
  summarise(auc_subj = mean(Amount, na.rm = TRUE), .groups = "drop") %>%
  group_by(SampleCat, bilenames) %>%
  summarise(
    n_subj    = n(),                 
    mean_auc = mean(auc_subj, na.rm = TRUE), 
    .groups  = "drop"
  ) %>%
  inner_join(row_anno, by = "bilenames") %>%
  mutate(
    bilenames = factor(bilenames, levels = intersect(row_order, bilenames))
  )

group_order <- c("Baseline","OnLen_control","OnLen_diarrhea")
present_groups <- group_order[group_order %in% levels(ba_group_mean$SampleCat)]
y_map <- c(Baseline = 0.30, OnLen_control = 0.00, OnLen_diarrhea = -0.30)
ba_group_mean <- ba_group_mean %>%
  mutate(y_row = unname(y_map[as.character(SampleCat)]))

cols_class <- c("primary/host" = "darkgoldenrod1",
                "secondary/microbe" = "#8BBF50")
cols_conj  <- c("unconjugated"       = "#87CEFA",
                "glycine conjugated" = "#9370DB",
                "taurine conjugated" = "#FFB6C1",
                "sulfated"           = "#A9A9A9")

strip1_y <- min(y_map[present_groups]) - 0.22  
strip2_y <- min(y_map[present_groups]) - 0.34  

y_labs <- c("Baseline","Len (control)","Len (diarrhea)")
names(y_labs) <- group_order
y_breaks <- y_map[present_groups]
y_labels <- y_labs[present_groups]

p_ba_rows <- ggplot(ba_group_mean, aes(x = bilenames)) +
  geom_hline(yintercept = y_breaks, colour = "grey92", linewidth = 0.4) +
  geom_tile(aes(y = strip1_y, fill = Classification), height = 0.085, width = 0.9) +
  scale_fill_manual(values = cols_class, name = "Classification") +
  ggnewscale::new_scale_fill() +
  geom_tile(aes(y = strip2_y, fill = Conjugation), height = 0.085, width = 0.9) +
  scale_fill_manual(values = cols_conj, name = "Conjugation") +
  geom_point(aes(y = y_row, size = mean_auc),
             shape = 21, fill = "#4C78A8", colour = "black", alpha = 0.9) +
  scale_size_continuous(
    trans  = "sqrt",
    range  = c(1.5, 10),
    name   = "Mean AUC",
    breaks = c(1e4, 1e5, 1e6, 1e7, 4e7),       
    labels = scales::label_scientific() 
  )+
  scale_y_continuous(
    breaks = y_breaks,      
    labels = y_labels,    
    limits = c(strip2_y - 0.06, max(y_breaks) + 0.18),
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  labs(
    x = NULL,
    y = NULL,
    title = "Relative mean abundance of bile acids by sample group"
  ) +
  coord_cartesian(clip = "off") +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x   = element_text(angle = 35, hjust = 1, size=8.5),
    axis.ticks.y  = element_blank(),
    legend.position = "right",
    panel.grid      = element_blank(),
    plot.margin = margin(t = 10, r = 20, b = 10, l = 20, unit = "mm")
  )

print(p_ba_rows)
