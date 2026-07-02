#####
# Figure 2, Supplemental Figure 2
#####

if (!exists("df_taxa_genus")) source("R/00_setup.R")

# ===== FIGURE 2A =====
# α-diversity (Inverse Simpson)

alpha_div_plot_lmer <- function(
    data,
    response         = "invsimpson",
    use_log10        = TRUE,
    font_base        = 14,
    bracket_top_pad  = 0.65,
    bracket_gap      = 0.08,
    y_limits         = NULL,
    y_breaks         = waiver(),
    axis_text_size   = 12,
    axis_title_size  = 12,
    title            = "Alpha-diversity",
    y_label_override = NULL
) {
  stopifnot(all(c("SubjectID","SampleCategory", response) %in% names(data)))
  
  df <- data %>%
    dplyr::filter(!is.na(.data[[response]]), !is.na(SubjectID), !is.na(SampleCategory)) %>%
    dplyr::transmute(
      SubjectID = as.factor(SubjectID),
      SampleCat = factor(SampleCategory, levels = c("Baseline","OnLen_control","OnLen_diarrhea")),
      val = as.numeric(.data[[response]])
    )
  if (isTRUE(use_log10)) {
    df <- df %>% dplyr::mutate(y = log10(val))
    y_lab <- "log10(Inv. Simpson)"
  } else {
    df <- df %>% dplyr::mutate(y = val)
    y_lab <- "Inverse Simpson"
  }
  
  m0 <- lme4::lmer(y ~ SampleCat + (1|SubjectID), data = df)
  p_global <- car::Anova(m0, type = 3, test.statistic = "F")["SampleCat","Pr(>F)"]
  
  group_summ <- df %>%
    dplyr::group_by(SampleCat) %>%
    dplyr::summarise(
      mean = mean(y, na.rm = TRUE),
      se   = stats::sd(y, na.rm = TRUE)/sqrt(dplyr::n()),
      .groups = "drop"
    ) %>%
    dplyr::mutate(ymin = mean - 1.96*se, ymax = mean + 1.96*se)
  
  emm <- emmeans::emmeans(m0, "SampleCat")
  pw  <- emmeans::contrast(emm, method = "pairwise", adjust = "tukey") %>%
    as.data.frame() %>%
    tidyr::separate(contrast, c("group1","group2"), " - ") %>%
    dplyr::mutate(label = cut(p.value, c(-Inf,.001,.01,.05,Inf),
                              c("***","**","*","ns"), right = FALSE))
  
  if (is.null(y_limits)) {
    yr <- diff(range(df$y, na.rm = TRUE))
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
    ggplot2::geom_line(data = df, ggplot2::aes(SampleCat, y, group = SubjectID),
                       colour = "grey70", linewidth = .3, alpha = .4) +
    gghalves::geom_half_violin(data = df, ggplot2::aes(SampleCat, y, fill = SampleCat),
                               side = "l", width = .9, alpha = .5, trim = FALSE) +
    ggplot2::geom_jitter(data = df, ggplot2::aes(SampleCat, y, colour = SampleCat),
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
                                expand = ggplot2::expansion(mult = c(0.02, if (is.null(y_limits)) 0.08 else 0))) +
    ggplot2::scale_fill_brewer(palette = "Dark2") +
    ggplot2::scale_colour_brewer(palette = "Dark2") +
    ggplot2::labs(x = NULL, y = y_lab, title = title) +
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
      size = font_base * 0.30
    ) +
    ggplot2::annotate(
      "text", x = Inf, y = ylim_up, hjust = 1.02, vjust = 0,
      label = sprintf("Overall p %s",
                      ifelse(is.na(p_global), "NA",
                             ifelse(p_global < .001, "<0.001", sprintf("%.3f", p_global)))),
      fontface = "italic", size = axis_title_size * 0.35
    ) +
    ggplot2::coord_cartesian(ylim = c(NA, ylim_up), clip = "off")
}

p_alpha_lmm <- alpha_div_plot_lmer(
  data      = metadata,
  response  = "invsimpson",
  use_log10 = TRUE,
  title     = "Alpha-diversity (Inverse Simpson)"
)
print(p_alpha_lmm)

# ===== FIGURE 2B =====
# Dominated samples and domination composition

THRESH <- 0.30
TOP_N  <- 10

RAW_LEVELS    <- c("Baseline", "OnLen_control", "OnLen_diarrhea")
PRETTY_LEVELS <- c("Baseline", "Len (control)", "Len (diarrhea)")

X_ids <- intersect(rownames(df_taxa_genus), metadata$Sample.ID)

X <- df_taxa_genus[X_ids, , drop = FALSE]

meta_use <- metadata |>
  dplyr::filter(Sample.ID %in% X_ids) |>
  dplyr::mutate(
    Group_raw = factor(SampleCategory, levels = RAW_LEVELS),
    Group     = forcats::fct_recode(
      Group_raw,
      "Baseline"       = "Baseline",
      "Len (control)"  = "OnLen_control",
      "Len (diarrhea)" = "OnLen_diarrhea"
    ) |>
      forcats::fct_relevel(PRETTY_LEVELS)
  )

row_tot <- rowSums(X, na.rm = TRUE)
X_rel   <- sweep(X, 1, ifelse(row_tot > 0, row_tot, 1), "/")

get_dom <- function(v) {
  if (sum(v, na.rm = TRUE) <= 0) {
    return(c(genus = NA_character_, prop = NA_real_))
  }
  j <- which.max(v)
  c(genus = colnames(X_rel)[j], prop = as.numeric(v[j]))
}

dom_mat <- t(apply(X_rel, 1, get_dom)) |>
  as.data.frame()

dom_df <- dom_mat |>
  dplyr::mutate(
    Sample.ID = rownames(X_rel),
    prop      = as.numeric(prop),
    genus     = ifelse(is.na(genus) | genus == "", "(unlabeled)", genus)
  ) |>
  dplyr::left_join(
    meta_use |>
      dplyr::select(Sample.ID, Group),
    by = "Sample.ID"
  ) |>
  dplyr::filter(!is.na(Group)) |>
  dplyr::mutate(dominated = is.finite(prop) & prop >= THRESH)

summary_dom <- dom_df |>
  dplyr::group_by(Group) |>
  dplyr::summarise(
    n_total     = dplyr::n(),
    n_dominated = sum(dominated, na.rm = TRUE),
    prop_dom    = n_dominated / n_total,
    .groups     = "drop"
  ) |>
  dplyr::arrange(forcats::fct_relevel(Group, PRETTY_LEVELS))

pairs_dom <- list(
  c("Baseline", "Len (control)"),
  c("Len (control)", "Len (diarrhea)"),
  c("Baseline", "Len (diarrhea)")
)

pair_test_dom <- function(summ, g1, g2) {
  s1 <- dplyr::filter(summ, Group == g1)
  s2 <- dplyr::filter(summ, Group == g2)
  if (nrow(s1) == 0 || nrow(s2) == 0 || s1$n_total < 1 || s2$n_total < 1) {
    return(tibble::tibble(group1 = g1, group2 = g2, p.value = NA_real_))
  }
  x <- c(s1$n_dominated, s2$n_dominated)
  n <- c(s1$n_total,     s2$n_total)
  extreme <- (sum(x) == 0) || (sum(x == n) == 2)
  
  p <- tryCatch(
    {
      if (extreme) {
        stats::fisher.test(
          matrix(c(x[1], n[1] - x[1], x[2], n[2] - x[2]), 2, byrow = TRUE)
        )$p.value
      } else {
        stats::prop.test(x = x, n = n, correct = TRUE)$p.value
      }
    },
    error = function(e) NA_real_
  )
  
  tibble::tibble(group1 = g1, group2 = g2, p.value = p)
}

pw_dom <- purrr::map_dfr(
  pairs_dom,
  ~pair_test_dom(summary_dom, .x[1], .x[2])
) |>
  dplyr::mutate(
    q.value = p.adjust(p.value, method = "BH"),
    label   = dplyr::case_when(
      is.na(q.value) ~ "NA",
      q.value < .001 ~ "***",
      q.value < .01  ~ "**",
      q.value < .05  ~ "*",
      TRUE           ~ "ns"
    )
  )

y_limits_dom <- c(0, 1)
yr_dom       <- diff(y_limits_dom)

ord_dom <- pw_dom |>
  dplyr::mutate(
    y.position = y_limits_dom[2] - 0.08 * yr_dom - 0.06 * yr_dom * (dplyr::row_number() - 1)
  )

tab_dom <- with(dom_df, table(dominated, Group))
overall_lab_dom <- if (ncol(tab_dom) < 2 || nrow(tab_dom) < 2) {
  "Overall Fisher p NA"
} else {
  ft <- tryCatch(
    stats::fisher.test(tab_dom, simulate.p.value = TRUE, B = 1e5),
    error = function(e) NULL
  )
  if (is.null(ft)) "Overall Fisher p NA"
  else if (ft$p.value < .001) "Overall Fisher p <0.001"
  else sprintf("Overall Fisher p %.3f", ft$p.value)
}

p_dom <- ggplot2::ggplot(summary_dom, ggplot2::aes(x = Group, y = prop_dom, fill = Group)) +
  ggplot2::geom_col(width = 0.7) +
  ggplot2::geom_text(
    ggplot2::aes(label = scales::percent(prop_dom, accuracy = 1)),
    vjust = -0.4,
    size  = 4
  ) +
  ggplot2::scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
    limits = c(0, 1),
    expand = ggplot2::expansion(mult = c(0.02, 0.10))
  ) +
  ggplot2::scale_fill_brewer(palette = "Dark2", guide = "none") +
  ggplot2::labs(
    x     = NULL,
    y     = "Dominated samples (%)",
    title = sprintf("Dominated samples (≥%s%% single genus)", THRESH * 100)
  ) +
  ggplot2::theme_classic(base_size = 13) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 20, hjust = 1)) +
  ggpubr::stat_pvalue_manual(
    ord_dom |>
      dplyr::filter(!is.na(label)),
    label         = "label",
    hide.ns       = FALSE,
    tip.length    = .01,
    bracket.size  = .6,
    step.increase = 0,
    size          = 4
  ) +
  ggplot2::annotate(
    "text",
    x        = Inf,
    y        = Inf,
    hjust    = 1.02,
    vjust    = 1.2,
    label    = overall_lab_dom,
    fontface = "italic",
    size     = 4
  ) +
  ggplot2::coord_cartesian(clip = "off")

print(p_dom)

top_genera <- dom_df |>
  dplyr::filter(dominated) |>
  dplyr::count(genus, sort = TRUE) |>
  dplyr::slice_head(n = TOP_N) |>
  dplyr::pull(genus)

comp_df <- dom_df |>
  dplyr::mutate(
    genus2 = dplyr::if_else(
      dominated,
      dplyr::if_else(genus %in% top_genera, genus, "Other"),
      NA_character_
    )
  ) |>
  dplyr::group_by(Group) |>
  dplyr::mutate(n_group = dplyr::n()) |>
  dplyr::ungroup() |>
  dplyr::filter(!is.na(genus2)) |>
  dplyr::count(Group, genus2, n_group, name = "n_dom_genus") |>
  dplyr::mutate(frac = n_dom_genus / n_group) |>
  tidyr::complete(Group, genus2, fill = list(n_dom_genus = 0, frac = 0)) |>
  dplyr::mutate(
    Group  = forcats::fct_relevel(Group, PRETTY_LEVELS),
    genus2 = forcats::fct_infreq(genus2)
  ) |>
  dplyr::arrange(Group, dplyr::desc(frac))

lvl       <- levels(comp_df$genus2)
pal_genus <- colorspace::qualitative_hcl(length(lvl), palette = "Set 2")
names(pal_genus) <- lvl

fig2_p_dom_comp <- ggplot2::ggplot(comp_df, ggplot2::aes(x = Group, y = frac, fill = genus2)) +
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

fig2_p_dom_comp_overlay <- fig2_p_dom_comp +
  ggpubr::stat_pvalue_manual(
    ord_dom |>
      dplyr::filter(!is.na(label)),
    label         = "label",
    hide.ns       = FALSE,
    tip.length    = .01,
    bracket.size  = .6,
    step.increase = 0,
    size          = 4
  ) +
  ggplot2::annotate(
    "text",
    x        = Inf,
    y        = Inf,
    hjust    = 1.02,
    vjust    = 1.2,
    label    = overall_lab_dom,
    fontface = "italic",
    size     = 4
  ) +
  ggplot2::coord_cartesian(ylim = c(0, 1.02), clip = "off")

print(fig2_p_dom_comp_overlay)

# ===== FIGURE 2C =====
# β-diversity (Bray–Curtis PCoA + PERMANOVA)

cat_levels <- c("Baseline", "OnLen_control", "OnLen_diarrhea")
cat_labels <- c(
  Baseline       = "Baseline",
  OnLen_control  = "Len (control)",
  OnLen_diarrhea = "Len (diarrhea)"
)

cat_cols  <- c(
  Baseline       = "#1B9E77",
  OnLen_control  = "#D95F02",
  OnLen_diarrhea = "#7570B3"
)

metadata <- metadata |>
  mutate(
    SampleCategory = factor(SampleCategory, levels = cat_levels),
    id             = factor(as.numeric(as.factor(SubjectID)))
  )

braydist_all <- phyloseq::distance(phy, method = "bray")
ord          <- phyloseq::ordinate(phy, method = "PCoA", distance = braydist_all)

ord_vectors <- ord$vectors[rownames(df_taxa_genus), , drop = FALSE]

pc_var <- ord$values$Relative_eig
vx <- if (!is.null(pc_var) && length(pc_var) >= 1) round(pc_var[1] * 100, 1) else NA_real_
vy <- if (!is.null(pc_var) && length(pc_var) >= 2) round(pc_var[2] * 100, 1) else NA_real_

x_lab <- if (is.finite(vx)) sprintf("PCoA 1 (%.1f%%)", vx) else "PCoA 1"
y_lab <- if (is.finite(vy)) sprintf("PCoA 2 (%.1f%%)", vy) else "PCoA 2"

braydist <- vegan::vegdist(df_taxa_genus, method = "bray")

permtab <- as.data.frame(
  vegan::adonis2(
    braydist ~ SampleCategory + id,
    by   = "terms",
    data = metadata
  )
)

f_col <- dplyr::case_when(
  "F"       %in% names(permtab) ~ "F",
  "F.Model" %in% names(permtab) ~ "F.Model",
  TRUE                          ~ names(permtab)[2]
)

F_stat <- as.numeric(permtab["SampleCategory", f_col])
P_val  <- as.numeric(permtab["SampleCategory", "Pr(>F)"])

ann_perm <- sprintf(
  "PERMANOVA F=%.2f, p=%s",
  F_stat,
  ifelse(P_val < .001, "<0.001", sprintf("%.3f", P_val))
)

df_ordination <- ord_vectors |>
  as.data.frame() |>
  rownames_to_column("Sample.ID") |>
  select(Sample.ID, Axis.1, Axis.2) |>
  left_join(metadata, by = "Sample.ID") |>
  mutate(SampleCategory = factor(SampleCategory, levels = cat_levels))

cat_alphas <- c(
  Baseline       = 0.03,
  OnLen_control  = 0.07,
  OnLen_diarrhea = 0.15
)

p_beta <- ggplot(df_ordination, aes(Axis.1, Axis.2, color = SampleCategory, fill = SampleCategory)) +
  stat_ellipse(
    aes(alpha = SampleCategory),
    geom        = "polygon",
    linewidth   = 1,
    level       = 0.95,
    show.legend = TRUE
  ) +
  geom_point(size = 2, alpha = 0.6) +
  scale_color_manual(
    values = cat_cols,
    labels = cat_labels,
    name   = ""
  ) +
  scale_fill_manual(
    values = cat_cols,
    labels = cat_labels,
    name   = ""
  ) +
  scale_alpha_manual(
    values = cat_alphas,
    guide  = "none" 
  ) +
  labs(
    x     = x_lab,
    y     = y_lab,
    title = "Beta-diversity (Bray-Curtis)"
  ) +
  scale_x_continuous(expand = expansion(mult = 0.12)) +
  scale_y_continuous(expand = expansion(mult = 0.12)) +
  coord_equal(expand = TRUE, clip = "off") +
  theme_classic(base_size = 14) +
  theme(
    plot.margin     = margin(6, 10, 6, 6),
    legend.position = "right",
    legend.direction = "vertical",
    panel.grid      = element_blank()
  ) +
  annotate(
    "text",
    x     = Inf,
    y     = -Inf,
    hjust = 1.02,
    vjust = -0.6,
    label = ann_perm,
    size  = 3.5
  )

print(p_beta)

# ===== FIGURE 2D =====
# Within-subject Bray–Curtis distance (baseline → on-treatment)

metadata <- metadata |>
  dplyr::mutate(
    SampleCategory = factor(
      SampleCategory,
      levels = c("Baseline", "OnLen_control", "OnLen_diarrhea")
    )
  )

keep_ids <- intersect(rownames(df_taxa_genus), metadata$Sample.ID)
abun <- df_taxa_genus[keep_ids, , drop = FALSE]
abun <- abun[, colSums(abun, na.rm = TRUE) > 0, drop = FALSE]

dm <- as.matrix(vegan::vegdist(abun, method = "bray"))

pairs_df <- metadata |>
  dplyr::filter(
    Sample.ID %in% rownames(abun),
    SampleCategory %in% c("Baseline", "OnLen_control", "OnLen_diarrhea")
  ) |>
  dplyr::group_by(SubjectID) |>
  dplyr::summarise(
    base_id  = dplyr::first(Sample.ID[SampleCategory == "Baseline"],       default = NA_character_),
    diarr_id = dplyr::first(Sample.ID[SampleCategory == "OnLen_diarrhea"], default = NA_character_),
    ctrl_id  = dplyr::first(Sample.ID[SampleCategory == "OnLen_control"],  default = NA_character_),
    .groups  = "drop"
  ) |>
  dplyr::mutate(
    onlen_id  = dplyr::if_else(!is.na(diarr_id), diarr_id, ctrl_id),
    onlen_cat = dplyr::if_else(
      !is.na(diarr_id), "OnLen_diarrhea",
      dplyr::if_else(!is.na(ctrl_id), "OnLen_control", NA_character_)
    )
  ) |>
  dplyr::filter(!is.na(base_id), !is.na(onlen_id))

get_pair_dist <- function(b, o) {
  if (is.na(b) || is.na(o) || !(b %in% rownames(dm)) || !(o %in% rownames(dm))) {
    return(NA_real_)
  }
  as.numeric(dm[b, o])
}

dist_df <- pairs_df |>
  dplyr::mutate(
    dist  = mapply(get_pair_dist, base_id, onlen_id),
    group = dplyr::recode(
      onlen_cat,
      "OnLen_control"  = "Len (control)",
      "OnLen_diarrhea" = "Len (diarrhea)"
    )
  ) |>
  dplyr::filter(!is.na(dist)) |>
  dplyr::mutate(
    group = factor(group, levels = c("Len (control)", "Len (diarrhea)"))
  )

sumtab <- dist_df |>
  dplyr::group_by(group) |>
  dplyr::summarise(
    n      = dplyr::n(),
    median = stats::median(dist),
    mean   = base::mean(dist),
    sd     = stats::sd(dist),
    .groups = "drop"
  )
print(sumtab)

wilx <- stats::wilcox.test(dist ~ group, data = dist_df, exact = FALSE)
p_val <- wilx$p.value

grp_cols  <- c(
  "Len (control)"  = "#D95F02",
  "Len (diarrhea)" = "#7570B3"
)

p_shift <- ggplot2::ggplot(dist_df, ggplot2::aes(x = group, y = dist, fill = group, colour = group)) +
  ggplot2::geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.30, linewidth = 0.7) +
  ggplot2::geom_jitter(width = 0.12, size = 1.6, alpha = 0.65) +
  ggplot2::scale_fill_manual(values = grp_cols, guide = "none") +
  ggplot2::scale_colour_manual(values = grp_cols, guide = "none") +
  ggplot2::labs(
    x     = NULL,
    y     = "Bray–Curtis distance",
    title = "Bray–Curtis distance (Baseline → Len)"
  ) +
  ggplot2::theme_classic(base_size = 13) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 20, hjust = 1)) +
  ggplot2::annotate(
    "text",
    x     = Inf,
    y     = Inf,
    hjust = 1.02,
    vjust = 1.2,
    label = sprintf(
      "Wilcoxon p%s",
      ifelse(p_val < .001, "<0.001", sprintf("= %.3f", p_val))
    ),
    fontface = "italic",
    size     = 4
  ) +
  ggplot2::coord_cartesian(clip = "off")

print(p_shift)

# ===== FIGURE 2E =====
# PICRUSt2 KO predictions (BSH, 7α-HSDH, bai operon)

make_ko_cpm <- function(ko_long, ko_ids, metadata, pseudo = 1) {
  df <- ko_long |>
    dplyr::mutate(Sample.ID = as.character(Sample.ID),
                SubjectID = sub("-.*$", "", Sample.ID)) |>
    dplyr::group_by(Sample.ID) |>
    dplyr::mutate(total_pred = sum(Abundance, na.rm = TRUE)) |>
    dplyr::ungroup() |>
    dplyr::filter(KO %in% ko_ids) |>
    dplyr::mutate(CPM = Abundance / total_pred * 1e6) |>
    dplyr::group_by(Sample.ID) |>
    dplyr::summarise(
      KO_CPM  = sum(CPM, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::left_join(
      metadata |>
        dplyr::select(Sample.ID, SubjectID, SampleCategory),
      by = "Sample.ID"
    ) |>
    dplyr::filter(!is.na(SubjectID), !is.na(SampleCategory)) |>
    dplyr::mutate(
      SubjectID       = factor(SubjectID),
      SampleCat = factor(
        SampleCategory,
        levels = c("Baseline", "OnLen_control", "OnLen_diarrhea")
      ),
      log_CPM   = log10(KO_CPM + pseudo)
    )
  
  if (nrow(df) == 0L) {
    stop(
      "make_ko_cpm: no rows after filtering for KO ids. ",
      "Check that picrust_ko_long$KO matches ko_ids (e.g. 'ko:K01442' vs 'K01442')."
    )
  }
  
  df
}

plot_ko_cpm <- function(
    df,
    title           = "Predicted KO abundance",
    y_label         = "log10(CPM)",
    font_base       = 14,
    axis_text_size  = 12,
    axis_title_size = 12,
    bracket_top_pad = 0.10,
    bracket_gap     = 0.06,
    y_limits        = NULL,
    y_breaks        = waiver(),
    p_override      = NULL  
) {
  m0 <- lme4::lmer(log_CPM ~ SampleCat + (1 | SubjectID), data = df)
  
  if (!is.null(p_override)) {
    p_global <- p_override
  } else {
    aov_tab <- tryCatch(
      car::Anova(m0, type = 3, test.statistic = "F"),
      error = function(e) NULL
    )
    if (!is.null(aov_tab) &&
        "SampleCat" %in% rownames(aov_tab) &&
        "Pr(>F)"   %in% colnames(aov_tab)) {
      p_global <- aov_tab["SampleCat", "Pr(>F)"]
    } else {
      p_global <- NA_real_
    }
  }
  
  group_summ <- df |>
    dplyr::group_by(SampleCat) |>
    dplyr::summarise(
      mean = mean(log_CPM),
      se   = sd(log_CPM) / sqrt(dplyr::n()),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      ymin = mean - 1.96 * se,
      ymax = mean + 1.96 * se
    )
  
  pw <- emmeans::emmeans(m0, "SampleCat") |>
    emmeans::contrast("pairwise", adjust = "tukey") |>
    as.data.frame() |>
    tidyr::separate(contrast, c("group1", "group2"), " - ") |>
    dplyr::mutate(
      label = cut(
        p.value,
        c(-Inf, .001, .01, .05, Inf),
        c("***", "**", "*", "ns"),
        right = FALSE
      )
    )
  
  ord <- tibble::tibble(
    group1 = c("Baseline", "OnLen_control", "Baseline"),
    group2 = c("OnLen_control", "OnLen_diarrhea", "OnLen_diarrhea")
  ) |>
    dplyr::left_join(pw, by = c("group1", "group2"))
  
  yr <- diff(range(df$log_CPM, na.rm = TRUE))
  if (!is.finite(yr) || yr == 0) yr <- 1
  
  if (is.null(y_limits)) {
    anchor <- max(group_summ$ymax, na.rm = TRUE) + bracket_top_pad * yr
    dir    <- +1
    ylim_up <- NA_real_
  } else {
    anchor <- y_limits[2] - bracket_top_pad * yr
    dir    <- -1
    ylim_up <- y_limits[2]
  }
  
  ord <- ord |>
    dplyr::mutate(
      y.position = anchor + dir * bracket_gap * yr * (dplyr::row_number() - 1)
    )
  
  ord_valid <- ord |>
    dplyr::filter(!is.na(label), !is.na(y.position), is.finite(y.position))
  
  if (is.na(ylim_up)) {
    if (nrow(ord_valid) > 0) {
      ylim_up <- max(ord_valid$y.position, na.rm = TRUE) + 0.06 * yr
    } else {
      ylim_up <- max(df$log_CPM, na.rm = TRUE) + 0.12 * yr
    }
  }
  
  ggplot2::ggplot() +
    ggplot2::geom_line(
      data   = df,
      ggplot2::aes(SampleCat, log_CPM, group = SubjectID),
      colour   = "grey70",
      linewidth = .3,
      alpha     = .4
    ) +
    gghalves::geom_half_violin(
      data  = df,
      ggplot2::aes(SampleCat, log_CPM, fill = SampleCat),
      side  = "l",
      width = .9,
      alpha = .5,
      trim  = FALSE
    ) +
    ggplot2::geom_jitter(
      data  = df,
      ggplot2::aes(SampleCat, log_CPM, colour = SampleCat),
      width = .1,
      size  = 1.4,
      alpha = .60
    ) +
    ggplot2::geom_ribbon(
      data        = group_summ,
      ggplot2::aes(SampleCat, ymin = ymin, ymax = ymax, group = 1),
      inherit.aes = FALSE,
      fill        = "#377eb8",
      alpha       = .20
    ) +
    ggplot2::geom_line(
      data   = group_summ,
      ggplot2::aes(SampleCat, mean, group = 1),
      colour   = "#377eb8",
      linewidth = 1.2
    ) +
    ggplot2::scale_x_discrete(
      labels = c(
        Baseline       = "Baseline",
        OnLen_control  = "Len (control)",
        OnLen_diarrhea = "Len (diarrhea)"
      )
    ) +
    ggplot2::scale_y_continuous(
      limits = y_limits,
      breaks = y_breaks,
      expand = ggplot2::expansion(mult = c(0.02, if (is.null(y_limits)) 0.08 else 0))
    ) +
    ggplot2::scale_fill_brewer(palette = "Dark2") +
    ggplot2::scale_colour_brewer(palette = "Dark2") +
    ggplot2::labs(
      x     = NULL,
      y     = y_label,
      title = title
    ) +
    ggplot2::theme_classic(base_size = font_base) +
    ggplot2::theme(
      legend.position = "none",
      axis.text.x  = ggplot2::element_text(size = axis_text_size, angle = 20, hjust = 1),
      axis.text.y  = ggplot2::element_text(size = axis_text_size),
      axis.title.y = ggplot2::element_text(
        size   = axis_title_size,
        margin = ggplot2::margin(r = 8)
      )
    ) +
    ggpubr::stat_pvalue_manual(
      ord_valid,
      label         = "label",
      hide.ns       = FALSE,
      tip.length    = .01,
      bracket.size  = .6,
      step.increase = 0,
      size          = font_base * 0.30
    ) +
    ggplot2::annotate(
      "text",
      x      = Inf,
      y      = ylim_up,
      hjust  = 1.02,
      vjust  = 1.2,
      label  = sprintf(
        "Overall p %s",
        ifelse(
          is.na(p_global), "NA",
          ifelse(p_global < .001, "<0.001", sprintf("%.3f", p_global))
        )
      ),
      fontface = "italic",
      size     = axis_title_size * 0.35
    ) +
    ggplot2::coord_cartesian(ylim = c(NA, ylim_up), clip = "off")
}

bsh_df <- make_ko_cpm(
  ko_long  = picrust_ko_long,
  ko_ids   = "ko:K01442",
  metadata = metadata,
  pseudo   = 1
)
p_bsh <- plot_ko_cpm(
  bsh_df,
  title    = "Bile salt hydrolase (BSH)",
  y_limits = c(1.8, 3.7),
  y_breaks = seq(2, 3.5, by = 0.5)
)
print(p_bsh)

hsdh_df <- make_ko_cpm(
  ko_long  = picrust_ko_long,
  ko_ids   = "ko:K00076",
  metadata = metadata,
  pseudo   = 1
)
p_hsdh <- plot_ko_cpm(
  hsdh_df,
  title    = "7α-HSDH",
  y_limits = c(0.5, 3.9),
  y_breaks = seq(0.5, 3.5, by = 0.5)
)
print(p_hsdh)

bai_kos <- c(
  "ko:K15868", "ko:K15869", "ko:K15870", "ko:K15871",
  "ko:K15872", "ko:K15873", "ko:K15874", "ko:K07007"
)
bai_df <- make_ko_cpm(
  ko_long  = picrust_ko_long,
  ko_ids   = bai_kos,
  metadata = metadata,
  pseudo   = 1
)
p_bai <- plot_ko_cpm(
  bai_df,
  title    = "bai operon",
  y_limits = c(2.6, 3.5),
  y_breaks = seq(2.75, 3.5, by = 0.25)
)
print(p_bai)

# ===== SUPPLEMENTAL FIGURE 2A =====
# Paired α-diversity (within-patient)

pts_both <- metadata %>%
  dplyr::filter(SampleCategory %in% c("OnLen_control", "OnLen_diarrhea"),
                !is.na(invsimpson)) %>%
  dplyr::group_by(SubjectID) %>%
  dplyr::filter(all(c("OnLen_control", "OnLen_diarrhea") %in% SampleCategory)) %>%
  dplyr::ungroup()

n_pts <- dplyr::n_distinct(pts_both$SubjectID)
cat(sprintf("\nPatients with both Len-control and Len-diarrhea samples: N=%d\n", n_pts))

meta_paired <- pts_both %>%
  dplyr::mutate(
    SubjectID = factor(SubjectID),
    SampleCat = factor(SampleCategory,
                       levels = c("OnLen_control", "OnLen_diarrhea")),
    y = log10(invsimpson)
  )

cat("\n--- Sample counts per patient ---\n")
print(table(meta_paired$SubjectID, meta_paired$SampleCat))

wide_paired <- meta_paired %>%
  dplyr::select(SubjectID, SampleCat, y) %>%
  tidyr::pivot_wider(names_from = SampleCat, values_from = y)

t_res <- t.test(wide_paired$OnLen_control, wide_paired$OnLen_diarrhea, paired = TRUE)
cat("\n--- Paired t-test ---\n")
print(t_res)

w_res <- wilcox.test(wide_paired$OnLen_control, wide_paired$OnLen_diarrhea, paired = TRUE)
cat("\n--- Wilcoxon signed-rank test ---\n")
print(w_res)

seg_df <- wide_paired %>%
  dplyr::transmute(
    SubjectID,
    x_num    = 1,
    xend_num = 2,
    y        = OnLen_control,
    yend     = OnLen_diarrhea,
    direction = ifelse(OnLen_diarrhea < OnLen_control, "Decreased", "Increased")
  )

cat("\n--- Direction of change (Len-control → Len-diarrhea) ---\n")
print(table(seg_df$direction))

dir_cols <- c("Decreased" = "#E31A1C", "Increased" = "#0072B2")

p_paired <- w_res$p.value
p_label <- sprintf("Wilcoxon signed-rank p %s",
                   ifelse(p_paired < .001, "< 0.001", sprintf("= %.3f", p_paired)))

p_paired_plot <- ggplot2::ggplot(meta_paired, ggplot2::aes(x = SampleCat, y = y)) +
  ggplot2::geom_boxplot(
    ggplot2::aes(fill = SampleCat),
    width = 0.4, alpha = 0.25, outlier.shape = NA,
    colour = "black", linewidth = 0.5
  ) +
  ggplot2::geom_segment(
    data = seg_df,
    ggplot2::aes(x = x_num, xend = xend_num, y = y, yend = yend,
                 linewidth = direction, alpha = direction),
    colour = "grey40"
  ) +
  ggplot2::scale_linewidth_manual(values = c("Decreased" = 0.8, "Increased" = 0.3), guide = "none") +
  ggplot2::scale_alpha_manual(values = c("Decreased" = 0.5, "Increased" = 0.3), guide = "none") +
  ggplot2::geom_point(
    ggplot2::aes(fill = SampleCat),
    shape = 21, size = 3, alpha = 0.85, colour = "black", stroke = 0.4
  ) +
  ggplot2::scale_x_discrete(
    labels = c(OnLen_control  = "Len (control)",
               OnLen_diarrhea = "Len (diarrhea)")
  ) +
  ggplot2::scale_fill_manual(
    values = c("OnLen_control"  = "#D95F02",
               "OnLen_diarrhea" = "#7570B3"),
    guide  = "none"
  ) +
  ggplot2::labs(
    x     = NULL,
    y     = "log10(Inv. Simpson)",
    title = bquote(alpha * "-diversity (within-patient, N=" * .(n_pts) * ")")
  ) +
  ggplot2::theme_classic(base_size = 14) +
  ggplot2::theme(
    legend.position  = "none",
    axis.text.x      = ggplot2::element_text(size = 12, angle = 20, hjust = 1),
    axis.text.y      = ggplot2::element_text(size = 12),
    axis.title.y     = ggplot2::element_text(size = 12, margin = ggplot2::margin(r = 8))
  ) +
  ggplot2::annotate(
    "text", x = 1.5, y = max(meta_paired$y, na.rm = TRUE) + 0.08,
    label = p_label,
    fontface = "italic", size = 3.5
  ) +
  ggplot2::coord_cartesian(
    ylim = c(min(meta_paired$y, na.rm = TRUE) - 0.05,
             max(meta_paired$y, na.rm = TRUE) + 0.15),
    clip = "off"
  )

print(p_paired_plot)

# ===== SUPPLEMENTAL FIGURE 2B =====
# Butyrate producers

butyrate_genera_vital2018 <- c(
  "Acetonema",
  "Acidaminococcus",
  "Acidipropionibacterium",  # formerly Propionibacterium acidifaciens
  "Agathobaculum",           # formerly Eubacterium desmolans
  "Agathobacter",            # formerly Eubacterium rectale
  "Alistipes",
  "Anaerococcus",
  "Anaerofustis",
  "Anaerostipes",
  "Anaerobutyricum",         # formerly Eubacterium hallii
  "Anaerotruncus",
  "Brachyspira",
  "Butyrivibrio",
  "Clostridioides",          # formerly Clostridium difficile
  "Clostridium",
  "Coprococcus",
  "Eubacterium",
  "Faecalibacterium",        
  "Faecalitalea",            # formerly Eubacterium cylindroides
  "Flavonifractor",          # formerly Clostridium orbiscindens
  "Fusobacterium",
  "Holdemanella",            # formerly Eubacterium biforme
  "Lachnoanaerobaculum",
  "Lachnoclostridium",       
  "Lacrimispora",            # formerly Clostridium saccharolyticum
  "Megasphaera",
  "Odoribacter",
  "Peptoniphilus",
  "Porphyromonas",
  "Propionibacterium",       
  "Pseudobutyrivibrio",     
  "Pseudoramibacter",
  "Roseburia",               
  "Shuttleworthia",
  "Subdoligranulum",
  "Treponema"
)

calc_butyrate_abundance <- function(genus_table, butyrate_genera) {
  rel_abun <- sweep(genus_table, 1, rowSums(genus_table), "/")
  
  bp_present <- intersect(colnames(rel_abun), butyrate_genera)
  
  cat("Butyrate-producing genera found in your data:\n")
  cat(paste(bp_present, collapse = ", "), "\n")
  cat("Number found:", length(bp_present), "out of", length(butyrate_genera), "\n\n")
  
  if (length(bp_present) == 0) {
    warning("No butyrate-producing genera found in dataset")
    return(data.frame(
      Sample.ID = rownames(genus_table),
      butyrate_rel_abun = 0
    ))
  }
  
  butyrate_abun <- rowSums(rel_abun[, bp_present, drop = FALSE])
  
  data.frame(
    Sample.ID = rownames(genus_table),
    butyrate_rel_abun = butyrate_abun,
    stringsAsFactors = FALSE
  )
}

bp_abun <- calc_butyrate_abundance(df_taxa_genus, butyrate_genera_vital2018)

bp_data <- metadata %>%
  left_join(bp_abun, by = "Sample.ID") %>%
  filter(!is.na(butyrate_rel_abun), !is.na(SubjectID), !is.na(SampleCategory)) %>%
  mutate(
    SubjectID = as.factor(SubjectID),
    SampleCat = factor(
      SampleCategory, 
      levels = c("Baseline", "OnLen_control", "OnLen_diarrhea")
    ),
    log_butyrate = log10(butyrate_rel_abun + 1e-6)
  )

bp_lmer <- lmer(log_butyrate ~ SampleCat + (1 | SubjectID), data = bp_data)

bp_anova <- car::Anova(bp_lmer, type = 3, test.statistic = "F")
p_global_bp <- bp_anova["SampleCat", "Pr(>F)"]

bp_group_summ <- bp_data %>%
  group_by(SampleCat) %>%
  summarise(
    mean = mean(log_butyrate, na.rm = TRUE),
    se   = sd(log_butyrate, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  ) %>%
  mutate(
    ymin = mean - 1.96 * se,
    ymax = mean + 1.96 * se
  )

bp_emm <- emmeans::emmeans(bp_lmer, "SampleCat")
bp_pw <- emmeans::contrast(bp_emm, method = "pairwise", adjust = "tukey") %>%
  as.data.frame() %>%
  tidyr::separate(contrast, c("group1", "group2"), " - ") %>%
  mutate(
    label = cut(
      p.value,
      c(-Inf, .001, .01, .05, Inf),
      c("***", "**", "*", "ns"),
      right = FALSE
    )
  )

bp_ord <- tibble(
  group1 = c("Baseline", "OnLen_control", "Baseline"),
  group2 = c("OnLen_control", "OnLen_diarrhea", "OnLen_diarrhea")
) %>%
  left_join(bp_pw, by = c("group1", "group2"))

yr_bp <- diff(range(bp_data$log_butyrate, na.rm = TRUE))
anchor_bp <- max(bp_group_summ$ymax, na.rm = TRUE) + 0.5 * yr_bp

bp_ord <- bp_ord %>%
  mutate(
    y.position = anchor_bp + 0.08 * yr_bp * (row_number() - 1)
  )

ylim_up_bp <- max(bp_ord$y.position, na.rm = TRUE) + 0.08 * yr_bp

p_butyrate <- ggplot() +
  geom_line(
    data = bp_data,
    aes(SampleCat, log_butyrate, group = SubjectID),
    colour = "grey70",
    linewidth = 0.3,
    alpha = 0.4
  ) +
  geom_half_violin(
    data = bp_data,
    aes(SampleCat, log_butyrate, fill = SampleCat),
    side = "l",
    width = 0.9,
    alpha = 0.5,
    trim = FALSE
  ) +
  geom_jitter(
    data = bp_data,
    aes(SampleCat, log_butyrate, colour = SampleCat),
    width = 0.1,
    size = 1.4,
    alpha = 0.60
  ) +
  geom_ribbon(
    data = bp_group_summ,
    aes(SampleCat, ymin = ymin, ymax = ymax, group = 1),
    inherit.aes = FALSE,
    fill = "#377eb8",
    alpha = 0.20
  ) +
  geom_line(
    data = bp_group_summ,
    aes(SampleCat, mean, group = 1),
    colour = "#377eb8",
    linewidth = 1.2
  ) +
  scale_x_discrete(
    labels = c(
      Baseline = "Baseline",
      OnLen_control = "Len (control)",
      OnLen_diarrhea = "Len (diarrhea)"
    )
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0.08, 0.12))
  ) +
  scale_fill_brewer(palette = "Dark2") +
  scale_colour_brewer(palette = "Dark2") +
  labs(
    x = NULL,
    y = "log10(Relative abundance)",
    title = "Butyrate-producing genera"
  ) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 12, angle = 20, hjust = 1),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 12, margin = margin(r = 8))
  ) +
  stat_pvalue_manual(
    bp_ord %>% filter(!is.na(label)),
    label = "label",
    hide.ns = FALSE,
    tip.length = 0.01,
    bracket.size = 0.6,
    step.increase = 0,
    size = 14 * 0.30
  ) +
  annotate(
    "text",
    x = Inf,
    y = ylim_up_bp,
    hjust = 1.02,
    vjust = 0,
    label = sprintf(
      "Overall p %s",
      ifelse(is.na(p_global_bp), "NA",
             ifelse(p_global_bp < .001, "<0.001", sprintf("%.3f", p_global_bp)))
    ),
    fontface = "italic",
    size = 12 * 0.35
  ) +
  coord_cartesian(ylim = c(NA, ylim_up_bp), clip = "off")

print(p_butyrate)

cat("\n=== Butyrate Producer Analysis Summary ===\n")
cat("Group means (log10 relative abundance):\n")
print(bp_group_summ)

cat("\nPairwise comparisons (Tukey-adjusted):\n")
print(bp_pw %>% select(group1, group2, estimate, p.value, label))

cat("\nOverall ANOVA p-value:", 
    ifelse(p_global_bp < .001, "<0.001", sprintf("%.3f", p_global_bp)), "\n")

# ===== SUPPLEMENTAL FIGURE 2C-2E =====
# FLORAL GEE

metadata <- readr::read_csv(
  "data/metadata.csv",
  show_col_types = FALSE
)

phy <- readRDS("data/len.rds")

df_taxa_genus <- left_join(
  tax_table(phy) %>% as.data.frame()  %>% rownames_to_column("key"),
  otu_table(phy) %>% as.data.frame() %>% rownames_to_column("key")
) %>%
  pivot_longer(cols = -c(key:species)) %>%
  dplyr::select(-c(key, kingdom, phylum, class, ordr, family, species)) %>%
  pivot_wider(names_from="genus", values_fn = sum, values_fill = 0) %>%
  column_to_rownames("name")

df_taxa_genus <- df_taxa_genus[rownames(df_taxa_genus) %in% metadata$Sample.ID,]

df_taxa_genus_filtered <- df_taxa_genus[,colMeans(df_taxa_genus > 0) > 0.1]

df_combined <- df_taxa_genus_filtered |>
  rownames_to_column(var="Sample.ID") |>
  left_join(metadata) |>
  mutate(SubjectID = sub("-.*$", "", Sample.ID),
         id        = as.numeric(factor(SubjectID)))

# part 1 - len diarrhea vs. baseline
df_diarrhea <- df_combined |> 
  filter(SampleCategory != "OnLen_control") |> 
  mutate(outcome = ifelse(SampleCategory == "OnLen_diarrhea",1,0))

x <- df_diarrhea %>% 
  select(Blautia:Oscillibacter) %>% 
  as.matrix()

floral.fit.diarrhea <- mcv.FLORAL(mcv=10,
                                  ncore=3,
                                  seed=123,
                                  x=x,
                                  y=df_diarrhea$outcome,
                                  ncov=0,
                                  pseudo=1,
                                  longitudinal = TRUE, 
                                  family="binomial",
                                  corstr = "exchangeable",
                                  intercept=TRUE,
                                  ncov.lambda.weight = 1,
                                  scalefix=FALSE,
                                  lambda.min.ratio = 1e-4,
                                  id=df_diarrhea$id,
                                  mu=1e5,
                                  pfilter=0.1,
                                  ncv=5,
                                  a=3.7,
                                  step2=TRUE,
                                  progress=TRUE,
                                  plot=TRUE)

floral.fit.diarrhea$p_min
floral.fit.diarrhea$p_min_ratio
floral.fit.diarrhea$p_1se
floral.fit.diarrhea$p_1se_ratio

# part 2 - len control vs. baseline
df_normal <- df_combined |> 
  filter(SampleCategory != "OnLen_diarrhea")|> 
  mutate(outcome = ifelse(SampleCategory == "OnLen_control",1,0))

x <- df_normal %>% 
  select(Blautia:Oscillibacter) %>% 
  as.matrix()

floral.fit.normal <- mcv.FLORAL(mcv=10,
                                ncore=5,
                                seed=123,
                                x=x,
                                y=df_normal$outcome,
                                ncov=0,
                                pseudo=1,
                                longitudinal = TRUE, 
                                family="binomial",
                                corstr = "exchangeable",
                                intercept=TRUE,
                                ncov.lambda.weight = 1,
                                scalefix=FALSE,
                                lambda.min.ratio = 1e-4,
                                id=df_normal$id,
                                mu=1e5,
                                pfilter=0.1,
                                ncv=5,
                                a=3.7,
                                step2=TRUE,
                                progress=TRUE,
                                plot=TRUE)

floral.fit.normal$p_min
floral.fit.normal$p_min_ratio
floral.fit.normal$p_1se
floral.fit.normal$p_1se_ratio

# part 3 - len diarrhea vs. len control
df_lenalidomide <- df_combined |> 
  filter(SampleCategory != "Baseline")|> 
  mutate(outcome = ifelse(SampleCategory == "OnLen_diarrhea",1,0))

x <- df_lenalidomide %>% 
  select(Blautia:Oscillibacter) %>% 
  as.matrix()

floral.fit.lenalidomide <- mcv.FLORAL(mcv=10,
                                ncore=5,
                                seed=123,
                                x=x,
                                y=df_lenalidomide$outcome,
                                ncov=0,
                                pseudo=1,
                                longitudinal = TRUE, 
                                family="binomial",
                                corstr = "exchangeable",
                                intercept=TRUE,
                                ncov.lambda.weight = 1,
                                scalefix=FALSE,
                                lambda.min.ratio = 1e-4,
                                id=df_lenalidomide$id,
                                mu=1e5,
                                pfilter=0.1,
                                ncv=5,
                                a=3.7,
                                step2=TRUE,
                                progress=TRUE,
                                plot=TRUE)

floral.fit.lenalidomide$p_min
floral.fit.lenalidomide$p_min_ratio
floral.fit.lenalidomide$p_1se
floral.fit.lenalidomide$p_1se_ratio
