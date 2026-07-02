#####
# 00_setup.R — shared setup sourced by every figure/table script.
#####

# ---- Packages (all packages used across the figure/table scripts) ----
suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(tibble)
  library(forcats); library(stringr); library(purrr)
  library(ggplot2); library(patchwork); library(scales)
  library(RColorBrewer); library(colorspace)
  library(gghalves); library(ggpubr); library(ggnewscale); library(ggrepel)
  library(vegan); library(permute); library(phyloseq)
  library(lme4); library(lmerTest); library(emmeans); library(car)
  library(geepack); library(gtsummary); library(gt)
  library(FLORAL)
  library(conflicted)
})

# ---- Resolve function-name conflicts ----
.prefer <- function(name, pkg)
  try(suppressWarnings(conflicted::conflict_prefer(name, pkg, quiet = TRUE)), silent = TRUE)
.prefer("filter","dplyr");    .prefer("select","dplyr");    .prefer("mutate","dplyr")
.prefer("rename","dplyr");    .prefer("summarise","dplyr"); .prefer("summarize","dplyr")
.prefer("slice","dplyr");     .prefer("count","dplyr");     .prefer("lag","dplyr")
.prefer("first","dplyr");     .prefer("expand","tidyr");    .prefer("lmer","lmerTest")
rm(.prefer)

# ---- Constants & colors ----
LOG10_OFFSET <- 1
cat_cols     <- c("No" = "#1B9E77", "Yes" = "#D55E00")

# ---- Input paths  ----
DATA_DIR    <- "data"
PATH_BILE   <- file.path(DATA_DIR, "bile_acids.csv")
PATH_META   <- file.path(DATA_DIR, "metadata.csv")
PATH_PHY    <- file.path(DATA_DIR, "len.rds")
PATH_PICRST <- file.path(DATA_DIR, "picrust_ko.csv")

# ---- Load inputs ----
bile_long <- readr::read_csv(PATH_BILE, show_col_types = FALSE) |>
  dplyr::mutate(Sample.ID = as.character(Sample.ID),
                SubjectID = sub("-.*$", "", Sample.ID)) 

metadata <- readr::read_csv(PATH_META, show_col_types = FALSE) |>
  dplyr::mutate(Sample.ID = as.character(Sample.ID),
                SubjectID = sub("-.*$", "", Sample.ID))

phy <- readRDS(PATH_PHY)

picrust_ko_long <- readr::read_csv(PATH_PICRST, show_col_types = FALSE)

# ---- Rebuild df_taxa_genus and alpha-diversity (invsimpson) ----
df_taxa_genus <- dplyr::left_join(
  phyloseq::tax_table(phy) %>% as.data.frame() %>% tibble::rownames_to_column("key"),
  phyloseq::otu_table(phy)  %>% as.data.frame() %>% tibble::rownames_to_column("key"),
  by = "key"
) |>
  tidyr::pivot_longer(cols = -c(key:species)) |>
  dplyr::select(-c(key, kingdom, phylum, class, ordr, family, species)) |>
  tidyr::pivot_wider(names_from = "genus", values_fn = sum, values_fill = 0) |>
  tibble::column_to_rownames("name")

metadata <- metadata |>
  dplyr::filter(Sample.ID %in% rownames(df_taxa_genus)) |>
  dplyr::arrange(Sample.ID)

df_taxa_genus <- df_taxa_genus[metadata$Sample.ID, , drop = FALSE]
rownames(metadata) <- metadata$Sample.ID

df_alpha <- data.frame(
  Sample.ID  = rownames(df_taxa_genus),
  invsimpson = apply(df_taxa_genus / rowSums(df_taxa_genus), 1,
                     function(x) 1 / sum(x^2))
)
metadata <- metadata |> dplyr::left_join(df_alpha, by = "Sample.ID")

# ---- Bile-acid summaries + ratios ----
div0    <- function(num, den) ifelse(is.finite(den) & den > 0, num / den, NA_real_)
logdiff <- function(a, b) log10(a + LOG10_OFFSET) - log10(b + LOG10_OFFSET)
safe_ratio <- div0

acid_variants <- list(
  CA   = c("cholate (CA)",            "cholic acid"),
  DCA  = c("deoxycholate (DCA)",      "deoxycholic acid"),
  CDCA = c("chenodeoxycholate (CDCA)","chenodeoxycholic acid"),
  LCA  = c("lithocholate (LCA)",      "lithocholic acid"),
  GCA  = c("glycocholate (GCA)",      "glycocholic acid"),
  TCA  = c("taurocholate (TCA)",      "taurocholic acid")
)

bile_focus <- bile_long |>
  dplyr::mutate(target = dplyr::case_when(
    bilenames %in% acid_variants$CA   ~ "CA",
    bilenames %in% acid_variants$DCA  ~ "DCA",
    bilenames %in% acid_variants$CDCA ~ "CDCA",
    bilenames %in% acid_variants$LCA  ~ "LCA",
    bilenames %in% acid_variants$GCA  ~ "GCA",
    bilenames %in% acid_variants$TCA  ~ "TCA",
    TRUE ~ NA_character_)) |>
  dplyr::filter(!is.na(target)) |>
  dplyr::group_by(SubjectID, Sample.ID, SampleCat, target) |>
  dplyr::summarise(value = sum(Amount, na.rm = TRUE), .groups = "drop") |>
  dplyr::mutate(value_log10 = log10(value + LOG10_OFFSET)) |>
  tidyr::pivot_wider(names_from = target, values_from = c(value, value_log10),
                     names_sep = "_") |>
  dplyr::rename_with(~ sub("^value_(.*)$", "\\1", .x), starts_with("value_"))

meta_dist <- metadata |> dplyr::distinct(Sample.ID, .keep_all = TRUE)

df_ratio <- bile_focus |>
  dplyr::mutate(
    CA_DCA_ratio   = safe_ratio(CA,   DCA),
    CA_DCA_log10   = logdiff(CA,      DCA),
    CDCA_LCA_ratio = safe_ratio(CDCA, LCA),
    CDCA_LCA_log10 = logdiff(CDCA,    LCA)
  ) |>
  dplyr::left_join(
    meta_dist |> dplyr::select(Sample.ID, LaterDiarrhea, SampleCategory),
    by = "Sample.ID"
  )

# ---- PICRUSt2 KO CPM tables (BSH, 7a-HSDH, bai operon) ----
make_ko_cpm <- function(ko_long, ko_ids, metadata, pseudo = 1) {
  df <- ko_long |>
    dplyr::mutate(Sample.ID = as.character(Sample.ID)) |>
    dplyr::group_by(Sample.ID) |>
    dplyr::mutate(total_pred = sum(Abundance, na.rm = TRUE)) |>
    dplyr::ungroup() |>
    dplyr::filter(KO %in% ko_ids) |>
    dplyr::mutate(CPM = Abundance / total_pred * 1e6) |>
    dplyr::group_by(Sample.ID) |>
    dplyr::summarise(KO_CPM = sum(CPM, na.rm = TRUE), .groups = "drop") |>
    dplyr::left_join(
      metadata |> dplyr::select(Sample.ID, SubjectID, SampleCategory),
      by = "Sample.ID"
    ) |>
    dplyr::filter(!is.na(SubjectID), !is.na(SampleCategory)) |>
    dplyr::mutate(
      SubjectID = factor(SubjectID),
      SampleCat = factor(SampleCategory,
                         levels = c("Baseline", "OnLen_control", "OnLen_diarrhea")),
      log_CPM   = log10(KO_CPM + pseudo)
    )
  if (nrow(df) == 0L) stop("make_ko_cpm: no rows after KO filter; check ko:Kxxxxx ids.")
  df
}

bsh_df  <- make_ko_cpm(picrust_ko_long, "ko:K01442", metadata)
hsdh_df <- make_ko_cpm(picrust_ko_long, "ko:K00076", metadata)
bai_kos <- c("ko:K15868","ko:K15869","ko:K15870","ko:K15871",
             "ko:K15872","ko:K15873","ko:K15874","ko:K07007")
bai_df  <- make_ko_cpm(picrust_ko_long, bai_kos, metadata)

# ---- Generic formatting helpers (one-liners, used across scripts) ----
fmt_p <- function(p) ifelse(!is.finite(p), "NA",
                     ifelse(p < .001, "<0.001", sprintf("=%.3f", p)))
fmt_p_simple <- fmt_p
wilcox_label <- function(p) paste0("Wilcoxon p", fmt_p(p))
lab_pq <- function(p, q = NULL)
  if (is.null(q)) paste0("p", fmt_p(p)) else paste0("p", fmt_p(p), ", q", fmt_p(q))

message("00_setup.R: loaded ", nrow(metadata), " samples, ",
        ncol(df_taxa_genus), " genera.")
