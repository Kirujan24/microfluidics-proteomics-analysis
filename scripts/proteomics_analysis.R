# scripts/proteomics_analysis.R
# LC–MS/MS proteomics (microfluidic urinary EV study)
# Marker check + optional overlap summaries + GO enrichment
# Put raw Excel export in: data/proteomics/
#
# Edit ONCE:
#   1) file_raw + description_col (match your sheet)
#   2) intensity_pattern (how your intensity columns are named)
#   3) marker gene lists (if your panel changes)
#   4) overlap_symbols (the overlap set you used for GO)
library(tidyverse)
library(readxl)
library(stringr)
library(clusterProfiler)
library(org.Rn.eg.db)
library(VennDiagram)
library(grid)
library(plotrix)

# paths
input_file <- "data/proteomics/Exosomes_lcms_raw.xlsx"
output_dir <- "outputs/proteomics"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# load data
proteins <- read_excel(input_file)

if (!"Description" %in% colnames(proteins)) {
  stop("Expected column 'Description' not found.")
}

# extract gene symbols from UniProt-style descriptions
proteins <- proteins %>%
  mutate(Gene = str_extract(Description, "(?<=GN=)[^ ]+"))

missing_genes <- proteins %>% filter(is.na(Gene) | Gene == "")
if (nrow(missing_genes) > 0) {
  write_csv(
    missing_genes,
    file.path(output_dir, "proteins_missing_gene_symbol.csv")
  )
}

# marker sets
markers_ev <- c("Cd9", "Pdcd6ip", "Mfge8", "Rab7a", "Ezr")
markers_liver <- c(
  "Aldob","Apoa1","Apoe","Fgb","Fgg",
  "Gsta1","Gsta3","Gstp1","Gclc","Gclm",
  "Gpx3","Hpx","Prdx2","Prdx5"
)
markers_housekeeping <- c("Gapdh","Actb")

markers <- c(markers_ev, markers_liver, markers_housekeeping)

# subset marker proteins
marker_df <- proteins %>%
  filter(Gene %in% markers) %>%
  mutate(Gene = factor(Gene, levels = markers))

# intensity columns (experiment-specific naming)
intensity_cols <- grep(
  "Pre|Pass|Control|Early|Intermediate|Late",
  colnames(marker_df),
  value = TRUE
)

if (length(intensity_cols) == 0) {
  stop("No intensity columns detected.")
}

# build matrix
intensity_mat <- marker_df %>%
  select(all_of(intensity_cols)) %>%
  mutate(across(everything(), as.numeric)) %>%
  as.matrix()

rownames(intensity_mat) <- marker_df$Gene
intensity_mat[intensity_mat == 0] <- NA

# log transform
log_mat <- log10(intensity_mat)
log_mat[!is.finite(log_mat)] <- NA
log_mat[is.na(log_mat)] <- min(log_mat, na.rm = TRUE) - 2

plot_df <- as.data.frame(log_mat) %>%
  rownames_to_column("Gene") %>%
  pivot_longer(-Gene, names_to = "Sample", values_to = "Intensity") %>%
  mutate(Gene = factor(Gene, levels = markers))

# heatmap
ggplot(plot_df, aes(Sample, Gene, fill = Intensity)) +
  geom_tile() +
  scale_fill_gradient2(
    midpoint = median(plot_df$Intensity, na.rm = TRUE),
    name = "log10 intensity"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    title = "Proteomics marker expression",
    x = NULL,
    y = NULL
  )

# GO enrichment on overlapping proteins
overlap_symbols <- c(
  "Nhlrc3","Cpne3","Atp6v1a","Smpd1","Egfr","Mdh1","Gsta1",
  "Ctsb","Hbb","Alb","Apoa1","Gsta3","Hp","Anxa5",
  "F2","Nme2","Acp2","Igh","Slc2a5","Lgals5","Anxa6",
  "Ctsc","Pgls","Ctbs","Hexa","Aadat","Scpep1","Slc9a3r1",
  "Sult1c2a"
)

entrez_map <- bitr(
  overlap_symbols,
  fromType = "SYMBOL",
  toType   = "ENTREZID",
  OrgDb    = org.Rn.eg.db
)

gene_ids <- unique(entrez_map$ENTREZID)

run_go <- function(ids, ontology) {
  enrichGO(
    gene          = ids,
    OrgDb         = org.Rn.eg.db,
    keyType       = "ENTREZID",
    ont           = ontology,
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.1,
    qvalueCutoff  = 0.2,
    readable      = TRUE
  ) %>%
    as.data.frame()
}

plot_go <- function(go_df, title_text) {
  if (nrow(go_df) == 0) return(NULL)

  go_df %>%
    arrange(p.adjust) %>%
    slice_head(n = 15) %>%
    mutate(score = -log10(p.adjust)) %>%
    ggplot(aes(score, reorder(Description, score))) +
    geom_col() +
    theme_bw() +
    labs(
      title = title_text,
      x = "-log10(adj p)",
      y = NULL
    )
}

plot_go(run_go(gene_ids, "BP"), "GO: Biological Process")
plot_go(run_go(gene_ids, "MF"), "GO: Molecular Function")
plot_go(run_go(gene_ids, "CC"), "GO: Cellular Component")

