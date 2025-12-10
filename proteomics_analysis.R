file_raw           <- "Exosomes lcms raw.xlsx"   # path to the raw Excel file
description_col    <- "Description"              # column with the protein header
missing_gene_output <- "proteins_needing_uniprot_lookup.csv"
# Marker gene sets
markers_ev <- c("Cd9", "Pdcd6ip", "Mfge8", "Rab7a", "Ezr")
markers_liver <- c(
  "Aldob","Apoa1","Apoe","Fgb","Fgg",
  "Gsta1","Gsta3","Gstp1","Gclc","Gclm",
  "Gpx3","Hpx","Prdx2","Prdx5"
)
markers_house <- c("Gapdh","Actb")
markers_all <- c(markers_ev, markers_liver, markers_house)
intensity_pattern <- "Control|Early|Intermediate|Late"
suppressPackageStartupMessages({
  library(tidyverse)
  library(readxl)
  library(stringr)
  library(clusterProfiler)
  library(org.Rn.eg.db)
  library(VennDiagram)
  library(grid)
  library(plotrix)
})
# extract gene symbol from UniProt-like description ("... GN=Thbs1 PE=1 SV=2")
extract_gene_from_header <- function(x) {
  g <- stringr::str_extract(x, "(?<=GN=)[^ ]+")
  return(g)   
}
run_go <- function(glist, ontology) {
  ego <- enrichGO(
    gene          = glist,
    OrgDb         = org.Rn.eg.db,
    keyType       = "ENTREZID",
    ont           = ontology,     
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.1,
    qvalueCutoff  = 0.2,
    readable      = TRUE
  )

  df <- as.data.frame(ego)
  if (nrow(df) == 0) return(df)

  if (nrow(df) > 15) df <- df[1:15, ]

  df$log10padj <- -log10(df$p.adjust)
  df$Label     <- paste0(df$Description, " (", df$ID, ")")
  df
}
plot_go <- function(df, title_text, bar_color) {
  if (nrow(df) == 0) return(NULL)

  ggplot(df, aes(x = log10padj, y = reorder(Label, log10padj))) +
    geom_col(fill = bar_color) +
    geom_text(
      aes(label = round(log10padj, 2)),
      hjust = -0.2,
      size  = 4
    ) +
    labs(
      x     = expression(-log[10]("adjusted p value")),
      y     = NULL,
      title = title_text
    ) +
    theme_bw() +
    theme(
      text            = element_text(size = 14),
      plot.title      = element_text(size = 16, face = "bold", hjust = 0.5),
      axis.text.y     = element_text(size = 10),
      axis.text.x     = element_text(size = 10),
      axis.title.x    = element_text(size = 12),
      panel.grid.major= element_blank(),
      panel.grid.minor= element_blank()
    ) +
    xlim(0, max(df$log10padj) + 0.5)
}
if (!exists("all_proteins")) {
  all_proteins  <- paste0("All", 1:50)
}
if (!exists("exo_proteins")) {
  exo_proteins  <- paste0("Exo", 1:30)
}
if (!exists("ref_proteome")) {
  ref_proteome  <- paste0("Ref", 1:40)
}

venn_plot <- venn.diagram(
  x = list(
    `All proteins` = all_proteins,
    Exosomes       = exo_proteins,
    Reference      = ref_proteome
  ),
  filename     = NULL,
  fill         = c("#F5F0F6", "#6667AB", "#CC2936"),
  alpha        = 0.6,
  cex          = 2,
  lwd          = 3,
  cat.cex      = 2,
  cat.fontface = "bold"
)

grid.newpage()
grid::grid.draw(venn_plot)

# pie chart
pie_values <- c(
  Cytosolic               = 3,
  `Liver-associated`      = 20,
  Membrane                = 3,
  `Unclassified proteins` = 73
)

pie_labels <- paste0(names(pie_values), " (", pie_values, "%)")

plotrix::pie3D(
  pie_values,
  labels   = pie_labels,
  explode  = 0.06,
  radius   = 1.3,
  height   = 0.08,
  theta    = pi / 6,
  labelcex = 0.9,
  border   = "grey40",
  col      = c("steelblue3", "chocolate1", "cadetblue3", "tan3"),
  main     = "Category distribution of overlapping proteins (n = 29)"
)
overlap_ids <- c(
  "Nhlrc3","Cpne3","Atp6v1a","Smpd1","Egfr","Mdh1","Gsta1",
  "Ctsb","Hbb","Alb","Apoa1","Gsta3","Hp","Anxa5",
  "F2","Nme2","Acp2","Igh","Slc2a5","Lgals5","Anxa6",
  "Ctsc","Pgls","Ctbs","Hexa","Aadat","Scpep1","Slc9a3r1",
  "Sult1c2a"
)
entrez_ids <- bitr(
  overlap_ids,
  fromType = "SYMBOL",
  toType   = "ENTREZID",
  OrgDb    = org.Rn.eg.db
)

gene_list <- unique(entrez_ids$ENTREZID)

go_bp <- run_go(gene_list, "BP")
go_mf <- run_go(gene_list, "MF")
go_cc <- run_go(gene_list, "CC")

cat("\n=== GO BP ===\n")
print(plot_go(go_bp, "Biological process", "#D62728"))

cat("\n=== GO MF ===\n")
print(plot_go(go_mf, "Molecular function", "#1F77B4"))

cat("\n=== GO CC ===\n")
print(plot_go(go_cc, "Cellular component", "#2CA02C"))
stopifnot(file.exists(file_raw))

raw <- readxl::read_excel(file_raw)
# check that the description column exists
if (!description_col %in% names(raw)) {
  stop("Column '", description_col, "' not found in the input file. ",
       "Please update 'description_col' at the top of the script.")
}

raw2 <- raw %>%
  dplyr::mutate(Gene = extract_gene_from_header(.data[[description_col]]))

# save proteins where GN= was not found so you can check them in UniProt
missing_genes <- raw2 %>%
  dplyr::filter(is.na(Gene))

if (nrow(missing_genes) > 0) {
  readr::write_csv(missing_genes, missing_gene_output)
  message(
    "There are ", nrow(missing_genes),
    " proteins without GN=. They were written to '",
    missing_gene_output,
    "'. Please look them up in UniProt and fill in Gene names manually if needed."
  )
}
df_markers <- raw2 %>%
  dplyr::filter(Gene %in% markers_all) %>%
  dplyr::mutate(Gene = factor(Gene, levels = markers_all))

# columns that contain intensities
intensity_cols <- grep(
  intensity_pattern,
  names(df_markers),
  value = TRUE
)

if (length(intensity_cols) == 0) {
  stop("No intensity columns found matching pattern '", intensity_pattern, "'. ",
       "Please adapt 'intensity_pattern' to your column names.")
}

expr <- df_markers %>%
  dplyr::select(dplyr::all_of(intensity_cols)) %>%
  dplyr::mutate(dplyr::across(everything(), as.numeric)) %>%
  as.matrix()

rownames(expr) <- df_markers$Gene

# log10 transform, handle zeros and missing values
expr[expr == 0] <- NA
expr_log <- log10(expr)
expr_log[!is.finite(expr_log)] <- NA

min_val <- min(expr_log, na.rm = TRUE)
expr_log[is.na(expr_log)] <- min_val - 2

plot_df <- expr_log %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Gene") %>%
  tidyr::pivot_longer(-Gene, names_to = "Sample", values_to = "Value") %>%
  dplyr::mutate(
    Gene = factor(Gene, levels = markers_all),
    Group = dplyr::case_when(
      Gene %in% markers_ev    ~ "EV markers",
      Gene %in% markers_liver ~ "Liver markers",
      Gene %in% markers_house ~ "Housekeeping",
      TRUE                    ~ "Other"
    ),
    Group = factor(Group, levels = c("EV markers", "Liver markers", "Housekeeping", "Other"))
  )

heatmap_plot <- ggplot(plot_df, aes(x = Sample, y = Gene, fill = Value)) +
  geom_tile() +
  scale_fill_gradient2(
    low      = "blue",
    mid      = "white",
    high     = "red",
    midpoint = median(plot_df$Value, na.rm = TRUE),
    name     = "log10 intensity"
  ) +
  facet_grid(
    rows   = vars(Group),
    scales = "free_y",
    space  = "free_y",
    switch = "y"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid        = element_blank(),
    strip.background  = element_blank(),
    strip.placement   = "outside",
    strip.text.y.left = element_text(angle = 0, face = "bold"),
    axis.title.x      = element_blank(),
    axis.title.y      = element_blank(),
    axis.text.x       = element_text(angle = 45, hjust = 1, vjust = 1)
  ) +
  labs(title = "Representative liver- and EV-associated protein intensities")

cat("\n=== Marker heatmap ===\n")
print(heatmap_plot)
