# scripts/nta_analysis.R
# NTA across sequential microfluidic passes: Pre-MF → Pass-1 → Pass-2 → Pass-3
# Put NanoSight CSV exports in: data/NTA/
#
# Edit ONCE:
#   1) column names (col_size, col_conc, optional col_pdi)
#   2) file_groups mapping (your exact filenames)

suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
  library(ggpubr)
})

# ----------------------------
# paths
# ----------------------------
nta_path <- "data/NTA"

# ----------------------------
# column names in your NanoSight export
# ----------------------------
col_size <- "Particle size (nm)"
col_conc <- "Concentration (particles/ml)"

# If your export includes PDI, put the exact column name here; otherwise keep NULL
col_pdi <- NULL
# col_pdi <- "PDI"

# ----------------------------
# file -> group mapping (EDIT THIS)
# ----------------------------
# Use exact filenames as they appear in data/NTA/
file_groups <- c(
  # "your_file_name_here.csv" = "Pre-MF",
  # "your_file_name_here.csv" = "Pass-1",
  # "your_file_name_here.csv" = "Pass-2",
  # "your_file_name_here.csv" = "Pass-3"

  "PreMF_rep1.csv"  = "Pre-MF",
  "PreMF_rep2.csv"  = "Pre-MF",
  "Pass1_rep1.csv"  = "Pass-1",
  "Pass1_rep2.csv"  = "Pass-1",
  "Pass2_rep1.csv"  = "Pass-2",
  "Pass2_rep2.csv"  = "Pass-2",
  "Pass3_rep1.csv"  = "Pass-3",
  "Pass3_rep2.csv"  = "Pass-3"
)

# ----------------------------
# read files
# ----------------------------
files <- list.files(nta_path, pattern = "\\.csv$", full.names = TRUE)
if (length(files) == 0) stop("No CSV files found in ", nta_path)

nta_all <- map_dfr(files, function(f) {
  df <- read_csv(f, show_col_types = FALSE)

  if (!col_size %in% names(df)) stop("Missing size col in ", basename(f), " (expected: ", col_size, ")")
  if (!col_conc %in% names(df)) stop("Missing conc col in ", basename(f), " (expected: ", col_conc, ")")

  grp <- unname(file_groups[basename(f)])

  out <- tibble(
    Size_nm       = as.numeric(df[[col_size]]),
    Concentration = as.numeric(df[[col_conc]]),
    File          = basename(f),
    Group         = grp
  )

  if (!is.null(col_pdi) && col_pdi %in% names(df)) {
    out <- mutate(out, PDI = as.numeric(df[[col_pdi]][1]))
  }

  out
}) %>%
  filter(!is.na(Group)) %>%
  mutate(Group = factor(Group, levels = c("Pre-MF", "Pass-1", "Pass-2", "Pass-3")))

if (nrow(nta_all) == 0) {
  stop("No rows after grouping. Check file_groups names match your CSV filenames exactly.")
}

# ----------------------------
# Panel A: size distribution (0–500 nm)
# ----------------------------
p_size <- ggplot(nta_all, aes(Size_nm, Concentration, color = Group)) +
  geom_line(alpha = 0.9, linewidth = 0.8) +
  coord_cartesian(xlim = c(0, 500)) +
  labs(
    title = "Size distribution",
    x = "Particle size (nm)",
    y = "Concentration (particles/mL)"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.title = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

print(p_size)

# ----------------------------
# Panel B: PDI (if present)
# ----------------------------
if ("PDI" %in% names(nta_all)) {
  pdi_df <- nta_all %>%
    distinct(File, Group, PDI) %>%
    filter(!is.na(PDI))

  if (nrow(pdi_df) > 0) {
    p_pdi <- ggplot(pdi_df, aes(Group, PDI)) +
      geom_boxplot(width = 0.55, outlier.shape = NA) +
      geom_point(size = 2, position = position_jitter(width = 0.08)) +
      labs(
        title = "PDI across microfluidic passes",
        x = NULL,
        y = "PDI"
      ) +
      theme_bw(base_size = 12) +
      theme(plot.title = element_text(face = "bold", hjust = 0.5))

    print(p_pdi)
  }
}

# ----------------------------
# Fractions <150 and >150 per file
# ----------------------------
cutoff_nm <- 150

frac_df <- nta_all %>%
  group_by(Group, File) %>%
  summarise(
    total = sum(Concentration, na.rm = TRUE),
    lt150 = sum(Concentration[Size_nm < cutoff_nm], na.rm = TRUE),
    gt150 = sum(Concentration[Size_nm > cutoff_nm], na.rm = TRUE),
    pct_lt150 = 100 * lt150 / total,
    pct_gt150 = 100 * gt150 / total,
    .groups = "drop"
  )

# ----------------------------
# stats + significance brackets (Tukey)
# ----------------------------
tukey_to_brackets <- function(aov_obj, y_start, step = 5) {
  tk <- as.data.frame(TukeyHSD(aov_obj)$Group)
  tk$comparison <- rownames(tk)

  tk %>%
    separate(comparison, into = c("group1", "group2"), sep = "-") %>%
    transmute(
      group1, group2,
      p = `p adj`,
      p.signif = case_when(
        p < 0.0001 ~ "****",
        p < 0.001  ~ "***",
        p < 0.01   ~ "**",
        p < 0.05   ~ "*",
        TRUE       ~ "ns"
      )
    ) %>%
    filter(p.signif != "ns") %>%
    arrange(p) %>%
    mutate(y.position = y_start + seq(0, by = step, length.out = n()))
}

aov_gt <- aov(pct_gt150 ~ Group, data = frac_df)
aov_lt <- aov(pct_lt150 ~ Group, data = frac_df)

# print results to console (handy for manuscript/reporting)
print(summary(aov_gt))
print(as.data.frame(TukeyHSD(aov_gt)$Group))

print(summary(aov_lt))
print(as.data.frame(TukeyHSD(aov_lt)$Group))

# ----------------------------
# Panel C: % >150 nm
# ----------------------------
ymax_gt <- max(frac_df$pct_gt150, na.rm = TRUE)
pvals_gt <- tukey_to_brackets(aov_gt, y_start = ymax_gt + 5, step = 5)

p_gt150 <- ggplot(frac_df, aes(Group, pct_gt150)) +
  stat_summary(fun = mean, geom = "bar", width = 0.65, fill = "grey60") +
  stat_summary(fun.data = mean_sd, geom = "errorbar", width = 0.2) +
  labs(
    title = "% of particles >150 nm",
    x = NULL,
    y = "Particles >150 nm (%)"
  ) +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))

if (nrow(pvals_gt) > 0) {
  p_gt150 <- p_gt150 +
    stat_pvalue_manual(
      pvals_gt,
      label = "p.signif",
      tip.length = 0.01,
      bracket.size = 0.4
    )
}
print(p_gt150)

# ----------------------------
# Panel D: % <150 nm
# ----------------------------
ymax_lt <- max(frac_df$pct_lt150, na.rm = TRUE)
pvals_lt <- tukey_to_brackets(aov_lt, y_start = ymax_lt + 5, step = 5)

p_lt150 <- ggplot(frac_df, aes(Group, pct_lt150)) +
  stat_summary(fun = mean, geom = "bar", width = 0.65, fill = "grey60") +
  stat_summary(fun.data = mean_sd, geom = "errorbar", width = 0.2) +
  labs(
    title = "% of particles <150 nm",
    x = NULL,
    y = "Particles <150 nm (%)"
  ) +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))

if (nrow(pvals_lt) > 0) {
  p_lt150 <- p_lt150 +
    stat_pvalue_manual(
      pvals_lt,
      label = "p.signif",
      tip.length = 0.01,
      bracket.size = 0.4
    )
}
print(p_lt150)
