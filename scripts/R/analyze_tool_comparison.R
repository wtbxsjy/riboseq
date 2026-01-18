#!/usr/bin/env Rscript
#############################################################################
# Multi-Tool ORF Prediction Comparison Analysis
#
# Purpose: Analyze concordance and specificity of different ORF prediction tools
# Input:
#   - GENCODE unified ORF annotation (.orfs.out from gencode-riboseqORFs)
# Output:
#   - Tool comparison statistics
#   - High-confidence ORF lists
#   - Tool-specific ORF lists
#   - Visualization plots (Venn diagrams, UpSet plots, heatmaps)
#
# Usage:
#   Rscript analyze_tool_comparison.R \
#     --input results/Mouse_AllTools.orfs.out \
#     --tools "ORFquant,RiboTISH,Ribotricer" \
#     --outdir tool_comparison_results
#############################################################################

suppressPackageStartupMessages({
  library(optparse)
  library(tidyverse)
  library(ggplot2)
  library(ggvenn)
  library(UpSetR)
  library(ComplexHeatmap)
  library(pheatmap)
  library(RColorBrewer)
})

# ========== Command Line Arguments ==========
option_list <- list(
  make_option(c("--input"), type = "character",
              help = "GENCODE ORF output file (.orfs.out)"),
  make_option(c("--tools"), type = "character",
              default = "ORFquant,RiboTISH,Ribotricer",
              help = "Comma-separated tool names [default: %default]"),
  make_option(c("--outdir"), type = "character",
              default = "tool_comparison_results",
              help = "Output directory [default: %default]"),
  make_option(c("--min-samples"), type = "integer", default = 10,
              help = "Minimum samples for high-confidence ORFs [default: %default]"),
  make_option(c("--min-tools"), type = "integer", default = 2,
              help = "Minimum tools for high-confidence ORFs [default: %default]")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Validate inputs
if (is.null(opt$input)) {
  print_help(opt_parser)
  stop("--input is required", call. = FALSE)
}

if (!file.exists(opt$input)) {
  stop("Input file not found: ", opt$input, call. = FALSE)
}

# Parse tool names
TOOLS <- str_split(opt$tools, ",")[[1]]

# Create output directory
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

cat("========================================\n")
cat("Multi-Tool ORF Prediction Comparison\n")
cat("========================================\n\n")

# ========== Step 1: Load Data ==========
cat("[Step 1/7] Loading GENCODE ORF annotation...\n")

df <- read_tsv(opt$input, show_col_types = FALSE)

cat(sprintf("  ✓ Loaded %d ORFs\n", nrow(df)))
cat(sprintf("  ✓ Columns: %d\n", ncol(df)))

# ========== Step 2: Extract Tool Detection Matrix ==========
cat("\n[Step 2/7] Extracting tool detection matrix...\n")

# For each tool, find columns matching the tool name
tool_detection <- list()

for (tool in TOOLS) {
  # Find columns containing tool name (case-insensitive)
  tool_cols <- grep(tool, colnames(df), ignore.case = TRUE, value = TRUE)

  if (length(tool_cols) == 0) {
    warning(sprintf("No columns found for tool: %s", tool))
    next
  }

  # Each ORF detected by tool if ANY sample has detection (max across samples)
  tool_detection[[tool]] <- apply(df[, tool_cols, drop = FALSE], 1, max, na.rm = TRUE)

  detected_count <- sum(tool_detection[[tool]] > 0, na.rm = TRUE)
  cat(sprintf("  ✓ %s: %d ORFs (%.1f%%) from %d sample columns\n",
              tool, detected_count, 100 * detected_count / nrow(df), length(tool_cols)))
}

# Convert to data frame
detection_df <- as_tibble(tool_detection)
detection_df$orf_id <- df$orf_id

cat(sprintf("  ✓ Detection matrix: %d ORFs × %d tools\n", nrow(detection_df), length(TOOLS)))

# ========== Step 3: Tool Overlap Analysis (Venn/UpSet) ==========
cat("\n[Step 3/7] Analyzing tool overlaps...\n")

# Calculate all possible tool combinations
library(gtools)  # for combinations

# Count ORFs detected by each tool combination
calc_tool_combinations <- function(detection_df, tools) {
  results <- list()

  # Single tools only
  for (tool in tools) {
    mask <- detection_df[[tool]] > 0
    for (other_tool in setdiff(tools, tool)) {
      mask <- mask & (detection_df[[other_tool]] == 0)
    }
    results[[paste0(tool, "_only")]] <- sum(mask, na.rm = TRUE)
  }

  # Pairwise intersections (excluding 3rd tool)
  if (length(tools) >= 2) {
    pairs <- combn(tools, 2, simplify = FALSE)
    for (pair in pairs) {
      mask <- (detection_df[[pair[1]]] > 0) & (detection_df[[pair[2]]] > 0)
      for (other_tool in setdiff(tools, pair)) {
        mask <- mask & (detection_df[[other_tool]] == 0)
      }
      results[[paste(pair, collapse = "_")]] <- sum(mask, na.rm = TRUE)
    }
  }

  # All tools intersection
  if (length(tools) >= 3) {
    mask <- Reduce(`&`, lapply(tools, function(t) detection_df[[t]] > 0))
    results[["all_tools"]] <- sum(mask, na.rm = TRUE)
  }

  return(results)
}

combinations <- calc_tool_combinations(detection_df, TOOLS)

cat("\n  Tool combination statistics:\n")
for (combo_name in names(combinations)) {
  cat(sprintf("    %-40s: %6d ORFs\n", combo_name, combinations[[combo_name]]))
}

# ========== Step 4: ORF Biotype Distribution per Tool ==========
cat("\n[Step 4/7] Analyzing ORF biotype distributions...\n")

orf_type_stats <- list()

for (tool in TOOLS) {
  if (!tool %in% names(detection_df)) next

  # ORFs detected by this tool
  tool_orfs <- df[detection_df[[tool]] > 0, ]

  # Count by biotype
  type_counts <- table(tool_orfs$orf_biotype)

  orf_type_stats[[tool]] <- tibble(
    tool = tool,
    orf_biotype = names(type_counts),
    count = as.integer(type_counts)
  )
}

orf_type_df <- bind_rows(orf_type_stats)

# Pivot for display
orf_type_wide <- orf_type_df %>%
  pivot_wider(names_from = tool, values_from = count, values_fill = 0)

cat("\n  ORF biotype distribution:\n")
print(orf_type_wide, n = Inf)

# ========== Step 5: Sample Support Analysis ==========
cat("\n[Step 5/7] Analyzing sample support per tool...\n")

for (tool in TOOLS) {
  tool_cols <- grep(tool, colnames(df), ignore.case = TRUE, value = TRUE)

  if (length(tool_cols) == 0) next

  # Count samples detecting each ORF
  sample_support <- rowSums(df[, tool_cols, drop = FALSE] > 0, na.rm = TRUE)

  cat(sprintf("\n  %s sample support distribution:\n", tool))
  cat(sprintf("    1-5 samples:   %6d ORFs (%.1f%%)\n",
              sum(sample_support <= 5), 100 * sum(sample_support <= 5) / nrow(df)))
  cat(sprintf("    6-10 samples:  %6d ORFs\n",
              sum(sample_support > 5 & sample_support <= 10)))
  cat(sprintf("    11-20 samples: %6d ORFs\n",
              sum(sample_support > 10 & sample_support <= 20)))
  cat(sprintf("    21-30 samples: %6d ORFs\n",
              sum(sample_support > 20 & sample_support <= 30)))
  cat(sprintf("    >30 samples:   %6d ORFs (%.1f%%)\n",
              sum(sample_support > 30), 100 * sum(sample_support > 30) / nrow(df)))
}

# ========== Step 6: High-Confidence ORF Filtering ==========
cat("\n[Step 6/7] Filtering high-confidence ORFs...\n")

# Number of tools detecting each ORF
n_tools_per_orf <- rowSums(detection_df[, TOOLS, drop = FALSE] > 0, na.rm = TRUE)

# Standard 1: All tools
all_tools_mask <- n_tools_per_orf == length(TOOLS)
high_conf_all_tools <- df[all_tools_mask, ]

cat(sprintf("\n  Standard 1 (all %d tools): %d ORFs\n",
            length(TOOLS), nrow(high_conf_all_tools)))

# Standard 2: At least N tools
at_least_n_mask <- n_tools_per_orf >= opt$`min-tools`
high_conf_min_tools <- df[at_least_n_mask, ]

cat(sprintf("  Standard 2 (≥%d tools):     %d ORFs\n",
            opt$`min-tools`, nrow(high_conf_min_tools)))

# Standard 3: At least N tools + M samples
# Count total samples across all tools
all_tool_cols <- unlist(lapply(TOOLS, function(t) {
  grep(t, colnames(df), ignore.case = TRUE, value = TRUE)
}))

total_sample_support <- rowSums(df[, all_tool_cols, drop = FALSE] > 0, na.rm = TRUE)

stringent_mask <- at_least_n_mask & (total_sample_support >= opt$`min-samples`)
high_conf_stringent <- df[stringent_mask, ]

cat(sprintf("  Standard 3 (≥%d tools + ≥%d samples): %d ORFs\n",
            opt$`min-tools`, opt$`min-samples`, nrow(high_conf_stringent)))

# ========== Step 7: Tool-Specific ORFs ==========
cat("\n[Step 7/7] Identifying tool-specific ORFs...\n")

tool_specific_orfs <- list()

for (tool in TOOLS) {
  if (!tool %in% names(detection_df)) next

  # Only detected by this tool
  mask <- detection_df[[tool]] > 0
  for (other_tool in setdiff(TOOLS, tool)) {
    if (other_tool %in% names(detection_df)) {
      mask <- mask & (detection_df[[other_tool]] == 0)
    }
  }

  specific_orfs <- df[mask, ]
  tool_specific_orfs[[tool]] <- specific_orfs

  cat(sprintf("\n  %s-specific ORFs: %d\n", tool, nrow(specific_orfs)))

  if (nrow(specific_orfs) > 0) {
    type_dist <- table(specific_orfs$orf_biotype)
    cat("    Top biotypes:\n")
    for (i in seq_len(min(5, length(type_dist)))) {
      cat(sprintf("      - %s: %d\n",
                  names(sort(type_dist, decreasing = TRUE))[i],
                  sort(type_dist, decreasing = TRUE)[i]))
    }
  }
}

# ========== Save Results ==========
cat("\n[Saving Results]\n")

# Save high-confidence ORFs
write_tsv(high_conf_stringent,
          file.path(opt$outdir, "high_confidence_orfs.tsv"))
cat(sprintf("  ✓ High-confidence ORFs: %s\n", "high_confidence_orfs.tsv"))

# Save tool-specific ORFs
for (tool in names(tool_specific_orfs)) {
  filename <- sprintf("%s_specific_orfs.tsv", tool)
  write_tsv(tool_specific_orfs[[tool]], file.path(opt$outdir, filename))
  cat(sprintf("  ✓ %s-specific ORFs: %s\n", tool, filename))
}

# Save ORF type distribution
write_tsv(orf_type_wide, file.path(opt$outdir, "orf_biotype_by_tool.tsv"))
cat(sprintf("  ✓ ORF biotype table: %s\n", "orf_biotype_by_tool.tsv"))

# Save summary report
report_file <- file.path(opt$outdir, "tool_comparison_summary.txt")
sink(report_file)
cat("=" %R% 80, "\n")
cat("Multi-Tool ORF Prediction Comparison Report\n")
cat("=" %R% 80, "\n\n")

cat(sprintf("Tools analyzed: %s\n", paste(TOOLS, collapse = ", ")))
cat(sprintf("Total unique ORFs: %d\n\n", nrow(df)))

cat("=" %R% 80, "\n")
cat("1. Tool Detection Statistics\n")
cat("=" %R% 80, "\n")
for (tool in TOOLS) {
  if (tool %in% names(detection_df)) {
    count <- sum(detection_df[[tool]] > 0, na.rm = TRUE)
    cat(sprintf("%-20s: %6d ORFs (%.1f%%)\n",
                tool, count, 100 * count / nrow(df)))
  }
}

cat("\n" %R% ("=" %R% 80), "\n")
cat("2. Tool Overlap Analysis\n")
cat("=" %R% 80, "\n")
for (combo_name in names(combinations)) {
  cat(sprintf("%-40s: %6d ORFs\n", combo_name, combinations[[combo_name]]))
}

cat("\n" %R% ("=" %R% 80), "\n")
cat("3. High-Confidence ORF Sets\n")
cat("=" %R% 80, "\n")
cat(sprintf("All %d tools:                    %6d ORFs\n",
            length(TOOLS), nrow(high_conf_all_tools)))
cat(sprintf("≥%d tools:                        %6d ORFs\n",
            opt$`min-tools`, nrow(high_conf_min_tools)))
cat(sprintf("≥%d tools + ≥%d samples:          %6d ORFs (RECOMMENDED)\n",
            opt$`min-tools`, opt$`min-samples`, nrow(high_conf_stringent)))

sink()
cat(sprintf("  ✓ Summary report: %s\n", "tool_comparison_summary.txt"))

# ========== Generate Visualizations ==========
cat("\n[Generating Visualizations]\n")

# 1. Venn diagram (if 2-3 tools)
if (length(TOOLS) %in% 2:3) {
  tryCatch({
    venn_data <- lapply(TOOLS, function(t) {
      which(detection_df[[t]] > 0)
    })
    names(venn_data) <- TOOLS

    p_venn <- ggvenn(venn_data,
                     fill_color = brewer.pal(length(TOOLS), "Set2"),
                     stroke_size = 0.5,
                     set_name_size = 4)

    ggsave(file.path(opt$outdir, "venn_diagram.pdf"),
           p_venn, width = 8, height = 6)
    cat("  ✓ Venn diagram: venn_diagram.pdf\n")
  }, error = function(e) {
    cat("  ✗ Venn diagram failed:", e$message, "\n")
  })
}

# 2. UpSet plot (for all tool combinations)
tryCatch({
  upset_data <- detection_df %>%
    select(all_of(TOOLS)) %>%
    mutate(across(everything(), ~ as.integer(. > 0)))

  pdf(file.path(opt$outdir, "upset_plot.pdf"), width = 10, height = 6)
  print(upset(as.data.frame(upset_data),
              sets = TOOLS,
              order.by = "freq",
              mainbar.y.label = "Number of ORFs",
              sets.x.label = "ORFs per Tool"))
  dev.off()
  cat("  ✓ UpSet plot: upset_plot.pdf\n")
}, error = function(e) {
  cat("  ✗ UpSet plot failed:", e$message, "\n")
})

# 3. ORF biotype distribution heatmap
tryCatch({
  biotype_matrix <- orf_type_wide %>%
    column_to_rownames("orf_biotype") %>%
    as.matrix()

  # Normalize by column (tool) for percentage
  biotype_pct <- sweep(biotype_matrix, 2, colSums(biotype_matrix), FUN = "/") * 100

  pdf(file.path(opt$outdir, "biotype_heatmap.pdf"), width = 8, height = 10)
  pheatmap(biotype_pct,
           cluster_rows = TRUE,
           cluster_cols = FALSE,
           display_numbers = TRUE,
           number_format = "%.1f",
           main = "ORF Biotype Distribution (%)",
           color = colorRampPalette(c("white", "steelblue"))(50))
  dev.off()
  cat("  ✓ Biotype heatmap: biotype_heatmap.pdf\n")
}, error = function(e) {
  cat("  ✗ Biotype heatmap failed:", e$message, "\n")
})

# 4. Tool concordance bar plot
tryCatch({
  concordance_data <- tibble(
    n_tools = 1:length(TOOLS),
    n_orfs = sapply(1:length(TOOLS), function(n) {
      sum(n_tools_per_orf >= n)
    })
  )

  p_concordance <- ggplot(concordance_data, aes(x = factor(n_tools), y = n_orfs)) +
    geom_bar(stat = "identity", fill = "steelblue", alpha = 0.8) +
    geom_text(aes(label = n_orfs), vjust = -0.5, size = 4) +
    labs(
      title = "ORF Detection Concordance Across Tools",
      x = "Minimum Number of Tools",
      y = "Number of ORFs"
    ) +
    theme_minimal(base_size = 12) +
    theme(panel.grid.minor = element_blank())

  ggsave(file.path(opt$outdir, "concordance_barplot.pdf"),
         p_concordance, width = 8, height = 6)
  cat("  ✓ Concordance plot: concordance_barplot.pdf\n")
}, error = function(e) {
  cat("  ✗ Concordance plot failed:", e$message, "\n")
})

# ========== Completion ==========
cat("\n========================================\n")
cat("✅ Tool Comparison Analysis Complete!\n")
cat("========================================\n\n")
cat("Output files:\n")
cat("  Data:\n")
cat("    - high_confidence_orfs.tsv\n")
cat("    - *_specific_orfs.tsv (per tool)\n")
cat("    - orf_biotype_by_tool.tsv\n")
cat("    - tool_comparison_summary.txt\n")
cat("\n  Figures:\n")
cat("    - venn_diagram.pdf\n")
cat("    - upset_plot.pdf\n")
cat("    - biotype_heatmap.pdf\n")
cat("    - concordance_barplot.pdf\n")
cat("\nNext steps:\n")
cat("  1. Review high-confidence ORFs for downstream validation\n")
cat("  2. Investigate tool-specific ORFs for method biases\n")
cat("  3. Integrate with quantification data for expression filtering\n")

# Save session info
sink(file.path(opt$outdir, "session_info.txt"))
sessionInfo()
sink()
