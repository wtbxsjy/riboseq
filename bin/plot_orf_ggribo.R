#!/usr/bin/env Rscript
# ─── ggRibo 覆盖度图 — Pipeline 集成版 ──────────────────────────────────
#
# 接受 pipeline 输出格式, 为指定 ORF 批量生成 ggRibo P-site 覆盖度图。
# 支持: --orf-ids (逗号分隔), --orf-ids-file (文件列表), --n-top-orfs (按 OCS)
#
# 用法:
#   plot_orf_ggribo.R \
#       --orf-meta unified_orfs.metadata.tsv \
#       --expression expression_summary.tsv \
#       --psites-dir riboseqc/ \
#       --gtf reference.gtf \
#       --orf-ids "ORF_1,ORF_2" \
#       --samples "Ribo_11,Ribo_14" \
#       --output-dir ggribo_plots/
#
# 依赖: ggRibo, ggplot2, txdbmaker, GenomicFeatures, Rsamtools
# 容器: 需预装 ggRibo (Bioconductor); 或使用 Singularity 容器

suppressPackageStartupMessages({
  library(ggRibo)
  library(ggplot2)
  library(optparse)
})

# ─── CLI ──────────────────────────────────────────────────────────────────

option_list <- list(
  make_option("--orf-meta", type="character", default=NULL,
              help="Unified ORF metadata TSV"),
  make_option("--expression", type="character", default=NULL,
              help="Expression summary TSV (可选, 用于读值标注)"),
  make_option("--psites-dir", type="character", default=NULL,
              help="RiboseQC P-site bedgraph 目录"),
  make_option("--gtf", type="character", default=NULL,
              help="参考 GTF (用于 exon 结构显示)"),
  make_option("--orf-ids", type="character", default=NULL,
              help="ORF ID 列表 (逗号分隔)"),
  make_option("--orf-ids-file", type="character", default=NULL,
              help="ORF ID 列表文件 (每行一个)"),
  make_option("--n-top-orfs", type="integer", default=0,
              help="按 OCS/reads 取 top N ORF (若未指定 --orf-ids)"),
  make_option("--samples", type="character", default=NULL,
              help="指定样本 (逗号分隔, 默认取 top 3 表达的样本)"),
  make_option("--output-dir", type="character", default="ggribo_plots",
              help="输出目录 [default: ggribo_plots]"),
  make_option("--extend", type="integer", default=200,
              help="ORF 两侧扩展 bp [default: 200]"),
  make_option("--n-samples-per-orf", type="integer", default=3,
              help="每个 ORF 最多显示样本数 [default: 3]"),
  make_option("--chrom-prefix", type="character", default="",
              help="染色体名前缀, 如 'chr' (空字符串则不添加)")
)
opt <- parse_args(OptionParser(option_list=option_list))

# ─── 验证参数 ─────────────────────────────────────────────────────────────

if (is.null(opt$`orf-meta`)) stop("必需: --orf-meta")
if (is.null(opt$`psites-dir`)) stop("必需: --psites-dir")

orf_meta_file  <- opt$`orf-meta`
expr_file      <- opt$`expression`
psites_dir     <- opt$`psites-dir`
gtf_file       <- opt$`gtf`
orf_ids_str    <- opt$`orf-ids`
orf_ids_file   <- opt$`orf-ids-file`
n_top          <- opt$`n-top-orfs`
samples_str    <- opt$`samples`
output_dir     <- opt$`output-dir`
extend_bp      <- opt$`extend`
n_samples_max  <- opt$`n-samples-per-orf`
chrom_prefix   <- opt$`chrom-prefix`

dir.create(output_dir, showWarnings=FALSE, recursive=TRUE)

cat("=== ggRibo ORF Coverage Plotter ===\n")
cat("ORF metadata:", orf_meta_file, "\n")
cat("P-sites dir:", psites_dir, "\n")
cat("Output dir:", output_dir, "\n")

# ─── 加载 ORF 数据 ────────────────────────────────────────────────────────

orf_meta <- read.delim(orf_meta_file, sep="\t", header=TRUE, stringsAsFactors=FALSE)
cat(sprintf("Loaded %d ORFs from metadata\n", nrow(orf_meta)))

# 验证必需列
required <- c("orf_id", "chrom", "start", "end", "strand")
missing <- setdiff(required, colnames(orf_meta))
if (length(missing) > 0) stop("Metadata 缺少列: ", paste(missing, collapse=", "))

# 染色体前缀
if (nzchar(chrom_prefix)) {
  orf_meta$chrom <- paste0(chrom_prefix, orf_meta$chrom)
}

# ─── 确定 ORF 列表 ────────────────────────────────────────────────────────

orf_ids <- NULL

if (!is.null(orf_ids_str)) {
  orf_ids <- strsplit(orf_ids_str, ",")[[1]]
  orf_ids <- trimws(orf_ids)
  cat(sprintf("指定 %d 个 ORF\n", length(orf_ids)))
} else if (!is.null(orf_ids_file) && file.exists(orf_ids_file)) {
  orf_ids <- readLines(orf_ids_file)
  orf_ids <- trimws(orf_ids[orf_ids != "" & !grepl("^#", orf_ids)])
  cat(sprintf("从文件读取 %d 个 ORF\n", length(orf_ids)))
} else if (n_top > 0) {
  # 按 total_reads 降序取 top N
  if (!is.null(expr_file) && file.exists(expr_file)) {
    expr <- read.delim(expr_file, sep="\t", header=TRUE, stringsAsFactors=FALSE)
    if ("total_reads" %in% colnames(expr)) {
      expr <- expr[order(-expr$total_reads), ]
      orf_ids <- head(expr$orf_id, n_top)
      cat(sprintf("按 total_reads 取 top %d ORF\n", length(orf_ids)))
    }
  }
  if (is.null(orf_ids)) {
    orf_ids <- head(orf_meta$orf_id, n_top)
    cat(sprintf("无表达量文件, 取前 %d ORF\n", length(orf_ids)))
  }
}

if (is.null(orf_ids) || length(orf_ids) == 0) {
  stop("未指定 ORF。请使用 --orf-ids, --orf-ids-file, 或 --n-top-orfs")
}

# ─── 加载表达量 (可选) ────────────────────────────────────────────────────

expr_data <- NULL
if (!is.null(expr_file) && file.exists(expr_file)) {
  expr_data <- read.delim(expr_file, sep="\t", header=TRUE, stringsAsFactors=FALSE)
  cat(sprintf("Loaded expression data: %d ORFs\n", nrow(expr_data)))
}

# ─── Custom Range_info + gtf_import (绕过 ggRibo 命名空间冲突) ───────────

Range_info <- R6::R6Class("Range_info",
  public = list(
    exonsByTx  = NULL,
    txByGene   = NULL,
    cdsByTx    = NULL,
    fiveUTR    = NULL,
    threeUTR   = NULL,
    tx_to_gene = NULL,
    initialize = function(exonsByTx, txByGene, cdsByTx, fiveUTR, threeUTR, tx_to_gene) {
      self$exonsByTx  <- exonsByTx
      self$txByGene   <- txByGene
      self$cdsByTx    <- cdsByTx
      self$fiveUTR    <- fiveUTR
      self$threeUTR   <- threeUTR
      self$tx_to_gene <- tx_to_gene
    }
  )
)

gtf_import_custom <- function(annotation, format="gtf", dataSource="", organism="") {
  txdb <- suppressWarnings(
    txdbmaker::makeTxDbFromGFF(file=annotation, format=format,
                               dataSource=dataSource, organism=organism)
  )
  exonsByTx  <- GenomicFeatures::exonsBy(txdb, by="tx", use.names=TRUE)
  txByGene   <- GenomicFeatures::transcriptsBy(txdb, by="gene")
  cdsByTx    <- GenomicFeatures::cdsBy(txdb, by="tx", use.names=TRUE)
  fiveUTR    <- GenomicFeatures::fiveUTRsByTranscript(txdb, use.names=TRUE)
  threeUTR   <- GenomicFeatures::threeUTRsByTranscript(txdb, use.names=TRUE)
  tx_to_gene <- AnnotationDbi::select(txdb,
    keys    = AnnotationDbi::keys(txdb, keytype="TXNAME"),
    columns = c("TXNAME", "GENEID"),
    keytype = "TXNAME"
  )
  colnames(tx_to_gene) <- c("tx_id", "gene_id")
  tx_to_gene <- tx_to_gene[order(tx_to_gene$tx_id), ]

  Txome_Range <- Range_info$new(
    exonsByTx  = exonsByTx,
    txByGene   = txByGene,
    cdsByTx    = cdsByTx,
    fiveUTR    = fiveUTR,
    threeUTR   = threeUTR,
    tx_to_gene = tx_to_gene
  )
  assign("Txome_Range", Txome_Range, envir = .GlobalEnv)
}

# ─── 辅助函数 ─────────────────────────────────────────────────────────────

# 创建单 exon ORF 的临时 GTF
create_orf_gtf <- function(orf_row, out_path) {
  gene_id  <- sprintf("ORF_%s_%s_%d_%d",
                      orf_row$orf_id, orf_row$chrom, orf_row$start, orf_row$end)
  tx_id    <- paste0(gene_id, "_T001")
  exon_id  <- paste0(gene_id, "_E001")
  prot_id  <- paste0(gene_id, "_P001")

  lines <- c(
    sprintf('%s\tORF\tgene\t%d\t%d\t.\t%s\t.\tgene_id "%s"; gene_source "ORF"; gene_biotype "protein_coding";',
            orf_row$chrom, orf_row$start, orf_row$end, orf_row$strand, gene_id),
    sprintf('%s\tORF\ttranscript\t%d\t%d\t.\t%s\t.\tgene_id "%s"; transcript_id "%s"; gene_source "ORF"; gene_biotype "protein_coding"; transcript_source "ORF"; transcript_biotype "protein_coding";',
            orf_row$chrom, orf_row$start, orf_row$end, orf_row$strand, gene_id, tx_id),
    sprintf('%s\tORF\tCDS\t%d\t%d\t.\t%s\t0\tgene_id "%s"; transcript_id "%s"; exon_number "1"; gene_source "ORF"; gene_biotype "protein_coding"; transcript_source "ORF"; transcript_biotype "protein_coding"; protein_id "%s";',
            orf_row$chrom, orf_row$start, orf_row$end, orf_row$strand, gene_id, tx_id, prot_id),
    sprintf('%s\tORF\texon\t%d\t%d\t.\t%s\t.\tgene_id "%s"; transcript_id "%s"; exon_number "1"; gene_source "ORF"; gene_biotype "protein_coding"; transcript_source "ORF"; transcript_biotype "protein_coding"; exon_id "%s";',
            orf_row$chrom, orf_row$start, orf_row$end, orf_row$strand, gene_id, tx_id, exon_id)
  )

  writeLines(lines, out_path)
  list(gene_id=gene_id, tx_id=tx_id)
}

# 为 ORF 查找 P-site bedgraph
find_bedgraph_files <- function(sample_name, psites_dir) {
  plus_file  <- file.path(psites_dir, paste0(sample_name, "_P_sites_plus.bedgraph"))
  minus_file <- file.path(psites_dir, paste0(sample_name, "_P_sites_minus.bedgraph"))

  if (file.exists(plus_file) && file.exists(minus_file)) {
    return(list(plus=plus_file, minus=minus_file))
  }
  return(NULL)
}

# 找到 ORF 表达量最高的 top K 样本
find_top_samples_for_orf <- function(orf_id, expr_df, n=3, specific_samples=NULL) {
  if (!is.null(specific_samples)) {
    return(specific_samples)
  }

  if (is.null(expr_df)) return(NULL)

  row <- expr_df[expr_df$orf_id == orf_id, ]
  if (nrow(row) == 0) return(NULL)

  read_cols <- grep("_reads$", names(row), value=TRUE)
  if (length(read_cols) == 0) return(NULL)

  reads <- as.numeric(row[1, read_cols])
  names(reads) <- read_cols
  reads <- sort(reads, decreasing=TRUE)
  reads <- reads[reads > 0]

  if (length(reads) > n) reads <- reads[1:n]
  if (length(reads) == 0) return(NULL)

  sub("_reads$", "", names(reads))
}

# ─── 绘图函数 ─────────────────────────────────────────────────────────────

plot_one_orf <- function(orf_row, samples_to_use, expr_df, psites_dir, output_dir, extend_bp) {
  orf_id <- orf_row$orf_id
  cat(sprintf("\n  Plotting %s ...\n", orf_id))

  # 创建临时 GTF
  tmp_gtf <- file.path(output_dir, paste0("tmp_", orf_id, ".gtf"))
  ids <- create_orf_gtf(orf_row, tmp_gtf)

  # 收集 bedgraph 文件
  bg_files <- list()
  sample_names <- c()
  for (s in samples_to_use) {
    bg <- find_bedgraph_files(s, psites_dir)
    if (!is.null(bg)) {
      bg_files[[length(bg_files) + 1]] <- bg
      sample_names <- c(sample_names, s)
    }
  }

  if (length(bg_files) == 0) {
    cat(sprintf("    WARNING: No bedgraph files found for samples: %s\n",
                paste(samples_to_use, collapse=", ")))
    return(invisible(NULL))
  }

  # 构建标题
  orf_len <- orf_row$end - orf_row$start
  title_parts <- sprintf("%s [%s:%d-%d:%s, %dnt]",
                         orf_id, orf_row$chrom, orf_row$start, orf_row$end,
                         orf_row$strand, orf_len)

  # 从表达量数据添加 pN/RPKM 信息
  if (!is.null(expr_df)) {
    erow <- expr_df[expr_df$orf_id == orf_id, ]
    if (nrow(erow) > 0) {
      total_reads <- erow$total_reads[1]
      max_pn <- 0
      pn_cols <- grep("_pN$", names(erow), value=TRUE)
      if (length(pn_cols) > 0) {
        max_pn <- max(as.numeric(erow[1, pn_cols]), na.rm=TRUE)
      }
      title_parts <- sprintf("%s | reads=%s, pN=%.1f",
                             title_parts,
                             format(total_reads, big.mark=","),
                             max_pn)
    }
  }

  n_samps <- length(bg_files)
  cat(sprintf("    %d samples: %s\n", n_samps, paste(sample_names, collapse=", ")))

  result <- tryCatch({
    # 导入自定义注释
    gtf_import_custom(tmp_gtf, format="gtf", dataSource="ORF", organism="Unknown")

    # 创建 seq 输入
    seq_result <- create_seq_input(
      ribo_files   = bg_files,
      sample_names = sample_names,
      include_rna  = FALSE
    )
    assign("inputs_full", seq_result, envir = .GlobalEnv)

    # 调用 ggRibo
    p <- ggRibo(
      gene_id    = ids$gene_id,
      tx_id      = ids$tx_id,
      Extend     = extend_bp,
      NAME       = title_parts,
      Riboseq    = seq_result$Riboseq,
      SampleNames = sample_names,
      GRangeInfo  = Txome_Range,
      data_types  = rep("Ribo-seq", n_samps),
      Y_scale     = "each",
      plot_genomic_direction = TRUE,
      show_seq    = FALSE,
      ribo_linewidth = 0.6
    )

    # 保存 —— 文件名不能太长
    safe_id <- gsub("[:+\\-]", "_", orf_id)
    if (nchar(safe_id) > 80) {
      safe_id <- paste0(substr(safe_id, 1, 60), "..._", digest::digest(orf_id))
    }
    plot_height <- max(5, n_samps * 2.5)
    out_file <- file.path(output_dir, sprintf("%s_ggribo.png", safe_id))
    ggsave(out_file, p, width=14, height=plot_height, dpi=150, limitsize=FALSE)
    cat(sprintf("    -> %s\n", out_file))

    p
  }, error = function(e) {
    cat(sprintf("    ERROR: %s\n", e$message))
    return(invisible(NULL))
  })

  # 清理临时 GTF
  unlink(tmp_gtf)

  return(invisible(result))
}

# ─── Main ──────────────────────────────────────────────────────────────────

cat(sprintf("\n=== Plotting %d ORFs ===\n", length(orf_ids)))

# 解析指定样本
specified_samples <- NULL
if (!is.null(samples_str)) {
  specified_samples <- trimws(strsplit(samples_str, ",")[[1]])
  cat(sprintf("指定样本: %s\n", paste(specified_samples, collapse=", ")))
}

success <- 0
for (oid in orf_ids) {
  orf_row <- orf_meta[orf_meta$orf_id == oid, ]
  if (nrow(orf_row) == 0) {
    cat(sprintf("\nWARNING: ORF '%s' not found in metadata\n", oid))
    next
  }

  # 确定样本
  samples_to_use <- specified_samples
  if (is.null(samples_to_use)) {
    samples_to_use <- find_top_samples_for_orf(oid, expr_data, n=n_samples_max)
  }
  if (is.null(samples_to_use) || length(samples_to_use) == 0) {
    cat(sprintf("\n  WARNING: No samples with data for %s\n", oid))
    next
  }

  res <- plot_one_orf(orf_row[1, ], samples_to_use, expr_data, psites_dir, output_dir, extend_bp)
  if (!is.null(res)) success <- success + 1
}

cat(sprintf("\n=== Done: %d/%d ORFs plotted ===\n", success, length(orf_ids)))
cat(sprintf("Output: %s\n", output_dir))
