#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(GenomicRanges)
  library(IRanges)
})

option_list <- list(
  make_option(c("--unified-bed"), type = "character",
              help = "Unified ORF BED12 file (e.g. unified_orfs.bed)"),
  make_option(c("--bedgraph-dir"), type = "character",
              help = "Directory containing RiboseQC bedgraph files"),
  make_option(c("--sample-pattern"), type = "character",
              default = "^(.+)_P_sites_(plus|minus)\\.bedgraph$",
              help = "Regex with 2 capture groups: sample, strand [default: %default]"),
  make_option(c("--output"), type = "character", default = "orf_sample_pn.long.tsv",
              help = "Output long-format TSV [default: %default]"),
  make_option(c("--drop-zero"), action = "store_true", default = FALSE,
              help = "Drop rows with p_sites == 0 [default: %default]")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$`unified-bed`) || is.null(opt$`bedgraph-dir`)) {
  print_help(opt_parser)
  stop("Both --unified-bed and --bedgraph-dir are required", call. = FALSE)
}
if (!file.exists(opt$`unified-bed`)) {
  stop("Unified BED not found: ", opt$`unified-bed`, call. = FALSE)
}
if (!dir.exists(opt$`bedgraph-dir`)) {
  stop("Bedgraph directory not found: ", opt$`bedgraph-dir`, call. = FALSE)
}

parse_int_list <- function(x) {
  if (is.na(x) || !nzchar(x)) {
    return(integer(0))
  }
  parts <- strsplit(x, ",", fixed = TRUE)[[1]]
  parts <- parts[nzchar(parts)]
  if (!length(parts)) {
    return(integer(0))
  }
  as.integer(parts)
}

load_orf_exons <- function(bed_path) {
  bed <- fread(
    bed_path,
    header = FALSE,
    sep = "\t",
    fill = TRUE,
    showProgress = FALSE
  )
  if (ncol(bed) < 12) {
    stop("Expected BED12 input, got ", ncol(bed), " columns.", call. = FALSE)
  }
  setnames(
    bed,
    c("chr", "start0", "end1", "orf_id", "score", "strand",
      "thickStart", "thickEnd", "itemRgb", "blockCount", "blockSizes", "blockStarts")
  )

  exon_list <- vector("list", nrow(bed))
  lengths <- integer(nrow(bed))
  ids <- as.character(bed$orf_id)

  for (i in seq_len(nrow(bed))) {
    sizes <- parse_int_list(bed$blockSizes[i])
    starts <- parse_int_list(bed$blockStarts[i])
    n <- min(length(sizes), length(starts))
    if (n == 0) {
      exon_list[[i]] <- GRanges()
      lengths[i] <- 0L
      next
    }
    sizes <- sizes[seq_len(n)]
    starts <- starts[seq_len(n)]
    exon_starts <- as.integer(bed$start0[i]) + starts + 1L
    exon_ends <- exon_starts + sizes - 1L

    exon_list[[i]] <- GRanges(
      seqnames = as.character(bed$chr[i]),
      ranges = IRanges(start = exon_starts, end = exon_ends),
      strand = as.character(bed$strand[i]),
      orf_id = ids[i]
    )
    lengths[i] <- sum(sizes)
  }

  grl <- GRangesList(exon_list)
  names(grl) <- ids
  exons <- unlist(grl, use.names = FALSE)
  mcols(exons)$orf_id <- rep(ids, lengths(grl))
  list(exons = exons, orf_ids = ids, orf_len = setNames(lengths, ids))
}

parse_sample_files <- function(bg_dir, sample_pattern) {
  files <- list.files(bg_dir, pattern = "\\.bedgraph$", full.names = TRUE)
  if (!length(files)) {
    stop("No *.bedgraph files found in ", bg_dir, call. = FALSE)
  }

  parsed <- lapply(files, function(fp) {
    bn <- basename(fp)
    m <- regexec(sample_pattern, bn, perl = TRUE)
    hit <- regmatches(bn, m)[[1]]
    if (length(hit) == 3) {
      data.table(path = fp, sample = hit[2], strand = hit[3])
    } else {
      NULL
    }
  })
  info <- rbindlist(parsed, use.names = TRUE, fill = TRUE)
  if (!nrow(info)) {
    stop("No bedgraph files matched --sample-pattern: ", sample_pattern, call. = FALSE)
  }

  info[, strand := tolower(strand)]
  info <- info[strand %in% c("plus", "minus")]
  if (!nrow(info)) {
    stop("Matched files but no plus/minus strand captured by pattern.", call. = FALSE)
  }
  info[]
}

load_bedgraph <- function(path, strand_symbol) {
  if (!file.exists(path)) {
    return(GRanges())
  }
  bg <- fread(
    path,
    header = FALSE,
    sep = "\t",
    fill = TRUE,
    showProgress = FALSE
  )
  if (ncol(bg) < 4 || !nrow(bg)) {
    return(GRanges())
  }
  bg <- bg[, .(chr = as.character(V1), start0 = suppressWarnings(as.integer(V2)),
               end0 = suppressWarnings(as.integer(V3)), score = suppressWarnings(as.numeric(V4)))]
  bg <- bg[!is.na(start0) & !is.na(end0) & !is.na(score) & end0 > start0]
  if (!nrow(bg)) {
    return(GRanges())
  }
  GRanges(
    seqnames = bg$chr,
    ranges = IRanges(start = bg$start0 + 1L, end = bg$end0),
    strand = strand_symbol,
    score = bg$score
  )
}

count_psites <- function(psite_gr, orf_exons) {
  if (!length(psite_gr) || !length(orf_exons)) {
    return(setNames(numeric(0), character(0)))
  }
  h <- findOverlaps(psite_gr, orf_exons, ignore.strand = FALSE)
  if (!length(h)) {
    return(setNames(numeric(0), character(0)))
  }

  q <- queryHits(h)
  s <- subjectHits(h)
  ov_bp <- pmax(
    0L,
    pmin(end(psite_gr)[q], end(orf_exons)[s]) - pmax(start(psite_gr)[q], start(orf_exons)[s]) + 1L
  )
  contrib <- mcols(psite_gr)$score[q] * ov_bp
  tapply(contrib, mcols(orf_exons)$orf_id[s], sum)
}

cat("Loading unified ORFs...\n")
orf <- load_orf_exons(opt$`unified-bed`)
orf_ids <- orf$orf_ids
orf_len <- orf$orf_len
orf_exons <- orf$exons
cat(sprintf("  Loaded %d ORFs\n", length(orf_ids)))

cat("Discovering sample bedgraphs...\n")
bg_info <- parse_sample_files(opt$`bedgraph-dir`, opt$`sample-pattern`)
samples <- sort(unique(bg_info$sample))
cat(sprintf("  Matched %d sample(s)\n", length(samples)))

res <- rbindlist(lapply(samples, function(sm) {
  counts <- setNames(numeric(length(orf_ids)), orf_ids)

  plus_file <- bg_info[sample == sm & strand == "plus", path]
  minus_file <- bg_info[sample == sm & strand == "minus", path]

  if (length(plus_file) >= 1) {
    plus_gr <- load_bedgraph(plus_file[1], "+")
    plus_cnt <- count_psites(plus_gr, orf_exons)
    if (length(plus_cnt)) {
      counts[names(plus_cnt)] <- counts[names(plus_cnt)] + as.numeric(plus_cnt)
    }
  }
  if (length(minus_file) >= 1) {
    minus_gr <- load_bedgraph(minus_file[1], "-")
    minus_cnt <- count_psites(minus_gr, orf_exons)
    if (length(minus_cnt)) {
      counts[names(minus_cnt)] <- counts[names(minus_cnt)] + as.numeric(minus_cnt)
    }
  }

  data.table(
    sample = sm,
    orf_id = orf_ids,
    length_nt = as.numeric(orf_len[orf_ids]),
    p_sites = as.numeric(counts[orf_ids]),
    pN = as.numeric(counts[orf_ids]) / as.numeric(orf_len[orf_ids])
  )
}), use.names = TRUE)

if (isTRUE(opt$`drop-zero`)) {
  res <- res[p_sites > 0]
}

fwrite(res, opt$output, sep = "\t")
cat("Done. Wrote:", opt$output, "\n")
