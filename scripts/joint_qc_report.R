#!/usr/bin/env Rscript
#
# joint_qc_report.R - Joint Ribo-seq QC report for nf-core/riboseq
#
# Integrates metrics from riboWaltz, RiboseQC, Ribo-TISH, and Ribotricer to
# produce a unified per-sample quality assessment with letter grades.
#
# Dependencies: data.table (installed at runtime if missing)
#
# Usage:
#   Rscript joint_qc_report.R \
#     --rw_offset_dir ribowaltz_psite/ --rw_region_dir ribowaltz_region/ \
#     --rq_psites_dir riboseqc_psites/ --rt_qual_dir ribotish_qual/ \
#     --rtr_summary_dir ribotricer_summary/ --output_prefix joint_riboseq_qc
#

# ── Ensure data.table ─────────────────────────────────────────────────────────
if (!requireNamespace("data.table", quietly = TRUE)) {
    tryCatch(
        install.packages("data.table", repos = "https://cloud.r-project.org", quiet = TRUE),
        error = function(e) stop("Cannot install data.table: ", e$message)
    )
}
suppressPackageStartupMessages(library(data.table))

# ── Parse command-line args ───────────────────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)
opt <- function(name, default = "") {
    idx <- which(args == name)
    if (length(idx) == 0) return(default)
    if (idx == length(args)) stop("Missing value for ", name)
    args[idx + 1]
}

rw_offset_dir   <- opt("--rw_offset_dir")
rw_region_dir   <- opt("--rw_region_dir")
rq_psites_dir   <- opt("--rq_psites_dir")
rt_qual_dir     <- opt("--rt_qual_dir")
rtr_summary_dir <- opt("--rtr_summary_dir")
output_prefix   <- opt("--output_prefix", "joint_riboseq_qc")

# ── Helpers ───────────────────────────────────────────────────────────────────

sample_name <- function(filepath, suffixes) {
    nm <- basename(filepath)
    for (s in suffixes) if (grepl(paste0(s, "$"), nm)) return(sub(paste0(s, "$"), "", nm))
    sub("\\.[^.]*$", "", nm)
}

list_files <- function(dir_path, pattern = NULL) {
    if (dir_path == "" || !dir.exists(dir_path)) return(character(0))
    fl <- list.files(dir_path, pattern = pattern, full.names = TRUE)
    fl[file.info(fl)$size > 0]
}

parse_py_dict <- function(s) {
    # Parse simple Python dict {key: val, key: val, ...} with numeric keys
    s <- trimws(s)
    s <- sub("^\\{", "", s); s <- sub("\\}$", "", s)
    if (nchar(s) == 0) return(list())
    pairs <- strsplit(s, ",\\s*")[[1]]
    result <- list()
    for (p in pairs) {
        kv <- strsplit(p, ":\\s*")[[1]]
        if (length(kv) == 2) result[[trimws(kv[1])]] <- trimws(kv[2])
    }
    result
}

# ── 1. Load riboWaltz P-site offsets ─────────────────────────────────────────

load_rw <- function(dir_path) {
    files <- list_files(dir_path, "_psite_offset\\.tsv$")
    if (length(files) == 0) return(data.table())
    rbindlist(lapply(files, function(f) {
        s <- sample_name(f, c("_psite_offset.tsv", "_psite_offset.txt"))
        dt <- tryCatch(fread(f), error = function(e) NULL)
        if (is.null(dt) || nrow(dt) == 0) return(NULL)
        dt[, sample := s]
        if ("length" %in% names(dt) && !"read_length" %in% names(dt))
            setnames(dt, "length", "read_length")
        if ("corrected_offset_from_5" %in% names(dt))
            dt[, rw_offset := corrected_offset_from_5]
        else if ("offset_from_5" %in% names(dt))
            dt[, rw_offset := offset_from_5]
        cols <- intersect(names(dt), c("sample", "read_length", "rw_offset", "total_percentage"))
        dt[, ..cols]
    }), fill = TRUE)
}

# ── 2. Load riboWaltz region distribution ──────────────────────────────────

load_rw_region <- function(dir_path) {
    files <- list_files(dir_path, "_region_distribution\\.tsv$")
    if (length(files) == 0) return(data.table())
    rbindlist(lapply(files, function(f) {
        s <- sample_name(f, "_region_distribution.tsv")
        dt <- tryCatch(fread(f), error = function(e) NULL)
        if (is.null(dt) || nrow(dt) == 0) return(NULL)
        dt[, sample := s]
        setnames(dt, old = names(dt), tolower, skip_absent = TRUE)
        dt
    }), fill = TRUE)
}

# ── 3. Load RiboseQC P_sites_calcs ──────────────────────────────────────────

load_rq <- function(dir_path) {
    files <- list_files(dir_path, "_P_sites_calcs$")
    if (length(files) == 0) return(data.table())
    rbindlist(lapply(files, function(f) {
        s <- sample_name(f, "_P_sites_calcs")
        fl <- readLines(f, n = 1, warn = FALSE)
        if (length(fl) == 0 || grepl("^#|placeholder|insufficient|failed", fl, ignore.case = TRUE))
            return(NULL)
        dt <- tryCatch(fread(f), error = function(e) {
            tryCatch({
                obj <- readRDS(f)
                if (is.list(obj) && !is.data.frame(obj)) {
                    for (nm in names(obj))
                        if (is.data.frame(obj[[nm]]) && "max_coverage" %in% names(obj[[nm]]))
                            { obj <- obj[[nm]]; break }
                }
                as.data.table(obj)
            }, error = function(e2) NULL)
        })
        if (is.null(dt) || nrow(dt) == 0) return(NULL)
        dt[, sample := s]
        if ("rl" %in% names(dt)) setnames(dt, "rl", "read_length")
        dt
    }), fill = TRUE)
}

# ── 4. Load Ribo-TISH quality ───────────────────────────────────────────────

load_rt <- function(dir_path) {
    files <- list_files(dir_path, "_qual\\.txt$")
    if (length(files) == 0) return(data.table())
    rbindlist(lapply(files, function(f) {
        s <- sample_name(f, "_qual.txt")
        lines <- readLines(f, n = 1, warn = FALSE)
        if (length(lines) == 0) return(NULL)
        rl <- tryCatch(parse_py_dict(lines[1]), error = function(e) NULL)
        if (is.null(rl) || length(rl) == 0) return(NULL)
        data.table(
            sample = s,
            read_length = as.integer(names(rl)),
            rt_reads = as.integer(unlist(rl))
        )
    }), fill = TRUE)
}

# ── 5. Load Ribotricer bam_summary ──────────────────────────────────────────

load_rtr <- function(dir_path) {
    files <- list_files(dir_path, "_bam_summary\\.txt$")
    if (length(files) == 0) return(data.table())
    rbindlist(lapply(files, function(f) {
        s <- sample_name(f, "_bam_summary.txt")
        lines <- readLines(f, warn = FALSE)
        in_lens <- FALSE; meta <- list(); rlens <- integer(0); rcounts <- integer(0)
        for (line in lines) {
            line <- trimws(line)
            if (line == "") next
            if (line == "length dist:") { in_lens <- TRUE; next }
            if (line == "summary:") next
            parts <- strsplit(line, ":\\s*")[[1]]
            if (length(parts) < 2) next
            if (!in_lens) {
                meta[[trimws(parts[1])]] <- trimws(parts[2])
            } else {
                rlens <- c(rlens, as.integer(trimws(parts[1])))
                rcounts <- c(rcounts, as.integer(trimws(parts[2])))
            }
        }
        if (length(rlens) == 0) return(NULL)
        dt <- data.table(sample = s, read_length = rlens, rtr_reads = rcounts)
        dt[, rtr_total := as.integer(meta[["total_reads"]])]
        dt
    }), fill = TRUE)
}

# ── Load data ─────────────────────────────────────────────────────────────────

cat("=== JOINT QC REPORT ===\nLoading data...\n")
rw     <- load_rw(rw_offset_dir)
rw_reg <- load_rw_region(rw_region_dir)
rq     <- load_rq(rq_psites_dir)
rt     <- load_rt(rt_qual_dir)
rtr    <- load_rtr(rtr_summary_dir)

cat(sprintf("  riboWaltz offsets: %d rows, %d samples\n",
    nrow(rw), if (nrow(rw) > 0) uniqueN(rw$sample) else 0))
cat(sprintf("  riboWaltz regions: %d rows, %d samples\n",
    nrow(rw_reg), if (nrow(rw_reg) > 0) uniqueN(rw_reg$sample) else 0))
cat(sprintf("  RiboseQC psites:  %d rows, %d samples\n",
    nrow(rq), if (nrow(rq) > 0) uniqueN(rq$sample) else 0))
cat(sprintf("  Ribo-TISH qual:   %d rows, %d samples\n",
    nrow(rt), if (nrow(rt) > 0) uniqueN(rt$sample) else 0))
cat(sprintf("  Ribotricer bam:   %d rows, %d samples\n\n",
    nrow(rtr), if (nrow(rtr) > 0) uniqueN(rtr$sample) else 0))

if (nrow(rw) == 0 && nrow(rq) == 0) stop("No data loaded. Check input directories.")

# ── Per-sample metrics ───────────────────────────────────────────────────────

all_samples <- sort(unique(c(rw$sample, rq$sample, rt$sample, rtr$sample)))

metrics <- rbindlist(lapply(all_samples, function(s) {
    rw_s <- if (nrow(rw) > 0) rw[sample == s] else NULL
    rq_s <- if (nrow(rq) > 0) rq[sample == s] else NULL

    # Offset agreement (riboWaltz vs RiboseQC)
    rw_off <- rq_off <- NULL
    if (!is.null(rw_s) && nrow(rw_s) > 0 && "rw_offset" %in% names(rw_s))
        rw_off <- rw_s[, .(read_length = as.integer(read_length), rw_offset)]
    if (!is.null(rq_s) && nrow(rq_s) > 0) {
        if ("max_coverage" %in% names(rq_s))
            rq_off <- rq_s[max_coverage == TRUE | max_coverage == "TRUE",
                .(read_length = as.integer(read_length), rq_offset = cutoff, frame_pref = frame_preference)]
        else
            rq_off <- rq_s[, .(read_length = as.integer(read_length), rq_offset = cutoff, frame_pref = frame_preference)]
    }

    agree_core_n <- agree_core_tot <- NA_integer_
    if (!is.null(rw_off) && !is.null(rq_off) && nrow(rw_off) > 0 && nrow(rq_off) > 0) {
        cmp <- merge(rw_off, rq_off, by = "read_length", all = TRUE)
        core <- cmp[read_length >= 25 & read_length <= 31]
        agree_core_n <- sum(core$rw_offset == core$rq_offset, na.rm = TRUE)
        agree_core_tot <- sum(!is.na(core$rw_offset) & !is.na(core$rq_offset))
    }

    # CDS enrichment
    cds_enr <- NA_real_
    rw_reg_s <- if (nrow(rw_reg) > 0) rw_reg[sample == s] else NULL
    if (!is.null(rw_reg_s) && nrow(rw_reg_s) > 0) {
        cds_n  <- sum(rw_reg_s[psite_region == "cds" | psite_region == "CDS"]$n, na.rm = TRUE)
        utr5_n <- sum(rw_reg_s[psite_region == "5utr" | psite_region == "5UTR"]$n, na.rm = TRUE)
        utr3_n <- sum(rw_reg_s[psite_region == "3utr" | psite_region == "3UTR"]$n, na.rm = TRUE)
        cds_enr <- round(cds_n / max(1, utr5_n + utr3_n), 1)
    }

    # Best periodicity
    best_per  <- NA_real_
    best_plen <- NA_integer_
    if (!is.null(rq_off) && nrow(rq_off) > 0 && "frame_pref" %in% names(rq_off)) {
        fp <- rq_off[!is.na(frame_pref) & frame_pref > 0]
        if (nrow(fp) > 0) {
            br <- fp[which.max(frame_pref)]
            best_per <- br$frame_pref[1]; best_plen <- br$read_length[1]
        }
    }

    # Primary footprint length
    pri_len <- NA_integer_; pri_pct <- NA_real_
    if (!is.null(rw_s) && nrow(rw_s) > 0 && "total_percentage" %in% names(rw_s)) {
        pos <- rw_s[!is.na(total_percentage) & total_percentage > 0]
        if (nrow(pos) > 0) {
            bp <- pos[which.max(total_percentage)]
            pri_len <- as.integer(bp$read_length[1]); pri_pct <- bp$total_percentage[1]
        }
    }

    # Total reads (prefer Ribotricer)
    tot_reads <- NA_integer_; n_rl <- NA_integer_
    if (!is.null(rtr) && nrow(rtr) > 0) {
        rtr_s <- rtr[sample == s]
        if (nrow(rtr_s) > 0) { tot_reads <- rtr_s$rtr_total[1]; n_rl <- uniqueN(rtr_s$read_length) }
    }
    if (is.na(tot_reads) && !is.null(rt) && nrow(rt) > 0) {
        rt_s <- rt[sample == s]
        if (nrow(rt_s) > 0) { tot_reads <- sum(rt_s$rt_reads, na.rm = TRUE); n_rl <- uniqueN(rt_s$read_length) }
    }

    # Grades
    grade <- function(val, cuts) {
        if (is.na(val) || val <= 0) return(NA_character_)
        if (val >= cuts[1]) "A" else if (val >= cuts[2]) "B" else if (val >= cuts[3]) "C" else "D"
    }
    o_grade <- if (!is.na(agree_core_tot) && agree_core_tot > 0)
        grade(agree_core_n / agree_core_tot, c(0.9, 0.7, 0.5)) else NA_character_
    p_grade <- grade(best_per, c(70, 50, 30))
    c_grade <- grade(cds_enr, c(5, 2, 1))

    scores <- c(
        if (!is.na(o_grade)) match(o_grade, LETTERS[1:4]) else NA_integer_,
        if (!is.na(p_grade)) match(p_grade, LETTERS[1:4]) else NA_integer_,
        if (!is.na(c_grade)) match(c_grade, LETTERS[1:4]) else NA_integer_
    )
    avg_score <- mean(scores, na.rm = TRUE)
    overall <- if (is.na(avg_score)) NA_character_
        else if (avg_score <= 1.5) "A" else if (avg_score <= 2.5) "B" else "C"

    data.table(
        Sample = s,
        TotalReads = tot_reads,
        PrimaryLength = pri_len, PrimaryPct = pri_pct,
        CDSEnrichment = cds_enr,
        BestPeriodicity = best_per, BestPeriodLen = best_plen,
        OffsetAgreeCore = agree_core_n, CoreLenCompared = agree_core_tot,
        NReadLengths = n_rl,
        OffsetGrade = o_grade, PeriodGrade = p_grade, CDSGrade = c_grade,
        OverallGrade = overall
    )
}))

metrics <- metrics[order(Sample)]

# ── Print summary table ──────────────────────────────────────────────────────

cat(paste(rep("=", 105), collapse = ""), "\n")
cat("  JOINT RIBO-SEQ QC SUMMARY\n")
cat(paste(rep("=", 105), collapse = ""), "\n\n")

cat(sprintf("%-18s %8s %5s %7s %6s %7s %7s %5s %5s %5s %6s\n",
    "Sample", "Reads", "Len", "Pct%", "CDSx", "Per%", "P-Len", "Agr%", "OG", "PG", "CG"))
cat(paste(rep("-", 105), collapse = ""), "\n")
for (i in seq_len(nrow(metrics))) {
    r <- metrics[i]
    ap <- if (!is.na(r$CoreLenCompared) && r$CoreLenCompared > 0)
        round(100 * r$OffsetAgreeCore / r$CoreLenCompared, 0) else NA
    cat(sprintf("%-18s %8s %5s %6.1f %5.1f %6.1f %6s %4s%% %5s %5s %5s %5s\n",
        r$Sample,
        if (is.na(r$TotalReads)) "-" else format(r$TotalReads, big.mark = ",", scientific = FALSE),
        if (is.na(r$PrimaryLength)) "-" else r$PrimaryLength,
        if (is.na(r$PrimaryPct)) 0 else r$PrimaryPct,
        if (is.na(r$CDSEnrichment)) 0 else r$CDSEnrichment,
        if (is.na(r$BestPeriodicity)) 0 else r$BestPeriodicity,
        if (is.na(r$BestPeriodLen)) "-" else r$BestPeriodLen,
        if (is.na(ap)) "-" else as.character(ap),
        if (is.na(r$OffsetGrade)) "-" else r$OffsetGrade,
        if (is.na(r$PeriodGrade)) "-" else r$PeriodGrade,
        if (is.na(r$CDSGrade)) "-" else r$CDSGrade,
        if (is.na(r$OverallGrade)) "-" else r$OverallGrade))
}

# ── Overall assessment ───────────────────────────────────────────────────────

cat("\n", paste(rep("=", 60), collapse = ""), "\n")
cat("  OVERALL ASSESSMENT\n")
cat(paste(rep("=", 60), collapse = ""), "\n\n")

ga <- if (sum(metrics$CoreLenCompared, na.rm = TRUE) > 0)
    round(100 * sum(metrics$OffsetAgreeCore, na.rm = TRUE) /
          sum(metrics$CoreLenCompared, na.rm = TRUE), 0) else NA
cds_m <- mean(metrics$CDSEnrichment, na.rm = TRUE)
per_m <- mean(metrics$BestPeriodicity[metrics$BestPeriodLen %in% c(29, 30)], na.rm = TRUE)

cat(sprintf("  Core offset agreement (25-31nt): %s%%\n",
    if (is.na(ga)) "-" else as.character(ga)))
cat(sprintf("  Mean CDS enrichment: %.1fx\n", if (is.na(cds_m)) 0 else cds_m))
cat(sprintf("  Best periodicity (29-30nt): %.1f%%\n", if (is.na(per_m)) 0 else per_m))

gc <- metrics[, .N, by = OverallGrade][order(OverallGrade)]
cat(sprintf("  Grade distribution: %s\n",
    paste(sprintf("%s=%d", gc$OverallGrade, gc$N), collapse = ", ")))

# ── Write outputs ────────────────────────────────────────────────────────────

tsv_file <- paste0(output_prefix, ".tsv")
fwrite(metrics, file = tsv_file, sep = "\t", quote = FALSE, na = "-")
cat(sprintf("\nWrote: %s\n", tsv_file))

# MultiQC custom content YAML
mqc_yaml <- paste0(output_prefix, "_mqc.yaml")
writeLines(c(
    "id: 'joint_riboseq_qc'",
    "section_name: 'Joint Ribo-seq QC'",
    "description: 'Integrated quality assessment combining riboWaltz, RiboseQC, Ribo-TISH, and Ribotricer.'",
    "plot_type: 'table'",
    "section_href: 'https://github.com/nf-core/riboseq'",
    "pconfig:",
    "  id: 'joint_riboseq_qc_table'",
    "  title: 'Joint Ribo-seq QC: Per-Sample Quality Metrics'",
    "  namespace: 'Joint Ribo-seq QC'",
    "  format: '{:.1f}'",
    "  col1_header: 'Sample'",
    "  only_defined_headers: false",
    "  table_columns_visible:",
    "    TotalReads: true",
    "    PrimaryLength: true",
    "    PrimaryPct: true",
    "    CDSEnrichment: true",
    "    BestPeriodicity: true",
    "    BestPeriodLen: true",
    "    OffsetAgreeCore: true",
    "    CoreLenCompared: true",
    "    OverallGrade: true",
    "  table_columns_placement:",
    "    OverallGrade: 1000",
    "    PrimaryLength: 500"
), mqc_yaml)
cat(sprintf("Wrote: %s\n", mqc_yaml))

# MultiQC data file
display_names <- c(
    TotalReads = "Total Reads",
    PrimaryLength = "Primary Length (nt)",
    PrimaryPct = "Primary Length %",
    CDSEnrichment = "CDS Enrichment (x)",
    BestPeriodicity = "Best Periodicity (%)",
    BestPeriodLen = "Period Length (nt)",
    OffsetAgreeCore = "Offset Agreement (core)",
    CoreLenCompared = "Core Lengths Compared",
    NReadLengths = "N Read Lengths",
    OffsetGrade = "Offset Grade",
    PeriodGrade = "Periodicity Grade",
    CDSGrade = "CDS Grade",
    OverallGrade = "Overall Grade"
)
mqc_df <- copy(metrics)
for (old in names(display_names)) {
    if (old %in% names(mqc_df)) setnames(mqc_df, old, display_names[old], skip_absent = TRUE)
}
mqc_data <- paste0(output_prefix, "_mqc.txt")
fwrite(mqc_df, file = mqc_data, sep = "\t", quote = FALSE, na = "-")
cat(sprintf("Wrote: %s\n", mqc_data))

cat("\n===== JOINT QC REPORT COMPLETE =====\n")
