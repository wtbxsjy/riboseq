#' Classify ORFs using ORFquant genomic categories
#'
#' This function annotates ORFs with ORFquant-compatible genomic categories
#' (e.g. exact_start_stop, Alt5_start, novel_Upstream) using CDS annotation.
#' It accepts ORF coordinates from data.frame, GRanges/GRangesList, or BED/GTF/GFF.
#'
#' @param orfs ORF coordinates. Supported: GRanges, GRangesList, data.frame, or a
#'             file path (BED/GTF/GFF) importable by rtracklayer::import.
#' @param annotation CDS annotation. Supported: list with $cds_genes and $cds_txs
#'                   (as in ORFquant), or a file path (GTF/GFF) containing CDS features.
#' @param gene_id_col Column name for gene id (for data.frame) or metadata column name
#'                    (for GRanges). Defaults to "gene_id".
#' @param transcript_id_col Column name for transcript id (for data.frame) or metadata
#'                    column name (for GRanges). Defaults to "transcript_id".
#' @param seqname_col Column name for seqname (for data.frame). Defaults to "seqnames".
#' @param start_col Column name for start (for data.frame). Defaults to "start".
#' @param end_col Column name for end (for data.frame). Defaults to "end".
#' @param strand_col Column name for strand (for data.frame). Defaults to "strand".
#' @param exons_col Optional column containing exon coordinates as
#'        "start-end,start-end". If provided, multi-exon ORFs are built from this column.
#' @param orf_id_col Optional ORF id column; if missing, row index is used.

#' Project genomic ORF blocks to 1-based transcript-space coordinates.
#'
#' Handles multi-exon ORFs correctly by walking through the exon chain.
#' For "+" strand the transcript runs 5'→3' from low to high genomic coordinate;
#' for "-" strand it runs from high to low.
#'
#' @param orf_gr  GRanges of ORF CDS exons (1-based genomic coordinates).
#' @param exon_gr GRanges of transcript exons (1-based genomic coordinates).
#' @param strand  Character "+" or "-".
#' @return Integer vector c(tx_start, tx_stop) 1-based (tx_start <= tx_stop),
#'         or NULL when no exon overlap is found.
#' @export
project_to_tx_coords <- function(orf_gr, exon_gr, strand) {
  if (is.null(exon_gr) || length(exon_gr) == 0 ||
      is.null(orf_gr)  || length(orf_gr)  == 0) return(NULL)

  ord  <- order(IRanges::start(exon_gr))
  ex_s <- IRanges::start(exon_gr)[ord]
  ex_e <- IRanges::end(exon_gr)[ord]
  n_ex <- length(ex_s)

  # For "-" strand the 5' end is at the highest genomic coordinate, so we
  # traverse exons from right to left.
  if (identical(strand, "-")) {
    ex_s <- rev(ex_s)
    ex_e <- rev(ex_e)
  }

  ex_len     <- ex_e - ex_s + 1L
  cum_before <- c(0L, cumsum(ex_len[-n_ex]))  # cumulative tx length before exon j

  # Map one genomic position to a 1-based transcript position.
  gpos_to_tx <- function(gpos) {
    # Vectorized: find which exon contains gpos
    in_exon <- which(gpos >= ex_s & gpos <= ex_e)
    if (length(in_exon) == 0L) return(NA_integer_)
    j <- in_exon[1L]
    offset <- if (identical(strand, "-")) (ex_e[j] - gpos) else (gpos - ex_s[j])
    cum_before[j] + offset + 1L
  }

  tx_pos <- vapply(c(IRanges::start(orf_gr), IRanges::end(orf_gr)), gpos_to_tx, integer(1))
  tx_pos <- tx_pos[!is.na(tx_pos)]
  if (length(tx_pos) == 0) return(NULL)
  c(min(tx_pos), max(tx_pos))
}

#' @return data.frame with ORF_category_Gen and (when possible) ORF_category_Tx / ORF_category_Tx_compatible.
#' @export

orfquant_classify_orfs <- function(orfs,
                                   annotation,
                                   gene_id_col = "gene_id",
                                   transcript_id_col = "transcript_id",
                                   seqname_col = "seqnames",
                                   start_col = "start",
                                   end_col = "end",
                                   strand_col = "strand",
                                   exons_col = NULL,
                                   orf_id_col = NULL) {
  if (!requireNamespace("GenomicRanges", quietly = TRUE)) {
    stop("GenomicRanges is required.")
  }
  if (!requireNamespace("IRanges", quietly = TRUE)) {
    stop("IRanges is required.")
  }
  if (!requireNamespace("rtracklayer", quietly = TRUE)) {
    stop("rtracklayer is required.")
  }

  normalize_orfs <- function(x) {
    if (is.character(x) && length(x) == 1 && file.exists(x)) {
      x <- rtracklayer::import(x)
    }

    if (inherits(x, "GRangesList")) {
      grl <- x
      # Propagate per-range metadata to list level if missing.
      # When split(gr, gr$orf_id) is used, only constant-per-group columns
      # (like orf_id itself) are promoted to list-level mcols; gene_id/
      # transcript_id may vary within a group and stay at range level only.
      inner_gr  <- unlist(grl, use.names = TRUE)
      inner_meta <- S4Vectors::mcols(inner_gr)
      el_idx    <- rep(seq_along(grl), lengths(grl))
      for (col in c("gene_id", "transcript_id")) {
        if (!col %in% colnames(S4Vectors::mcols(grl))) {
          if (col %in% colnames(inner_meta)) {
            vals <- as.character(inner_meta[[col]])
            S4Vectors::mcols(grl)[[col]] <- vapply(seq_along(grl), function(i) {
              v <- vals[el_idx == i]
              if (length(v) == 0 || all(is.na(v))) NA_character_ else v[!is.na(v)][1]
            }, character(1))
          } else {
            S4Vectors::mcols(grl)[[col]] <- NA_character_
          }
        }
      }
      # Propagate orf_id from GRangesList names (= split key) before falling back to seq index.
      if (!"orf_id" %in% colnames(S4Vectors::mcols(grl)) && !is.null(names(grl))) {
        S4Vectors::mcols(grl)$orf_id <- names(grl)
      }
    } else if (inherits(x, "GRanges")) {
      grl <- GenomicRanges::GRangesList(split(x, seq_len(length(x))))
    } else if (is.data.frame(x)) {
      if (!is.null(exons_col) && exons_col %in% colnames(x)) {
        grl <- GenomicRanges::GRangesList(lapply(seq_len(nrow(x)), function(i) {
          exon_str <- x[[exons_col]][i]
          if (is.na(exon_str) || exon_str == "") {
            return(GenomicRanges::GRanges())
          }
          parts <- unlist(strsplit(exon_str, ","))
          coords <- lapply(parts, function(p) {
            se <- unlist(strsplit(p, "-"))
            if (length(se) != 2) return(NULL)
            c(as.integer(se[1]), as.integer(se[2]))
          })
          coords <- coords[!vapply(coords, is.null, logical(1))]
          if (length(coords) == 0) {
            return(GenomicRanges::GRanges())
          }
          starts <- vapply(coords, function(v) v[1], numeric(1))
          ends <- vapply(coords, function(v) v[2], numeric(1))
          GenomicRanges::GRanges(
            seqnames = x[[seqname_col]][i],
            ranges = IRanges::IRanges(start = starts, end = ends),
            strand = x[[strand_col]][i]
          )
        }))
      } else {
        gr <- GenomicRanges::GRanges(
          seqnames = x[[seqname_col]],
          ranges = IRanges::IRanges(start = x[[start_col]], end = x[[end_col]]),
          strand = x[[strand_col]]
        )
        grl <- GenomicRanges::GRangesList(split(gr, seq_len(length(gr))))
      }

      mcols_grl <- S4Vectors::DataFrame(
        gene_id = if (gene_id_col %in% colnames(x)) x[[gene_id_col]] else NA_character_,
        transcript_id = if (transcript_id_col %in% colnames(x)) x[[transcript_id_col]] else NA_character_,
        orf_id = if (!is.null(orf_id_col) && orf_id_col %in% colnames(x)) x[[orf_id_col]] else as.character(seq_len(nrow(x)))
      )
      S4Vectors::mcols(grl) <- mcols_grl
    } else {
      stop("Unsupported ORF input. Use GRanges/GRangesList, data.frame, or BED/GTF/GFF path.")
    }

    # Ensure metadata on GRangesList
    if (is.null(S4Vectors::mcols(grl)) || nrow(S4Vectors::mcols(grl)) == 0) {
      S4Vectors::mcols(grl) <- S4Vectors::DataFrame(
        gene_id = NA_character_,
        transcript_id = NA_character_,
        orf_id = as.character(seq_len(length(grl)))
      )
    } else if (!"orf_id" %in% colnames(S4Vectors::mcols(grl))) {
      S4Vectors::mcols(grl)$orf_id <- as.character(seq_len(length(grl)))
    }

    return(grl)
  }

  normalize_annotation <- function(a) {
    # Allow a pre-built annotation list only if it already has all required fields.
    if (is.list(a) && !is.null(a$cds_genes) && !is.null(a$cds_txs) &&
        !is.null(a$exon_txs)) {
      return(a)
    }
    if (is.character(a) && length(a) == 1 && file.exists(a)) {
      # Check for cached RDS annotation to avoid repeated slow GTF imports
      rds_cache <- paste0(a, ".ann_cache.rds")
      if (file.exists(rds_cache) &&
          file.mtime(rds_cache) >= file.mtime(a)) {
        message("  [normalize_annotation] Loading cached annotation from: ", rds_cache)
        return(readRDS(rds_cache))
      }
      message("  [normalize_annotation] Importing GTF: ", a, " (", format(file.size(a), big.mark=","), " bytes)")
      ann <- rtracklayer::import(a)
      message("  [normalize_annotation] Building CDS gene/transcript structures...")
      cds <- ann[ann$type == "CDS"]
      if (length(cds) == 0) stop("No CDS features found in annotation.")
      if (is.null(cds$gene_id))    stop("CDS annotation must contain gene_id.")
      if (is.null(cds$transcript_id)) cds$transcript_id <- cds$gene_id

      cds_genes <- GenomicRanges::reduce(split(cds, cds$gene_id))
      cds_txs   <- split(cds, cds$transcript_id)

      # Build per-transcript exon structure for coordinate projection.
      exon <- ann[ann$type == "exon"]
      exon_txs <- if (!is.null(exon$transcript_id) && length(exon) > 0) {
        split(exon, exon$transcript_id)
      } else {
        cds_txs  # fall back to CDS exons as proxy
      }

      # Build lncRNA gene set from 'gene' features in the annotation.
      # Any gene with a biotype in the LNCRNA_BIOTYPES set is recorded so
      # that ORFs whose gene_id is a non-coding gene can be classified as
      # "lncRNA" rather than falling back to a positional category.
      LNCRNA_BIOTYPES <- c(
        # Human / mouse (Ensembl / GENCODE)
        "lncRNA", "lincRNA", "antisense", "processed_transcript",
        "sense_intronic", "sense_overlapping", "non_coding",
        "3prime_overlapping_ncrna", "bidirectional_promoter_lncrna",
        # Plant species (Ensembl Plants)
        "ncRNA",            # Oryza sativa (rice): general long non-coding RNA
        "antisense_RNA",    # Oryza sativa: antisense long non-coding RNA
        "misc_non_coding"   # Zea mays (maize): catch-all non-coding gene category
      )
      gene_feat <- ann[ann$type == "gene"]
      lncrna_genes <- character(0)
      if (length(gene_feat) > 0 && !is.null(gene_feat$gene_biotype)) {
        lnc_mask <- gene_feat$gene_biotype %in% LNCRNA_BIOTYPES
        if (!is.null(gene_feat$gene_id)) {
          lncrna_genes <- unique(gene_feat$gene_id[lnc_mask])
        }
        message("  [normalize_annotation] ", length(lncrna_genes),
                " lncRNA genes identified.")
      }

      # Tx-coord projection is computed ON-DEMAND during classification
      # rather than pre-computed for all 66K transcripts (which takes 10+ min).
      # Only transcripts that have ORFs mapped to them need projection.
      # The classify_* functions check for NULL cds_txs_tx_coords and
      # call project_to_tx_coords() on-the-fly as needed.
      message("  [normalize_annotation] Skipping tx-coord pre-projection (on-demand)")
      cds_txs_tx_coords <- list()  # placeholder, filled on-demand

      # ── Pre-compute derived structures (saves ~minutes per run) ─────────
      message("  [normalize_annotation] Building derived lookup tables...")
      t_derived <- Sys.time()

      # gene → transcript_id map (split-based, O(n log n) not O(n²))
      # gene → transcript_id map (unlisted, O(n) single pass)
      flat_cds <- unlist(cds_txs)
      part_idx <- rep(seq_along(cds_txs), lengths(cds_txs))
      first_occ <- !duplicated(part_idx)
      tx_gene_v <- setNames(
        as.character(S4Vectors::mcols(flat_cds)$gene_id)[first_occ],
        names(cds_txs))
      gene_to_txids <- split(names(tx_gene_v), tx_gene_v)

      # Max CDS per gene (pre-compute widths, O(n) single pass)
      idx_by_gene <- split(seq_along(cds_txs), tx_gene_v)
      # Compute CDS widths for ALL transcripts in one pass (unlisted = fast)
      tx_widths <- tapply(width(flat_cds), part_idx, sum)
      max_cds_by_gene <- lapply(names(idx_by_gene), function(g) {
        idxs <- idx_by_gene[[g]]
        if (length(idxs) == 0) return(NULL)
        best <- idxs[which.max(tx_widths[idxs])]
        if (length(best) == 0) return(NULL)
        cds_txs[[best]]
      })
      names(max_cds_by_gene) <- names(idx_by_gene)
      max_cds_by_gene <- max_cds_by_gene[!vapply(max_cds_by_gene, is.null, logical(1))]

      # cds_genes is already reduced (created by GenomicRanges::reduce)
      # No need to re-reduce — just reference it directly
      cds_genes_reduced <- cds_genes  # already reduced

      # Merged all-CDS (for ORFs without a matching gene)
      all_cds_merged <- IRanges::reduce(unlist(cds_genes))

      message("  [normalize_annotation] Derived structures built in ",
              round(difftime(Sys.time(), t_derived, units = "secs"), 1), "s")

      result <- list(cds_genes = cds_genes, cds_txs = cds_txs,
                     exon_txs = exon_txs, cds_txs_tx_coords = cds_txs_tx_coords,
                     lncrna_genes = lncrna_genes,
                     gene_to_txids = gene_to_txids,
                     max_cds_by_gene = max_cds_by_gene,
                     cds_genes_reduced = cds_genes_reduced,
                     all_cds_merged = all_cds_merged)

      # Save cache for future use
      tryCatch({
        saveRDS(result, rds_cache)
        message("  [normalize_annotation] Cached annotation to: ", rds_cache)
      }, error = function(e) {
        message("  [normalize_annotation] Warning: could not save cache: ", e$message)
      })

      return(result)
    }
  }

  get_max_cds_by_gene <- function(cds_txs) {
    tx_gene <- vapply(cds_txs, function(x) {
      gid <- unique(x$gene_id)
      if (length(gid) == 0) NA_character_ else gid[1]
    }, character(1))

    gene_ids <- unique(tx_gene[!is.na(tx_gene)])
    max_cds <- list()
    for (g in gene_ids) {
      txs <- cds_txs[tx_gene == g]
      if (length(txs) == 0) next
      tx_len <- vapply(txs, function(x) sum(IRanges::width(IRanges::reduce(x))), numeric(1))
      max_cds[[g]] <- txs[[which.max(tx_len)]]
    }
    return(max_cds)
  }

  classify_genomic <- function(orf_gen, cds_gene, max_cdsok) {
    if (length(cds_gene) == 0 || length(orf_gen) == 0) {
      return("novel")
    }

    overl <- orf_gen %over% cds_gene

    if (sum(overl) == 0) {
      category <- "novel"

      nearest_cds <- cds_gene[GenomicRanges::nearest(orf_gen, cds_gene)[1]]
      overl_whole <- orf_gen@ranges %over% IRanges::IRanges(
        start = min(IRanges::start(nearest_cds)),
        end = max(IRanges::end(nearest_cds))
      )

      strd <- as.vector(GenomicRanges::strand(orf_gen[1]))
      if (strd == "+") {
        if (sum(overl_whole) == 0) {
          if (min(IRanges::start(orf_gen)) < min(IRanges::start(nearest_cds))) category <- "novel_Upstream"
          if (min(IRanges::start(orf_gen)) > max(IRanges::end(nearest_cds))) category <- "novel_Downstream"
        }
        if (sum(overl_whole) > 0) category <- "novel_Internal"
      }
      if (strd == "-") {
        if (sum(overl_whole) == 0) {
          if (max(IRanges::end(orf_gen)) > max(IRanges::end(nearest_cds))) category <- "novel_Upstream"
          if (min(IRanges::start(orf_gen)) < min(IRanges::start(nearest_cds))) category <- "novel_Downstream"
        }
        if (sum(overl_whole) > 0) category <- "novel_Internal"
      }
      return(category)
    }

    # overlaps CDS
    if (length(max_cdsok) == 0) {
      return("overlaps_CDS")
    }

    strd <- as.vector(GenomicRanges::strand(orf_gen[1]))
    if (strd == "+") {
      gen_sta <- min(IRanges::start(max_cdsok))
      gen_sto <- max(IRanges::end(max_cdsok))
      sta_or <- min(IRanges::start(orf_gen))
      sto_or <- max(IRanges::end(orf_gen)) + 3

      if (sto_or == gen_sto) {
        if (sta_or == gen_sta) return("exact_start_stop")
        if (sta_or < gen_sta) return("Alt5_start")
        if (sta_or > gen_sta) return("Alt3_start")
      }

      if (sto_or != gen_sto) {
        if (sta_or < gen_sta && sto_or < gen_sto) return("Alt5_start_Alt5_stop")
        if (sta_or < gen_sta && sto_or > gen_sto) return("Alt5_start_Alt3_stop")
        if (sta_or > gen_sta && sto_or > gen_sto) return("Alt3_start_Alt3_stop")
        if (sta_or > gen_sta && sto_or < gen_sto) return("Alt3_start_Alt5_stop")
        if (sta_or == gen_sta && sto_or < gen_sto) return("Alt5_stop")
        if (sta_or == gen_sta && sto_or > gen_sto) return("Alt3_stop")
      }
    }

    if (strd == "-") {
      gen_sta <- max(IRanges::end(max_cdsok))
      gen_sto <- min(IRanges::start(max_cdsok))
      sta_or <- max(IRanges::end(orf_gen))
      sto_or <- min(IRanges::start(orf_gen)) - 3

      if (sto_or == gen_sto) {
        if (sta_or == gen_sta) return("exact_start_stop")
        if (sta_or < gen_sta) return("Alt3_start")
        if (sta_or > gen_sta) return("Alt5_start")
      }

      if (sto_or != gen_sto) {
        if (sta_or > gen_sta && sto_or > gen_sto) return("Alt5_start_Alt5_stop")
        if (sta_or > gen_sta && sto_or < gen_sto) return("Alt5_start_Alt3_stop")
        if (sta_or < gen_sta && sto_or < gen_sto) return("Alt3_start_Alt3_stop")
        if (sta_or < gen_sta && sto_or > gen_sto) return("Alt3_start_Alt5_stop")
        if (sta_or == gen_sta && sto_or > gen_sto) return("Alt5_stop")
        if (sta_or == gen_sta && sto_or < gen_sto) return("Alt3_stop")
      }
    }

    return("overlaps_CDS")
  }

  classify_transcript <- function(orf_sta, orf_sto, ann_sta, ann_sto) {
    if (is.na(ann_sta) || is.na(ann_sto)) return(NA_character_)

    if (orf_sto == ann_sto) {
      if (orf_sta == ann_sta) return("ORF_annotated")
      if (orf_sta < ann_sta) return("N_extension")
      if (orf_sta > ann_sta) return("N_truncation")
    }

    if (orf_sto != ann_sto) {
      # Check more-specific (non-overlapping) cases BEFORE the overlapping cases.
      # uORF: entirely upstream of CDS start; overl_uORF: overlaps CDS 5' end.
      # dORF: entirely downstream of CDS stop; overl_dORF: overlaps CDS 3' end.
      if (orf_sta < ann_sta && orf_sto < ann_sta) return("uORF")
      if (orf_sta < ann_sta && orf_sto < ann_sto) return("overl_uORF")
      if (orf_sta < ann_sta && orf_sto > ann_sto) return("NC_extension")
      if (orf_sta > ann_sto && orf_sto > ann_sto) return("dORF")
      if (orf_sta > ann_sta && orf_sto > ann_sto) return("overl_dORF")
      if (orf_sta > ann_sta && orf_sto < ann_sto) return("nested_ORF")
      if (orf_sta == ann_sta && orf_sto < ann_sto) return("C_truncation")
      if (orf_sta == ann_sta && orf_sto > ann_sto) return("C_extension")
    }

    return(NA_character_)
  }

  orf_grl <- normalize_orfs(orfs)
  ann <- normalize_annotation(annotation)

  cds_genes <- ann$cds_genes
  cds_txs <- ann$cds_txs

  # Use pre-computed derived structures when available (from RDS cache),
  # falling back to on-the-fly computation for ad-hoc annotation objects.
  if (!is.null(ann$max_cds_by_gene)) {
    max_cds_by_gene <- ann$max_cds_by_gene
  } else {
    max_cds_by_gene <- get_max_cds_by_gene(cds_txs)
  }

  if (!is.null(ann$cds_genes_reduced)) {
    cds_genes_reduced <- ann$cds_genes_reduced
  } else {
    # cds_genes is already GenomicRanges::reduce()'d
    cds_genes_reduced <- cds_genes
  }

  if (!is.null(ann$all_cds_merged)) {
    all_cds <- ann$all_cds_merged
  } else {
    all_cds <- IRanges::reduce(unlist(cds_genes))
  }

  # lncRNA gene set from annotation (may be empty if annotation or GTF lacks gene features)
  lncrna_genes <- if (!is.null(ann$lncrna_genes)) ann$lncrna_genes else character(0)

  # On-demand tx-coord projection: compute and cache for a single transcript.
  # Avoids pre-computing for all 66K transcripts (saves 10+ min at startup).
  get_tx_coords <- function(tx_id) {
    if (is.null(ann$cds_txs_tx_coords)) {
      ann$cds_txs_tx_coords <<- list()
    }
    if (is.null(ann$cds_txs_tx_coords[[tx_id]])) {
      cds_gr  <- cds_txs[[tx_id]]
      exon_gr <- if (tx_id %in% names(ann$exon_txs)) ann$exon_txs[[tx_id]] else cds_gr
      if (length(cds_gr) == 0 || length(exon_gr) == 0) {
        ann$cds_txs_tx_coords[[tx_id]] <<- IRanges::IRanges()
      } else {
        strd <- as.character(GenomicRanges::strand(cds_gr)[1])
        tx_c <- project_to_tx_coords(cds_gr, exon_gr, strd)
        if (is.null(tx_c)) {
          ann$cds_txs_tx_coords[[tx_id]] <<- IRanges::IRanges()
        } else {
          ann$cds_txs_tx_coords[[tx_id]] <<- IRanges::IRanges(start = tx_c[1], end = tx_c[2])
        }
      }
    }
    ann$cds_txs_tx_coords[[tx_id]]
  }

  # Priority order for ORF_category_Tx_compatible: lower number = better match.
  TX_CLASS_PRIORITY <- c(
    ORF_annotated    = 1L,
    N_extension      = 2L,
    C_extension      = 3L,
    N_truncation     = 4L,
    C_truncation     = 5L,
    NC_extension     = 6L,
    overl_uORF       = 7L,
    overl_dORF       = 8L,
    nested_ORF       = 9L,
    uORF             = 10L,
    dORF             = 11L,
    novel            = 12L,
    novel_antisense  = 13L
  )

  # Build gene → transcript_id lookup (used for ORF_category_Tx_compatible)
  # Use pre-computed version when available, otherwise build on the fly.
  gene_to_txids <- if (!is.null(ann$gene_to_txids)) {
    ann$gene_to_txids
  } else if (!is.null(ann$cds_txs_tx_coords)) {
    tx_gene_map <- vapply(cds_txs, function(tx) {
      gids <- unique(tx$gene_id)
      if (length(gids) > 0) as.character(gids[1]) else NA_character_
    }, character(1))
    split(names(tx_gene_map), tx_gene_map)
  } else {
    list()
  }

  seqnames_vec <- vapply(orf_grl, function(x) {
    if (length(x) == 0) return(NA_character_)
    as.character(GenomicRanges::seqnames(x)[1])
  }, character(1))
  start_vec <- vapply(orf_grl, function(x) {
    if (length(x) == 0) return(NA_integer_)
    min(IRanges::start(x))
  }, integer(1))
  end_vec <- vapply(orf_grl, function(x) {
    if (length(x) == 0) return(NA_integer_)
    max(IRanges::end(x))
  }, integer(1))
  strand_vec <- vapply(orf_grl, function(x) {
    if (length(x) == 0) return(NA_character_)
    as.character(GenomicRanges::strand(x)[1])
  }, character(1))

  res <- data.frame(
    orf_id = S4Vectors::mcols(orf_grl)$orf_id,
    gene_id = S4Vectors::mcols(orf_grl)$gene_id,
    seqnames = seqnames_vec,
    start = start_vec,
    end = end_vec,
    strand = strand_vec,
    ORF_category_Gen = NA_character_,
    ORF_category_Tx = NA_character_,
    ORF_category_Tx_compatible = NA_character_,
    stringsAsFactors = FALSE
  )

  for (i in seq_along(orf_grl)) {
    orf_gen <- orf_grl[[i]]
    gene_id <- S4Vectors::mcols(orf_grl)$gene_id[i]

    # If the ORF's gene is a known lncRNA gene, assign "lncRNA" immediately
    # and skip the CDS-based positional classification (which is meaningless
    # for non-coding genes).
    if (!is.na(gene_id) && length(lncrna_genes) > 0 && gene_id %in% lncrna_genes) {
      res$ORF_category_Gen[i] <- "lncRNA"
      res$ORF_category_Tx[i]  <- "lncRNA"
      res$ORF_category_Tx_compatible[i] <- "lncRNA"
      next
    }

    if (!is.na(gene_id) && gene_id %in% names(cds_genes)) {
      cds_gene <- cds_genes_reduced[[gene_id]]
      max_cdsok <- max_cds_by_gene[[gene_id]]
    } else {
      cds_gene <- all_cds
      max_cdsok <- GenomicRanges::GRanges()
    }

    res$ORF_category_Gen[i] <- classify_genomic(orf_gen, cds_gene, max_cdsok)

    # Transcript-level categories using proper transcript-space projection.
    # project_to_tx_coords() maps the ORF's genomic blocks through the
    # transcript exon chain to get 1-based transcript coordinates.  This is
    # correct for both single-exon and multi-exon ORFs and matches the
    # coordinate system used natively by ORFquant (ORF_id_tr format:
    # transcript_id_txstart_txstop).
    if (!is.null(ann$cds_txs_tx_coords) && !is.null(ann$exon_txs)) {
      tx_id <- S4Vectors::mcols(orf_grl)$transcript_id[i]

      # ── ORF_category_Tx: classify against the ORF's own transcript ──────────
      if (!is.na(tx_id) && tx_id %in% names(ann$exon_txs)) {
        orf_tx <- project_to_tx_coords(orf_grl[[i]], ann$exon_txs[[tx_id]], strand_vec[i])
        if (!is.null(orf_tx) && tx_id %in% names(ann$cds_txs_tx_coords)) {
          ref_tx <- get_tx_coords(tx_id)
          if (length(ref_tx) > 0) {
            res$ORF_category_Tx[i] <- classify_transcript(
              orf_tx[1], orf_tx[2],
              IRanges::start(ref_tx), IRanges::end(ref_tx))
          } else {
            res$ORF_category_Tx[i] <- "novel"
          }
        } else if (!is.null(orf_tx)) {
          # transcript exists in annotation but has no CDS (non-coding) → novel
          res$ORF_category_Tx[i] <- "novel"
        }
      }

      # ── ORF_category_Tx_compatible: best across all transcripts of the gene ─
      # For each annotated transcript of the gene, project the ORF into that
      # transcript's coordinate system and keep the highest-priority category.
      if (!is.na(gene_id) && gene_id %in% names(gene_to_txids)) {
        best_class    <- "novel"
        best_priority <- TX_CLASS_PRIORITY["novel"]
        for (tx_id_c in gene_to_txids[[gene_id]]) {
          if (!tx_id_c %in% names(ann$exon_txs) ||
              !tx_id_c %in% names(ann$cds_txs_tx_coords)) next
          ref_tx_c <- get_tx_coords(tx_id_c)
          if (length(ref_tx_c) == 0) next
          orf_tx_c <- project_to_tx_coords(orf_grl[[i]], ann$exon_txs[[tx_id_c]], strand_vec[i])
          if (is.null(orf_tx_c)) next
          cls <- classify_transcript(orf_tx_c[1], orf_tx_c[2],
                                     IRanges::start(ref_tx_c), IRanges::end(ref_tx_c))
          if (!is.na(cls)) {
            prio <- TX_CLASS_PRIORITY[cls]
            if (is.na(prio)) prio <- TX_CLASS_PRIORITY["novel"]
            if (prio < best_priority) {
              best_priority <- prio
              best_class    <- cls
            }
          }
          if (best_priority <= 1L) break  # can't improve on ORF_annotated
        }
        res$ORF_category_Tx_compatible[i] <- best_class
      } else if (!is.na(gene_id)) {
        # Gene has no CDS-containing transcripts (non-coding gene) → novel
        res$ORF_category_Tx_compatible[i] <- "novel"
      }
    }
  }

  return(res)
}

#' Classify ORFs using mirai-based parallelization
#'
#' Splits ORFs into chunks and dispatches each chunk to a mirai daemon.
#' Annotation objects are passed via RDS file path to avoid serialization
#' overhead — each daemon loads its own copy from disk.
#'
#' @inheritParams orfquant_classify_orfs
#' @param n_cores Number of worker processes (default: detected cores - 1).
#' @param use_parallel Logical; if FALSE, falls back to serial orfquant_classify_orfs.
#' @return data.frame with ORF_category_Gen / ORF_category_Tx / ORF_category_Tx_compatible.
#' @export
orfquant_classify_orfs_parallel <- function(orfs,
                                            annotation,
                                            gene_id_col = "gene_id",
                                            transcript_id_col = "transcript_id",
                                            seqname_col = "seqnames",
                                            start_col = "start",
                                            end_col = "end",
                                            strand_col = "strand",
                                            exons_col = NULL,
                                            orf_id_col = NULL,
                                            n_cores = NULL,
                                            use_parallel = TRUE) {
  if (!use_parallel) {
    return(orfquant_classify_orfs(
      orfs = orfs, annotation = annotation,
      gene_id_col = gene_id_col, transcript_id_col = transcript_id_col,
      seqname_col = seqname_col, start_col = start_col, end_col = end_col,
      strand_col = strand_col, exons_col = exons_col, orf_id_col = orf_id_col
    ))
  }

  if (!requireNamespace("mirai", quietly = TRUE)) {
    message("[orfquant_classify_orfs_parallel] mirai not installed; falling back to serial")
    return(orfquant_classify_orfs(
      orfs = orfs, annotation = annotation,
      gene_id_col = gene_id_col, transcript_id_col = transcript_id_col,
      seqname_col = seqname_col, start_col = start_col, end_col = end_col,
      strand_col = strand_col, exons_col = exons_col, orf_id_col = orf_id_col
    ))
  }

  message("[orfquant_classify_orfs_parallel] Loading & normalizing inputs...")
  # Use the serial function to normalize inputs (reuses its normalize_* helpers)
  # We call orfquant_classify_orfs with use_cache_only=TRUE to get normalized objects
  # without running the full classification loop.
  # Actually, just call normalize_* via the serial function's environment.
  # The simplest fix: use orfquant_classify_orfs for normalization by
  # passing the ORFs and annotation through it.
  orf_grl <- orfquant_classify_orfs(orfs = orfs, annotation = annotation,
    gene_id_col = gene_id_col, transcript_id_col = transcript_id_col,
    seqname_col = seqname_col, start_col = start_col, end_col = end_col,
    strand_col = strand_col, exons_col = exons_col, orf_id_col = orf_id_col)
  # But this returns results, not normalized objects...
  # Better approach: just use the same code as the serial function
  # Import and normalize directly
  if (is.character(orfs) && length(orfs) == 1 && file.exists(orfs)) {
    orfs <- rtracklayer::import(orfs)
  }
  if (inherits(orfs, "GRanges")) {
    orfs <- orfs[orfs$type == "CDS"]
    orfs <- split(orfs, orfs$orf_id)
  } else if (inherits(orfs, "GRangesList")) {
    # Already normalized
    orf_grl <- orfs
  } else {
    stop("Unsupported ORF input type: ", class(orfs))
  }
  
  # Normalize annotation
  if (is.character(annotation) && length(annotation) == 1 && file.exists(annotation)) {
    ann <- readRDS(paste0(annotation, ".ann_cache.rds"))
  } else if (is.list(annotation)) {
    ann <- annotation
  } else {
    stop("Unsupported annotation input type: ", class(annotation))
  }
  
  # Ensure gene_id and transcript_id are propagated to GRangesList metadata
  inner_gr <- unlist(orf_grl, use.names = TRUE)
  inner_meta <- S4Vectors::mcols(inner_gr)
  el_idx <- rep(seq_along(orf_grl), lengths(orf_grl))
  for (col in c("gene_id", "transcript_id")) {
    if (!col %in% colnames(S4Vectors::mcols(orf_grl)) && col %in% colnames(inner_meta)) {
      vals <- as.character(inner_meta[[col]])
      S4Vectors::mcols(orf_grl)[[col]] <- vapply(seq_along(orf_grl), function(i) {
        v <- vals[el_idx == i]
        if (length(v) == 0 || all(is.na(v))) NA_character_ else v[!is.na(v)][1]
      }, character(1))
    } else if (!col %in% colnames(S4Vectors::mcols(orf_grl))) {
      S4Vectors::mcols(orf_grl)[[col]] <- NA_character_
    }
  }
  if (!"orf_id" %in% colnames(S4Vectors::mcols(orf_grl)) && !is.null(names(orf_grl))) {
    S4Vectors::mcols(orf_grl)$orf_id <- names(orf_grl)
  }

  n_orfs <- length(orf_grl)
  message("[orfquant_classify_orfs_parallel] ", n_orfs, " ORFs to classify")

  # ── Determine parallelism ──────────────────────────────────────────────
  if (is.null(n_cores)) {
    n_cores <- 1L  # serial is faster than fork-based mclapply for this workload
  }
  n_cores <- min(n_cores, n_orfs)  # don't create more workers than ORFs

  if (n_cores <= 1L) {
    message("[orfquant_classify_orfs_parallel] Only 1 core; using serial path")
    return(orfquant_classify_orfs(
      orfs = orf_grl, annotation = ann,
      gene_id_col = gene_id_col, transcript_id_col = transcript_id_col,
      seqname_col = seqname_col, start_col = start_col, end_col = end_col,
      strand_col = strand_col, exons_col = exons_col, orf_id_col = orf_id_col
    ))
  }

  # ── Save annotation + ORFs to RDS for daemon-side loading ────────────
  ann_rds <- tempfile(pattern = "orfquant_ann_", fileext = ".rds")
  saveRDS(ann, ann_rds, compress = FALSE)  # no compression = faster load
  on.exit(try(unlink(ann_rds), silent = TRUE), add = TRUE)

  orf_rds <- tempfile(pattern = "orfquant_orf_", fileext = ".rds")
  saveRDS(orf_grl, orf_rds, compress = FALSE)
  on.exit(try(unlink(orf_rds), silent = TRUE), add = TRUE)

  # ── Start mirai daemons ───────────────────────────────────────────────
  message("[orfquant_classify_orfs_parallel] Starting ", n_cores, " mirai daemons...")
  mirai::daemons(n_cores, dispatcher = TRUE)
  on.exit(mirai::daemons(0), add = TRUE)

  # ── Split ORFs into balanced chunks ───────────────────────────────────
  chunk_size <- ceiling(n_orfs / n_cores)
  idx_chunks <- split(seq_len(n_orfs),
                       ceiling(seq_len(n_orfs) / chunk_size))
  message("[orfquant_classify_orfs_parallel] ",
          length(idx_chunks), " chunks, ~", chunk_size, " ORFs each")

  # ── Dispatch chunks in parallel ───────────────────────────────────────
  message("[orfquant_classify_orfs_parallel] Dispatching to daemons...")
  futures <- lapply(names(idx_chunks), function(chunk_name) {
    idxs <- idx_chunks[[chunk_name]]
    mirai::mirai(
      .expr = {
        # ─── In-daemon: load packages and data from RDS ──────────────────
        suppressPackageStartupMessages({
          library(GenomicRanges)
          library(IRanges)
          library(GenomeInfoDb)
          library(S4Vectors)
        })
        cat(sprintf("[daemon %s] Loading data from RDS...\n", chunk_name))
        ann <- readRDS(ann_rds)
        orf_grl <- readRDS(orf_rds)

        # ─── Helper functions (re-declared inside daemon) ───────────────
        project_to_tx_coords <- function(orf_gr, exon_gr, strand) {
          if (is.null(exon_gr) || length(exon_gr) == 0 ||
              is.null(orf_gr)  || length(orf_gr)  == 0) return(NULL)
          ord  <- order(IRanges::start(exon_gr))
          ex_s <- IRanges::start(exon_gr)[ord]
          ex_e <- IRanges::end(exon_gr)[ord]
          n_ex <- length(ex_s)
          if (identical(strand, "-")) {
            ex_s <- rev(ex_s)
            ex_e <- rev(ex_e)
          }
          ex_len     <- ex_e - ex_s + 1L
          cum_before <- c(0L, cumsum(ex_len[-n_ex]))
          gpos_to_tx <- function(gpos) {
            in_exon <- which(gpos >= ex_s & gpos <= ex_e)
            if (length(in_exon) == 0L) return(NA_integer_)
            j <- in_exon[1L]
            offset <- if (identical(strand, "-")) (ex_e[j] - gpos) else (gpos - ex_s[j])
            cum_before[j] + offset + 1L
          }
          tx_pos <- vapply(c(IRanges::start(orf_gr), IRanges::end(orf_gr)),
                            gpos_to_tx, integer(1))
          tx_pos <- tx_pos[!is.na(tx_pos)]
          if (length(tx_pos) == 0) return(NULL)
          c(min(tx_pos), max(tx_pos))
        }

        classify_genomic <- function(orf_gen, cds_gene, max_cdsok) {
          if (length(cds_gene) == 0 || length(orf_gen) == 0) return("novel")
          overl <- orf_gen %over% cds_gene
          if (sum(overl) == 0) {
            category <- "novel"
            nearest_cds <- cds_gene[GenomicRanges::nearest(orf_gen, cds_gene)[1]]
            overl_whole <- orf_gen@ranges %over% IRanges::IRanges(
              start = min(IRanges::start(nearest_cds)),
              end = max(IRanges::end(nearest_cds)))
            strd <- as.vector(GenomicRanges::strand(orf_gen[1]))
            if (strd == "+") {
              if (sum(overl_whole) == 0) {
                if (min(IRanges::start(orf_gen)) < min(IRanges::start(nearest_cds))) category <- "novel_Upstream"
                if (min(IRanges::start(orf_gen)) > max(IRanges::end(nearest_cds))) category <- "novel_Downstream"
              }
              if (sum(overl_whole) > 0) category <- "novel_Internal"
            }
            if (strd == "-") {
              if (sum(overl_whole) == 0) {
                if (max(IRanges::end(orf_gen)) > max(IRanges::end(nearest_cds))) category <- "novel_Upstream"
                if (min(IRanges::start(orf_gen)) < min(IRanges::start(nearest_cds))) category <- "novel_Downstream"
              }
              if (sum(overl_whole) > 0) category <- "novel_Internal"
            }
            return(category)
          }
          if (length(max_cdsok) == 0) return("overlaps_CDS")
          strd <- as.vector(GenomicRanges::strand(orf_gen[1]))
          if (strd == "+") {
            gen_sta <- min(IRanges::start(max_cdsok))
            gen_sto <- max(IRanges::end(max_cdsok))
            sta_or <- min(IRanges::start(orf_gen))
            sto_or <- max(IRanges::end(orf_gen)) + 3
            if (sto_or == gen_sto) {
              if (sta_or == gen_sta) return("exact_start_stop")
              if (sta_or < gen_sta) return("Alt5_start")
              if (sta_or > gen_sta) return("Alt3_start")
            }
            if (sto_or != gen_sto) {
              if (sta_or < gen_sta && sto_or < gen_sto) return("Alt5_start_Alt5_stop")
              if (sta_or < gen_sta && sto_or > gen_sto) return("Alt5_start_Alt3_stop")
              if (sta_or > gen_sta && sto_or > gen_sto) return("Alt3_start_Alt3_stop")
              if (sta_or > gen_sta && sto_or < gen_sto) return("Alt3_start_Alt5_stop")
              if (sta_or == gen_sta && sto_or < gen_sto) return("Alt5_stop")
              if (sta_or == gen_sta && sto_or > gen_sto) return("Alt3_stop")
            }
          }
          if (strd == "-") {
            gen_sta <- max(IRanges::end(max_cdsok))
            gen_sto <- min(IRanges::start(max_cdsok))
            sta_or <- max(IRanges::end(orf_gen))
            sto_or <- min(IRanges::start(orf_gen)) - 3
            if (sto_or == gen_sto) {
              if (sta_or == gen_sta) return("exact_start_stop")
              if (sta_or < gen_sta) return("Alt3_start")
              if (sta_or > gen_sta) return("Alt5_start")
            }
            if (sto_or != gen_sto) {
              if (sta_or > gen_sta && sto_or > gen_sto) return("Alt5_start_Alt5_stop")
              if (sta_or > gen_sta && sto_or < gen_sto) return("Alt5_start_Alt3_stop")
              if (sta_or < gen_sta && sto_or < gen_sto) return("Alt3_start_Alt3_stop")
              if (sta_or < gen_sta && sto_or > gen_sto) return("Alt3_start_Alt5_stop")
              if (sta_or == gen_sta && sto_or > gen_sto) return("Alt5_stop")
              if (sta_or == gen_sta && sto_or < gen_sto) return("Alt3_stop")
            }
          }
          return("overlaps_CDS")
        }

        classify_transcript <- function(orf_sta, orf_sto, ann_sta, ann_sto) {
          if (is.na(ann_sta) || is.na(ann_sto)) return(NA_character_)
          if (orf_sto == ann_sto) {
            if (orf_sta == ann_sta) return("ORF_annotated")
            if (orf_sta < ann_sta) return("N_extension")
            if (orf_sta > ann_sta) return("N_truncation")
          }
          if (orf_sto != ann_sto) {
            if (orf_sta < ann_sta && orf_sto < ann_sta) return("uORF")
            if (orf_sta < ann_sta && orf_sto < ann_sto) return("overl_uORF")
            if (orf_sta < ann_sta && orf_sto > ann_sto) return("NC_extension")
            if (orf_sta > ann_sto && orf_sto > ann_sto) return("dORF")
            if (orf_sta > ann_sta && orf_sto > ann_sto) return("overl_dORF")
            if (orf_sta > ann_sta && orf_sto < ann_sto) return("nested_ORF")
            if (orf_sta == ann_sta && orf_sto < ann_sto) return("C_truncation")
            if (orf_sta == ann_sta && orf_sto > ann_sto) return("C_extension")
          }
          return(NA_character_)
        }

        get_max_cds_by_gene <- function(cds_txs) {
          tx_gene <- vapply(cds_txs, function(x) {
            gid <- unique(x$gene_id)
            if (length(gid) == 0) NA_character_ else gid[1]
          }, character(1))
          gene_ids <- unique(tx_gene[!is.na(tx_gene)])
          max_cds <- list()
          for (g in gene_ids) {
            txs <- cds_txs[tx_gene == g]
            if (length(txs) == 0) next
            tx_len <- vapply(txs, function(x) sum(IRanges::width(IRanges::reduce(x))), numeric(1))
            max_cds[[g]] <- txs[[which.max(tx_len)]]
          }
          return(max_cds)
        }

        # ─── Unpack annotation ──────────────────────────────────────────
        cds_genes <- ann$cds_genes
        cds_txs <- ann$cds_txs
        exon_txs <- ann$exon_txs
        cds_txs_tx_coords <- ann$cds_txs_tx_coords
        lncrna_genes <- if (!is.null(ann$lncrna_genes)) ann$lncrna_genes else character(0)

        # Use pre-computed derived structures from annotation cache
        if (!is.null(ann$max_cds_by_gene)) {
          max_cds_by_gene <- ann$max_cds_by_gene
        } else {
          max_cds_by_gene <- get_max_cds_by_gene(cds_txs)
        }

        if (!is.null(ann$cds_genes_reduced)) {
          cds_genes_reduced <- ann$cds_genes_reduced
        } else {
          cds_genes_reduced <- cds_genes  # already reduced
        }

        if (!is.null(ann$gene_to_txids)) {
          gene_to_txids <- ann$gene_to_txids
        } else if (!is.null(cds_txs_tx_coords)) {
          tx_gene_map <- vapply(cds_txs, function(tx) {
            gids <- unique(tx$gene_id)
            if (length(gids) > 0) as.character(gids[1]) else NA_character_
          }, character(1))
          gene_to_txids <- split(names(tx_gene_map), tx_gene_map)
        } else { gene_to_txids <- list() }

        if (!is.null(ann$all_cds_merged)) {
          all_cds <- ann$all_cds_merged
        } else {
          all_cds <- IRanges::reduce(unlist(cds_genes))
        }

        TX_CLASS_PRIORITY <- c(
          ORF_annotated = 1L, N_extension = 2L, C_extension = 3L,
          N_truncation = 4L, C_truncation = 5L, NC_extension = 6L,
          overl_uORF = 7L, overl_dORF = 8L, nested_ORF = 9L,
          uORF = 10L, dORF = 11L, novel = 12L, novel_antisense = 13L
        )

        # ─── Process this chunk of ORFs ─────────────────────────────────
        orf_subset <- orf_grl[idxs]
        n_chunk <- length(orf_subset)

        seqnames_v <- vapply(orf_subset, function(x) {
          if (length(x) == 0) NA_character_ else as.character(GenomicRanges::seqnames(x)[1])
        }, character(1))
        start_v <- vapply(orf_subset, function(x) {
          if (length(x) == 0) NA_integer_ else min(IRanges::start(x))
        }, integer(1))
        end_v <- vapply(orf_subset, function(x) {
          if (length(x) == 0) NA_integer_ else max(IRanges::end(x))
        }, integer(1))
        strand_v <- vapply(orf_subset, function(x) {
          if (length(x) == 0) NA_character_ else as.character(GenomicRanges::strand(x)[1])
        }, character(1))

        res <- data.frame(
          orf_id = S4Vectors::mcols(orf_subset)$orf_id,
          gene_id = S4Vectors::mcols(orf_subset)$gene_id,
          seqnames = seqnames_v,
          start = start_v,
          end = end_v,
          strand = strand_v,
          ORF_category_Gen = NA_character_,
          ORF_category_Tx = NA_character_,
          ORF_category_Tx_compatible = NA_character_,
          stringsAsFactors = FALSE
        )

        for (j in seq_len(n_chunk)) {
          orf_gen <- orf_subset[[j]]
          gene_id <- S4Vectors::mcols(orf_subset)$gene_id[j]

          if (!is.na(gene_id) && length(lncrna_genes) > 0 && gene_id %in% lncrna_genes) {
            res$ORF_category_Gen[j] <- "lncRNA"
            res$ORF_category_Tx[j] <- "lncRNA"
            res$ORF_category_Tx_compatible[j] <- "lncRNA"
            next
          }

          if (!is.na(gene_id) && gene_id %in% names(cds_genes)) {
            cds_gene <- cds_genes_reduced[[gene_id]]
            max_cdsok <- max_cds_by_gene[[gene_id]]
          } else {
            cds_gene <- all_cds
            max_cdsok <- GenomicRanges::GRanges()
          }

          res$ORF_category_Gen[j] <- classify_genomic(orf_gen, cds_gene, max_cdsok)

          if (!is.null(cds_txs_tx_coords) && !is.null(exon_txs)) {
            tx_id <- S4Vectors::mcols(orf_subset)$transcript_id[j]
            strand_j <- strand_v[j]

            if (!is.na(tx_id) && tx_id %in% names(exon_txs)) {
              orf_tx <- project_to_tx_coords(orf_subset[[j]], exon_txs[[tx_id]], strand_j)
              if (!is.null(orf_tx) && tx_id %in% names(cds_txs_tx_coords)) {
                ref_tx <- get_tx_coords(tx_id)
                if (length(ref_tx) > 0) {
                  res$ORF_category_Tx[j] <- classify_transcript(
                    orf_tx[1], orf_tx[2],
                    IRanges::start(ref_tx), IRanges::end(ref_tx))
                } else {
                  res$ORF_category_Tx[j] <- "novel"
                }
              } else if (!is.null(orf_tx)) {
                res$ORF_category_Tx[j] <- "novel"
              }
            }

            if (!is.na(gene_id) && gene_id %in% names(gene_to_txids)) {
              best_class <- "novel"
              best_priority <- TX_CLASS_PRIORITY["novel"]
              for (tx_id_c in gene_to_txids[[gene_id]]) {
                if (!tx_id_c %in% names(exon_txs) ||
                    !tx_id_c %in% names(cds_txs_tx_coords)) next
                ref_tx_c <- get_tx_coords(tx_id_c)
                if (length(ref_tx_c) == 0) next
                orf_tx_c <- project_to_tx_coords(orf_subset[[j]], exon_txs[[tx_id_c]], strand_j)
                if (is.null(orf_tx_c)) next
                cls <- classify_transcript(orf_tx_c[1], orf_tx_c[2],
                                           IRanges::start(ref_tx_c), IRanges::end(ref_tx_c))
                if (!is.na(cls)) {
                  prio <- TX_CLASS_PRIORITY[cls]
                  if (is.na(prio)) prio <- TX_CLASS_PRIORITY["novel"]
                  if (prio < best_priority) {
                    best_priority <- prio
                    best_class <- cls
                  }
                }
                if (best_priority <= 1L) break
              }
              res$ORF_category_Tx_compatible[j] <- best_class
            } else if (!is.na(gene_id)) {
              res$ORF_category_Tx_compatible[j] <- "novel"
            }
          }
        }

        cat(sprintf("[daemon %s] Done: %d ORFs classified\n", chunk_name, n_chunk))
        res
      },
      .args = list(
        idxs = idxs,
        orf_rds = orf_rds,
        ann_rds = ann_rds,
        chunk_name = chunk_name
      )
    )
  })

  # ── Collect results ───────────────────────────────────────────────────
  message("[orfquant_classify_orfs_parallel] Collecting results from daemons...")
  chunk_results <- lapply(futures, function(f) mirai::call_mirai(f)$data)

  errors <- vapply(chunk_results, inherits, logical(1), "miraiError")
  if (any(errors)) {
    err_msg <- paste(vapply(chunk_results[errors], function(e) e$message, character(1)),
                     collapse = "; ")
    message("[orfquant_classify_orfs_parallel] ERROR in daemon(s): ", err_msg)
    message("[orfquant_classify_orfs_parallel] Falling back to serial")
    return(orfquant_classify_orfs(
      orfs = orf_grl, annotation = ann,
      gene_id_col = gene_id_col, transcript_id_col = transcript_id_col,
      seqname_col = seqname_col, start_col = start_col, end_col = end_col,
      strand_col = strand_col, exons_col = exons_col, orf_id_col = orf_id_col
    ))
  }

  final <- do.call(rbind, chunk_results)
  final <- final[order(match(final$orf_id, S4Vectors::mcols(orf_grl)$orf_id)), ]
  rownames(final) <- NULL

  message("[orfquant_classify_orfs_parallel] Done: ", nrow(final), " ORFs classified")
  return(final)
}

#' Classify ORFs using class_ORFtype.py logic (gene-level only)
#'
#' This function mirrors the classification logic in class_ORFtype.py.
#' It returns gene-level categories and leaves transcript-level categories empty.
#'
#' @inheritParams orfquant_classify_orfs
#' @return data.frame with ORF_type_py and basic ORF metadata. Transcript-level categories are NA.
#' @export

orfquant_classify_orfs_py <- function(orfs,
                                      annotation,
                                      gene_id_col = "gene_id",
                                      transcript_id_col = "transcript_id",
                                      seqname_col = "seqnames",
                                      start_col = "start",
                                      end_col = "end",
                                      strand_col = "strand",
                                      exons_col = NULL,
                                      orf_id_col = NULL) {
  if (!requireNamespace("GenomicRanges", quietly = TRUE)) {
    stop("GenomicRanges is required.")
  }
  if (!requireNamespace("IRanges", quietly = TRUE)) {
    stop("IRanges is required.")
  }
  if (!requireNamespace("rtracklayer", quietly = TRUE)) {
    stop("rtracklayer is required.")
  }

  normalize_orfs <- function(x) {
    if (is.character(x) && length(x) == 1 && file.exists(x)) {
      x <- rtracklayer::import(x)
    }

    if (inherits(x, "GRangesList")) {
      grl <- x
      # Propagate per-range metadata to list level if missing.
      inner_gr  <- unlist(grl, use.names = TRUE)
      inner_meta <- S4Vectors::mcols(inner_gr)
      el_idx    <- rep(seq_along(grl), lengths(grl))
      for (col in c("gene_id", "transcript_id")) {
        if (!col %in% colnames(S4Vectors::mcols(grl))) {
          if (col %in% colnames(inner_meta)) {
            vals <- as.character(inner_meta[[col]])
            S4Vectors::mcols(grl)[[col]] <- vapply(seq_along(grl), function(i) {
              v <- vals[el_idx == i]
              if (length(v) == 0 || all(is.na(v))) NA_character_ else v[!is.na(v)][1]
            }, character(1))
          } else {
            S4Vectors::mcols(grl)[[col]] <- NA_character_
          }
        }
      }
      # Propagate orf_id from GRangesList names (= split key) before falling back to seq index.
      if (!"orf_id" %in% colnames(S4Vectors::mcols(grl)) && !is.null(names(grl))) {
        S4Vectors::mcols(grl)$orf_id <- names(grl)
      }
    } else if (inherits(x, "GRanges")) {
      grl <- GenomicRanges::GRangesList(split(x, seq_len(length(x))))
    } else if (is.data.frame(x)) {
      if (!is.null(exons_col) && exons_col %in% colnames(x)) {
        grl <- GenomicRanges::GRangesList(lapply(seq_len(nrow(x)), function(i) {
          exon_str <- x[[exons_col]][i]
          if (is.na(exon_str) || exon_str == "") {
            return(GenomicRanges::GRanges())
          }
          parts <- unlist(strsplit(exon_str, ","))
          coords <- lapply(parts, function(p) {
            se <- unlist(strsplit(p, "-"))
            if (length(se) != 2) return(NULL)
            c(as.integer(se[1]), as.integer(se[2]))
          })
          coords <- coords[!vapply(coords, is.null, logical(1))]
          if (length(coords) == 0) {
            return(GenomicRanges::GRanges())
          }
          starts <- vapply(coords, function(v) v[1], numeric(1))
          ends <- vapply(coords, function(v) v[2], numeric(1))
          GenomicRanges::GRanges(
            seqnames = x[[seqname_col]][i],
            ranges = IRanges::IRanges(start = starts, end = ends),
            strand = x[[strand_col]][i]
          )
        }))
      } else {
        gr <- GenomicRanges::GRanges(
          seqnames = x[[seqname_col]],
          ranges = IRanges::IRanges(start = x[[start_col]], end = x[[end_col]]),
          strand = x[[strand_col]]
        )
        grl <- GenomicRanges::GRangesList(split(gr, seq_len(length(gr))))
      }

      mcols_grl <- S4Vectors::DataFrame(
        gene_id = if (gene_id_col %in% colnames(x)) x[[gene_id_col]] else NA_character_,
        transcript_id = if (transcript_id_col %in% colnames(x)) x[[transcript_id_col]] else NA_character_,
        orf_id = if (!is.null(orf_id_col) && orf_id_col %in% colnames(x)) x[[orf_id_col]] else as.character(seq_len(nrow(x)))
      )
      S4Vectors::mcols(grl) <- mcols_grl
    } else {
      stop("Unsupported ORF input. Use GRanges/GRangesList, data.frame, or BED/GTF/GFF path.")
    }

    if (is.null(S4Vectors::mcols(grl)) || nrow(S4Vectors::mcols(grl)) == 0) {
      S4Vectors::mcols(grl) <- S4Vectors::DataFrame(
        gene_id = NA_character_,
        transcript_id = NA_character_,
        orf_id = as.character(seq_len(length(grl)))
      )
    } else if (!"orf_id" %in% colnames(S4Vectors::mcols(grl))) {
      S4Vectors::mcols(grl)$orf_id <- as.character(seq_len(length(grl)))
    }

    return(grl)
  }

  normalize_annotation <- function(a) {
    if (is.list(a) && !is.null(a$cds_genes) && !is.null(a$cds_txs)) {
      return(a)
    }
    if (is.character(a) && length(a) == 1 && file.exists(a)) {
      ann <- rtracklayer::import(a)
      cds <- ann[ann$type == "CDS"]
      if (length(cds) == 0) {
        stop("No CDS features found in annotation.")
      }
      if (is.null(cds$gene_id)) {
        stop("CDS annotation must contain gene_id.")
      }
      if (is.null(cds$transcript_id)) {
        cds$transcript_id <- cds$gene_id
      }
      cds_genes <- GenomicRanges::reduce(split(cds, cds$gene_id))
      cds_txs <- split(cds, cds$transcript_id)
      # Build per-transcript CDS bounds in strand-normalised coordinates.
      # For + strand: IRanges(start=min_genomic, end=max_genomic)
      # For - strand: IRanges(start=-max_genomic, end=-min_genomic)
      # Both ORF and reference include the stop codon in genomic coords,
      # so no -3 adjustment is needed when comparing to orf coordinates.
      cds_txs_coords <- lapply(cds_txs, function(tx_cds) {
        if (length(tx_cds) == 0) return(IRanges::IRanges())
        strd   <- as.character(GenomicRanges::strand(tx_cds)[1])
        cds_lo <- min(IRanges::start(tx_cds))
        cds_hi <- max(IRanges::end(tx_cds))
        if (strd == "-") IRanges::IRanges(start = -cds_hi, end = -cds_lo)
        else             IRanges::IRanges(start = cds_lo,  end = cds_hi)
      })
      return(list(cds_genes = cds_genes, cds_txs = cds_txs, cds_txs_coords = cds_txs_coords))
    }
    stop("Unsupported annotation input. Use ORFquant Annotation or GTF/GFF path.")
  }

  check_overlap <- function(exons1, exons2) {
    if (length(exons1) == 0 || length(exons2) == 0) return(FALSE)
    ov <- GenomicRanges::findOverlaps(exons1, exons2)
    return(length(ov) > 0)
  }

  classify_py <- function(orf_exons, orf_strand, cds_exons, cds_strand) {
    if (length(cds_exons) == 0) return("novel")
    if (orf_strand != cds_strand) return("novel_antisenese")

    if (check_overlap(orf_exons, cds_exons)) {
      strand <- orf_strand
      cds_min_coord <- min(IRanges::start(cds_exons))
      cds_max_coord <- max(IRanges::end(cds_exons))
      gen_sta <- if (strand == "+") cds_min_coord else cds_max_coord
      gen_sto <- if (strand == "+") cds_max_coord else cds_min_coord

      orf_min_coord <- min(IRanges::start(orf_exons))
      orf_max_coord <- max(IRanges::end(orf_exons))
      sta_or <- if (strand == "+") orf_min_coord else orf_max_coord
      sto_or <- if (strand == "+") orf_max_coord else orf_min_coord

      if (sto_or == gen_sto) {
        if (sta_or == gen_sta) return("exact_start_stop")
        if ((strand == "+" && sta_or < gen_sta) || (strand == "-" && sta_or > gen_sta)) return("Alt5_start")
        if ((strand == "+" && sta_or > gen_sta) || (strand == "-" && sta_or < gen_sta)) return("Alt3_start")
      } else {
        sta_is_upstream <- (strand == "+" && sta_or < gen_sta) || (strand == "-" && sta_or > gen_sta)
        sta_is_downstream <- (strand == "+" && sta_or > gen_sta) || (strand == "-" && sta_or < gen_sta)
        sto_is_upstream <- (strand == "+" && sto_or < gen_sto) || (strand == "-" && sto_or > gen_sto)
        sto_is_downstream <- (strand == "+" && sto_or > gen_sto) || (strand == "-" && sto_or < gen_sto)

        if (sta_or == gen_sta) {
          return(if (sto_is_upstream) "Alt5_stop" else "Alt3_stop")
        } else if (sta_is_upstream) {
          if (sto_is_upstream) return("Alt5_start_Alt5_stop")
          if (sto_is_downstream) return("Alt5_start_Alt3_stop")
          return("Alt5_start")
        } else if (sta_is_downstream) {
          if (sto_is_downstream) return("Alt3_start_Alt3_stop")
          if (sto_is_upstream) return("Alt3_start_Alt5_stop")
          return("Alt3_start")
        }
      }

      return("complex_variant_unhandled")
    }

    cds_min_coord <- min(IRanges::start(cds_exons))
    cds_max_coord <- max(IRanges::end(cds_exons))
    orf_min_coord <- min(IRanges::start(orf_exons))
    orf_max_coord <- max(IRanges::end(orf_exons))

    if (max(orf_min_coord, cds_min_coord) <= min(orf_max_coord, cds_max_coord)) {
      return("novel_Internal")
    }

    if (orf_strand == "+") {
      if (orf_max_coord < cds_min_coord) return("novel_Upstream")
      if (orf_min_coord > cds_max_coord) return("novel_Downstream")
    } else {
      if (orf_min_coord > cds_max_coord) return("novel_Upstream")
      if (orf_max_coord < cds_min_coord) return("novel_Downstream")
    }

    return("novel_other")
  }

  orf_grl <- normalize_orfs(orfs)
  ann <- normalize_annotation(annotation)

  cds_genes <- ann$cds_genes
  cds_txs <- ann$cds_txs

  # Pre-compute IRanges::reduce() per gene ONCE (same pattern as orfquant_classify_orfs).
  cds_genes_reduced2 <- if (length(cds_genes) > 0) {
    setNames(
      lapply(names(cds_genes), function(g) GenomicRanges::reduce(unlist(cds_genes[g]))),
      names(cds_genes)
    )
  } else {
    list()
  }

  # Pre-compute gene → CDS strand map so the per-ORF loop does O(1) lookup
  # instead of iterating ALL CDS transcripts via sapply() for each ORF.
  gene_to_cds_strand2 <- if (length(cds_txs) > 0) {
    gene_ids_per_tx <- vapply(cds_txs, function(x) {
      gid <- unique(S4Vectors::mcols(x)$gene_id)
      if (length(gid) > 0 && !is.na(gid[1])) gid[1] else NA_character_
    }, character(1))
    tapply(seq_along(cds_txs), gene_ids_per_tx, function(idxs) {
      as.vector(GenomicRanges::strand(cds_txs[[idxs[1]]][1]))
    })
  } else {
    list()
  }

  seqnames_vec <- vapply(orf_grl, function(x) {
    if (length(x) == 0) return(NA_character_)
    as.character(GenomicRanges::seqnames(x)[1])
  }, character(1))
  start_vec <- vapply(orf_grl, function(x) {
    if (length(x) == 0) return(NA_integer_)
    min(IRanges::start(x))
  }, integer(1))
  end_vec <- vapply(orf_grl, function(x) {
    if (length(x) == 0) return(NA_integer_)
    max(IRanges::end(x))
  }, integer(1))
  strand_vec <- vapply(orf_grl, function(x) {
    if (length(x) == 0) return(NA_character_)
    as.character(GenomicRanges::strand(x)[1])
  }, character(1))

  res <- data.frame(
    orf_id = S4Vectors::mcols(orf_grl)$orf_id,
    gene_id = S4Vectors::mcols(orf_grl)$gene_id,
    seqnames = seqnames_vec,
    start = start_vec,
    end = end_vec,
    strand = strand_vec,
    ORF_type_py = NA_character_,
    ORF_category_Tx = NA_character_,
    ORF_category_Tx_compatible = NA_character_,
    stringsAsFactors = FALSE
  )

  for (i in seq_along(orf_grl)) {
    orf_gen <- orf_grl[[i]]
    gene_id <- S4Vectors::mcols(orf_grl)$gene_id[i]
    orf_strand <- as.vector(GenomicRanges::strand(orf_gen[1]))

    cds_gene <- GenomicRanges::GRanges()
    cds_strand <- orf_strand

    if (!is.na(gene_id) && gene_id %in% names(cds_genes)) {
      cds_gene <- cds_genes_reduced2[[gene_id]]
      strand_lookup <- gene_to_cds_strand2[[gene_id]]
      if (!is.null(strand_lookup)) cds_strand <- strand_lookup
    }

    res$ORF_type_py[i] <- classify_py(orf_gen, orf_strand, cds_gene, cds_strand)
  }

  return(res)
}
