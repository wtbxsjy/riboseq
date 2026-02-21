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
      if (orf_sta < ann_sta && orf_sto < ann_sto) return("overl_uORF")
      if (orf_sta < ann_sta && orf_sto < ann_sta) return("uORF")
      if (orf_sta < ann_sta && orf_sto > ann_sto) return("NC_extension")
      if (orf_sta > ann_sta && orf_sto > ann_sto) return("overl_dORF")
      if (orf_sta > ann_sto && orf_sto > ann_sto) return("dORF")
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
  max_cds_by_gene <- get_max_cds_by_gene(cds_txs)

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

  all_cds <- IRanges::reduce(unlist(cds_genes))

  for (i in seq_along(orf_grl)) {
    orf_gen <- orf_grl[[i]]
    gene_id <- S4Vectors::mcols(orf_grl)$gene_id[i]

    if (!is.na(gene_id) && gene_id %in% names(cds_genes)) {
      cds_gene <- IRanges::reduce(unlist(cds_genes[gene_id]))
      max_cdsok <- max_cds_by_gene[[gene_id]]
    } else {
      cds_gene <- all_cds
      max_cdsok <- GenomicRanges::GRanges()
    }

    res$ORF_category_Gen[i] <- classify_genomic(orf_gen, cds_gene, max_cdsok)

    # Transcript-level categories (requires transcript coordinates)
    if (!is.null(ann$cds_txs_coords)) {
      tx_id <- S4Vectors::mcols(orf_grl)$transcript_id[i]
      if (!is.na(tx_id) && tx_id %in% names(ann$cds_txs_coords)) {
        annotated_ORF <- ann$cds_txs_coords[[tx_id]]
        if (length(annotated_ORF) > 0) {
          # cds_txs_coords stores strand-normalised coords (negated for - strand)
          # so smaller value is always the 5' end; no -3 needed (genomic coords
          # include stop codon in both ORF and reference).
          ann_sta <- IRanges::start(annotated_ORF)
          ann_sto <- IRanges::end(annotated_ORF)
          strd_orf <- strand_vec[i]
          orf_lo  <- min(IRanges::start(orf_grl[[i]]))
          orf_hi  <- max(IRanges::end(orf_grl[[i]]))
          if (!is.na(strd_orf) && strd_orf == "-") {
            orf_sta <- -orf_hi
            orf_sto <- -orf_lo
          } else {
            orf_sta <- orf_lo
            orf_sto <- orf_hi
          }
          res$ORF_category_Tx[i] <- classify_transcript(orf_sta, orf_sto, ann_sta, ann_sto)
        } else {
          res$ORF_category_Tx[i] <- "novel"
        }
      }

      if ("compatible_ORF_id_tr" %in% colnames(S4Vectors::mcols(orf_grl))) {
        comp_id <- S4Vectors::mcols(orf_grl)$compatible_ORF_id_tr[i]
        if (!is.na(comp_id) && grepl("_", comp_id)) {
          comp_tx <- paste(strsplit(comp_id, "_")[[1]][-((length(strsplit(comp_id, "_")[[1]])-1):length(strsplit(comp_id, "_")[[1]]))], collapse = "_")
          comp_sta <- as.numeric(strsplit(comp_id, "_")[[1]][length(strsplit(comp_id, "_")[[1]])-1])
          comp_sto <- as.numeric(strsplit(comp_id, "_")[[1]][length(strsplit(comp_id, "_")[[1]])])
          if (!is.na(comp_tx) && comp_tx %in% names(ann$cds_txs_coords)) {
            annotated_ORF_comp <- ann$cds_txs_coords[[comp_tx]]
            if (length(annotated_ORF_comp) > 0) {
              ann_sta <- IRanges::start(annotated_ORF_comp)
              ann_sto <- IRanges::end(annotated_ORF_comp)
              res$ORF_category_Tx_compatible[i] <- classify_transcript(comp_sta, comp_sto, ann_sta, ann_sto)
            } else {
              res$ORF_category_Tx_compatible[i] <- "novel"
            }
          }
        }
      }
    }
  }

  return(res)
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
      cds_gene <- GenomicRanges::reduce(unlist(cds_genes[gene_id]))
      if (length(cds_txs) > 0) {
        txs_gene <- cds_txs[sapply(cds_txs, function(x) { gid <- unique(x$gene_id); length(gid) > 0 && gid[1] == gene_id })]
        if (length(txs_gene) > 0) {
          cds_strand <- as.vector(GenomicRanges::strand(txs_gene[[1]][1]))
        }
      }
    }

    res$ORF_type_py[i] <- classify_py(orf_gen, orf_strand, cds_gene, cds_strand)
  }

  return(res)
}
