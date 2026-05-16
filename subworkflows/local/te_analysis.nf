/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    TE ANALYSIS SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Translational Efficiency analysis using DESeq2 interaction model (deltaTE method).

Inputs:
  - ch_te_bams: All BAMs for TE quantification (RNA-seq + Ribo-seq)
                Format: [meta, bam, bai] where meta includes sample_type, treatment
  - ch_unified_orfs_bed: Unified ORF annotation in BED12 format (single file)
  - ch_contrasts: Contrast specifications channel with [meta, contrast_variable, reference, target]
  - ch_gtf: Reference GTF (value channel)

Steps:
  1. Quantify ORF expression for all samples (featureCounts)
  2. Build combined count matrix (genes x samples) via MERGE_COUNTS
  3. Create sample metadata sheet for DESeq2
  4. Run DESeq2 deltaTE per contrast
  5. Emit results and plots
*/

include { QUANTIFY_ORFS  } from '../../modules/local/quantify_orfs/main'
include { MERGE_COUNTS   } from '../../modules/local/merge_counts/main'
include { DESEQ2_DELTATE } from '../../modules/local/deseq2_deltate/main'

workflow TE_ANALYSIS {
    take:
    ch_te_bams           // [meta, bam, bai] for all samples (RNaseq + Riboseq)
    ch_unified_orfs_bed  // single BED12 file (value channel)
    ch_contrasts         // [meta, contrast_variable, reference, target]
    ch_gtf               // value channel: GTF file

    main:
    ch_versions = Channel.empty()

    // 1. Quantify ORF expression for all samples
    QUANTIFY_ORFS(
        ch_te_bams,
        ch_unified_orfs_bed
    )
    ch_versions = ch_versions.mix(QUANTIFY_ORFS.out.versions)

    // 2. Build sample sheet from metadata
    // The metadata from QUANTIFY_ORFS output carries the original meta
    ch_sample_sheet = QUANTIFY_ORFS.out.counts
        .map { meta, f -> [ meta.id, meta.sample_type ?: 'unknown', meta.treatment ?: '', f ] }
        .collect()
        .map { all_rows ->
            def lines = ["sample,type,treatment"]
            def count_paths = []
            all_rows.each { sample_id, sample_type, treatment, count_file ->
                lines << "${sample_id},${sample_type},${treatment}"
                count_paths << count_file
            }
            def ssheet = File.createTempFile('sample_sheet', '.csv')
            ssheet.text = lines.join('\\n')
            [ ssheet, count_paths ]
        }

    // 3. Merge all count files into matrix
    MERGE_COUNTS(
        ch_sample_sheet.map { ssheet, count_paths -> count_paths }.flatMap(),
        ch_sample_sheet.map { ssheet, _ -> ssheet }
    )
    ch_versions = ch_versions.mix(MERGE_COUNTS.out.versions)

    // 4. Run DESeq2 deltaTE per contrast
    DESEQ2_DELTATE(
        ch_contrasts,
        MERGE_COUNTS.out.counts
            .combine(MERGE_COUNTS.out.sample_sheet)
            .map { counts, ssheet -> [ [id:'deltate'], ssheet, counts ] }
    )
    ch_versions = ch_versions.mix(DESEQ2_DELTATE.out.versions)

    emit:
    te_results       = DESEQ2_DELTATE.out.translation
    te_translated    = DESEQ2_DELTATE.out.translated_mrna
    te_total_mrna    = DESEQ2_DELTATE.out.total_mrna
    te_plots         = DESEQ2_DELTATE.out.fold_change_plot
    te_genes         = DESEQ2_DELTATE.out.dtegs
    versions         = ch_versions
}
