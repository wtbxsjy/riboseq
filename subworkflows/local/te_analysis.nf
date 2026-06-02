//
// Translational Efficiency (TE) Analysis Subworkflow
//
// Orchestrates three steps:
//   1. QUANTIFY_ORFS  — featureCounts on each sample BAM against ORF annotation
//   2. MERGE_COUNTS   — combines per-sample counts into a gene-by-sample matrix
//   3. DESEQ2_DELTATE — deltaTE analysis with DESeq2 interaction model
//
// Input:
//   ch_te_bams    : tuple val(meta), path(bam), path(bai)
//                   meta must contain: id, sample_type ("riboseq" or "rnaseq"), group
//   ch_orfs_bed   : path to ORF annotation BED12
//   ch_contrasts  : tuple val(meta), val(variable), val(reference), val(target)
//                   e.g. [id: "treatment_vs_control"], "group", "control", "treatment"
//   ch_gtf        : path to reference GTF (passed through to downstream)
//

include { QUANTIFY_ORFS  } from '../../modules/local/quantify_orfs'
include { MERGE_COUNTS   } from '../../modules/local/merge_counts'
include { DESEQ2_DELTATE } from '../../modules/local/deseq2_deltate'

workflow TE_ANALYSIS {
    take:
        ch_te_bams          // channel: tuple val(meta), path(bam), path(bai)
        ch_orfs_bed         // channel: path to ORF BED12 annotation
        ch_contrasts        // channel: tuple val(meta), val(variable), val(reference), val(target)
        ch_gtf              // channel: path to reference GTF (for downstream)

    main:
        ch_versions = Channel.empty()

        // Step 1: Quantify ORFs with featureCounts for each sample BAM
        // Input: (meta, bam, bai) + ORF BED → Output: (meta, *_counts.tsv)
        QUANTIFY_ORFS(
            ch_te_bams,
            ch_orfs_bed
        )
        ch_versions = ch_versions.mix(QUANTIFY_ORFS.out.versions)

        // Collect all count files
        ch_all_counts = QUANTIFY_ORFS.out.counts
            .map { meta, counts -> counts }
            .collect()

        // Step 2: Merge per-sample counts into a single matrix
        // Input: all *_counts.tsv + sample sheet → Output: merged_counts.tsv
        // The sample sheet is forwarded from the pipeline input
        ch_sample_sheet = Channel.fromPath(params.input, checkIfExists: true)

        MERGE_COUNTS(
            ch_all_counts,
            ch_sample_sheet
        )
        ch_versions = ch_versions.mix(MERGE_COUNTS.out.versions)

        // Step 3: DeltaTE analysis for each contrast
        // Input: (contrast meta, variable, reference, target) + (counts, samplesheet)
        // The MERGE_COUNTS output is combined into a meta+data tuple for DESEQ2_DELTATE
        ch_merged_data = MERGE_COUNTS.out.sample_sheet
            .combine(MERGE_COUNTS.out.counts)
            .map { samplesheet, counts ->
                def meta = [id: 'merged']
                [meta, samplesheet, counts]
            }

        DESEQ2_DELTATE(
            ch_contrasts,
            ch_merged_data
        )
        ch_versions = ch_versions.mix(DESEQ2_DELTATE.out.versions)

    emit:
        te_results         = DESEQ2_DELTATE.out.translation         // *.translation.deltate.results.tsv
        te_genes           = DESEQ2_DELTATE.out.dtegs               // *.dtegs.deltate.genes.tsv
        te_all_results     = Channel.empty()                         // placeholder for future MultiQC
        versions           = ch_versions

        // Additional outputs for reporting
        translated_mrna    = DESEQ2_DELTATE.out.translated_mrna
        total_mrna         = DESEQ2_DELTATE.out.total_mrna
        mrna_abundance     = DESEQ2_DELTATE.out.mrna_abundance
        translation_genes  = DESEQ2_DELTATE.out.translation_genes
        intensified        = DESEQ2_DELTATE.out.intensified
        buffering          = DESEQ2_DELTATE.out.buffering
        fold_change_plot   = DESEQ2_DELTATE.out.fold_change_plot
        rdata              = DESEQ2_DELTATE.out.rdata
}
