/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS as BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS_GENOME        } from '../../subworkflows/nf-core/bam_dedup_stats_samtools_umitools/main'
include { BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS as BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS_TRANSCRIPTOME } from '../../subworkflows/nf-core/bam_dedup_stats_samtools_umitools/main'
include { BAM_SORT_STATS_SAMTOOLS                                                            } from '../../subworkflows/nf-core/bam_sort_stats_samtools/main'
include { FASTQ_QC_TRIM_FILTER_SETSTRANDEDNESS                                                 } from '../../subworkflows/nf-core/fastq_qc_trim_filter_setstrandedness/main'
include { BAM_DEDUP_UMI      } from '../../subworkflows/nf-core/bam_dedup_umi'
include { FASTQ_ALIGN_STAR   } from '../../subworkflows/nf-core/fastq_align_star'
include { FASTQ_ALIGN_HISAT2 } from '../../subworkflows/local/fastq_align_hisat2'
include { RPBP               } from '../../subworkflows/local/rpbp'
include { RIBOCODE           } from '../../subworkflows/local/ribocode'
include { RIBOSEQC as RIBOSEQC_PREFILTER  } from '../../subworkflows/local/riboseqc'
include { RIBOSEQC as RIBOSEQC_POSTFILTER } from '../../subworkflows/local/riboseqc'
include { ORFQUANT           } from '../../subworkflows/local/orfquant'

// Local module: sORF BAM filtering (unique mapping + contig exclusion + read length)
include { SORF_BAM_FILTER } from '../../modules/local/sorf_bam_filter'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { SAMTOOLS_INDEX                                       } from '../../modules/nf-core/samtools/index'
include { MULTIQC                                              } from '../../modules/nf-core/multiqc/main'
include { SAMTOOLS_SORT                                        } from '../../modules/nf-core/samtools/sort'
include { UMITOOLS_PREPAREFORRSEM as UMITOOLS_PREPAREFORSALMON } from '../../modules/nf-core/umitools/prepareforrsem'
include { RIBOTISH_QUALITY as RIBOTISH_QUALITY_RIBOSEQ         } from '../../modules/nf-core/ribotish/quality'
include { RIBOTISH_QUALITY as RIBOTISH_QUALITY_TISEQ           } from '../../modules/nf-core/ribotish/quality'
include { RIBOTISH_PREDICT as RIBOTISH_PREDICT_PREFILTER       } from '../../modules/nf-core/ribotish/predict'
include { RIBOTISH_PREDICT as RIBOTISH_PREDICT_POSTFILTER      } from '../../modules/nf-core/ribotish/predict'
include { RIBOTISH_PREDICT as RIBOTISH_PREDICT_ALL             } from '../../modules/nf-core/ribotish/predict'
include { RIBOTRICER_PREPAREORFS                               } from '../../modules/nf-core/ribotricer/prepareorfs'
include { RIBOTRICER_DETECTORFS as RIBOTRICER_DETECTORFS_PREFILTER } from '../../modules/nf-core/ribotricer/detectorfs'
include { RIBOTRICER_DETECTORFS as RIBOTRICER_DETECTORFS_POSTFILTER } from '../../modules/nf-core/ribotricer/detectorfs'
include { HISAT2_EXTRACTSPLICESITES                            } from '../../modules/nf-core/hisat2/extractsplicesites/main'
include { UNIFY_ORF_PREDICTIONS                                } from '../../modules/local/unify_orf_predictions/main'
include { CLASSIFY_ORFS_GENCODE; CLASSIFY_ORFS_ORFQUANT; CLASSIFY_ORFS_ORF_TYPE } from '../../modules/local/classify_orfs/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { paramsSummaryMap         } from 'plugin/nf-schema'
include { samplesheetToList        } from 'plugin/nf-schema'
include { paramsSummaryMultiqc     } from '../../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML   } from '../../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText   } from '../../subworkflows/local/utils_nfcore_riboseq_pipeline'
include { validateInputSamplesheet } from '../../subworkflows/local/utils_nfcore_riboseq_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow RIBOSEQ {

    take:
    ch_samplesheet      // channel: path(sample_sheet.csv) OR [ meta, bam, bai ] for BAM input
    ch_versions         // channel: [ path(versions.yml) ]
    ch_fasta            // channel: path(genome.fasta)
    ch_gtf              // channel: path(genome.gtf)
    ch_fai              // channel: path(genome.fai)
    ch_chrom_sizes      // channel: path(genome.sizes)
    ch_transcript_fasta // channel: path(transcript.fasta)
    ch_star_index       // channel: path(star/index/)
    ch_hisat2_index     // channel: path(hisat2/index/)
    ch_hisat2_transcriptome_index // channel: path(hisat2/transcriptome_index/)
    ch_salmon_index     // channel: path(salmon/index/)
    ch_contaminant_index // channel: path(contaminant/index/)

    main:

    // Get BAM input mode from global params (set by PIPELINE_INITIALISATION)
    def is_bam_input = params.is_bam_input ?: false

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        VALIDATE INPUTS
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    // Determine whether to filter the GTF or not
    def filterGtf =
        ((
            // Condition 1: Alignment is required and aligner is set
            !params.skip_alignment && params.aligner
        ) ||
        (
            // Condition 2: Transcript FASTA file is not provided
            !params.transcript_fasta
        )) &&
        (
            // Condition 3: --skip_gtf_filter is not provided
            !params.skip_gtf_filter
        )

    ch_multiqc_files = Channel.empty()
    ch_splicesites   = Channel.empty()

    // Initialize BAM channels
    ch_genome_bam        = Channel.empty()
    ch_genome_bam_index  = Channel.empty()
    ch_transcriptome_bam = Channel.empty()
    ch_transcriptome_bai = Channel.empty()
    ch_fastq             = Channel.empty()

    //
    // BAM INPUT MODE: Skip preprocessing and alignment
    //
    if (is_bam_input) {
        log.info "=".multiply(60)
        log.info "BAM INPUT MODE ACTIVATED"
        log.info "Skipping: preprocessing, alignment, UMI deduplication, RiboCode"
        log.info "=".multiply(60)

        // ch_samplesheet contains [ meta, bam, bai ] tuples.
        // For downstream tools (e.g. Ribo-TISH), we require a coordinate-sorted BAM plus an index.
        // We therefore always sort + index BAM inputs here and collect samtools QC stats for MultiQC.

        BAM_SORT_STATS_SAMTOOLS(
            ch_samplesheet.map { meta, bam, bai -> [ meta, bam ] },
            ch_fasta.map { [ [:], it ] }
        )

        ch_versions = ch_versions.mix(BAM_SORT_STATS_SAMTOOLS.out.versions)

        // Set genome BAM channels
        ch_genome_bam       = BAM_SORT_STATS_SAMTOOLS.out.bam
        ch_genome_bam_index = BAM_SORT_STATS_SAMTOOLS.out.bai.mix(BAM_SORT_STATS_SAMTOOLS.out.csi)

        // MultiQC inputs (samtools stats/flagstat/idxstats)
        ch_multiqc_files = ch_multiqc_files
            .mix(BAM_SORT_STATS_SAMTOOLS.out.stats.collect{ it[1] })
            .mix(BAM_SORT_STATS_SAMTOOLS.out.flagstat.collect{ it[1] })
            .mix(BAM_SORT_STATS_SAMTOOLS.out.idxstats.collect{ it[1] })

    } else {
        //
        // FASTQ INPUT MODE: Standard preprocessing and alignment
        //

        //
        // Create input channel from input file provided through params.input
        //
        Channel
            .fromList(samplesheetToList(params.input, "${projectDir}/assets/schema_input.json"))
            .map {
                meta, fastq_1, fastq_2, _bam, _bam_index ->
                    if (!fastq_2) {
                        return [ meta.id, meta + [ single_end:true ], [ fastq_1 ] ]
                    } else {
                        return [ meta.id, meta + [ single_end:false ], [ fastq_1, fastq_2 ] ]
                    }
            }
            .groupTuple()
            .map {
                validateInputSamplesheet(it)
            }
            .set { ch_fastq }

        //
        // SUBWORKFLOW: preprocess reads for RNA-seq. Includes trimming,
        // contaminant removal, strandedness inference
        //

        // The subworkflow only has to do Salmon indexing if it discovers 'auto'
        // samples, and if we haven't already made one elsewhere
        salmon_index_available = params.salmon_index || (!params.skip_pseudo_alignment && params.pseudo_aligner == 'salmon')

        FASTQ_QC_TRIM_FILTER_SETSTRANDEDNESS (
            ch_fastq,
            ch_fasta,
            ch_transcript_fasta,
            ch_gtf,
            ch_salmon_index,
            ch_contaminant_index,
            params.skip_contaminant_filter,
            params.filter_aligner,
            params.save_contaminant_reads,
            params.skip_fastqc || params.skip_qc,
            params.skip_trimming,
            params.skip_umi_extract,
            !salmon_index_available,
            params.trimmer,
            params.min_trimmed_reads,
            params.save_trimmed,
            params.with_umi,
            params.umi_discard_read,
            params.stranded_threshold,
            params.unstranded_threshold,
            params.skip_linting
        )

        ch_multiqc_files = ch_multiqc_files.mix(FASTQ_QC_TRIM_FILTER_SETSTRANDEDNESS.out.multiqc_files)
        ch_versions      = ch_versions.mix(FASTQ_QC_TRIM_FILTER_SETSTRANDEDNESS.out.versions)

        //
        // SUBWORKFLOW: align with STAR or HISAT2
        //

        if (params.aligner == 'star') {
            FASTQ_ALIGN_STAR(
                FASTQ_QC_TRIM_FILTER_SETSTRANDEDNESS.out.reads,
                ch_star_index.map { [ [:], it ] },
                ch_gtf.map { [ [:], it ] },
                params.star_ignore_sjdbgtf,
                '',
                params.seq_center ?: '',
                ch_fasta.map { [ [:], it ] },
                ch_transcript_fasta.map { [ [:], it ] }
            )

            ch_genome_bam              = FASTQ_ALIGN_STAR.out.bam
            ch_genome_bam_index        = FASTQ_ALIGN_STAR.out.bai
            ch_transcriptome_bam       = FASTQ_ALIGN_STAR.out.orig_bam_transcript
            ch_transcriptome_bai       = FASTQ_ALIGN_STAR.out.bai_transcript
            ch_versions                = ch_versions.mix(FASTQ_ALIGN_STAR.out.versions)

            ch_multiqc_files = ch_multiqc_files
                .mix(FASTQ_ALIGN_STAR.out.stats.collect{it[1]})
                .mix(FASTQ_ALIGN_STAR.out.flagstat.collect{it[1]})
                .mix(FASTQ_ALIGN_STAR.out.idxstats.collect{it[1]})
                .mix(FASTQ_ALIGN_STAR.out.log_final.collect{it[1]})
        } else if (params.aligner == 'hisat2') {

            // Extract splice sites for HISAT2
            HISAT2_EXTRACTSPLICESITES ( ch_gtf.map { [ [:], it ] } )
            ch_splicesites = HISAT2_EXTRACTSPLICESITES.out.txt
            ch_versions = ch_versions.mix(HISAT2_EXTRACTSPLICESITES.out.versions)

            FASTQ_ALIGN_HISAT2(
                FASTQ_QC_TRIM_FILTER_SETSTRANDEDNESS.out.reads,
                ch_hisat2_index.map { [ [:], it ] },
                ch_hisat2_transcriptome_index.map { [ [:], it ] },
                ch_splicesites,
                ch_fasta.map { [ [:], it ] }
            )

            ch_genome_bam        = FASTQ_ALIGN_HISAT2.out.bam
            ch_transcriptome_bam = FASTQ_ALIGN_HISAT2.out.transcriptome_bam
            ch_transcriptome_bai = FASTQ_ALIGN_HISAT2.out.transcriptome_bai
            ch_genome_bam_index  = FASTQ_ALIGN_HISAT2.out.bai
            ch_versions          = ch_versions.mix(FASTQ_ALIGN_HISAT2.out.versions)

            ch_multiqc_files = ch_multiqc_files
                .mix(FASTQ_ALIGN_HISAT2.out.stats.collect{it[1]})
                .mix(FASTQ_ALIGN_HISAT2.out.flagstat.collect{it[1]})
                .mix(FASTQ_ALIGN_HISAT2.out.idxstats.collect{it[1]})
        }

        //
        // SUBWORKFLOW: Remove duplicate reads from BAM file based on UMIs
        //

        if (params.with_umi) {

            BAM_DEDUP_UMI(
                ch_genome_bam.join(ch_genome_bam_index, by: [0]),
                ch_fasta.map { [ [:], it ] },
                params.umi_dedup_tool,
                params.umitools_dedup_stats,
                params.bam_csi_index,
                ch_transcriptome_bam,
                ch_transcript_fasta.map { [ [:], it ] }
            )

            ch_genome_bam        = BAM_DEDUP_UMI.out.bam
            ch_transcriptome_bam = BAM_DEDUP_UMI.out.transcriptome_bam
            ch_genome_bam_index  = BAM_DEDUP_UMI.out.bai
            ch_versions          = ch_versions.mix(BAM_DEDUP_UMI.out.versions)

            ch_multiqc_files = ch_multiqc_files
                .mix(BAM_DEDUP_UMI.out.multiqc_files)
        }
    }  // End of FASTQ input mode

    //
    // Take the riboseq samples and route to ribotish
    //

    ch_genome_bam
        .branch { meta, bam ->
            riboseq: meta.sample_type == 'riboseq'
                return [ meta, bam ]
            tiseq: meta.sample_type == 'tiseq'
                return [ meta, bam ]
            rnaseq: meta.sample_type == 'rnaseq'
                return [ meta, bam ]
        }
        .set{
            ch_genome_bam_by_type
        }

    ch_bams_for_analysis = ch_genome_bam_by_type.riboseq.join(ch_genome_bam_index)
    ch_fasta_gtf = ch_fasta.combine(ch_gtf).map{ fasta, gtf -> [ [:], fasta, gtf ] }

    // Pre-filter: using unfiltered BAMs (with MT reads) - suffix 'prefilter' for output organization
    ch_bams_for_prefilter = ch_bams_for_analysis.map { meta, bam, bai -> [ meta + [ filter_status: 'prefilter' ], bam, bai ] }

    // ORF prediction outputs for unified post-processing
    ch_ribotish_predictions = Channel.empty()
    ch_ribotricer_orfs      = Channel.empty()
    ch_orfquant_gtf         = Channel.empty()

    if (params.sorf_filter) {
        SORF_BAM_FILTER(
            ch_bams_for_analysis,
            ch_fai,
            params.sorf_unique_mode,
            params.sorf_unique_mapq,
            params.sorf_read_len_min,
            params.sorf_read_len_max,
            params.sorf_exclude_contigs_regex
        )
        ch_versions = ch_versions.mix(SORF_BAM_FILTER.out.versions)
        ch_multiqc_files = ch_multiqc_files.mix(SORF_BAM_FILTER.out.stats)

        SAMTOOLS_INDEX(
            SORF_BAM_FILTER.out.bam
        )
        ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions)

        def ch_filtered_index = params.bam_csi_index ? SAMTOOLS_INDEX.out.csi : SAMTOOLS_INDEX.out.bai
        ch_bams_for_sorf_prediction = SORF_BAM_FILTER.out.bam.join(ch_filtered_index)
        
        // Post-filter: using filtered BAMs (MT reads removed) for RiboseQC/ORFquant
        // When sorf_filter is enabled, RiboseQC and ORFquant should use the same filtered BAMs
        ch_bams_for_postfilter = ch_bams_for_sorf_prediction.map { meta, bam, bai -> [ meta + [ filter_status: 'postfilter' ], bam, bai ] }
    } else {
        // If sorf_filter is disabled, use original BAMs for all downstream analysis
        ch_bams_for_sorf_prediction = ch_bams_for_analysis
        ch_bams_for_postfilter = ch_bams_for_analysis.map { meta, bam, bai -> [ meta + [ filter_status: 'postfilter' ], bam, bai ] }
    }

    if (!params.skip_ribotish){
        // Ribotish quality analysis for P-site offset calculation
        RIBOTISH_QUALITY_RIBOSEQ(
            ch_bams_for_sorf_prediction,
            ch_gtf.map { [ [:], it ] }
        )
        ch_versions      = ch_versions.mix(RIBOTISH_QUALITY_RIBOSEQ.out.versions)

        // Prefilter: using unfiltered BAMs (optional, for QC comparison)
        if (params.run_prefilter_qc) {
            ribotish_prefilter_inputs = ch_bams_for_prefilter
                .join(RIBOTISH_QUALITY_RIBOSEQ.out.offset)
                .multiMap{ meta, bam, bai, offset ->
                    bam: [ meta, bam, bai ]
                    offset: [ meta, offset ]
                }

            RIBOTISH_PREDICT_PREFILTER(
                ribotish_prefilter_inputs.bam,
                [[:],[],[]],
                ch_fasta_gtf,
                [[:],[]],
                ribotish_prefilter_inputs.offset,
                [[:],[]]
            )
            ch_versions = ch_versions.mix(RIBOTISH_PREDICT_PREFILTER.out.versions)
        }

        // Postfilter: using filtered BAMs (MT removed)
        // ch_bams_for_postfilter already has filter_status set
        ribotish_postfilter_inputs = ch_bams_for_postfilter
            .join(RIBOTISH_QUALITY_RIBOSEQ.out.offset)
            .multiMap{ meta, bam, bai, offset ->
                bam: [ meta, bam, bai ]
                offset: [ meta, offset ]
            }

        RIBOTISH_PREDICT_POSTFILTER(
            ribotish_postfilter_inputs.bam,
            [[:],[],[]],
            ch_fasta_gtf,
            [[:],[]],
            ribotish_postfilter_inputs.offset,
            [[:],[]]
        )
        ch_versions = ch_versions.mix(RIBOTISH_PREDICT_POSTFILTER.out.versions)
        ch_ribotish_predictions = RIBOTISH_PREDICT_POSTFILTER.out.predictions

        if (params.sorf_predict_pooled) {
            RIBOTISH_PREDICT_ALL(
                ribotish_postfilter_inputs.bam.map{meta, bam, bai -> [[id:'allsamples', filter_status:'postfilter'], bam, bai]}.groupTuple(),
                [[:],[],[]],
                ch_fasta_gtf,
                [[:],[]],
                ribotish_postfilter_inputs.offset.map{meta, offset -> [[id:'allsamples'], offset]}.groupTuple(),
                [[:],[]]
            )
            ch_versions = ch_versions.mix(RIBOTISH_PREDICT_ALL.out.versions)
        } else {
            log.info "Pooled(all-samples) RiboTISH prediction is disabled (set --sorf_predict_pooled to enable)."
        }
    }

    if (!params.skip_ribotricer){
        RIBOTRICER_PREPAREORFS(
            ch_fasta_gtf
        )
        ch_versions = ch_versions.mix(RIBOTRICER_PREPAREORFS.out.versions)

        // Prefilter: using unfiltered BAMs (with MT reads) - optional for QC comparison
        if (params.run_prefilter_qc) {
            RIBOTRICER_DETECTORFS_PREFILTER(
                ch_bams_for_prefilter,
                RIBOTRICER_PREPAREORFS.out.candidate_orfs
            )
            ch_versions = ch_versions.mix(RIBOTRICER_DETECTORFS_PREFILTER.out.versions)
        }

        // Postfilter: using filtered BAMs (MT removed)
        // ch_bams_for_postfilter already has filter_status set
        RIBOTRICER_DETECTORFS_POSTFILTER(
            ch_bams_for_postfilter,
            RIBOTRICER_PREPAREORFS.out.candidate_orfs
        )
        ch_versions = ch_versions.mix(RIBOTRICER_DETECTORFS_POSTFILTER.out.versions)
        ch_ribotricer_orfs = RIBOTRICER_DETECTORFS_POSTFILTER.out.orfs
    }

    if (!params.skip_rpbp){
        def ribosomal_fasta = params.contaminant_fasta ? file(params.contaminant_fasta) : []
        if (!ribosomal_fasta) {
             error "RPBP requires a contaminant FASTA file. Please specify --contaminant_fasta."
        }

        RPBP(
            ch_bams_for_sorf_prediction,
            ch_fasta,
            ch_gtf,
            ribosomal_fasta
        )
        ch_versions = ch_versions.mix(RPBP.out.versions)
    }

    if (!params.skip_ribocode && !is_bam_input) {
        if (params.aligner != 'star' && params.aligner != 'hisat2') {
            log.warn "RiboCode requires STAR or HISAT2 alignment to generate transcriptome BAMs. Skipping RiboCode."
        } else {
             ch_transcriptome_bam
                .join(ch_transcriptome_bai)
                .filter { meta, bam, bai -> meta.sample_type == 'riboseq' }
                .set { ch_riboseq_transcriptome_bam }

             RIBOCODE(
                 ch_riboseq_transcriptome_bam,
                 ch_gtf,
                 ch_fasta
             )
             ch_versions = ch_versions.mix(RIBOCODE.out.versions)
        }
    } else if (!params.skip_ribocode && is_bam_input) {
        log.warn "RiboCode requires transcriptome BAM which is not available in BAM input mode. Skipping RiboCode."
    }

    //
    // RiboseQC: Comprehensive quality control for Ribo-seq data
    //
    ch_riboseqc_annotation = Channel.empty()
    ch_riboseqc_orfquant   = Channel.empty()

    if (!params.skip_riboseqc) {
        // Prefilter: using unfiltered BAMs (with MT reads) - optional for QC comparison
        if (params.run_prefilter_qc) {
            RIBOSEQC_PREFILTER(
                ch_bams_for_prefilter,
                ch_gtf,
                ch_fasta
            )
            ch_versions = ch_versions.mix(RIBOSEQC_PREFILTER.out.versions)
        }

        // Postfilter: RiboseQC uses filtered BAMs when sorf_filter is enabled
        // Data consistency: RiboseQC and ORFquant must use the same BAM source
        RIBOSEQC_POSTFILTER(
            ch_bams_for_postfilter,
            ch_gtf,
            ch_fasta
        )
        ch_versions = ch_versions.mix(RIBOSEQC_POSTFILTER.out.versions)

        // Store RiboseQC postfilter outputs for ORFquant
        ch_riboseqc_annotation = RIBOSEQC_POSTFILTER.out.annotation
        ch_riboseqc_orfquant   = RIBOSEQC_POSTFILTER.out.orfquant
    }

    //
    // ORFquant: ORF detection and quantification using RiboseQC output
    // Requires RiboseQC to generate the *_for_ORFquant input files
    //
    if (!params.skip_orfquant && !params.skip_riboseqc) {
        // Prepare ORFquant package channel
        ch_orfquant_pkg = params.orfquant_pkg ? Channel.value(file(params.orfquant_pkg, checkIfExists: true)) : Channel.value(file('NO_FILE'))

        ORFQUANT(
            ch_riboseqc_orfquant,
            ch_riboseqc_annotation,
            ch_fasta,
            ch_orfquant_pkg
        )
        ch_versions = ch_versions.mix(ORFQUANT.out.versions)
        ch_orfquant_gtf = ORFQUANT.out.gtf
    } else if (!params.skip_orfquant && params.skip_riboseqc) {
        log.warn "ORFquant requires RiboseQC output. Skipping ORFquant because RiboseQC is skipped."
    }

    //
    // Unified ORF predictions (scripts/unify_orf_predictions.py)
    //
    ch_unify_metadata = Channel.empty()
    ch_unify_bed      = Channel.empty()
    ch_unify_gtf      = Channel.empty()

    def has_unify_inputs = (!params.skip_ribotish) || (!params.skip_ribotricer) || (!params.skip_orfquant && !params.skip_riboseqc)

    if (!params.skip_unify_orf_predictions) {
        if (!has_unify_inputs) {
            log.warn "Unified ORF prediction is enabled but no ORF prediction tool ran; skipping."
        } else {
            def unify_prefix = (params.unify_orf_predictions_prefix ?: 'unified_orfs').tokenize('/').last()

            ch_ribotish_list = (params.skip_ribotish ? Channel.value([]) : ch_ribotish_predictions.map { meta, file -> file }.collect())
                .map { it ?: [] }
            ch_ribotricer_list = (params.skip_ribotricer ? Channel.value([]) : ch_ribotricer_orfs.map { meta, file -> file }.collect())
                .map { it ?: [] }
            ch_orfquant_list = ((params.skip_orfquant || params.skip_riboseqc) ? Channel.value([]) : ch_orfquant_gtf.map { meta, file -> file }.collect())
                .map { it ?: [] }

            // Combine three channels using chained combine operators
            ch_unify_inputs = ch_ribotish_list
                .combine(ch_ribotricer_list)
                .combine(ch_orfquant_list)
                .map { ribotish_files, ribotricer_files, orfquant_files ->
                    def all_files = []
                    // Extract file objects for staging
                    if (ribotish_files) { all_files.addAll(ribotish_files instanceof List ? ribotish_files : [ribotish_files]) }
                    if (ribotricer_files) { all_files.addAll(ribotricer_files instanceof List ? ribotricer_files : [ribotricer_files]) }
                    if (orfquant_files) { all_files.addAll(orfquant_files instanceof List ? orfquant_files : [orfquant_files]) }
                    // Extract filenames as strings for command line arguments
                    def ribotish_names = ribotish_files ? (ribotish_files instanceof List ? ribotish_files.collect{ it.getName() } : [ribotish_files.getName()]) : []
                    def ribotricer_names = ribotricer_files ? (ribotricer_files instanceof List ? ribotricer_files.collect{ it.getName() } : [ribotricer_files.getName()]) : []
                    def orfquant_names = orfquant_files ? (orfquant_files instanceof List ? orfquant_files.collect{ it.getName() } : [orfquant_files.getName()]) : []
                    [ ribotish_names, ribotricer_names, orfquant_names, all_files ]
                }

            UNIFY_ORF_PREDICTIONS(
                ch_unify_inputs,
                ch_gtf,
                ch_fasta,
                file("${workflow.projectDir}/scripts/unify_orf_predictions.py", checkIfExists: true)
            )
            ch_versions = ch_versions.mix(UNIFY_ORF_PREDICTIONS.out.versions)
            ch_unify_metadata = UNIFY_ORF_PREDICTIONS.out.metadata
            ch_unify_bed      = UNIFY_ORF_PREDICTIONS.out.bed
            ch_unify_gtf      = UNIFY_ORF_PREDICTIONS.out.gtf
        }
    }

    //
    // ORF classification (scripts/classify_orfs_wrapper.py)
    //
    if (!params.skip_orf_classification) {
        if (params.skip_unify_orf_predictions) {
            error "ORF classification requires unified ORF predictions. Please disable --skip_unify_orf_predictions."
        }

        def classify_mode = (params.orf_classify_mode ?: 'orf_type').toLowerCase()
        def classify_prefix = (params.unify_orf_predictions_prefix ?: 'unified_orfs').tokenize('/').last()
        def classify_wrapper = file("${workflow.projectDir}/scripts/classify_orfs_wrapper.py", checkIfExists: true)
        def class_orf_dir = file("${workflow.projectDir}/scripts/class_orf", checkIfExists: true)
        def gencode_orf_dir = file("${workflow.projectDir}/scripts/gencode-riboseqORFs", checkIfExists: true)

        if (classify_mode == 'gencode') {
            if (!params.orf_classify_ensembl_dir) {
                error "ORF classification mode 'gencode' requires --orf_classify_ensembl_dir."
            }
            CLASSIFY_ORFS_GENCODE(
                ch_unify_bed,
                ch_unify_metadata,
                classify_prefix,
                classify_wrapper,
                class_orf_dir,
                gencode_orf_dir,
                file(params.orf_classify_ensembl_dir, checkIfExists: true)
            )
            ch_versions = ch_versions.mix(CLASSIFY_ORFS_GENCODE.out.versions)
        } else if (classify_mode == 'orfquant') {
            CLASSIFY_ORFS_ORFQUANT(
                ch_unify_gtf,
                classify_prefix,
                classify_wrapper,
                class_orf_dir,
                ch_gtf
            )
            ch_versions = ch_versions.mix(CLASSIFY_ORFS_ORFQUANT.out.versions)
        } else if (classify_mode == 'orf_type') {
            CLASSIFY_ORFS_ORF_TYPE(
                ch_unify_metadata,
                classify_prefix,
                classify_wrapper,
                class_orf_dir,
                ch_gtf
            )
            ch_versions = ch_versions.mix(CLASSIFY_ORFS_ORF_TYPE.out.versions)
        } else {
            error "Unsupported orf_classify_mode: ${params.orf_classify_mode}. Use gencode, orfquant, or orf_type."
        }
    }

    //
    // Collate and save software versions
    //
    ch_versions = ch_versions.filter{it != null}

    softwareVersionsToYAML(ch_versions)
        .collectFile(storeDir: "${params.outdir}/pipeline_info", name: 'nf_core_pipeline_software_mqc_versions.yml', sort: true, newLine: true)
        .set { ch_collated_versions }

    //
    // MODULE: MultiQC
    //
    if (!params.skip_multiqc) {
        ch_multiqc_config                     = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
        ch_multiqc_custom_config              = params.multiqc_config ? Channel.fromPath(params.multiqc_config, checkIfExists: true) : Channel.empty()
        ch_multiqc_logo                       = params.multiqc_logo ? Channel.fromPath(params.multiqc_logo, checkIfExists: true) : Channel.empty()
        summary_params                        = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
        ch_workflow_summary                   = Channel.value(paramsSummaryMultiqc(summary_params))
        ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
        ch_methods_description                = Channel.value(methodsDescriptionText(ch_multiqc_custom_methods_description))
        ch_multiqc_files                      = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
        ch_multiqc_files                      = ch_multiqc_files.mix(ch_collated_versions)
        ch_multiqc_files                      = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml', sort: true))

        ch_name_replacements = is_bam_input ?
            Channel.value(file("$projectDir/assets/name_replacement.empty.txt", checkIfExists: true)) :
            ch_fastq
                .map{ meta, reads ->
                    def name1 = file(reads[0][0]).simpleName + "\t" + meta.id + '_1'
                    def fastqcnames = meta.id + "_raw\t" + meta.id + "\n" + meta.id + "_trimmed\t" + meta.id
                    if (reads[0][1] ){
                        def name2 = file(reads[0][1]).simpleName + "\t" + meta.id + '_2'
                        def fastqcnames1 = meta.id + "_raw_1\t" + meta.id + "_1\n" + meta.id + "_trimmed_1\t" + meta.id + "_1"
                        def fastqcnames2 = meta.id + "_raw_2\t" + meta.id + "_2\n" + meta.id + "_trimmed_2\t" + meta.id + "_2"
                        return [ name1, name2, fastqcnames1, fastqcnames2 ]
                    } else{
                        return [ name1, fastqcnames ]
                    }
                }
                .flatten()
                .collectFile(name: 'name_replacement.txt', newLine: true)

        MULTIQC (
            ch_multiqc_files.collect(),
            ch_multiqc_config.toList(),
            ch_multiqc_custom_config.toList(),
            ch_multiqc_logo.toList(),
            ch_name_replacements,
            []
        )
    ch_multiqc_report = MULTIQC.out.report.toList()
    } else {
        ch_multiqc_report = Channel.empty()
    }

    emit:
    multiqc_report = ch_multiqc_report   // channel: /path/to/multiqc_report.html
    versions       = ch_versions         // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
