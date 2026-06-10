/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
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
include { SORF_BAM_FILTER as SORF_BAM_FILTER_PATHOGEN } from '../../modules/local/sorf_bam_filter'
include { SORF_BAM_FILTER as SORF_BAM_FILTER_TRANSCRIPTOME } from '../../modules/local/sorf_bam_filter'
// Local module: split combined BAM by pathogen contig names for dual-genome analysis
include { SPLIT_BAM_BY_CONTIG } from '../../modules/local/split_bam_by_contig/main'
// Local module: merge replicate BAMs (same-group samples) before ORF prediction
include { SAMTOOLS_MERGE_REPLICATES; SAMTOOLS_MERGE_REPLICATES as SAMTOOLS_MERGE_REPLICATES_TRANSCRIPTOME } from '../../modules/local/samtools_merge/main'
// Local module: riboWaltz P-site analysis and QC
include { RIBOWALTZ                                         } from '../../modules/local/ribowaltz/main'
// Local subworkflow: TE analysis (RNA-seq + Ribo-seq integration)
include { TE_ANALYSIS                                       } from '../../subworkflows/local/te_analysis'
include { TE_ANALYSIS as TE_ANALYSIS_PATHOGEN              } from '../../subworkflows/local/te_analysis'
include { TE_ANALYSIS as TE_ANALYSIS_LNCRNA               } from '../../subworkflows/local/te_analysis'
// Local module: PRICE ORF detection (GEDI platform)
include { PRICE                                             } from '../../modules/local/price/main'
include { GTF2BED as GTF2BED_PATHOGEN                        } from '../../modules/local/gtf2bed'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { SAMTOOLS_INDEX                                       } from '../../modules/nf-core/samtools/index'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_HOST                 } from '../../modules/nf-core/samtools/index'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_PATHOGEN              } from '../../modules/nf-core/samtools/index'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_PATHOGEN_FILTERED     } from '../../modules/nf-core/samtools/index'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_TRANSCRIPTOME_FILTERED } from '../../modules/nf-core/samtools/index'
include { MULTIQC                                              } from '../../modules/nf-core/multiqc/main'
include { RIBOTISH_QUALITY as RIBOTISH_QUALITY_RIBOSEQ         } from '../../modules/nf-core/ribotish/quality'
include { RIBOTISH_PREDICT as RIBOTISH_PREDICT_PREFILTER       } from '../../modules/nf-core/ribotish/predict'
include { RIBOTISH_PREDICT as RIBOTISH_PREDICT_POSTFILTER      } from '../../modules/nf-core/ribotish/predict'
include { RIBOTISH_PREDICT as RIBOTISH_PREDICT_ALL             } from '../../modules/nf-core/ribotish/predict'
include { RIBOTRICER_PREPAREORFS                               } from '../../modules/nf-core/ribotricer/prepareorfs'
include { RIBOTRICER_DETECTORFS as RIBOTRICER_DETECTORFS_PREFILTER } from '../../modules/nf-core/ribotricer/detectorfs'
include { RIBOTRICER_DETECTORFS as RIBOTRICER_DETECTORFS_POSTFILTER } from '../../modules/nf-core/ribotricer/detectorfs'
include { HISAT2_EXTRACTSPLICESITES                            } from '../../modules/nf-core/hisat2/extractsplicesites/main'
include { UNIFY_ORF_PREDICTIONS; UNIFY_ORF_PREDICTIONS_PER_TOOL         } from '../../modules/local/unify_orf_predictions/main'
include { CLASSIFY_ORFS_GENCODE; CLASSIFY_ORFS_ORFQUANT; CLASSIFY_ORFS_ORF_TYPE } from '../../modules/local/classify_orfs/main'
include { CLASSIFY_ORFS_GENCODE   as CLASSIFY_ORFS_GENCODE_RIBOTISH   } from '../../modules/local/classify_orfs/main'
include { CLASSIFY_ORFS_GENCODE   as CLASSIFY_ORFS_GENCODE_RIBOTRICER } from '../../modules/local/classify_orfs/main'
include { CLASSIFY_ORFS_GENCODE   as CLASSIFY_ORFS_GENCODE_RIBOCODE   } from '../../modules/local/classify_orfs/main'
include { CLASSIFY_ORFS_GENCODE   as CLASSIFY_ORFS_GENCODE_ORFQUANT   } from '../../modules/local/classify_orfs/main'
include { CLASSIFY_ORFS_ORFQUANT  as CLASSIFY_ORFS_ORFQUANT_RIBOTISH  } from '../../modules/local/classify_orfs/main'
include { CLASSIFY_ORFS_ORFQUANT  as CLASSIFY_ORFS_ORFQUANT_RIBOTRICER} from '../../modules/local/classify_orfs/main'
include { CLASSIFY_ORFS_ORFQUANT  as CLASSIFY_ORFS_ORFQUANT_RIBOCODE  } from '../../modules/local/classify_orfs/main'
include { CLASSIFY_ORFS_ORFQUANT  as CLASSIFY_ORFS_ORFQUANT_ORFQUANT  } from '../../modules/local/classify_orfs/main'
include { CLASSIFY_ORFS_ORF_TYPE  as CLASSIFY_ORFS_ORF_TYPE_RIBOTISH  } from '../../modules/local/classify_orfs/main'
include { CLASSIFY_ORFS_ORF_TYPE  as CLASSIFY_ORFS_ORF_TYPE_RIBOTRICER} from '../../modules/local/classify_orfs/main'
include { CLASSIFY_ORFS_ORF_TYPE  as CLASSIFY_ORFS_ORF_TYPE_RIBOCODE  } from '../../modules/local/classify_orfs/main'
include { CLASSIFY_ORFS_ORF_TYPE  as CLASSIFY_ORFS_ORF_TYPE_ORFQUANT  } from '../../modules/local/classify_orfs/main'
include { COLLECT_QC_STATS                                     } from '../../modules/local/collect_qc_stats/main'
include { ORF_QC                                               } from '../../modules/local/orf_qc/main'
include { EXPRESSION_QUANT                                     } from '../../modules/local/expression_quant/main'

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
    def can_run_ribocode = !params.skip_ribocode && !is_bam_input && ['star', 'hisat2'].contains(params.aligner)

    ch_multiqc_files     = Channel.empty()
    ch_splicesites       = Channel.empty()
    ch_ribowaltz_psite   = Channel.empty()

    // Initialize BAM channels
    ch_genome_bam        = Channel.empty()
    ch_genome_bam_index  = Channel.empty()
    ch_transcriptome_bam = Channel.empty()
    ch_transcriptome_bai = Channel.empty()
    ch_fastq             = Channel.empty()

    // Initialize pathogen-specific channels
    ch_pathogen_bam       = Channel.empty()
    ch_pathogen_bam_index = Channel.empty()

    // Initialize QC stats channels (populated later, used by COLLECT_QC_STATS)
    ch_qc_star_logs      = Channel.empty()
    ch_qc_sorf_stats     = Channel.empty()
    ch_qc_psites_calcs   = Channel.empty()
    ch_qc_ribotish_all   = Channel.empty()
    ch_qc_ribotricer_orfs = Channel.empty()
    ch_qc_orfquant_results = Channel.empty()
    ch_qc_ribocode_txt    = Channel.empty()
    ch_qc_rt_qual         = Channel.empty()
    ch_qc_rtr_bam_summary = Channel.empty()
    ch_qc_rw_region       = Channel.empty()

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
            ch_qc_star_logs            = FASTQ_ALIGN_STAR.out.log_final.map { meta, log -> log }

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
    // BAM SPLITTING: If pathogen_contig_pattern is set, split combined BAM into host and pathogen
    //
    if (params.pathogen_contig_pattern && !params.skip_pathogen_analysis) {
        SPLIT_BAM_BY_CONTIG(
            ch_genome_bam.join(ch_genome_bam_index),
            ch_fai,
            params.pathogen_contig_pattern
        )
        ch_versions = ch_versions.mix(SPLIT_BAM_BY_CONTIG.out.versions)

        ch_host_bam       = SPLIT_BAM_BY_CONTIG.out.host_bam
        ch_pathogen_bam   = SPLIT_BAM_BY_CONTIG.out.pathogen_bam

        // Index host BAMs
        SAMTOOLS_INDEX_HOST( ch_host_bam )
        ch_host_bam_index = params.bam_csi_index ? SAMTOOLS_INDEX_HOST.out.csi : SAMTOOLS_INDEX_HOST.out.bai
        ch_versions = ch_versions.mix(SAMTOOLS_INDEX_HOST.out.versions)

        // Index pathogen BAMs
        SAMTOOLS_INDEX_PATHOGEN( ch_pathogen_bam )
        ch_pathogen_bam_index = params.bam_csi_index ? SAMTOOLS_INDEX_PATHOGEN.out.csi : SAMTOOLS_INDEX_PATHOGEN.out.bai
        ch_versions = ch_versions.mix(SAMTOOLS_INDEX_PATHOGEN.out.versions)
    } else {
        // No pathogen splitting: entire BAM is treated as host
        ch_host_bam       = ch_genome_bam
        ch_host_bam_index = ch_genome_bam_index
    }

    //
    // Take the riboseq samples and route to ribotish
    //

    ch_host_bam
        .branch { meta, bam ->
            riboseq: meta.sample_type == 'riboseq'
                return [ meta, bam ]
            tiseq: meta.sample_type == 'tiseq'
                return [ meta, bam ]
            rnaseq: meta.sample_type == 'rnaseq'
                return [ meta, bam ]
            lncrna: meta.sample_type == 'lncrna'
                return [ meta, bam ]
        }
        .set{
            ch_genome_bam_by_type
        }

    ch_bams_for_analysis = ch_genome_bam_by_type.riboseq.join(ch_host_bam_index)

    // Prepare channels for TE analysis (RNA-seq + Ribo-seq BAMs)
    ch_rnaseq_bam_bai = ch_genome_bam_by_type.rnaseq.join(ch_host_bam_index)
    ch_lncrna_bam_bai = ch_genome_bam_by_type.lncrna.join(ch_host_bam_index)
    ch_te_bams = Channel.empty()
    ch_lncrna_te_bams = Channel.empty()
    // Detect lncRNA samples for conditional TE_ANALYSIS_LNCRNA guard
    has_lncrna_samples = file(params.input, checkIfExists: true).readLines().any { it.split(',')[4]?.trim() == 'lncrna' }

    //
    // Pathogen BAM routing (dual-genome mode)
    //
    if (params.pathogen_contig_pattern && !params.skip_pathogen_analysis) {
        ch_pathogen_bam
            .branch { meta, bam ->
                riboseq: meta.sample_type == 'riboseq'
                    return [ meta, bam ]
                tiseq: meta.sample_type == 'tiseq'
                    return [ meta, bam ]
                rnaseq: meta.sample_type == 'rnaseq'
                    return [ meta, bam ]
            }
            .set{ ch_pathogen_bam_by_type }

        ch_pathogen_bams_for_analysis = ch_pathogen_bam_by_type.riboseq.join(ch_pathogen_bam_index)
        ch_pathogen_rnaseq_bam_bai = ch_pathogen_bam_by_type.rnaseq.join(ch_pathogen_bam_index)

        // Build pathogen FASTA+GTF channel for ORF prediction tools
        ch_pathogen_fasta_file = file(params.pathogen_fasta, checkIfExists: true)
        ch_pathogen_gtf_file   = file(params.pathogen_gtf, checkIfExists: true)
        ch_pathogen_fasta_gtf  = Channel.value([
            [id: 'pathogen'],
            ch_pathogen_fasta_file,
            ch_pathogen_gtf_file
        ])
    }

    // Create fasta+gtf tuple with a meaningful meta id for tools that need it (e.g., ribotricer)
    // Use .first() to convert to a value channel so it broadcasts to all per-sample invocations
    ch_fasta_gtf = ch_fasta.combine(ch_gtf).map{ fasta, gtf ->
        def genome_name = fasta.simpleName.replaceAll(/\.(genome|transcripts|dna|cdna).*/, '')
        [ [id: genome_name], fasta, gtf ]
    }.first()

    // Pre-filter: using unfiltered BAMs (with MT reads) - suffix 'prefilter' for output organization
    ch_bams_for_prefilter = ch_bams_for_analysis.map { meta, bam, bai -> [ meta + [ filter_status: 'prefilter' ], bam, bai ] }

    // ORF prediction outputs for unified post-processing
    ch_ribotish_predictions = Channel.empty()
    ch_ribotricer_orfs      = Channel.empty()
    ch_orfquant_gtf         = Channel.empty()
    ch_ribocode_gtf         = Channel.empty()

    if (params.sorf_filter) {
        SORF_BAM_FILTER(
            ch_bams_for_analysis,
            ch_fai,
            params.sorf_unique_mode,
            params.sorf_unique_mapq,
            params.sorf_read_len_min,
            params.sorf_read_len_max,
            params.sorf_exclude_contigs_regex,
            ''    // No GTF
        )
        ch_versions = ch_versions.mix(SORF_BAM_FILTER.out.versions)
        ch_multiqc_files = ch_multiqc_files.mix(SORF_BAM_FILTER.out.stats)
        ch_qc_sorf_stats = SORF_BAM_FILTER.out.stats.map { it instanceof List ? it[1] : it }

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

    //
    // Pathogen sORF BAM filter (dual-genome mode)
    //
    if (params.pathogen_contig_pattern && !params.skip_pathogen_analysis) {
        if (params.sorf_filter) {
            SORF_BAM_FILTER_PATHOGEN(
                ch_pathogen_bams_for_analysis,
                ch_fai,
                params.sorf_unique_mode,
                params.sorf_unique_mapq,
                params.sorf_read_len_min,
                params.sorf_read_len_max,
                '',
                ''
            )
            ch_versions = ch_versions.mix(SORF_BAM_FILTER_PATHOGEN.out.versions)

            SAMTOOLS_INDEX_PATHOGEN_FILTERED(
                SORF_BAM_FILTER_PATHOGEN.out.bam
            )
            ch_versions = ch_versions.mix(SAMTOOLS_INDEX_PATHOGEN_FILTERED.out.versions)

            ch_pathogen_bams_for_prediction = SORF_BAM_FILTER_PATHOGEN.out.bam
                .join(params.bam_csi_index ? SAMTOOLS_INDEX_PATHOGEN_FILTERED.out.csi : SAMTOOLS_INDEX_PATHOGEN_FILTERED.out.bai)
        } else {
            ch_pathogen_bams_for_prediction = ch_pathogen_bams_for_analysis
        }
    }

    //
    // REPLICATE BAM MERGING: Merge same-group replicates after sORF filtering
    // Produces additional merged BAMs that run through ALL ORF prediction tools
    // alongside the individual replicate BAMs.
    //
    if (params.merge_replicates) {
        // Only consider samples that have a group assigned
        ch_bams_with_group = ch_bams_for_sorf_prediction
            .filter { meta, bam, bai -> meta.group != null && meta.group != '' }

        // Group by (group, sample_type, strandedness) — all replicates in a group
        // must share the same sample_type and strandedness.
        // Build a merged meta: id="{group}_merged", preserve sample_type & strandedness.
        ch_grouped_bams = ch_bams_with_group
            .map { meta, bam, bai -> [ meta.group, meta, bam ] }
            .groupTuple(by: 0)
            .map { group, metas, bams ->
                def merged_meta = [
                    id          : "${group}_merged",
                    group       : group,
                    sample_type : metas[0].sample_type,
                    strandedness: metas[0].strandedness,
                    single_end  : metas[0].single_end,
                    is_merged   : true
                ]
                [ merged_meta, bams ]
            }

        SAMTOOLS_MERGE_REPLICATES( ch_grouped_bams )
        ch_versions = ch_versions.mix(SAMTOOLS_MERGE_REPLICATES.out.versions)

        // Build merged BAM channel: join bam + bai outputs
        ch_merged_bams = SAMTOOLS_MERGE_REPLICATES.out.bam
            .join(SAMTOOLS_MERGE_REPLICATES.out.bai)

        // Extend prediction channels with merged BAMs
        ch_bams_for_sorf_prediction = ch_bams_for_sorf_prediction.mix(ch_merged_bams)

        // Postfilter channel: also add merged BAMs (with filter_status tag)
        ch_bams_for_postfilter = ch_bams_for_postfilter.mix(
            ch_merged_bams.map { meta, bam, bai -> [ meta + [ filter_status: 'postfilter' ], bam, bai ] }
        )

        log.info "Replicate merging enabled: merged BAMs will be created for each 'group' and run through all ORF prediction tools."
    } else {
        log.info "Replicate merging disabled (set --merge_replicates to enable)."
    }

    //
    // Prepare TE analysis input channel: combine RNA-seq + Ribo-seq BAMs
    // RNA-seq BAMs use unfiltered data; Ribo-seq BAMs use filtered data (post-sORF filter)
    //
    if (!params.skip_te_analysis && params.contrasts) {
        // Ribo-seq BAMs for TE: use postfilter BAMs (stripping filter_status for consistent meta)
        ch_riboseq_te_bams = ch_bams_for_postfilter
            .map { meta, bam, bai ->
                def te_meta = meta.clone()
                te_meta.sample_type = 'riboseq'
                [ te_meta, bam, bai ]
            }
        ch_te_bams = ch_rnaseq_bam_bai.mix(ch_riboseq_te_bams)

        // lncRNA TE BAMs: mix lncRNA-seq + Ribo-seq for separate TE analysis
        if (has_lncrna_samples) {
            ch_lncrna_te_bams = ch_lncrna_bam_bai.mix(ch_riboseq_te_bams)
        }

        // Pathogen TE BAMs: assemble from pathogen BAM channels (dual-genome mode)
        if (params.pathogen_contig_pattern && !params.skip_pathogen_analysis) {
            ch_pathogen_riboseq_te_bams = ch_pathogen_bams_for_prediction
                .map { meta, bam, bai ->
                    def te_meta = meta.clone()
                    te_meta.sample_type = 'riboseq'
                    [ te_meta, bam, bai ]
                }
            ch_pathogen_te_bams = ch_pathogen_rnaseq_bam_bai.mix(ch_pathogen_riboseq_te_bams)

            // Generate CDS BED from pathogen GTF for featureCounts annotation
            GTF2BED_PATHOGEN( ch_pathogen_gtf_file )
        }
    }

    if (!params.skip_ribotish){
        // Ribotish quality analysis for P-site offset calculation
        RIBOTISH_QUALITY_RIBOSEQ(
            ch_bams_for_sorf_prediction,
            ch_gtf.map { [ [:], it ] }
        )
        ch_versions      = ch_versions.mix(RIBOTISH_QUALITY_RIBOSEQ.out.versions)
        ch_qc_rt_qual    = RIBOTISH_QUALITY_RIBOSEQ.out.distribution

        // Prefilter: using unfiltered BAMs (optional, for QC comparison)
        if (params.run_prefilter_qc) {
            // Join by sample id only, not by full meta
            ch_bams_prefilter_keyed = ch_bams_for_prefilter
                .map { meta, bam, bai -> [ meta.id, meta, bam, bai ] }
            ch_offset_keyed = RIBOTISH_QUALITY_RIBOSEQ.out.offset
                .map { meta, offset -> [ meta.id, offset ] }
            
            ribotish_prefilter_inputs = ch_bams_prefilter_keyed
                .join(ch_offset_keyed, by: 0)
                .map { id, meta, bam, bai, offset -> [ meta, bam, bai, offset ] }
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
        // Join by sample id only, not by full meta (since meta may differ in filter_status)
        // First, prepare channels with id as the join key
        ch_bams_for_ribotish = ch_bams_for_postfilter
            .map { meta, bam, bai -> [ meta.id, meta, bam, bai ] }
        ch_offset_for_ribotish = RIBOTISH_QUALITY_RIBOSEQ.out.offset
            .map { meta, offset -> [ meta.id, offset ] }
        
        ribotish_postfilter_inputs = ch_bams_for_ribotish
            .join(ch_offset_for_ribotish, by: 0)
            .map { id, meta, bam, bai, offset -> [ meta, bam, bai, offset ] }
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
        ch_qc_ribotish_all = RIBOTISH_PREDICT_POSTFILTER.out.all.map { meta, f -> f }

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
                RIBOTRICER_PREPAREORFS.out.candidate_orfs.first()
            )
            ch_versions = ch_versions.mix(RIBOTRICER_DETECTORFS_PREFILTER.out.versions)
        }

        // Postfilter: using filtered BAMs (MT removed)
        // ch_bams_for_postfilter already has filter_status set
        // Use .first() on candidate_orfs to broadcast the single genome-level file to all per-sample BAMs
        RIBOTRICER_DETECTORFS_POSTFILTER(
            ch_bams_for_postfilter,
            RIBOTRICER_PREPAREORFS.out.candidate_orfs.first()
        )
        ch_versions = ch_versions.mix(RIBOTRICER_DETECTORFS_POSTFILTER.out.versions)
        ch_ribotricer_orfs = RIBOTRICER_DETECTORFS_POSTFILTER.out.orfs
        ch_qc_ribotricer_orfs = RIBOTRICER_DETECTORFS_POSTFILTER.out.orfs.map { meta, f -> f }
        ch_qc_rtr_bam_summary = RIBOTRICER_DETECTORFS_POSTFILTER.out.bam_summary
    }

    ch_rpbp_bayes = Channel.empty()
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
        ch_rpbp_bayes = RPBP.out.bayes_factors
    }

    // Create shared transcriptome BAM channel for riboseq samples
    // Used by RiboCode and riboWaltz
    if (!is_bam_input) {
         ch_transcriptome_bam
            .join(ch_transcriptome_bai)
            .filter { meta, bam, bai -> meta.sample_type == 'riboseq' }
            .set { ch_riboseq_transcriptome_bam }
    } else {
        ch_riboseq_transcriptome_bam = Channel.empty()
    }

    //
    // Transcriptome BAM sORF filtering: Apply the same read-length, unique-mapping,
    // and contig-exclusion filters to transcriptome BAMs so RiboCode and riboWaltz
    // receive data that is consistent with the genome-BAM path used by other tools.
    // Without this step, RiboCode processes unfiltered transcriptome alignments
    // (multi-mappers, all read lengths, no contig exclusion) while Ribo-TISH,
    // Ribotricer, PRICE, and ORFquant all receive sORF-filtered genome BAMs.
    //
    if (!is_bam_input) {
        if (params.sorf_filter) {
            SORF_BAM_FILTER_TRANSCRIPTOME(
                ch_riboseq_transcriptome_bam,
                ch_fai,
                params.sorf_unique_mode,
                params.sorf_unique_mapq,
                params.sorf_read_len_min,
                params.sorf_read_len_max,
                params.sorf_exclude_contigs_regex,   // Also used to extract MT/Pt transcript IDs from GTF below
                ch_gtf.first()                        // GTF → extract transcript IDs for MT/Pt chromosomes
            )
            ch_versions = ch_versions.mix(SORF_BAM_FILTER_TRANSCRIPTOME.out.versions)

            SAMTOOLS_INDEX_TRANSCRIPTOME_FILTERED(
                SORF_BAM_FILTER_TRANSCRIPTOME.out.bam
            )
            ch_versions = ch_versions.mix(SAMTOOLS_INDEX_TRANSCRIPTOME_FILTERED.out.versions)

            def ch_tx_filtered_index = params.bam_csi_index ?
                SAMTOOLS_INDEX_TRANSCRIPTOME_FILTERED.out.csi :
                SAMTOOLS_INDEX_TRANSCRIPTOME_FILTERED.out.bai
            ch_riboseq_transcriptome_bam_filtered = SORF_BAM_FILTER_TRANSCRIPTOME.out.bam
                .join(ch_tx_filtered_index)
        } else {
            ch_riboseq_transcriptome_bam_filtered = ch_riboseq_transcriptome_bam
        }
    } else {
        ch_riboseq_transcriptome_bam_filtered = Channel.empty()
    }

    //
    // Merge transcriptome BAMs for the same replicate groups.
    // Must run AFTER the sORF filter so ch_riboseq_transcriptome_bam_filtered is available.
    //
    if (params.merge_replicates && !is_bam_input) {
        ch_tx_bams_with_group = ch_riboseq_transcriptome_bam_filtered
            .filter { meta, bam, bai -> meta.group != null && meta.group != '' }

        ch_grouped_tx_bams = ch_tx_bams_with_group
            .map { meta, bam, bai -> [ meta.group, meta, bam ] }
            .groupTuple(by: 0)
            .map { group, metas, bams ->
                def merged_meta = [
                    id          : "${group}_merged",
                    group       : group,
                    sample_type : metas[0].sample_type,
                    strandedness: metas[0].strandedness,
                    single_end  : metas[0].single_end,
                    is_merged   : true
                ]
                [ merged_meta, bams ]
            }

        SAMTOOLS_MERGE_REPLICATES_TRANSCRIPTOME( ch_grouped_tx_bams )
        ch_versions = ch_versions.mix(SAMTOOLS_MERGE_REPLICATES_TRANSCRIPTOME.out.versions)

        ch_merged_tx_bams = SAMTOOLS_MERGE_REPLICATES_TRANSCRIPTOME.out.bam
            .join(SAMTOOLS_MERGE_REPLICATES_TRANSCRIPTOME.out.bai)

        ch_riboseq_transcriptome_bam_filtered = ch_riboseq_transcriptome_bam_filtered
            .mix(ch_merged_tx_bams)

        log.info "Merged transcriptome BAMs created for each group — riboWaltz and RiboCode will process them."
    }

    if (!params.skip_ribocode && !is_bam_input) {
        if (!can_run_ribocode) {
            log.warn "RiboCode requires STAR or HISAT2 alignment to generate transcriptome BAMs. Skipping RiboCode."
        } else {
             RIBOCODE(
                 ch_riboseq_transcriptome_bam_filtered,
                 ch_gtf,
                 ch_fasta
             )
             ch_versions = ch_versions.mix(RIBOCODE.out.versions)
             ch_ribocode_gtf = RIBOCODE.out.gtf
             ch_qc_ribocode_txt = RIBOCODE.out.collapsed.map { meta, f -> f }
        }
    } else if (!params.skip_ribocode && is_bam_input) {
        log.warn "RiboCode requires transcriptome BAM which is not available in BAM input mode. Skipping RiboCode."
    }

    //
    // PRICE: GEDI-based ORF detection pipeline (Erhard Lab)
    //
    ch_price_gtf = Channel.empty()
    if (!params.skip_price) {
        PRICE(
            ch_bams_for_sorf_prediction,
            ch_fasta,
            ch_gtf
        )
        ch_versions = ch_versions.mix(PRICE.out.versions)
        ch_price_gtf = PRICE.out.orfs_tsv
        // Fall back to GTF if TSV is empty (stub mode)
        ch_price_gtf = ch_price_gtf.ifEmpty( PRICE.out.gtf )
    }

    //
    // riboWaltz: P-site offset calculation and QC analysis (runs before RiboseQC
    // so its offsets can serve as fallback when RiboseQC P_sites_calcs is empty)
    //
    if (!params.skip_ribowaltz) {
        // Use transcriptome BAM in alignment mode, genome BAM in BAM-input mode
        def ch_ribowaltz_bam = is_bam_input ? ch_bams_for_postfilter : ch_riboseq_transcriptome_bam_filtered

        RIBOWALTZ(
            ch_ribowaltz_bam,
            ch_gtf,
            ch_fasta
        )
        ch_versions = ch_versions.mix(RIBOWALTZ.out.versions)
        ch_ribowaltz_psite = RIBOWALTZ.out.psite_offset
        ch_qc_rw_region = RIBOWALTZ.out.region_distribution
        // Feed QC outputs to MultiQC-compatible collection
        ch_multiqc_files = ch_multiqc_files.mix(
            RIBOWALTZ.out.psite_offset.map { meta, f -> f },
            RIBOWALTZ.out.cds_coverage.map { meta, f -> f }.ifEmpty(Channel.empty()),
            RIBOWALTZ.out.codon_usage.map { meta, f -> f }.ifEmpty(Channel.empty()),
            RIBOWALTZ.out.frame_distribution.map { meta, f -> f }.ifEmpty(Channel.empty())
        )
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
                ch_fasta,
                Channel.empty()
            )
            ch_versions = ch_versions.mix(RIBOSEQC_PREFILTER.out.versions)
        }

        // Postfilter: RiboseQC uses filtered BAMs when sorf_filter is enabled
        // Data consistency: RiboseQC and ORFquant must use the same BAM source
        RIBOSEQC_POSTFILTER(
            ch_bams_for_postfilter,
            ch_gtf,
            ch_fasta,
            ch_ribowaltz_psite
        )
        ch_versions = ch_versions.mix(RIBOSEQC_POSTFILTER.out.versions)

        // Store RiboseQC postfilter outputs for ORFquant
        ch_riboseqc_annotation = RIBOSEQC_POSTFILTER.out.annotation
        ch_riboseqc_orfquant   = RIBOSEQC_POSTFILTER.out.orfquant
        ch_qc_psites_calcs     = RIBOSEQC_POSTFILTER.out.psites_calcs.map { meta, f -> f }
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
        ch_qc_orfquant_results = ORFQUANT.out.results.map { meta, f -> f }
    } else if (!params.skip_orfquant && params.skip_riboseqc) {
        log.warn "ORFquant requires RiboseQC output. Skipping ORFquant because RiboseQC is skipped."
    }

    //
    // Unified ORF predictions (scripts/unify_orf_predictions.py)
    //
    ch_unify_metadata = Channel.empty()
    ch_unify_bed      = Channel.empty()
    ch_unify_gtf      = Channel.empty()

    def has_unify_inputs = (!params.skip_ribotish) || (!params.skip_ribotricer) || can_run_ribocode || (!params.skip_orfquant && !params.skip_riboseqc) || (!params.skip_price)

    if (!params.skip_unify_orf_predictions) {
        if (!has_unify_inputs) {
            log.warn "Unified ORF prediction is enabled but no ORF prediction tool ran; skipping."
        } else {
            def unify_prefix = (params.unify_orf_predictions_prefix ?: 'unified_orfs').tokenize('/').last()

            ch_ribotish_list = (params.skip_ribotish ? Channel.value([]) : ch_ribotish_predictions.map { meta, file -> file }.collect())
                .map { it ?: [] }
            ch_ribotricer_list = (params.skip_ribotricer ? Channel.value([]) : ch_ribotricer_orfs.map { meta, file -> file }.collect())
                .map { it ?: [] }
            ch_ribocode_list = (can_run_ribocode ? ch_ribocode_gtf.map { meta, file -> file }.collect() : Channel.value([]))
                .map { it ?: [] }
            ch_orfquant_list = ((params.skip_orfquant || params.skip_riboseqc) ? Channel.value([]) : ch_orfquant_gtf.map { meta, file -> file }.collect())
                .map { it ?: [] }
            ch_price_list = (params.skip_price ? Channel.value([]) : ch_price_gtf.map { meta, file -> file }.collect())
                .map { it ?: [] }

            // Collect RiboseQC P-site bedgraph files for unified P-site statistics
            // Each sample produces two bedgraph files: *_P_sites_plus.bedgraph and *_P_sites_minus.bedgraph
            ch_psites_bedgraph = params.skip_riboseqc ? 
                Channel.value([]) :
                RIBOSEQC_POSTFILTER.out.psites_bedgraph
                    .map { meta, files -> files }
                    .flatten()
                    .collect()
                    .map { it ?: [] }
            
            // Collect sample names from RiboseQC output
            ch_sample_list = params.skip_riboseqc ?
                Channel.value([]) :
                RIBOSEQC_POSTFILTER.out.psites_bedgraph
                    .map { meta, files -> meta.id }
                    .collect()
                    .map { it ?: [] }

            // Combine five channels with robust handling for nested or flat structures
            ch_unify_inputs = ch_ribotish_list
                .combine(ch_ribotricer_list)
                .combine(ch_ribocode_list)
                .combine(ch_orfquant_list)
                .combine(ch_price_list)
                .map { combined ->
                    def ribotish_files = []
                    def ribotricer_files = []
                    def ribocode_files = []
                    def orfquant_files = []
                    def price_files = []

                    if (combined instanceof List && combined.size() == 2 && combined[0] instanceof List && combined[0].size() == 2 && combined[0][0] instanceof List && combined[0][0].size() == 2 && combined[0][0][0] instanceof List) {
                        // 5-level nested: [[[[ribotish, ribotricer], ribocode], orfquant], price]
                        ribotish_files = combined[0][0][0][0]
                        ribotricer_files = combined[0][0][0][1]
                        ribocode_files = combined[0][0][1]
                        orfquant_files = combined[0][1]
                        price_files = combined[1]
                    } else if (combined instanceof List && combined.size() == 5) {
                        // Flat 5-element structure
                        ribotish_files = combined[0]
                        ribotricer_files = combined[1]
                        ribocode_files = combined[2]
                        orfquant_files = combined[3]
                        price_files = combined[4]
                    } else if (combined instanceof List && combined.size() == 4) {
                        // Flat 4-element (no price)
                        ribotish_files = combined[0]
                        ribotricer_files = combined[1]
                        ribocode_files = combined[2]
                        orfquant_files = combined[3]
                    } else if (combined instanceof List && combined[0] instanceof List && combined[0].size() == 2 && combined[0][0] instanceof List) {
                        // 4-level nested: [[[ribotish, ribotricer], ribocode], orfquant]
                        ribotish_files = combined[0][0][0]
                        ribotricer_files = combined[0][0][1]
                        ribocode_files = combined[0][1]
                        orfquant_files = combined[1]
                    } else if (combined instanceof List) {
                        // Flat list of files - split by tool-specific suffix
                        ribotish_files = combined.findAll { it.getName().endsWith('_pred.txt') }
                        ribotricer_files = combined.findAll { it.getName().endsWith('_translating_ORFs.tsv') }
                        ribocode_files = combined.findAll { it.getName().endsWith('.gtf.gz') && !it.getName().endsWith('_Detected_ORFs.gtf.gz') }
                        orfquant_files = combined.findAll { it.getName().endsWith('_Detected_ORFs.gtf.gz') && !it.getName().contains('PRICE') }
                        price_files = combined.findAll { it.getName().endsWith('_Detected_ORFs.gtf.gz') && it.getName().contains('PRICE') || it.getName().endsWith('.orfs.tsv') }
                    } else {
                        ribotish_files = combined
                    }

                    def ribotish_list = ribotish_files instanceof List ? ribotish_files : (ribotish_files ? [ribotish_files] : [])
                    def ribotricer_list = ribotricer_files instanceof List ? ribotricer_files : (ribotricer_files ? [ribotricer_files] : [])
                    def ribocode_list = ribocode_files instanceof List ? ribocode_files : (ribocode_files ? [ribocode_files] : [])
                    def orfquant_list = orfquant_files instanceof List ? orfquant_files : (orfquant_files ? [orfquant_files] : [])
                    def price_list = price_files instanceof List ? price_files : (price_files ? [price_files] : [])

                    def all_files = []
                    all_files.addAll(ribotish_list)
                    all_files.addAll(ribotricer_list)
                    all_files.addAll(ribocode_list)
                    all_files.addAll(orfquant_list)
                    all_files.addAll(price_list)
                    def ribotish_names = ribotish_list.collect{ it.getName() }
                    def ribotricer_names = ribotricer_list.collect{ it.getName() }
                    def ribocode_names = ribocode_list.collect{ it.getName() }
                    def orfquant_names = orfquant_list.collect{ it.getName() }
                    def price_names = price_list.collect{ it.getName() }
                    [ ribotish_names, ribotricer_names, ribocode_names, orfquant_names, price_names, all_files ]
                }

            UNIFY_ORF_PREDICTIONS(
                ch_unify_inputs,
                ch_gtf,
                ch_fasta,
                file("${workflow.projectDir}/scripts/unify_orf_predictions.py", checkIfExists: true),
                file("${workflow.projectDir}/scripts/run_orf.py"),
                ch_psites_bedgraph,
                ch_sample_list
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
        def classify_prefix = (params.unify_orf_predictions_prefix ?: 'unified_orfs').tokenize('/').last()
        def classify_wrapper = file("${workflow.projectDir}/scripts/classify_orfs_wrapper.py", checkIfExists: true)
        def class_orf_dir = file("${workflow.projectDir}/scripts/class_orf", checkIfExists: true)
        def gencode_orf_dir = file("${workflow.projectDir}/scripts/gencode-riboseqORFs", checkIfExists: true)
        def base_classify_dir = params.orf_classify_output_dir ?: 'orf_classification'

        if (params.orf_classify_mode) {
            log.warn "orf_classify_mode is ignored; running all ORF classification modes."
        }

        if (params.skip_unify_orf_predictions) {
            // Per-tool mode: exact-match dedup per tool, no cross-tool merging, classify each independently
            if (!has_unify_inputs) {
                log.warn "Per-tool ORF classification is enabled but no ORF prediction tool ran; skipping."
            } else {
                // Re-use the same input channel construction as the unified path
                def ch_ribotish_list_pt = (params.skip_ribotish ? Channel.value([]) : ch_ribotish_predictions.map { meta, file -> file }.collect())
                    .map { it ?: [] }
                def ch_ribotricer_list_pt = (params.skip_ribotricer ? Channel.value([]) : ch_ribotricer_orfs.map { meta, file -> file }.collect())
                    .map { it ?: [] }
                def ch_ribocode_list_pt = (can_run_ribocode ? ch_ribocode_gtf.map { meta, file -> file }.collect() : Channel.value([]))
                    .map { it ?: [] }
                def ch_orfquant_list_pt = ((params.skip_orfquant || params.skip_riboseqc) ? Channel.value([]) : ch_orfquant_gtf.map { meta, file -> file }.collect())
                    .map { it ?: [] }
                def ch_psites_bedgraph_pt = params.skip_riboseqc ?
                    Channel.value([]) :
                    RIBOSEQC_POSTFILTER.out.psites_bedgraph.map { meta, files -> files }.flatten().collect().map { it ?: [] }
                def ch_sample_list_pt = params.skip_riboseqc ?
                    Channel.value([]) :
                    RIBOSEQC_POSTFILTER.out.psites_bedgraph.map { meta, files -> meta.id }.collect().map { it ?: [] }

                def ch_unify_inputs_pt = ch_ribotish_list_pt
                    .combine(ch_ribotricer_list_pt)
                    .combine(ch_ribocode_list_pt)
                    .combine(ch_orfquant_list_pt)
                    .map { combined ->
                        def ribotish_files = []; def ribotricer_files = []; def ribocode_files = []; def orfquant_files = []
                        if (combined instanceof List && combined.size() == 2 && combined[0] instanceof List && combined[0].size() == 2 && combined[0][0] instanceof List) {
                            ribotish_files = combined[0][0][0]; ribotricer_files = combined[0][0][1]; ribocode_files = combined[0][1]; orfquant_files = combined[1]
                        } else if (combined instanceof List && combined.size() == 4) {
                            ribotish_files = combined[0]; ribotricer_files = combined[1]; ribocode_files = combined[2]; orfquant_files = combined[3]
                        } else if (combined instanceof List) {
                            ribotish_files  = combined.findAll { it.getName().endsWith('_pred.txt') }
                            ribotricer_files = combined.findAll { it.getName().endsWith('_translating_ORFs.tsv') }
                            ribocode_files  = combined.findAll { it.getName().endsWith('.gtf.gz') && !it.getName().endsWith('_Detected_ORFs.gtf.gz') }
                            orfquant_files  = combined.findAll { it.getName().endsWith('_Detected_ORFs.gtf.gz') }
                        } else { ribotish_files = combined }
                        def rt  = ribotish_files  instanceof List ? ribotish_files  : (ribotish_files  ? [ribotish_files]  : [])
                        def rtr = ribotricer_files instanceof List ? ribotricer_files : (ribotricer_files ? [ribotricer_files] : [])
                        def rc  = ribocode_files   instanceof List ? ribocode_files   : (ribocode_files   ? [ribocode_files]   : [])
                        def oq  = orfquant_files   instanceof List ? orfquant_files   : (orfquant_files   ? [orfquant_files]   : [])
                        def all = []; all.addAll(rt); all.addAll(rtr); all.addAll(rc); all.addAll(oq)
                        [ rt.collect{ it.getName() }, rtr.collect{ it.getName() }, rc.collect{ it.getName() }, oq.collect{ it.getName() }, all ]
                    }

                UNIFY_ORF_PREDICTIONS_PER_TOOL(
                    ch_unify_inputs_pt,
                    ch_gtf,
                    ch_fasta,
                    file("${workflow.projectDir}/scripts/unify_orf_predictions.py", checkIfExists: true),
                    file("${workflow.projectDir}/scripts/run_orf.py"),
                    ch_psites_bedgraph_pt,
                    ch_sample_list_pt
                )
                ch_versions = ch_versions.mix(UNIFY_ORF_PREDICTIONS_PER_TOOL.out.versions)

                // Run all three classifiers for each tool's exact-deduplicated outputs.
                def ch_ensembl_dir = params.orf_classify_ensembl_dir ?
                    Channel.value(file(params.orf_classify_ensembl_dir, checkIfExists: true)) :
                    Channel.empty()

                if (!params.skip_ribotish) {
                    def rt_prefix = "${classify_prefix}_ribotish"
                    if (params.orf_classify_ensembl_dir) {
                        CLASSIFY_ORFS_GENCODE_RIBOTISH(
                            UNIFY_ORF_PREDICTIONS_PER_TOOL.out.ribotish_bed,
                            UNIFY_ORF_PREDICTIONS_PER_TOOL.out.ribotish_metadata,
                            rt_prefix,
                            classify_wrapper,
                            class_orf_dir,
                            gencode_orf_dir,
                            file(params.orf_classify_ensembl_dir, checkIfExists: true),
                            "${base_classify_dir}/per_tool/ribotish/gencode"
                        )
                        ch_versions = ch_versions.mix(CLASSIFY_ORFS_GENCODE_RIBOTISH.out.versions)
                    }
                    if (!params.skip_orf_classify_orfquant) {
                        CLASSIFY_ORFS_ORFQUANT_RIBOTISH(
                            UNIFY_ORF_PREDICTIONS_PER_TOOL.out.ribotish_gtf,
                            UNIFY_ORF_PREDICTIONS_PER_TOOL.out.ribotish_metadata,
                            rt_prefix,
                            classify_wrapper,
                            class_orf_dir,
                            ch_gtf,
                            "${base_classify_dir}/per_tool/ribotish/orfquant"
                        )
                        ch_versions = ch_versions.mix(CLASSIFY_ORFS_ORFQUANT_RIBOTISH.out.versions)
                    }
                    CLASSIFY_ORFS_ORF_TYPE_RIBOTISH(
                        UNIFY_ORF_PREDICTIONS_PER_TOOL.out.ribotish_metadata,
                        rt_prefix,
                        classify_wrapper,
                        class_orf_dir,
                        ch_gtf,
                        "${base_classify_dir}/per_tool/ribotish/orf_type"
                    )
                    ch_versions = ch_versions.mix(CLASSIFY_ORFS_ORF_TYPE_RIBOTISH.out.versions)
                }

                if (!params.skip_ribotricer) {
                    def rtr_prefix = "${classify_prefix}_ribotricer"
                    if (params.orf_classify_ensembl_dir) {
                        CLASSIFY_ORFS_GENCODE_RIBOTRICER(
                            UNIFY_ORF_PREDICTIONS_PER_TOOL.out.ribotricer_bed,
                            UNIFY_ORF_PREDICTIONS_PER_TOOL.out.ribotricer_metadata,
                            rtr_prefix,
                            classify_wrapper,
                            class_orf_dir,
                            gencode_orf_dir,
                            file(params.orf_classify_ensembl_dir, checkIfExists: true),
                            "${base_classify_dir}/per_tool/ribotricer/gencode"
                        )
                        ch_versions = ch_versions.mix(CLASSIFY_ORFS_GENCODE_RIBOTRICER.out.versions)
                    }
                    if (!params.skip_orf_classify_orfquant) {
                        CLASSIFY_ORFS_ORFQUANT_RIBOTRICER(
                            UNIFY_ORF_PREDICTIONS_PER_TOOL.out.ribotricer_gtf,
                            UNIFY_ORF_PREDICTIONS_PER_TOOL.out.ribotricer_metadata,
                            rtr_prefix,
                            classify_wrapper,
                            class_orf_dir,
                            ch_gtf,
                            "${base_classify_dir}/per_tool/ribotricer/orfquant"
                        )
                        ch_versions = ch_versions.mix(CLASSIFY_ORFS_ORFQUANT_RIBOTRICER.out.versions)
                    }
                    CLASSIFY_ORFS_ORF_TYPE_RIBOTRICER(
                        UNIFY_ORF_PREDICTIONS_PER_TOOL.out.ribotricer_metadata,
                        rtr_prefix,
                        classify_wrapper,
                        class_orf_dir,
                        ch_gtf,
                        "${base_classify_dir}/per_tool/ribotricer/orf_type"
                    )
                    ch_versions = ch_versions.mix(CLASSIFY_ORFS_ORF_TYPE_RIBOTRICER.out.versions)
                }

                if (can_run_ribocode) {
                    def rc_prefix = "${classify_prefix}_ribocode"
                    if (params.orf_classify_ensembl_dir) {
                        CLASSIFY_ORFS_GENCODE_RIBOCODE(
                            UNIFY_ORF_PREDICTIONS_PER_TOOL.out.ribocode_bed,
                            UNIFY_ORF_PREDICTIONS_PER_TOOL.out.ribocode_metadata,
                            rc_prefix,
                            classify_wrapper,
                            class_orf_dir,
                            gencode_orf_dir,
                            file(params.orf_classify_ensembl_dir, checkIfExists: true),
                            "${base_classify_dir}/per_tool/ribocode/gencode"
                        )
                        ch_versions = ch_versions.mix(CLASSIFY_ORFS_GENCODE_RIBOCODE.out.versions)
                    }
                    if (!params.skip_orf_classify_orfquant) {
                        CLASSIFY_ORFS_ORFQUANT_RIBOCODE(
                            UNIFY_ORF_PREDICTIONS_PER_TOOL.out.ribocode_gtf,
                            UNIFY_ORF_PREDICTIONS_PER_TOOL.out.ribocode_metadata,
                            rc_prefix,
                            classify_wrapper,
                            class_orf_dir,
                            ch_gtf,
                            "${base_classify_dir}/per_tool/ribocode/orfquant"
                        )
                        ch_versions = ch_versions.mix(CLASSIFY_ORFS_ORFQUANT_RIBOCODE.out.versions)
                    }
                    CLASSIFY_ORFS_ORF_TYPE_RIBOCODE(
                        UNIFY_ORF_PREDICTIONS_PER_TOOL.out.ribocode_metadata,
                        rc_prefix,
                        classify_wrapper,
                        class_orf_dir,
                        ch_gtf,
                        "${base_classify_dir}/per_tool/ribocode/orf_type"
                    )
                    ch_versions = ch_versions.mix(CLASSIFY_ORFS_ORF_TYPE_RIBOCODE.out.versions)
                }

                if (!params.skip_orfquant && !params.skip_riboseqc) {
                    def oq_prefix = "${classify_prefix}_orfquant"
                    if (params.orf_classify_ensembl_dir) {
                        CLASSIFY_ORFS_GENCODE_ORFQUANT(
                            UNIFY_ORF_PREDICTIONS_PER_TOOL.out.orfquant_bed,
                            UNIFY_ORF_PREDICTIONS_PER_TOOL.out.orfquant_metadata,
                            oq_prefix,
                            classify_wrapper,
                            class_orf_dir,
                            gencode_orf_dir,
                            file(params.orf_classify_ensembl_dir, checkIfExists: true),
                            "${base_classify_dir}/per_tool/orfquant/gencode"
                        )
                        ch_versions = ch_versions.mix(CLASSIFY_ORFS_GENCODE_ORFQUANT.out.versions)
                    }
                    if (!params.skip_orf_classify_orfquant) {
                        CLASSIFY_ORFS_ORFQUANT_ORFQUANT(
                            UNIFY_ORF_PREDICTIONS_PER_TOOL.out.orfquant_gtf,
                            UNIFY_ORF_PREDICTIONS_PER_TOOL.out.orfquant_metadata,
                            oq_prefix,
                            classify_wrapper,
                            class_orf_dir,
                            ch_gtf,
                            "${base_classify_dir}/per_tool/orfquant/orfquant"
                        )
                        ch_versions = ch_versions.mix(CLASSIFY_ORFS_ORFQUANT_ORFQUANT.out.versions)
                    }
                    CLASSIFY_ORFS_ORF_TYPE_ORFQUANT(
                        UNIFY_ORF_PREDICTIONS_PER_TOOL.out.orfquant_metadata,
                        oq_prefix,
                        classify_wrapper,
                        class_orf_dir,
                        ch_gtf,
                        "${base_classify_dir}/per_tool/orfquant/orf_type"
                    )
                    ch_versions = ch_versions.mix(CLASSIFY_ORFS_ORF_TYPE_ORFQUANT.out.versions)
                }
            }
        } else {
            // Unified classification path (default: skip_unify_orf_predictions = false)
            def gencode_outdir = "${base_classify_dir}/gencode"
            def orfquant_outdir = "${base_classify_dir}/orfquant"
            def orftype_outdir = "${base_classify_dir}/orf_type"

            if (!params.orf_classify_ensembl_dir) {
                error "ORF classification requires --orf_classify_ensembl_dir to run gencode classification."
            }
            CLASSIFY_ORFS_GENCODE(
                ch_unify_bed,
                ch_unify_metadata,
                classify_prefix,
                classify_wrapper,
                class_orf_dir,
                gencode_orf_dir,
                file(params.orf_classify_ensembl_dir, checkIfExists: true),
                gencode_outdir
            )
            ch_versions = ch_versions.mix(CLASSIFY_ORFS_GENCODE.out.versions)

            if (!params.skip_orf_classify_orfquant) {
                CLASSIFY_ORFS_ORFQUANT(
                    ch_unify_gtf,
                    ch_unify_metadata,
                    classify_prefix,
                    classify_wrapper,
                    class_orf_dir,
                    ch_gtf,
                    orfquant_outdir
                )
                ch_versions = ch_versions.mix(CLASSIFY_ORFS_ORFQUANT.out.versions)
            }

            CLASSIFY_ORFS_ORF_TYPE(
                ch_unify_metadata,
                classify_prefix,
                classify_wrapper,
                class_orf_dir,
                ch_gtf,
                orftype_outdir
            )
            ch_versions = ch_versions.mix(CLASSIFY_ORFS_ORF_TYPE.out.versions)
        }
    }

    //
    // Translational Efficiency (TE) Analysis
    // Integrates RNA-seq and Ribo-seq data to detect differentially translated genes
    //
    if (!params.skip_te_analysis && params.contrasts) {
        // Parse contrasts CSV file
        def contrasts_file = file(params.contrasts, checkIfExists: true)
        ch_contrasts = Channel.fromPath(params.contrasts, checkIfExists: true)
            .splitCsv(header: true)
            .map { row ->
                def contrast_meta = [id: row.id]
                [ contrast_meta, row.variable, row.reference, row.target ]
            }

        // Get unified ORFs BED for quantification annotation
        def ch_unified_bed_for_te = params.skip_unify_orf_predictions ?
            Channel.empty() :
            ch_unify_bed.first()

        if (ch_unified_bed_for_te) {
            TE_ANALYSIS(
                ch_te_bams,
                ch_unified_bed_for_te,
                ch_contrasts,
                ch_gtf
            )
            ch_versions = ch_versions.mix(TE_ANALYSIS.out.versions)

            // Feed TE outputs to MultiQC
            ch_multiqc_files = ch_multiqc_files.mix(
                TE_ANALYSIS.out.te_results.map { meta, f -> f },
                TE_ANALYSIS.out.te_genes.map { meta, f -> f }
            )
        } else {
            log.warn "TE analysis requires unified ORF BED annotation. Run with --skip_unify_orf_predictions=false."
        }

        // lncRNA TE analysis: compare Ribo-seq vs lncRNA-seq (separate from RNA-seq)
        if (ch_unified_bed_for_te && has_lncrna_samples) {
            TE_ANALYSIS_LNCRNA(
                ch_lncrna_te_bams,
                ch_unified_bed_for_te,
                ch_contrasts,
                ch_gtf
            )
            ch_versions = ch_versions.mix(TE_ANALYSIS_LNCRNA.out.versions)

            ch_multiqc_files = ch_multiqc_files.mix(
                TE_ANALYSIS_LNCRNA.out.te_results.map { meta, f -> f },
                TE_ANALYSIS_LNCRNA.out.te_genes.map { meta, f -> f }
            )
        }

        // Pathogen TE analysis (dual-genome mode)
        if (params.pathogen_contig_pattern && !params.skip_pathogen_analysis && !params.skip_te_analysis_pathogen) {
            def ch_pathogen_cds_bed = GTF2BED_PATHOGEN.out.bed
            if (ch_pathogen_cds_bed) {
                TE_ANALYSIS_PATHOGEN(
                    ch_pathogen_te_bams,
                    ch_pathogen_cds_bed,
                    ch_contrasts,
                    ch_pathogen_gtf_file
                )
                ch_versions = ch_versions.mix(TE_ANALYSIS_PATHOGEN.out.versions)
            }
        }
    }

    //
    // COLLECT QC STATS: Aggregate per-sample QC metrics into CSVs for downstream plotting
    //
    if (!params.skip_collect_qc_stats) {
        def ch_collect_script = Channel.value(file("${projectDir}/scripts/collect_qc_stats.py", checkIfExists: true))

        COLLECT_QC_STATS(
            ch_qc_star_logs.collect().ifEmpty([]),
            ch_qc_sorf_stats.collect(),
            ch_qc_psites_calcs.collect().ifEmpty([]),
            ch_qc_ribotish_all.collect().ifEmpty([]),
            ch_qc_ribotricer_orfs.collect().ifEmpty([]),
            ch_qc_ribocode_txt.collect().ifEmpty([]),
            ch_qc_orfquant_results.collect().ifEmpty([]),
            ch_collect_script
        )
        ch_versions = ch_versions.mix(COLLECT_QC_STATS.out.versions)
    }

    //
    // ORF QC Module: Unified quality control across all ORF prediction tools.
    // Replaces the previous JOINT_QC_REPORT — provides read-level QC, per-ORF
    // confidence scores, cross-tool comparison, and MultiQC integration.
    // Runs post-unification to harmonize metrics, compute cross-tool agreement,
    // and assign per-ORF confidence scores (OCS).
    //
    // NOTE: On -resume with a newly-added ORF_QC, PoisonPill from unresolved
    // input channels may propagate to MULTIQC. If using -resume after a code
    // update, run a clean start or manually invoke scripts in bin/ instead.
    // See docs/orf_qc_usage.md for manual execution.
    if (!params.skip_orf_qc) {
        def ch_orf_qc_unified = ch_unify_bed.join(ch_unify_metadata)

        if (ch_orf_qc_unified) {
            ORF_QC(
                ch_orf_qc_unified,
                ch_qc_ribocode_txt.collect().ifEmpty([]),
                ch_qc_psites_calcs.collect().ifEmpty([]),
                ch_ribowaltz_psite.map { meta, f -> f }.collect().ifEmpty([]),
                ch_qc_rw_region.map { meta, f -> f }.collect().ifEmpty([]),
                ch_qc_ribotricer_orfs.collect().ifEmpty([]),
                ch_ribotish_predictions.map { meta, f -> f }.collect().ifEmpty([]),
                ch_offset_for_ribotish.map { meta, f -> f }.collect().ifEmpty([]),
                ch_price_gtf.map { meta, f -> f }.collect().ifEmpty([]),
                ch_rpbp_bayes.map { meta, f -> f }.collect().ifEmpty([]),
                ch_orfquant_gtf.map { meta, f -> f }.collect().ifEmpty([])
            )
        } else {
            log.warn "ORF QC module enabled but no unified ORFs available — skipping."
        }
    }

    //
    // ORF Expression Quantification: per-ORF per-sample P-site reads/pN + RPKM/TPM
    // Runs after ORF_QC (needs confidence scores). Queries RiboseQC bedgraphs.
    //
    if (!params.skip_expression_quant) {
        // ORF confidence: from ORF_QC if available, otherwise placeholder
        ch_orf_conf = !params.skip_orf_qc ?
            ORF_QC.out.confidence.map { meta, f -> f }.collect().ifEmpty([]) :
            Channel.value([])

        if (ch_unify_metadata && ch_unify_bed) {
            // Collect P-site + coverage bedgraphs from postfilter RiboseQC
            ch_psites_bg = RIBOSEQC_POSTFILTER.out.psites_bedgraph
                .map { meta, f -> f }.collect().ifEmpty([])

            ch_coverage_bg = RIBOSEQC_POSTFILTER.out.coverage
                .map { meta, f -> f }.collect().ifEmpty([])

            EXPRESSION_QUANT(
                ch_unify_metadata.first(),
                ch_unify_bed.first(),
                ch_orf_conf.first(),
                ch_psites_bg,
                ch_coverage_bg
            )
            ch_versions = ch_versions.mix(EXPRESSION_QUANT.out.versions)

            // Feed expression outputs to MultiQC
            ch_multiqc_files = ch_multiqc_files.mix(
                EXPRESSION_QUANT.out.expression,
                EXPRESSION_QUANT.out.rpkm_tpm
            )
        } else {
            log.warn "EXPRESSION_QUANT requires unified ORFs — skipping."
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
