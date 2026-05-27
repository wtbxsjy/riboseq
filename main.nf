#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf-core/riboseq
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/nf-core/riboseq
    Website: https://nf-co.re/riboseq
    Slack  : https://nfcore.slack.com/channels/riboseq
----------------------------------------------------------------------------------------
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { RIBOSEQ                 } from './workflows/riboseq'
include { PREPARE_GENOME          } from './subworkflows/local/prepare_genome'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_riboseq_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_riboseq_pipeline'
include { getGenomeAttribute      } from './subworkflows/local/utils_nfcore_riboseq_pipeline'
include { checkMaxContigSize      } from './subworkflows/local/utils_nfcore_riboseq_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GENOME PARAMETER VALUES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

params.fasta             = getGenomeAttribute('fasta')
params.transcript_fasta  = getGenomeAttribute('transcript_fasta')
params.additional_fasta  = getGenomeAttribute('additional_fasta')
params.gtf               = getGenomeAttribute('gtf')
params.gff               = getGenomeAttribute('gff')
params.contaminant_fasta = getGenomeAttribute('contaminant_fasta')
params.bowtie_index      = getGenomeAttribute('bowtie')
params.bowtie2_index     = getGenomeAttribute('bowtie2')
params.star_index        = getGenomeAttribute('star')
params.hisat2_index      = getGenomeAttribute('hisat2')
params.salmon_index      = getGenomeAttribute('salmon')

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Run main analysis pipeline depending on type of input
//
workflow NFCORE_RIBOSEQ {

    take:
    ch_samplesheet  // channel: [ meta, bam, bai ] for BAM input OR [ meta, [ fastqs ] ] for FASTQ input

    main:

    // Get BAM input mode from global params (set by PIPELINE_INITIALISATION)
    def is_bam_input = params.is_bam_input ?: false

    ch_versions = channel.empty()

    // Validate contaminant filtering parameters
    if (!params.skip_contaminant_filter && !params.contaminant_fasta) {
        error "Contaminant filtering is enabled but no contaminant FASTA file provided. Please specify --contaminant_fasta or --skip_contaminant_filter."
    }

    // Validate pathogen dual-genome analysis parameters
    if (params.pathogen_contig_pattern && !params.skip_pathogen_analysis) {
        if (!params.pathogen_fasta) {
            error "--pathogen_contig_pattern requires --pathogen_fasta to be set."
        }
        if (!params.pathogen_gtf) {
            error "--pathogen_contig_pattern requires --pathogen_gtf to be set."
        }
    }

    //
    // SUBWORKFLOW: Prepare reference genome files
    //
    PREPARE_GENOME (
        params.fasta,
        params.gtf,
        params.gff,
        params.additional_fasta,
        params.transcript_fasta,
        params.contaminant_fasta,
        params.star_index,
        params.salmon_index,
        params.hisat2_index,
        params.bowtie_index,
        params.bowtie2_index,
        params.gencode,
        params.aligner,
        params.skip_gtf_filter,
        params.skip_contaminant_filter,
        params.filter_aligner,
        params.skip_alignment
    )
    ch_versions = ch_versions.mix(PREPARE_GENOME.out.versions)

    // Check if contigs in genome fasta file > 512 Mbp
    if (!params.skip_alignment && !params.bam_csi_index) {
        PREPARE_GENOME
            .out
            .fai
            .map { checkMaxContigSize(it) }
    }

    //
    // WORKFLOW: Run nf-core/riboseq workflow
    //
    //
    // WORKFLOW: Run nf-core/riboseq workflow
    //

    RIBOSEQ (
        ch_samplesheet,
        ch_versions,
        PREPARE_GENOME.out.fasta,
        PREPARE_GENOME.out.gtf,
        PREPARE_GENOME.out.fai,
        PREPARE_GENOME.out.chrom_sizes,
        PREPARE_GENOME.out.transcript_fasta,
        PREPARE_GENOME.out.star_index,
        PREPARE_GENOME.out.hisat2_index,
        PREPARE_GENOME.out.hisat2_transcriptome_index,
        PREPARE_GENOME.out.salmon_index,
        PREPARE_GENOME.out.contaminant_index
    )
    ch_versions = ch_versions.mix(RIBOSEQ.out.versions)

    emit:
    multiqc_report = RIBOSEQ.out.multiqc_report // channel: /path/to/multiqc_report.html
    versions       = ch_versions                // channel: [version1, version2, ...]


}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    main:

    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    PIPELINE_INITIALISATION (
        params.version,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        params.input
    )

    //
    // WORKFLOW: Run main workflow
    //
    NFCORE_RIBOSEQ (
        PIPELINE_INITIALISATION.out.samplesheet
    )

    //
    // SUBWORKFLOW: Run completion tasks
    //
    PIPELINE_COMPLETION (
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.monochrome_logs,
        params.hook_url,
        NFCORE_RIBOSEQ.out.multiqc_report
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
