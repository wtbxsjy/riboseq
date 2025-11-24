//
// Check input samplesheet and get read channels
//

include { RPBP_PREPARE_GENOME } from '../../modules/local/rpbp/prepare_genome/main'
include { RPBP_PREDICT        } from '../../modules/local/rpbp/predict/main'

workflow RPBP {
    take:
    bam_bai          // channel: [ val(meta), [ bam ], [ bai ] ]
    fasta            // channel: /path/to/genome.fasta
    gtf              // channel: /path/to/genome.gtf
    ribosomal_fasta  // channel: /path/to/ribosomal.fasta

    main:

    ch_versions = channel.empty()

    //
    // Prepare Genome
    //
    RPBP_PREPARE_GENOME (
        fasta,
        gtf,
        ribosomal_fasta
    )
    ch_versions = ch_versions.mix(RPBP_PREPARE_GENOME.out.versions)

    //
    // Predict
    //
    RPBP_PREDICT (
        bam_bai.map{ meta, bam, bai -> [ meta, bam, bai ] }, // Ensure structure
        RPBP_PREPARE_GENOME.out.orfs_genomic,
        RPBP_PREPARE_GENOME.out.orfs_exons
    )
    ch_versions = ch_versions.mix(RPBP_PREDICT.out.versions)

    emit:
    predictions   = RPBP_PREDICT.out.predictions
    bayes_factors = RPBP_PREDICT.out.bayes_factors
    versions      = ch_versions
}
