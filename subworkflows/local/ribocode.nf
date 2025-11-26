//
// Subworkflow to run RiboCode for ORF prediction
//

include { RIBOCODE_DETECT  } from '../../modules/local/ribocode/detect/main'

workflow RIBOCODE {
    take:
    ch_bam       // channel: [ val(meta), path(bam) ] (Transcriptome BAM)
    ch_gtf       // channel: path(gtf)
    ch_fasta     // channel: path(fasta)

    main:
    ch_versions = Channel.empty()

    //
    // Run RiboCode detection
    //
    RIBOCODE_DETECT (
        ch_bam,
        ch_gtf,
        ch_fasta
    )
    ch_versions = ch_versions.mix(RIBOCODE_DETECT.out.versions)

    emit:
    results  = RIBOCODE_DETECT.out.results
    versions = ch_versions
}
