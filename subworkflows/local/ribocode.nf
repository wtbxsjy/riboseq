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
    txt       = RIBOCODE_DETECT.out.txt
    collapsed = RIBOCODE_DETECT.out.collapsed
    gtf       = RIBOCODE_DETECT.out.gtf
    bed       = RIBOCODE_DETECT.out.bed
    results   = RIBOCODE_DETECT.out.results
    versions  = ch_versions
}
