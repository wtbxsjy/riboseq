//
// Subworkflow to run RiboseQC for quality control analysis
//

include { RIBOSEQC_PREPAREANNOTATION } from '../../modules/local/riboseqc/prepareannotation/main'
include { RIBOSEQC_ANALYSIS          } from '../../modules/local/riboseqc/analysis/main'

workflow RIBOSEQC {
    take:
    ch_bam       // channel: [ val(meta), path(bam), path(bai) ] (Genome BAM)
    ch_gtf       // channel: path(gtf)
    ch_fasta     // channel: path(fasta)

    main:
    ch_versions = Channel.empty()

    //
    // Prepare RiboseQC annotation (once per GTF/FASTA combination)
    //
    RIBOSEQC_PREPAREANNOTATION (
        ch_gtf,
        ch_fasta
    )
    ch_versions = ch_versions.mix(RIBOSEQC_PREPAREANNOTATION.out.versions)

    //
    // Run RiboseQC analysis on each BAM
    //
    RIBOSEQC_ANALYSIS (
        ch_bam,
        RIBOSEQC_PREPAREANNOTATION.out.annotation,
        ch_fasta
    )
    ch_versions = ch_versions.mix(RIBOSEQC_ANALYSIS.out.versions)

    emit:
    results  = RIBOSEQC_ANALYSIS.out.results   // channel: [ val(meta), path(results) ]
    html     = RIBOSEQC_ANALYSIS.out.html      // channel: [ val(meta), path(html) ]
    orfquant = RIBOSEQC_ANALYSIS.out.orfquant  // channel: [ val(meta), path(orfquant) ]
    versions = ch_versions                      // channel: [ path(versions.yml) ]
}
