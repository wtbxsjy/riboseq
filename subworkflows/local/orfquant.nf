//
// Subworkflow to run ORFquant for ORF detection and quantification
// ORFquant uses output from RiboseQC for P-site analysis
//

include { ORFQUANT_RUN } from '../../modules/local/orfquant/main'

workflow ORFQUANT {
    take:
    ch_for_orfquant   // channel: [ val(meta), path(for_orfquant) ] - from RiboseQC
    ch_annotation     // channel: path(annotation) - *_Rannot file from RiboseQC
    ch_fasta          // channel: path(fasta) - Genome fasta file

    main:
    ch_versions = Channel.empty()

    //
    // Run ORFquant analysis
    //
    ORFQUANT_RUN (
        ch_for_orfquant,
        ch_annotation,
        ch_fasta
    )
    ch_versions = ch_versions.mix(ORFQUANT_RUN.out.versions)

    emit:
    results      = ORFQUANT_RUN.out.results       // channel: [ val(meta), path(results) ]
    gtf          = ORFQUANT_RUN.out.gtf           // channel: [ val(meta), path(gtf) ]
    proteins     = ORFQUANT_RUN.out.proteins      // channel: [ val(meta), path(fasta) ]
    tmp_results  = ORFQUANT_RUN.out.tmp_results   // channel: [ val(meta), path(tmp_results) ]
    plots_data   = ORFQUANT_RUN.out.plots_data    // channel: [ val(meta), path(plots_data) ]
    plots_dir    = ORFQUANT_RUN.out.plots_dir     // channel: [ val(meta), path(plots_dir) ]
    versions     = ch_versions                     // channel: [ path(versions.yml) ]
}
