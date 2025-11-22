//
// Alignment with HISAT2
//

include { HISAT2_ALIGN } from '../../modules/nf-core/hisat2/align/main'
include { BAM_SORT_STATS_SAMTOOLS } from '../nf-core/bam_sort_stats_samtools/main'

workflow FASTQ_ALIGN_HISAT2 {
    take:
    ch_reads          // channel: [ val(meta), [ path(reads) ] ]
    ch_index          // channel: [ val(meta), path(index) ]
    ch_splicesites    // channel: [ val(meta), path(splicesites) ]
    ch_fasta          // channel: [ val(meta), path(fasta) ]

    main:

    ch_versions = Channel.empty()

    //
    // Map reads with HISAT2
    //
    HISAT2_ALIGN ( ch_reads, ch_index, ch_splicesites )
    ch_versions = ch_versions.mix(HISAT2_ALIGN.out.versions)

    //
    // Sort, index BAM file and run samtools stats, flagstat and idxstats
    //
    BAM_SORT_STATS_SAMTOOLS ( HISAT2_ALIGN.out.bam, ch_fasta )
    ch_versions = ch_versions.mix(BAM_SORT_STATS_SAMTOOLS.out.versions)

    emit:
    bam      = BAM_SORT_STATS_SAMTOOLS.out.bam      // channel: [ val(meta), path(bam) ]
    bai      = BAM_SORT_STATS_SAMTOOLS.out.bai      // channel: [ val(meta), path(bai) ]
    stats    = BAM_SORT_STATS_SAMTOOLS.out.stats    // channel: [ val(meta), path(stats) ]
    flagstat = BAM_SORT_STATS_SAMTOOLS.out.flagstat // channel: [ val(meta), path(flagstat) ]
    idxstats = BAM_SORT_STATS_SAMTOOLS.out.idxstats // channel: [ val(meta), path(idxstats) ]
    versions = ch_versions                          // channel: [ path(versions.yml) ]
}
