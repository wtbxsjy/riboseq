//
// Alignment with HISAT2
//

include { HISAT2_ALIGN } from '../../modules/nf-core/hisat2/align/main'
include { HISAT2_ALIGN as HISAT2_ALIGN_TRANSCRIPTOME } from '../../modules/nf-core/hisat2/align/main'
include { BAM_SORT_STATS_SAMTOOLS } from '../nf-core/bam_sort_stats_samtools/main'
include { SAMTOOLS_SORT } from '../../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_INDEX } from '../../modules/nf-core/samtools/index/main'

workflow FASTQ_ALIGN_HISAT2 {
    take:
    ch_reads          // channel: [ val(meta), [ path(reads) ] ]
    ch_index          // channel: [ val(meta), path(index) ]
    ch_transcriptome_index // channel: [ val(meta), path(index) ]
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
    // Map reads with HISAT2 to transcriptome
    //
    HISAT2_ALIGN_TRANSCRIPTOME ( ch_reads, ch_transcriptome_index, Channel.value([[:], []]) )
    ch_versions = ch_versions.mix(HISAT2_ALIGN_TRANSCRIPTOME.out.versions)

    //
    // Sort transcriptome BAM
    //
    SAMTOOLS_SORT ( HISAT2_ALIGN_TRANSCRIPTOME.out.bam, channel.value([[:], []]) )
    ch_versions = ch_versions.mix(SAMTOOLS_SORT.out.versions)

    //
    // Index transcriptome BAM
    //
    SAMTOOLS_INDEX ( SAMTOOLS_SORT.out.bam )
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions)

    //
    // Sort, index BAM file and run samtools stats, flagstat and idxstats
    //
    BAM_SORT_STATS_SAMTOOLS ( HISAT2_ALIGN.out.bam, ch_fasta )
    ch_versions = ch_versions.mix(BAM_SORT_STATS_SAMTOOLS.out.versions)

    emit:
    bam      = BAM_SORT_STATS_SAMTOOLS.out.bam      // channel: [ val(meta), path(bam) ]
    transcriptome_bam = SAMTOOLS_SORT.out.bam       // channel: [ val(meta), path(bam) ]
    transcriptome_bai = SAMTOOLS_INDEX.out.bai      // channel: [ val(meta), path(bai) ]
    bai      = BAM_SORT_STATS_SAMTOOLS.out.bai      // channel: [ val(meta), path(bai) ]
    stats    = BAM_SORT_STATS_SAMTOOLS.out.stats    // channel: [ val(meta), path(stats) ]
    flagstat = BAM_SORT_STATS_SAMTOOLS.out.flagstat // channel: [ val(meta), path(flagstat) ]
    idxstats = BAM_SORT_STATS_SAMTOOLS.out.idxstats // channel: [ val(meta), path(idxstats) ]
    versions = ch_versions                          // channel: [ path(versions.yml) ]
}
