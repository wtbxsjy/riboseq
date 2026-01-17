//
// Subworkflow: Unified ORF annotation using gencode-riboseqORFs
// Integrates ORF predictions from multiple tools and maps to GENCODE/Ensembl
//

include { DOWNLOAD_ENSEMBL_FILES        } from '../../modules/local/prepare_ensembl_annotation/main'
include { CALCULATE_PSITE_BED           } from '../../modules/local/prepare_ensembl_annotation/main'
include { CONVERT_RIBOTISH_TO_GENCODE   } from '../../modules/local/convert_ribotish_to_gencode/main'
include { CONVERT_RIBOTRICER_TO_GENCODE } from '../../modules/local/convert_ribotricer_to_gencode/main'
include { GENCODE_ORF_MAPPER            } from '../../modules/local/gencode_orf_mapper/main'

workflow GENCODE_ORF_ANNOTATION {
    take:
    ch_ribotish_orfs      // channel: [ meta, predict, quality ]
    ch_ribotricer_orfs    // channel: [ meta, orfs ]
    ch_fasta              // channel: path(genome.fasta)
    ch_gtf                // channel: path(genome.gtf)
    ensembl_release       // val: Ensembl release number
    genome_assembly       // val: Genome assembly (e.g., GRCh38)
    species               // val: Species name (e.g., homo_sapiens)
    project_id            // val: Project identifier

    main:
    ch_versions = Channel.empty()

    //
    // Prepare Ensembl annotation (run once per pipeline execution)
    //
    DOWNLOAD_ENSEMBL_FILES(
        ensembl_release,
        genome_assembly,
        species
    )
    ch_versions = ch_versions.mix(DOWNLOAD_ENSEMBL_FILES.out.versions)

    CALCULATE_PSITE_BED(
        DOWNLOAD_ENSEMBL_FILES.out.ensembl_dir,
        ensembl_release,
        genome_assembly,
        species
    )
    ch_versions = ch_versions.mix(CALCULATE_PSITE_BED.out.versions)

    //
    // Convert ORF formats to gencode-compatible format
    //
    ch_converted_orfs_fasta = Channel.empty()
    ch_converted_orfs_bed = Channel.empty()

    // Convert Ribo-TISH predictions
    if (ch_ribotish_orfs) {
        CONVERT_RIBOTISH_TO_GENCODE(
            ch_ribotish_orfs,
            ch_fasta,
            ch_gtf
        )
        ch_converted_orfs_fasta = ch_converted_orfs_fasta.mix(
            CONVERT_RIBOTISH_TO_GENCODE.out.fasta.map { meta, fa -> fa }
        )
        ch_converted_orfs_bed = ch_converted_orfs_bed.mix(
            CONVERT_RIBOTISH_TO_GENCODE.out.bed.map { meta, bed -> bed }
        )
        ch_versions = ch_versions.mix(CONVERT_RIBOTISH_TO_GENCODE.out.versions.first())
    }

    // Convert Ribotricer predictions
    if (ch_ribotricer_orfs) {
        CONVERT_RIBOTRICER_TO_GENCODE(
            ch_ribotricer_orfs,
            ch_fasta,
            ch_gtf
        )
        ch_converted_orfs_fasta = ch_converted_orfs_fasta.mix(
            CONVERT_RIBOTRICER_TO_GENCODE.out.fasta.map { meta, fa -> fa }
        )
        ch_converted_orfs_bed = ch_converted_orfs_bed.mix(
            CONVERT_RIBOTRICER_TO_GENCODE.out.bed.map { meta, bed -> bed }
        )
        ch_versions = ch_versions.mix(CONVERT_RIBOTRICER_TO_GENCODE.out.versions.first())
    }

    //
    // Merge all ORF files
    //
    ch_merged_fasta = ch_converted_orfs_fasta
        .collectFile(
            name: "${project_id}_all_orfs.fa",
            newLine: false,
            storeDir: params.outdir ? "${params.outdir}/gencode_annotation/merged" : null
        )

    ch_merged_bed = ch_converted_orfs_bed
        .collectFile(
            name: "${project_id}_all_orfs.bed",
            newLine: true,
            storeDir: params.outdir ? "${params.outdir}/gencode_annotation/merged" : null
        )

    //
    // Run GENCODE ORF mapper for unified annotation
    //
    GENCODE_ORF_MAPPER(
        ch_merged_fasta,
        ch_merged_bed,
        CALCULATE_PSITE_BED.out.ensembl_dir_complete,
        project_id
    )
    ch_versions = ch_versions.mix(GENCODE_ORF_MAPPER.out.versions)

    emit:
    unified_fasta = GENCODE_ORF_MAPPER.out.unified_fasta       // channel: path(*.orfs.fa)
    unified_bed   = GENCODE_ORF_MAPPER.out.unified_bed         // channel: path(*.orfs.bed)
    unified_gtf   = GENCODE_ORF_MAPPER.out.unified_gtf         // channel: path(*.orfs.gtf)
    unified_table = GENCODE_ORF_MAPPER.out.unified_table       // channel: path(*.orfs.out)
    altmapped     = GENCODE_ORF_MAPPER.out.altmapped           // channel: path(*.altmapped)
    unmapped      = GENCODE_ORF_MAPPER.out.unmapped            // channel: path(*.unmapped)
    merged_fasta  = ch_merged_fasta                             // channel: path(*_all_orfs.fa)
    merged_bed    = ch_merged_bed                               // channel: path(*_all_orfs.bed)
    versions      = ch_versions                                 // channel: path(versions.yml)
}
