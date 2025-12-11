//
// Uncompress and prepare reference genome files
//

include { GUNZIP as GUNZIP_FASTA              } from '../../../modules/nf-core/gunzip'
include { GUNZIP as GUNZIP_GTF                } from '../../../modules/nf-core/gunzip'
include { GUNZIP as GUNZIP_GFF                } from '../../../modules/nf-core/gunzip'
include { GUNZIP as GUNZIP_GENE_BED           } from '../../../modules/nf-core/gunzip'
include { GUNZIP as GUNZIP_TRANSCRIPT_FASTA   } from '../../../modules/nf-core/gunzip'
include { GUNZIP as GUNZIP_ADDITIONAL_FASTA   } from '../../../modules/nf-core/gunzip'
include { GUNZIP as GUNZIP_CONTAM_FASTA       } from '../../../modules/nf-core/gunzip'

include { UNTAR as UNTAR_BOWTIE_INDEX        } from '../../../modules/nf-core/untar'
include { UNTAR as UNTAR_BOWTIE2_INDEX       } from '../../../modules/nf-core/untar'
include { UNTAR as UNTAR_STAR_INDEX          } from '../../../modules/nf-core/untar'
include { UNTAR as UNTAR_SALMON_INDEX        } from '../../../modules/nf-core/untar'
include { UNTAR as UNTAR_HISAT2_INDEX        } from '../../../modules/nf-core/untar'

include { CUSTOM_CATADDITIONALFASTA         } from '../../../modules/nf-core/custom/catadditionalfasta'
include { CUSTOM_GETCHROMSIZES              } from '../../../modules/nf-core/custom/getchromsizes'
include { GFFREAD                           } from '../../../modules/nf-core/gffread'
include { BOWTIE_BUILD                      } from '../../../modules/local/bowtie/build'
include { BOWTIE2_BUILD                     } from '../../../modules/local/bowtie2/build'
include { STAR_GENOMEGENERATE               } from '../../../modules/nf-core/star/genomegenerate'
include { SALMON_INDEX                      } from '../../../modules/nf-core/salmon/index'
include { HISAT2_BUILD                      } from '../../../modules/nf-core/hisat2/build'
include { HISAT2_BUILD as HISAT2_BUILD_TRANSCRIPTOME } from '../../../modules/nf-core/hisat2/build'
include { HISAT2_EXTRACTSPLICESITES         } from '../../../modules/nf-core/hisat2/extractsplicesites'
include { RSEM_PREPAREREFERENCE as RSEM_PREPAREREFERENCE_GENOME } from '../../../modules/nf-core/rsem/preparereference'
include { RSEM_PREPAREREFERENCE as MAKE_TRANSCRIPTS_FASTA       } from '../../../modules/nf-core/rsem/preparereference'

include { PREPROCESS_TRANSCRIPTS_FASTA_GENCODE } from '../../../modules/local/preprocess_transcripts_fasta_gencode'
include { GTF2BED                              } from '../../../modules/local/gtf2bed'
include { GTF_FILTER                           } from '../../../modules/local/gtf_filter'
include { STAR_GENOMEGENERATE_IGENOMES         } from '../../../modules/local/star_genomegenerate_igenomes'

workflow PREPARE_GENOME {
    take:
    fasta                    //      file: /path/to/genome.fasta
    gtf                      //      file: /path/to/genome.gtf
    gff                      //      file: /path/to/genome.gff
    additional_fasta         //      file: /path/to/additional.fasta
    transcript_fasta         //      file: /path/to/transcript.fasta
    contaminant_fasta        //      file: /path/to/contaminants.fasta
    star_index               // directory: /path/to/star/index/
    salmon_index             // directory: /path/to/salmon/index/
    hisat2_index             // directory: /path/to/hisat2/index/
    bowtie_index             // directory: /path/to/bowtie/index/
    bowtie2_index            // directory: /path/to/bowtie2/index/
    gencode                  //   boolean: whether the genome is from GENCODE
    aligner                  //    string: Specifies the alignment algorithm to use - available options are 'star'
    skip_gtf_filter          //   boolean: Skip filtering of GTF for valid scaffolds and/ or transcript IDs
    skip_contaminant_filter  //   boolean: Skip Bowtie/Bowtie2-based contaminant filtering
    filter_aligner           //    string: Specifies contaminant aligner - available options are 'bowtie' or 'bowtie2'
    skip_alignment           //   boolean: Skip all of the alignment-based processes within the pipeline

    main:
    ch_versions = Channel.empty()

    //
    // Uncompress genome fasta file if required
    //
    if (fasta.endsWith('.gz')) {
        ch_fasta    = GUNZIP_FASTA ( [ [:], fasta ] ).gunzip.map { it[1] }
        ch_versions = ch_versions.mix(GUNZIP_FASTA.out.versions)
    } else {
        ch_fasta = Channel.value(file(fasta))
    }

    //
    // Uncompress GTF annotation file or create from GFF3 if required
    //
    if (gtf || gff) {
        if (gtf) {
            if (gtf.endsWith('.gz')) {
                ch_gtf      = GUNZIP_GTF ( [ [:], gtf ] ).gunzip.map { it[1] }
                ch_versions = ch_versions.mix(GUNZIP_GTF.out.versions)
            } else {
                ch_gtf = Channel.value(file(gtf))
            }
        } else if (gff) {
            if (gff.endsWith('.gz')) {
                ch_gff      = GUNZIP_GFF ( [ [:], gff ] ).gunzip.map { it[1] }
                ch_versions = ch_versions.mix(GUNZIP_GFF.out.versions)
            } else {
                ch_gff = Channel.value(file(gff))
            }
            ch_gtf      = GFFREAD ( ch_gff ).gtf
            ch_versions = ch_versions.mix(GFFREAD.out.versions)
        }

        // Determine whether to filter the GTF or not
        def filter_gtf =
            ((
                // Condition 1: Alignment is required and aligner is set
                !skip_alignment && aligner
            ) ||
            (
                // Condition 2: Transcript FASTA file is not provided
                !transcript_fasta
            )) &&
            (
                // Condition 4: --skip_gtf_filter is not provided
                !skip_gtf_filter
            )
        if (filter_gtf) {
            GTF_FILTER ( ch_fasta, ch_gtf )
            ch_gtf = GTF_FILTER.out.genome_gtf
            ch_versions = ch_versions.mix(GTF_FILTER.out.versions)
        }
    }

    //
    // Uncompress additional fasta file and concatenate with reference fasta and gtf files
    //
    def biotype = gencode ? "gene_type" : "gene_biotype"
    if (additional_fasta) {
        if (additional_fasta.endsWith('.gz')) {
            ch_add_fasta = GUNZIP_ADDITIONAL_FASTA ( [ [:], additional_fasta ] ).gunzip.map { it[1] }
            ch_versions  = ch_versions.mix(GUNZIP_ADDITIONAL_FASTA.out.versions)
        } else {
            ch_add_fasta = Channel.value(file(additional_fasta))
        }

        CUSTOM_CATADDITIONALFASTA(
            ch_fasta.combine(ch_gtf).map{fasta, gtf -> [[:], fasta, gtf]},
            ch_add_fasta.map{[[:], it]},
            biotype
        )
        ch_fasta    = CUSTOM_CATADDITIONALFASTA.out.fasta.map{it[1]}.first()
        ch_gtf      = CUSTOM_CATADDITIONALFASTA.out.gtf.map{it[1]}.first()
        ch_versions = ch_versions.mix(CUSTOM_CATADDITIONALFASTA.out.versions)
    }

    //
    // Uncompress transcript fasta file / create if required
    //
    if (transcript_fasta) {
        if (transcript_fasta.endsWith('.gz')) {
            ch_transcript_fasta = GUNZIP_TRANSCRIPT_FASTA ( [ [:], transcript_fasta ] ).gunzip.map { it[1] }
            ch_versions         = ch_versions.mix(GUNZIP_TRANSCRIPT_FASTA.out.versions)
        } else {
            ch_transcript_fasta = Channel.value(file(transcript_fasta))
        }
        if (gencode) {
            PREPROCESS_TRANSCRIPTS_FASTA_GENCODE ( ch_transcript_fasta )
            ch_transcript_fasta = PREPROCESS_TRANSCRIPTS_FASTA_GENCODE.out.fasta
            ch_versions         = ch_versions.mix(PREPROCESS_TRANSCRIPTS_FASTA_GENCODE.out.versions)
        }
    } else {
        ch_transcript_fasta = MAKE_TRANSCRIPTS_FASTA ( ch_fasta, ch_gtf ).transcript_fasta
        ch_versions         = ch_versions.mix(MAKE_TRANSCRIPTS_FASTA.out.versions)
    }

    //
    // Create chromosome sizes file
    //
    CUSTOM_GETCHROMSIZES ( ch_fasta.map { [ [:], it ] } )
    ch_fai         = CUSTOM_GETCHROMSIZES.out.fai.map { it[1] }
    ch_chrom_sizes = CUSTOM_GETCHROMSIZES.out.sizes.map { it[1] }
    ch_versions    = ch_versions.mix(CUSTOM_GETCHROMSIZES.out.versions)

    //
    // Get list of indices that need to be created
    //
    def prepare_tool_indices = []
    if (!skip_alignment && aligner) { prepare_tool_indices << aligner }

    //
    // Prepare contaminant index for Bowtie/Bowtie2 filtering if required
    //
    ch_contaminant_index = Channel.empty()
    if (!skip_contaminant_filter) {
        def contaminantAligner = (filter_aligner ?: 'bowtie').toLowerCase()
        if (!(contaminantAligner in ['bowtie', 'bowtie2'])) {
            exit 1, "Unsupported --filter_aligner '${filter_aligner}'. Choose 'bowtie' or 'bowtie2'."
        }

        def provided_index = contaminantAligner == 'bowtie2' ? bowtie2_index : bowtie_index
        def provided_index_path = provided_index ? provided_index.toString() : null
        if (!provided_index_path && !contaminant_fasta) {
            exit 1, 'Contaminant filtering requested but neither --contaminant_fasta nor a pre-built index was provided.'
        }

        if (provided_index_path) {
            if (provided_index_path.endsWith('.tar.gz')) {
                if (contaminantAligner == 'bowtie2') {
                    ch_contaminant_index = UNTAR_BOWTIE2_INDEX ( [ [:], provided_index_path ] ).untar.map { it[1] }
                    ch_versions          = ch_versions.mix(UNTAR_BOWTIE2_INDEX.out.versions)
                } else {
                    ch_contaminant_index = UNTAR_BOWTIE_INDEX ( [ [:], provided_index_path ] ).untar.map { it[1] }
                    ch_versions          = ch_versions.mix(UNTAR_BOWTIE_INDEX.out.versions)
                }
            } else {
                ch_contaminant_index = Channel.value(file(provided_index_path, checkIfExists: true))
            }
        } else {
            if (!contaminant_fasta) {
                exit 1, 'Contaminant FASTA is required to build Bowtie/Bowtie2 index.'
            }

            // Work on a local copy of the path to avoid shadowing the workflow input name
            def contaminant_path = contaminant_fasta
            def contaminant_file = file(contaminant_path, checkIfExists: true)
            def contaminant_channel
            if (contaminant_path.endsWith('.gz')) {
                contaminant_channel = GUNZIP_CONTAM_FASTA ( [ [:], contaminant_file ] ).gunzip.map { it[1] }
                ch_versions         = ch_versions.mix(GUNZIP_CONTAM_FASTA.out.versions)
            } else {
                contaminant_channel = Channel.value(contaminant_file)
            }

            def meta = [ id: 'contaminants' ]
            if (contaminantAligner == 'bowtie2') {
                BOWTIE2_BUILD (
                    contaminant_channel.map { [ meta, it ] }
                )
                ch_contaminant_index = BOWTIE2_BUILD.out.index.map { it[1] }
                ch_versions          = ch_versions.mix(BOWTIE2_BUILD.out.versions)
            } else {
                BOWTIE_BUILD (
                    contaminant_channel.map { [ meta, it ] }
                )
                ch_contaminant_index = BOWTIE_BUILD.out.index
                ch_versions          = ch_versions.mix(BOWTIE_BUILD.out.versions)
            }
        }
    }

    //
    // Uncompress STAR index or generate from scratch if required
    //
    ch_star_index = Channel.empty()
    if ('star' in prepare_tool_indices) {
        if (star_index) {
            if (star_index.endsWith('.tar.gz')) {
                ch_star_index = UNTAR_STAR_INDEX ( [ [:], star_index ] ).untar.map { it[1] }
                ch_versions   = ch_versions.mix(UNTAR_STAR_INDEX.out.versions)
            } else {
                ch_star_index = Channel.value(file(star_index))
            }
        } else {
            // Check if an AWS iGenome has been provided to use the appropriate version of STAR
            def is_aws_igenome = false
            if (fasta && gtf) {
                if ((file(fasta).getName() - '.gz' == 'genome.fa') && (file(gtf).getName() - '.gz' == 'genes.gtf')) {
                    is_aws_igenome = true
                }
            }
            if (is_aws_igenome) {
                ch_star_index = STAR_GENOMEGENERATE_IGENOMES ( ch_fasta, ch_gtf ).index
                ch_versions   = ch_versions.mix(STAR_GENOMEGENERATE_IGENOMES.out.versions)
            } else {
                ch_star_index = STAR_GENOMEGENERATE ( ch_fasta.map { [ [:], it ] }, ch_gtf.map { [ [:], it ] } ).index.map { it[1] }
                ch_versions   = ch_versions.mix(STAR_GENOMEGENERATE.out.versions)
            }
        }
    }

    //
    // Uncompress HISAT2 index or generate from scratch if required
    //
    ch_hisat2_index = Channel.empty()
    ch_hisat2_transcriptome_index = Channel.empty()
    if ('hisat2' in prepare_tool_indices) {
        // Handle genome index: use provided or build from scratch
        if (hisat2_index) {
            if (hisat2_index.endsWith('.tar.gz')) {
                ch_hisat2_index = UNTAR_HISAT2_INDEX ( [ [:], hisat2_index ] ).untar.map { it[1] }
                ch_versions     = ch_versions.mix(UNTAR_HISAT2_INDEX.out.versions)
            } else {
                ch_hisat2_index = Channel.value(file(hisat2_index))
            }
        } else {
            HISAT2_EXTRACTSPLICESITES ( ch_gtf.map { [ [:], it ] } )
            ch_splicesites  = HISAT2_EXTRACTSPLICESITES.out.txt.map { it[1] }
            ch_versions     = ch_versions.mix(HISAT2_EXTRACTSPLICESITES.out.versions)

            HISAT2_BUILD ( ch_fasta.map { [ [:], it ] }, ch_gtf.map { [ [:], it ] }, ch_splicesites.map { [ [:], it ] } )
            ch_hisat2_index = HISAT2_BUILD.out.index.map { it[1] }
            ch_versions     = ch_versions.mix(HISAT2_BUILD.out.versions)
        }

        // Always build transcriptome index for HISAT2 (regardless of whether genome index was provided)
        // Users rarely provide transcriptome index from command line, so we always build it
        HISAT2_BUILD_TRANSCRIPTOME ( ch_transcript_fasta.map { [ [:], it ] }, ch_gtf.map { [ [:], it ] }, Channel.value([[:], []]) )
        ch_hisat2_transcriptome_index = HISAT2_BUILD_TRANSCRIPTOME.out.index.map { it[1] }
        ch_versions     = ch_versions.mix(HISAT2_BUILD_TRANSCRIPTOME.out.versions)
    }

    //
    // Uncompress Salmon index or generate from scratch if required
    //
    ch_salmon_index = Channel.empty()
    if (salmon_index) {
        if (salmon_index.endsWith('.tar.gz')) {
            ch_salmon_index = UNTAR_SALMON_INDEX ( [ [:], salmon_index ] ).untar.map { it[1] }
            ch_versions     = ch_versions.mix(UNTAR_SALMON_INDEX.out.versions)
        } else {
            ch_salmon_index = Channel.value(file(salmon_index))
        }
    } else {
        if ('salmon' in prepare_tool_indices) {
            ch_salmon_index = SALMON_INDEX ( ch_fasta, ch_transcript_fasta ).index
            ch_versions     = ch_versions.mix(SALMON_INDEX.out.versions)
        }
    }

    emit:
    fasta            = ch_fasta                  // channel: path(genome.fasta)
    gtf              = ch_gtf                    // channel: path(genome.gtf)
    fai              = ch_fai                    // channel: path(genome.fai)
    transcript_fasta = ch_transcript_fasta       // channel: path(transcript.fasta)
    chrom_sizes      = ch_chrom_sizes            // channel: path(genome.sizes)
    contaminant_index = ch_contaminant_index.first() // channel: path(contaminant/index/)
    star_index       = ch_star_index             // channel: path(star/index/)
    hisat2_index     = ch_hisat2_index           // channel: path(hisat2/index/)
    hisat2_transcriptome_index = ch_hisat2_transcriptome_index // channel: path(hisat2/transcriptome_index/)
    salmon_index     = ch_salmon_index           // channel: path(salmon/index/)
    versions         = ch_versions.ifEmpty(null) // channel: [ versions.yml ]
}
