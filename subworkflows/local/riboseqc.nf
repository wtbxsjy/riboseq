//
// Subworkflow to run RiboseQC for quality control analysis
// Includes optional P-site offset correction for improved ORFquant accuracy
//

include { RIBOSEQC_PREPAREANNOTATION     } from '../../modules/local/riboseqc/prepareannotation/main'
include { RIBOSEQC_ANALYSIS              } from '../../modules/local/riboseqc/analysis/main'
include { EXTRACT_RL_CUTOFF              } from '../../modules/local/extract_rl_cutoff/main'
include { PREPARE_FOR_ORFQUANT_CORRECTED } from '../../modules/local/prepare_for_orfquant/main'

workflow RIBOSEQC {
    take:
    ch_bam               // channel: [ val(meta), path(bam), path(bai) ] (Genome BAM)
    ch_gtf               // channel: path(gtf)
    ch_fasta             // channel: path(fasta)
    ch_ribowaltz_psite   // channel: [ val(meta), path(psite_offset.tsv) ] - optional fallback for offset correction

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

    // Initialize output channels
    ch_orfquant_final = Channel.empty()
    ch_rl_cutoff      = Channel.empty()

    //
    // P-site offset correction for ORFquant (enabled by default via params.orfquant_psite_correction)
    // Extracts rows where max_coverage is TRUE from P_sites_calcs and regenerates for_ORFquant file
    // When RiboseQC P_sites_calcs is empty/invalid, falls back to riboWaltz-derived offsets
    //
    if (params.orfquant_psite_correction) {
        // Placeholder file used when riboWaltz offsets are unavailable.
        // Must exist (Nextflow validates input files) but is never read because use_rw=false.
        def placeholder = file("${projectDir}/assets/samplesheet.csv")

        if (params.skip_ribowaltz) {
            // No riboWaltz: attach dummy fallback to every sample (no join needed)
            ch_rl_inputs = RIBOSEQC_ANALYSIS.out.psites_calcs
                .map { meta, psites ->
                    tuple(meta, psites, [id: '_NO_RW_'], placeholder, false)
                }
        } else {
            // riboWaltz available: join by sample ID for per-sample offset correction.
            // Unmatched samples (no riboWaltz data) get placeholder fallback via remainder.
            def ch_rw_keyed = ch_ribowaltz_psite
                .map { meta, f -> tuple(meta.id, meta, f) }

            def ch_psites_keyed = RIBOSEQC_ANALYSIS.out.psites_calcs
                .map { meta, f -> tuple(meta.id, meta, f) }

            ch_rl_inputs = ch_psites_keyed
                .join(ch_rw_keyed, remainder: true, by: 0)
                .map { id, meta_qc, psites, meta_rw, rw_file ->
                    def has_rw = (meta_rw != null && meta_rw.id != '_NO_RW_' && rw_file != null)
                    if (!has_rw) {
                        return tuple(meta_qc, psites, [id: '_NO_RW_'], placeholder, false)
                    }
                    return tuple(meta_qc, psites, meta_rw, rw_file, true)
                }
        }

        //
        // Extract read length to P-site offset (cutoff) mapping from RiboseQC results
        // Uses riboWaltz offsets as fallback when RiboseQC data is invalid
        //
        EXTRACT_RL_CUTOFF (
            ch_rl_inputs.map { meta, psites, meta_rw, rw_file, use_rw ->
                [ meta, psites ]
            },
            ch_rl_inputs.map { meta, psites, meta_rw, rw_file, use_rw ->
                [ meta_rw, rw_file, use_rw ]
            }
        )
        ch_versions = ch_versions.mix(EXTRACT_RL_CUTOFF.out.versions)
        ch_rl_cutoff = EXTRACT_RL_CUTOFF.out.rl_cutoff

        //
        // Regenerate for_ORFquant file with corrected P-site offsets
        //
        // Join BAM channel with rl_cutoff by sample id
        ch_bam_for_correction = ch_bam
            .map { meta, bam, bai -> [ meta.id, meta, bam, bai ] }
        ch_rl_cutoff_keyed = EXTRACT_RL_CUTOFF.out.rl_cutoff
            .map { meta, rl_cutoff -> [ meta.id, meta, rl_cutoff ] }

        ch_prepare_inputs = ch_bam_for_correction
            .join(ch_rl_cutoff_keyed, by: 0)
            .map { id, meta_bam, bam, bai, meta_rl, rl_cutoff ->
                [ meta_bam, bam, bai, meta_rl, rl_cutoff ]
            }

        PREPARE_FOR_ORFQUANT_CORRECTED (
            ch_prepare_inputs.map { meta, bam, bai, meta_rl, rl_cutoff -> [ meta, bam, bai ] },
            RIBOSEQC_PREPAREANNOTATION.out.annotation,
            ch_prepare_inputs.map { meta, bam, bai, meta_rl, rl_cutoff -> [ meta_rl, rl_cutoff ] }
        )
        ch_versions = ch_versions.mix(PREPARE_FOR_ORFQUANT_CORRECTED.out.versions)

        // Use corrected for_ORFquant file
        ch_orfquant_final = PREPARE_FOR_ORFQUANT_CORRECTED.out.for_orfquant
    } else {
        // Use original for_ORFquant file from RiboseQC
        ch_orfquant_final = RIBOSEQC_ANALYSIS.out.orfquant
    }

    emit:
    annotation      = RIBOSEQC_PREPAREANNOTATION.out.annotation  // channel: path(annotation) - *_Rannot file
    results         = RIBOSEQC_ANALYSIS.out.results              // channel: [ val(meta), path(results) ]
    results_all     = RIBOSEQC_ANALYSIS.out.results_all          // channel: [ val(meta), path(results_all) ]
    orfquant        = ch_orfquant_final                          // channel: [ val(meta), path(for_orfquant) ] - Corrected if enabled
    orfquant_orig   = RIBOSEQC_ANALYSIS.out.orfquant             // channel: [ val(meta), path(orfquant) ] - Original
    coverage        = RIBOSEQC_ANALYSIS.out.coverage             // channel: [ val(meta), path(bedgraph) ]
    psites_bedgraph = RIBOSEQC_ANALYSIS.out.psites_bedgraph      // channel: [ val(meta), path(bedgraph) ]
    psites_calcs    = RIBOSEQC_ANALYSIS.out.psites_calcs         // channel: [ val(meta), path(psites_calcs) ]
    rl_cutoff       = ch_rl_cutoff                               // channel: [ val(meta), path(rl_cutoff) ] - Only if correction enabled
    junctions       = RIBOSEQC_ANALYSIS.out.junctions            // channel: [ val(meta), path(junctions) ]
    versions        = ch_versions                                 // channel: [ path(versions.yml) ]
}
