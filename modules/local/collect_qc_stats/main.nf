process COLLECT_QC_STATS {
    tag "qc_stats"
    label 'process_single'

    publishDir "${params.outdir}/qc_stats", mode: params.publish_dir_mode

    conda "conda-forge::python=3.9 conda-forge::biopython=1.79"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.79' :
        'quay.io/biocontainers/biopython:1.79' }"

    input:
    path star_logs,        stageAs: 'star_logs/*'        // STAR *Log.final.out per sample (optional)
    path sorf_stats,       stageAs: 'sorf/*'             // *.sorf.filter_stats.tsv per sample
    path psites_calcs,     stageAs: 'psites/*'           // *_P_sites_calcs per sample (optional)
    path ribotish_all,     stageAs: 'ribotish/*'         // *_all.txt per sample (optional)
    path ribotricer_orfs,  stageAs: 'ribotricer/*'       // *_translating_ORFs.tsv per sample (optional)
    path orfquant_results, stageAs: 'orfquant/*'         // *_final_ORFquant_results per sample (optional)
    path collect_script

    output:
    path "alignment_stats.csv",         emit: alignment_stats,  optional: true
    path "sorf_stats.csv",              emit: sorf_stats
    path "psite_periodicity_stats.csv", emit: psite_stats,      optional: true
    path "orf_counts.csv",              emit: orf_counts,        optional: true
    path "qc_summary.csv",              emit: qc_summary
    path "versions.yml",                emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def star_arg      = star_logs instanceof List && star_logs.size() > 0 ?
                            "--star_logs ${star_logs.collect{ "star_logs/${it.name}" }.join(' ')}" : ''
    def sorf_arg      = sorf_stats instanceof List && sorf_stats.size() > 0 ?
                            "--sorf_stats ${sorf_stats.collect{ "sorf/${it.name}" }.join(' ')}" :
                            sorf_stats && sorf_stats.name != 'NO_FILE' ?
                            "--sorf_stats sorf/${sorf_stats.name}" : ''
    def psite_arg     = psites_calcs instanceof List && psites_calcs.size() > 0 ?
                            "--psites_calcs ${psites_calcs.collect{ "psites/${it.name}" }.join(' ')}" : ''
    def rtish_arg     = ribotish_all instanceof List && ribotish_all.size() > 0 ?
                            "--ribotish_all ${ribotish_all.collect{ "ribotish/${it.name}" }.join(' ')}" : ''
    def rtricer_arg   = ribotricer_orfs instanceof List && ribotricer_orfs.size() > 0 ?
                            "--ribotricer_orfs ${ribotricer_orfs.collect{ "ribotricer/${it.name}" }.join(' ')}" : ''
    def orfquant_arg  = orfquant_results instanceof List && orfquant_results.size() > 0 ?
                            "--orfquant ${orfquant_results.collect{ "orfquant/${it.name}" }.join(' ')}" : ''
    """
    python3 ${collect_script} \\
        ${star_arg} \\
        ${sorf_arg} \\
        ${psite_arg} \\
        ${rtish_arg} \\
        ${rtricer_arg} \\
        ${orfquant_arg} \\
        --output_dir .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    touch alignment_stats.csv sorf_stats.csv psite_periodicity_stats.csv orf_counts.csv qc_summary.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: "3.9.0"
    END_VERSIONS
    """
}
