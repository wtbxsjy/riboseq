/*
 * Joint Ribo-seq QC Report
 *
 * Integrates quality metrics from multiple tools:
 *   - riboWaltz: P-site offsets, region distribution
 *   - RiboseQC: P-site periodicity (frame_preference)
 *   - Ribo-TISH: read length distribution, metagene peaks
 *   - Ribotricer: read length distribution, total reads
 *
 * Generates a comprehensive per-sample quality assessment with letter grades
 * and MultiQC-compatible output for integration into the final report.
 */

process JOINT_QC_REPORT {
    tag 'joint_qc_report'
    label 'process_single'

    conda "conda-forge::r-base=4.1 conda-forge::r-data.table=1.14"
    // Reuse ribowaltz container (already has data.table + jsonlite).
    // Falls back to Docker rocker/r-base if ribowaltz_container is not set.
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        (params.ribowaltz_container ?: 'oras://community.wave.seqera.io/library/ribowaltz_bioconductor-genomicfeatures_r-data.table:latest') :
        (params.ribowaltz_container ?: 'community.wave.seqera.io/library/ribowaltz_bioconductor-genomicfeatures_r-data.table:latest') }"

    input:
    path rw_offsets         // riboWaltz *_psite_offset.tsv files
    path rw_regions         // riboWaltz *_region_distribution.tsv files
    path rq_psites          // RiboseQC *_P_sites_calcs files
    path rt_qual            // Ribo-TISH *_qual.txt files
    path rtr_summaries      // Ribotricer *_bam_summary.txt files
    path report_script

    output:
    path "joint_riboseq_qc.tsv"          , emit: summary
    path "joint_riboseq_qc_mqc.yaml"     , emit: mqc_yaml
    path "joint_riboseq_qc_mqc.txt"      , emit: mqc_data
    path "versions.yml"                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def rw_off_dir   = rw_offsets instanceof List && rw_offsets.size() > 0 ? 'rw_offset' : ''
    def rw_reg_dir   = rw_regions instanceof List && rw_regions.size() > 0 ? 'rw_region' : ''
    def rq_dir       = rq_psites instanceof List && rq_psites.size() > 0 ? 'rq_psites' : ''
    def rt_dir       = rt_qual instanceof List && rt_qual.size() > 0 ? 'rt_qual' : ''
    def rtr_dir      = rtr_summaries instanceof List && rtr_summaries.size() > 0 ? 'rtr_summary' : ''

    def rw_off_arg   = rw_off_dir ? "--rw_offset_dir ${rw_off_dir}" : ''
    def rw_reg_arg   = rw_reg_dir ? "--rw_region_dir ${rw_reg_dir}" : ''
    def rq_arg        = rq_dir ? "--rq_psites_dir ${rq_dir}" : ''
    def rt_arg        = rt_dir ? "--rt_qual_dir ${rt_dir}" : ''
    def rtr_arg       = rtr_dir ? "--rtr_summary_dir ${rtr_dir}" : ''
    """
    # Stage files into per-tool directories
    mkdir -p ${rw_off_dir} ${rw_reg_dir} ${rq_dir} ${rt_dir} ${rtr_dir}

    # riboWaltz P-site offset files: symlink, falling back to copy
    for f in ${rw_offsets instanceof List ? rw_offsets.collect{it.name}.join(' ') : rw_offsets ? rw_offsets.name : ''}; do
        [ -z "\$f" ] && continue
        if [ -L "\$f" ]; then cp -rL "\$f" ${rw_off_dir}/ 2>/dev/null || cp -r "\$f" ${rw_off_dir}/
        elif [ -f "\$f" ]; then cp "\$f" ${rw_off_dir}/
        fi
    done

    for f in ${rw_regions instanceof List ? rw_regions.collect{it.name}.join(' ') : rw_regions ? rw_regions.name : ''}; do
        [ -z "\$f" ] && continue
        if [ -L "\$f" ]; then cp -rL "\$f" ${rw_reg_dir}/ 2>/dev/null || cp -r "\$f" ${rw_reg_dir}/
        elif [ -f "\$f" ]; then cp "\$f" ${rw_reg_dir}/
        fi
    done

    for f in ${rq_psites instanceof List ? rq_psites.collect{it.name}.join(' ') : rq_psites ? rq_psites.name : ''}; do
        [ -z "\$f" ] && continue
        if [ -L "\$f" ]; then cp -rL "\$f" ${rq_dir}/ 2>/dev/null || cp -r "\$f" ${rq_dir}/
        elif [ -f "\$f" ]; then cp "\$f" ${rq_dir}/
        fi
    done

    for f in ${rt_qual instanceof List ? rt_qual.collect{it.name}.join(' ') : rt_qual ? rt_qual.name : ''}; do
        [ -z "\$f" ] && continue
        if [ -L "\$f" ]; then cp -rL "\$f" ${rt_dir}/ 2>/dev/null || cp -r "\$f" ${rt_dir}/
        elif [ -f "\$f" ]; then cp "\$f" ${rt_dir}/
        fi
    done

    for f in ${rtr_summaries instanceof List ? rtr_summaries.collect{it.name}.join(' ') : rtr_summaries ? rtr_summaries.name : ''}; do
        [ -z "\$f" ] && continue
        if [ -L "\$f" ]; then cp -rL "\$f" ${rtr_dir}/ 2>/dev/null || cp -r "\$f" ${rtr_dir}/
        elif [ -f "\$f" ]; then cp "\$f" ${rtr_dir}/
        fi
    done

    # Run the joint QC report
    Rscript ${report_script} \\
        ${rw_off_arg} \\
        ${rw_reg_arg} \\
        ${rq_arg} \\
        ${rt_arg} \\
        ${rtr_arg} \\
        --output_prefix joint_riboseq_qc

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(Rscript -e 'cat(paste(R.version[["major"]], R.version[["minor"]], sep="."))')
    END_VERSIONS
    """

    stub:
    """
    touch joint_riboseq_qc.tsv joint_riboseq_qc_mqc.yaml joint_riboseq_qc_mqc.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: "4.1.0"
    END_VERSIONS
    """
}
