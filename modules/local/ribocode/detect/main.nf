process RIBOCODE_DETECT {
    tag "$meta.id"
    label 'process_medium'

    // Allow process to complete even if RiboCode fails due to low periodicity
    errorStrategy 'ignore'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ribocode:1.2.11--pyh145b6a8_1' :
        'quay.io/biocontainers/ribocode:1.2.11--pyh145b6a8_1' }"

    input:
    tuple val(meta), path(bam), path(bai)
    path gtf
    path fasta

    output:
    tuple val(meta), path("${meta.id}.txt")          , emit: txt
    tuple val(meta), path("${meta.id}_collapsed.txt"), emit: collapsed
    tuple val(meta), path("${meta.id}.gtf.gz")       , emit: gtf
    tuple val(meta), path("${meta.id}.bed.gz")       , emit: bed
    tuple val(meta), path("${meta.id}*")             , emit: results
    path "versions.yml"                              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def strandedness = 'yes'
    if (meta.strandedness == 'reverse') {
        strandedness = 'reverse'
    } else if (meta.strandedness == 'unstranded') {
        strandedness = 'no'
    }
    // RiboCode strandedness options: yes, reverse, no.
    // nf-core/riboseq strandedness: forward, reverse, unstranded.
    // forward -> yes
    // reverse -> reverse
    // unstranded -> no

    """
    # Ensure temporary/config directories are writable inside the container
    export TMPDIR="\$PWD/tmp"
    mkdir -p "\$TMPDIR"
    export MPLCONFIGDIR="\$PWD/mplconfig"
    mkdir -p "\$MPLCONFIGDIR"

    # RiboCode can crash when transcript "level" is missing for some features.
    # Normalize the GTF by adding a default level to any record lacking it.
    GTF_IN="${gtf}"
    if [[ "\$GTF_IN" == *.gz ]]; then
        gunzip -c "\$GTF_IN" > input.gtf
        GTF_IN="input.gtf"
    fi

    python - "\$GTF_IN" ribocode.gtf <<'PY'
import re
import sys

gtf_in = sys.argv[1]
gtf_out = sys.argv[2]

level_re = re.compile(r'(?:^|;\\s*)level\\s+"[^"]*"\\s*;')

with open(gtf_in, 'rt', encoding='utf-8', errors='replace') as fin, open(gtf_out, 'wt', encoding='utf-8') as fout:
    for line in fin:
        if line.startswith('#') or not line.strip():
            fout.write(line)
            continue
        parts = line.rstrip('\\n').split('\\t')
        if len(parts) < 9:
            fout.write(line)
            continue
        attrs = parts[8].strip()
        if not level_re.search(attrs):
            if attrs and not attrs.endswith(';'):
                attrs += ';'
            attrs += ' level "NA";'
            parts[8] = attrs
        fout.write('\\t'.join(parts) + '\\n')
PY

    # Run RiboCode, capture exit status
    RiboCode_onestep \\
        -g ribocode.gtf \\
        -f $fasta \\
        -r $bam \\
        --stranded $strandedness \\
        -o ${meta.id} \\
        -outgtf \\
        -outbed \\
        $args || {
        if [ "\${RIBOCODE_FAIL_ON_EMPTY:-false}" = "true" ]; then
            echo "FATAL: RiboCode produced no ORF predictions for ${meta.id} and --ribocode_fail_on_empty is true"
            exit 1
        fi
        echo "RiboCode failed for ${meta.id} - likely due to insufficient periodicity in data"
        echo "This is common for low-depth or test datasets"
        # Create marker and placeholder files so downstream channels stay stable.
        echo "FAILED: No periodicity detected" > ${meta.id}.ribocode_failed.txt
        echo -e "ORF_ID\\tORF_type\\ttranscript_id\\tgene_id\\tchrom\\tstrand\\tORF_length\\tORF_gstart\\tORF_gstop\\tpval_combined\\tadjusted_pval" > ${meta.id}.txt
        cp ${meta.id}.txt ${meta.id}_collapsed.txt
        touch ${meta.id}.gtf ${meta.id}_collapsed.gtf ${meta.id}.bed ${meta.id}_collapsed.bed
        gzip -f ${meta.id}.gtf ${meta.id}_collapsed.gtf ${meta.id}.bed ${meta.id}_collapsed.bed
    }

    # Older RiboCode failures can leave only partial output; normalize to an empty
    # tabular/interval set rather than letting downstream workflow wiring fail.
    if [ ! -f ${meta.id}.txt ]; then
        echo -e "ORF_ID\\tORF_type\\ttranscript_id\\tgene_id\\tchrom\\tstrand\\tORF_length\\tORF_gstart\\tORF_gstop\\tpval_combined\\tadjusted_pval" > ${meta.id}.txt
    fi
    if [ ! -f ${meta.id}_collapsed.txt ]; then
        cp ${meta.id}.txt ${meta.id}_collapsed.txt
    fi
    touch ${meta.id}.gtf ${meta.id}_collapsed.gtf ${meta.id}.bed ${meta.id}_collapsed.bed

    # Compress text outputs to save disk space
    for f in ${meta.id}.gtf ${meta.id}_collapsed.gtf ${meta.id}.bed ${meta.id}_collapsed.bed; do
        [ -f "\$f" ] && gzip -f "\$f" || true
    done

    ribocode_ver="\$(RiboCode_onestep --version 2>&1 | sed 's/^.* //')"
    printf '"%s":\\n    ribocode: "%s"\\n' "${task.process}" "\$ribocode_ver" > versions.yml
    """.stripIndent()
}
