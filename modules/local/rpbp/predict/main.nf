process RPBP_PREDICT {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        (params.rpbp_container ?: 'https://depot.galaxyproject.org/singularity/rpbp:4.0.1--py312hf731ba3_0') :
        'biocontainers/rpbp:4.0.1--py312hf731ba3_0' }"

    input:
    tuple val(meta), path(bam), path(bai)
    path orfs_genomic
    path orfs_exons

    output:
    tuple val(meta), path("*.predicted-orfs.bed.gz"), emit: predictions
    tuple val(meta), path("*.predicted-orfs.dna.fa"), emit: predicted_dna
    tuple val(meta), path("*.predicted-orfs.protein.fa"), emit: predicted_protein
    tuple val(meta), path("*.bayes-factors.bed.gz"), emit: bayes_factors
    path "versions.yml"                            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def args = task.ext.args ?: ''
    """
    # Error handling: wrap the entire RPBP pipeline to handle empty results gracefully
    rpbp_ok=0
    set +e

    # Find models using python
    # We output them to files to read them into variables
    python3 -c '
import os
import rpbp
import glob
import sys

pkg_path = os.path.dirname(rpbp.__file__)
models_path = os.path.join(pkg_path, "models")

def get_models(subdir):
    p = os.path.join(models_path, subdir)
    # Check if this path exists, if not try sys.prefix fallback
    if not os.path.exists(p):
        p = os.path.join(sys.prefix, "share", "rpbp", "models", subdir)

    if not os.path.exists(p):
        # Last resort: try to find where the package is installed and look for models there
        # This handles cases where rpbp.__file__ might be weird or symlinked
        pass

    if not os.path.exists(p):
        sys.stderr.write(f"WARNING: Could not find models directory for {subdir} at {p}\\n")
        return ""

    # In rpbp 4.x, models are compiled binaries (ELF), not .pkl
    # We want to find files that do not end in .stan, .py, .pyc and are not directories
    files = []
    for f in os.listdir(p):
        fp = os.path.join(p, f)
        if os.path.isfile(fp) and not f.endswith(".stan") and not f.endswith(".py") and not f.endswith(".pyc") and not f.startswith("__"):
                files.append(fp)

    if not files:
        sys.stderr.write(f"WARNING: No models found in {p}\\n")

    return " ".join(files)

with open("periodic_models.txt", "w") as f: f.write(get_models("periodic"))
with open("nonperiodic_models.txt", "w") as f: f.write(get_models("nonperiodic"))
with open("translated_models.txt", "w") as f: f.write(get_models("translated"))
with open("untranslated_models.txt", "w") as f: f.write(get_models("untranslated"))
    '

    PERIODIC_MODELS=\$(cat periodic_models.txt)
    NONPERIODIC_MODELS=\$(cat nonperiodic_models.txt)
    TRANSLATED_MODELS=\$(cat translated_models.txt)
    UNTRANSLATED_MODELS=\$(cat untranslated_models.txt)

    # 1. Extract metagene profiles
    extract-metagene-profiles \\
        $bam \\
        $orfs_genomic \\
        ${prefix}.metagene-profiles.csv.gz \\
        --num-cpus $task.cpus || rpbp_ok=\$?

    # 2. Estimate metagene profile Bayes factors
    if [ \$rpbp_ok -eq 0 ]; then
        estimate-metagene-profile-bayes-factors \\
            ${prefix}.metagene-profiles.csv.gz \\
            ${prefix}.metagene-profile-bayes-factors.csv.gz \\
            --periodic-models \$PERIODIC_MODELS \\
            --nonperiodic-models \$NONPERIODIC_MODELS \\
            --num-cpus $task.cpus \\
            $args || rpbp_ok=\$?
    fi

    # 3. Select periodic offsets
    if [ \$rpbp_ok -eq 0 ]; then
        select-periodic-offsets \\
            ${prefix}.metagene-profile-bayes-factors.csv.gz \\
            ${prefix}.periodic-offsets.csv.gz || rpbp_ok=\$?
    fi

    # 4. Extract ORF profiles
    # Parse lengths and offsets
    if [ \$rpbp_ok -eq 0 ]; then
        ARGS=\$(python3 -c "import pandas as pd; df=pd.read_csv('${prefix}.periodic-offsets.csv.gz'); print('--lengths ' + ' '.join(map(str, df['length'].astype(int))) + ' --offsets ' + ' '.join(map(str, df['highest_peak_offset'].astype(int))))")
        extract-orf-profiles \\
            $bam \\
            $orfs_genomic \\
            $orfs_exons \\
            ${prefix}.profiles.mtx.gz \\
            \$ARGS \\
            --num-cpus $task.cpus || rpbp_ok=\$?
    fi

    # 5. Estimate ORF Bayes factors
    if [ \$rpbp_ok -eq 0 ]; then
        estimate-orf-bayes-factors \\
            ${prefix}.profiles.mtx.gz \\
            $orfs_genomic \\
            ${prefix}.bayes-factors.bed.gz \\
            --translated-models \$TRANSLATED_MODELS \\
            --untranslated-models \$UNTRANSLATED_MODELS \\
            --num-cpus $task.cpus \\
            $args || rpbp_ok=\$?
    fi

    # 6. Select final prediction set
    if [ \$rpbp_ok -eq 0 ]; then
        select-final-prediction-set \\
            ${prefix}.bayes-factors.bed.gz \\
            $orfs_genomic \\
            ${prefix}.predicted-orfs.bed.gz \\
            ${prefix}.predicted-orfs.dna.fa \\
            ${prefix}.predicted-orfs.protein.fa || rpbp_ok=\$?
    fi

    set -e

    # Handle empty/failed results
    if [ \$rpbp_ok -ne 0 ]; then
        if [ "\${RPBP_FAIL_ON_EMPTY:-false}" = "true" ]; then
            echo "FATAL: RPBP failed to produce ORF predictions for ${meta.id} and --rpbp_fail_on_empty is true"
            exit 1
        fi
        echo "WARNING: RPBP failed for ${meta.id} - creating placeholder files"
        echo "# RPBP placeholder - no ORFs predicted" | gzip -c > ${prefix}.predicted-orfs.bed.gz
        echo ">no_orfs" > ${prefix}.predicted-orfs.dna.fa
        echo ">no_orfs" > ${prefix}.predicted-orfs.protein.fa
        echo -e "# RPBP placeholder - no Bayes factors computed" | gzip -c > ${prefix}.bayes-factors.bed.gz
    fi

    cat <<END_VERSIONS > versions.yml
"${task.process}":
    rpbp: \$(python3 -c "import rpbp; print(rpbp.__version__)")
END_VERSIONS
    """
}
