process SPLIT_BAM_BY_CONTIG {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.21--h50ea8bc_0' :
        'biocontainers/samtools:1.21--h50ea8bc_0' }"

    input:
    tuple val(meta), path(bam), path(bai)
    path fai
    val pathogen_regex

    output:
    tuple val(meta), path("*${meta.id}.host.bam"),   emit: host_bam
    tuple val(meta), path("*${meta.id}.pathogen.bam"), emit: pathogen_bam
    path "*_split_stats.tsv",                        emit: stats
    path "versions.yml",                             emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = meta.id
    """
    set -euo pipefail

    prefix='${prefix}'
    re='${pathogen_regex}'
    fai='${fai}'
    bam='${bam}'

    # Identify pathogen contigs from FAI
    awk -v re="\$re" '\$1 ~ re {print \$1}' "\$fai" > "\${prefix}.pathogen_contigs.txt" || true

    pathogen_count=\$(wc -l < "\${prefix}.pathogen_contigs.txt" 2>/dev/null || echo 0)
    echo "Found \${pathogen_count} pathogen contig(s) matching pattern: \${re}"

    if [ "\${pathogen_count}" -eq 0 ]; then
        echo "WARNING: No pathogen contigs found matching pattern '\${re}'"
        echo "Creating empty pathogen BAM and passing through all reads as host"
        samtools view -b "\$bam" -o "\${prefix}.host.bam"
        samtools view -H "\$bam" | samtools view -b -o "\${prefix}.pathogen.bam" -
    else
        # Extract pathogen-mapped reads
        samtools view -b "\$bam" \$(cat "\${prefix}.pathogen_contigs.txt") -o "\${prefix}.pathogen.bam"

        # Extract host-mapped reads: exclude pathogen contigs AND their mates
        # Strategy: write reads on host contigs only, including mates that may be on pathogen contigs
        # First approach: use samtools with inverted contig filter
        samtools view -h "\$bam" \\
            | awk -v contigs_file="\${prefix}.pathogen_contigs.txt" '
              BEGIN {
                  while ((getline < contigs_file) > 0) { pathogen[\$1] = 1 }
                  close(contigs_file)
              }
              /^@/ { print; next }
              {
                  rna = \$3
                  if (rna == "*") { print; next }
                  if (!(rna in pathogen)) { print }
              }
            ' \\
            | samtools view -b -o "\${prefix}.host.bam" -
    fi

    # Compute split statistics
    host_reads=\$(samtools view -c "\${prefix}.host.bam" 2>/dev/null || echo 0)
    pathogen_reads=\$(samtools view -c "\${prefix}.pathogen.bam" 2>/dev/null || echo 0)
    total_reads=\$((host_reads + pathogen_reads))

    printf "species\\ttotal_reads\\tpct\\n" > "\${prefix}_split_stats.tsv"
    if [ "\$total_reads" -gt 0 ]; then
        host_pct=\$(awk -v h="\$host_reads" -v t="\$total_reads" 'BEGIN { printf "%.2f", 100.0 * h / t }')
        pathogen_pct=\$(awk -v p="\$pathogen_reads" -v t="\$total_reads" 'BEGIN { printf "%.2f", 100.0 * p / t }')
    else
        host_pct="0.00"
        pathogen_pct="0.00"
    fi
    printf "host\\t%d\\t%s\\n" "\$host_reads" "\$host_pct" >> "\${prefix}_split_stats.tsv"
    printf "pathogen\\t%d\\t%s\\n" "\$pathogen_reads" "\$pathogen_pct" >> "\${prefix}_split_stats.tsv"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | head -n 1 | sed 's/samtools //')
        awk: \$(awk --version 2>/dev/null | head -n 1 || echo "unknown")
    END_VERSIONS
    """
}
