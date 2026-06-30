process SORF_BAM_FILTER {
    tag "${meta.id}"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.21--h50ea8bc_0' :
        'biocontainers/samtools:1.21--h50ea8bc_0' }"

    input:
    tuple val(meta), path(bam), path(bai)
    path fai
    val unique_mode
    val unique_mapq
    val read_len_min
    val read_len_max
    val exclude_contigs_regex
    val gtf                     // GTF path (string) — used to extract MT/Pt transcript IDs

    output:
    tuple val(meta), path("*.sorf.filtered.bam"), emit: bam
    path "*.sorf.filter_stats.tsv", emit: stats
    path "*.sorf.excluded_contigs.txt", optional: true, emit: excluded_contigs
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = meta.id
    def re = exclude_contigs_regex ?: ''
    def rlmin = read_len_min ?: 0
    def rlmax = read_len_max ?: 0
    def umode = unique_mode ?: 'auto'
    def mapq = unique_mapq ?: 0
    def has_gtf = (gtf && gtf != '' && !(gtf instanceof List && gtf.isEmpty())) ? ( gtf instanceof List ? gtf[0] : gtf ) : false

    """
    set -euo pipefail

    prefix='${prefix}'
    mode='${umode}'
    mapq=${mapq}
    rlmin=${rlmin}
    rlmax=${rlmax}
    re='${re}'
    gtf='${has_gtf}'

    # Build exclusion list from FAI contig names matching the regex
    if [[ -n "\$re" ]]; then
      awk -v re="\$re" 'BEGIN{FS="\\t"; OFS="\\t"} \$1 ~ re {print \$1}' ${fai} | sort -u > ${prefix}.sorf.excluded_contigs.txt || true
    fi

    # If GTF is provided, also extract MT/Pt transcript IDs (for transcriptome BAM filtering)
    # Transcriptome BAMs use transcript IDs as RNAME, not chromosome names
    if [[ -n "\$gtf" && -f "\$gtf" && -n "\$re" ]]; then
      # Extract chromosome names matching the regex from GTF col 1,
      # then get the transcript_id values from those same lines
      awk -v re="\$re" -F'\\t' '
        BEGIN { OFS="\\t" }
        \$1 ~ re {
          # Extract transcript_id from the attributes column (col 9)
          s=\$9
          while (match(s, /transcript_id "[^"]+"/)) {
            tx = substr(s, RSTART+16, RLENGTH-17)
            print tx
            s = substr(s, RSTART+RLENGTH)
          }
        }
      ' "\$gtf" | sort -u > ${prefix}.sorf.excluded_transcript_ids.txt || true
      # Combine FAI contigs and transcript IDs; pass to awk as a combined exclusion set
      cat ${prefix}.sorf.excluded_contigs.txt ${prefix}.sorf.excluded_transcript_ids.txt 2>/dev/null \\\
        | sort -u > ${prefix}.sorf.excluded_combined.txt || true
    fi

    # Count primary mapped reads (unfiltered baseline)
    # Exclude: unmapped(0x4), secondary(0x100), duplicate(0x400), supplementary(0x800)
    set +e
    total_primary_mapped=`samtools view -c -F 0xD04 ${bam} 2> samtools_total.err`
    total_status=\$?
    set -e
    if [[ \$total_status -ne 0 ]]; then
      echo "WARNING: samtools failed counting reads; treating as 0 (see samtools_total.err)" >&2
      total_primary_mapped=0
    fi

    # Build combined exclusion set from FAI contigs + transcript IDs (if GTF provided)
    # Used for transcriptome BAM filtering where MT/Pt genes have transcript ID names
    if [[ -f ${prefix}.sorf.excluded_combined.txt ]]; then
      awk '{print \$1}' ${prefix}.sorf.excluded_combined.txt > ${prefix}.sorf.excluded_list.txt
      use_exclude_file="${prefix}.sorf.excluded_list.txt"
    else
      use_exclude_file=""
    fi

    # Filter: keep header lines; drop contigs matching regex; enforce read length; enforce unique mapping.
    set +e
    samtools view -h -F 0xD04 ${bam} \
        | awk -v mode="\$mode" -v mapq="\$mapq" -v rlmin="\$rlmin" -v rlmax="\$rlmax" -v re="\$re" -v excl_file="\$use_exclude_file" '
          BEGIN {
            OFS="\\t"
            # Load exclusion file if provided
            if (excl_file != "") {
              while ((getline line < excl_file) > 0) {
                excluded[line] = 1
              }
              close(excl_file)
            }
          }
          /^@/ {print; next}
          {
            rname=\$3
            if (excl_file != "" && rname in excluded) next
            if (re != "" && rname ~ re) next

            # Read length from SEQ column
            seqlen=length(\$10)
            if (rlmin > 0 && seqlen < rlmin) next
            if (rlmax > 0 && seqlen > rlmax) next

            if (mode == "mapq") {
              if (\$5 < mapq) next
            } else if (mode == "nh" || mode == "auto") {
              nh=""
              for (i=12; i<=NF; i++) {
                if (\$i ~ /^NH:i:/) { split(\$i,a,":"); nh=a[3]; break }
              }
              if (mode == "nh") {
                if (nh == "" || nh != 1) next
              } else {
                # auto: prefer NH when present; fall back to MAPQ when NH is absent
                if (nh != "") {
                  if (nh != 1) next
                } else {
                  if (\$5 < mapq) next
                }
              }
            }

            print
          }
        ' \
      | samtools view -b | samtools sort -o ${prefix}.sorf.filtered.bam -
    # CRITICAL: Capture PIPESTATUS in one line before it resets (required with set -u)
    pipe_statuses=("\${PIPESTATUS[@]}")
    status_view=\${pipe_statuses[0]:-1}
    status_awk=\${pipe_statuses[1]:-1}
    status_bam=\${pipe_statuses[2]:-1}
    status_sort=\${pipe_statuses[3]:-1}
    set -e

    if [[ \$status_view -ne 0 || \$status_awk -ne 0 || \$status_bam -ne 0 || \$status_sort -ne 0 ]]; then
      echo "WARNING: samtools/awk pipeline failed; creating empty filtered BAM (see stderr)" >&2
      set +e
      samtools view -H ${bam} | samtools view -b | samtools sort -o ${prefix}.sorf.filtered.bam -
      header_status=\$?
      set -e
      if [[ \$header_status -ne 0 ]]; then
        echo "WARNING: failed to write header-only BAM; creating empty placeholder" >&2
        rm -f ${prefix}.sorf.filtered.bam
        touch ${prefix}.sorf.filtered.bam
      fi
      kept_primary_mapped=0
    else
      kept_primary_mapped=`samtools view -c -F 0xD04 ${prefix}.sorf.filtered.bam`
    fi

    printf "sample\\ttotal_primary_mapped\\tkept_primary_mapped\\tpct_kept\\n" > ${prefix}.sorf.filter_stats.tsv
    awk -v s="${prefix}" -v t="\$total_primary_mapped" -v k="\$kept_primary_mapped" 'BEGIN{pct=(t>0)?(100.0*k/t):0; printf "%s\\t%d\\t%d\\t%.2f\\n", s, t, k, pct}' >> ${prefix}.sorf.filter_stats.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
      samtools: \$(samtools --version | head -n 1 | sed 's/samtools //')
      awk: \$(awk --version 2>/dev/null | head -n 1 || echo "unknown")
    END_VERSIONS
    """

    stub:
    def prefix_stub = meta.id
    """
    touch ${prefix_stub}.sorf.filtered.bam
    printf "sample\ttotal_primary_mapped\tkept_primary_mapped\tpct_kept\n" > ${prefix_stub}.sorf.filter_stats.tsv
    printf "${prefix_stub}\t0\t0\t0.00\n" >> ${prefix_stub}.sorf.filter_stats.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: "1.21"
        awk: "gawk"
    END_VERSIONS
    """
}
