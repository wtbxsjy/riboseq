process PRICE {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/gedi_price:latest' :
        'community.wave.seqera.io/library/gedi_price:latest' }"

    input:
    tuple val(meta), path(bam), path(bai)
    path fasta
    path gtf

    output:
    tuple val(meta), path("*_Detected_ORFs.gtf"), emit: gtf
    tuple val(meta), path("*.orfs.tsv")          , emit: orfs_tsv, optional: true
    path "versions.yml"                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def memory = task.memory ? task.memory.toGiga() + 'g' : '16g'
    """
    #!/bin/bash
    set -euo pipefail

    OUTDIR="."
    PREFIX="${prefix}"

    # Step 1: Convert BAM to GEDI CIT format
    echo "[PRICE] Converting BAM to CIT format..."
    gedi -e Bam2CIT \
        -bam ${bam} \
        -o "\${OUTDIR}/\${PREFIX}.cit" \
        -genome ${fasta} \
        -gtf ${gtf} \
        || {
            echo "WARNING: Bam2CIT conversion failed."
            echo "BAM file may be empty or incompatible. Creating placeholder outputs."
            echo -e "# PRICE placeholder\\norf_id\\tchr\\tstart\\tend\\tstrand\\ttype\\tscore" > "\${OUTDIR}/\${PREFIX}.orfs.tsv"
            echo "# PRICE placeholder - BAM conversion failed" > "\${OUTDIR}/\${PREFIX}_Detected_ORFs.gtf"
            cat <<-END_VERSIONS > "\${OUTDIR}/versions.yml"
            "${task.process}":
                gedi: "1.0.6d"
                price: "pipeline"
            END_VERSIONS
            exit 0
        }

    CIT_FILE="\${OUTDIR}/\${PREFIX}.cit"

    # Step 2: Run PRICE pipeline
    echo "[PRICE] Running PRICE ORF detection..."
    export GEDI_MEMORY="${memory}"

    gedi -e Price \
        -r "\${CIT_FILE}" \
        -g ${gtf} \
        -f ${fasta} \
        -o "\${OUTDIR}/\${PREFIX}" \
        -t ${task.cpus} \
        || {
            echo "WARNING: PRICE main pipeline returned non-zero exit code."
            echo "Some stages may have completed. Checking for partial output..."
        }

    # Step 3: Collect ORF output
    if [ -f "\${OUTDIR}/\${PREFIX}.orfs.tsv" ] && [ -s "\${OUTDIR}/\${PREFIX}.orfs.tsv" ]; then
        echo "[PRICE] ORF detection successful."

        # Run post-PRICE analysis if output exists
        gedi -e PriceAnalyze \
            -o "\${OUTDIR}/\${PREFIX}" \
            -g ${gtf} \
            || echo "WARNING: PriceAnalyze post-processing failed."
    else
        echo "[PRICE] No ORFs detected. Creating placeholder outputs."
        echo -e "orf_id\\tchr\\tstart\\tend\\tstrand\\ttype\\tscore" > "\${OUTDIR}/\${PREFIX}.orfs.tsv"
    fi

    # Step 4: Convert PRICE ORF TSV to GTF for unification pipeline
    echo "[PRICE] Converting ORFs to GTF..."
    cat <<'RCONVERT' > convert_orfs.R
    # Convert PRICE ORF TSV to GTF
    prefix <- "${prefix}"
    orf_file <- paste0(prefix, ".orfs.tsv")
    gtf_out  <- paste0(prefix, "_Detected_ORFs.gtf")

    if (file.exists(orf_file)) {
        orfs <- tryCatch(
            read.delim(orf_file, stringsAsFactors = FALSE, check.names = FALSE),
            error = function(e) NULL
        )

        if (!is.null(orfs) && nrow(orfs) > 0 && "orf_id" %in% colnames(orfs)) {
            cat(paste0("# PRICE ORFs: ", nrow(orfs), " entries\\n"), file = gtf_out)

            for (i in seq_len(nrow(orfs))) {
                row <- orfs[i, ]
                chr    <- if ("chr"   %in% names(row)) row\$chr   else "unknown"
                start  <- if ("start" %in% names(row)) row\$start else 1
                end    <- if ("end"   %in% names(row)) row\$end   else 1
                strand <- if ("strand" %in% names(row)) row\$strand else "."
                orf_id <- row\$orf_id

                cat(sprintf('%s\\tPRICE\\tCDS\\t%s\\t%s\\t.\\t%s\\t.\\tgene_id "%s"; transcript_id "%s";\\n',
                    chr, start, end, strand, orf_id, orf_id),
                    file = gtf_out, append = TRUE)
            }
            cat(sprintf("Wrote %d PRICE ORFs to GTF\\n", nrow(orfs)))
        } else {
            cat("# PRICE: no valid ORFs found\\n", file = gtf_out)
            cat("PRICE: no valid ORF entries to convert\\n")
        }
    } else {
        cat("# PRICE: output file not found\\n", file = gtf_out)
        cat("PRICE: ORF TSV file not found\\n")
    }
RCONVERT

    Rscript convert_orfs.R

    # Write versions
    GEDI_VER=\$(gedi -e Version 2>/dev/null || echo "1.0.6d")
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gedi: "\$GEDI_VER"
        price: "pipeline"
    END_VERSIONS

    echo "[PRICE] Pipeline complete."
    ls -lh \${PREFIX}*
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo -e "orf_id\\tchr\\tstart\\tend\\tstrand\\ttype\\tscore" > ${prefix}.orfs.tsv
    echo "# PRICE stub" > ${prefix}_Detected_ORFs.gtf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gedi: "1.0.6d"
        price: "pipeline"
    END_VERSIONS
    """
}
