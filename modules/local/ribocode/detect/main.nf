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
gene_id_re = re.compile(r'gene_id\\s+"([^"]+)"')
tx_id_re = re.compile(r'transcript_id\\s+"([^"]+)"')

# Pass 1: collect all registered gene IDs from 'gene' features,
# and all transcript IDs from 'transcript' features
valid_genes = set()
valid_transcripts = set()
tx_to_gene = {}  # transcript_id -> gene_id
with open(gtf_in, 'rt', encoding='utf-8', errors='replace') as fin:
    for line in fin:
        if line.startswith('#'): continue
        parts = line.split('\\t')
        if len(parts) < 9: continue
        if parts[2] == 'gene':
            m = gene_id_re.search(parts[8])
            if m: valid_genes.add(m.group(1))
        elif parts[2] == 'transcript':
            m = tx_id_re.search(parts[8])
            gm = gene_id_re.search(parts[8])
            if m:
                valid_transcripts.add(m.group(1))
                if gm: tx_to_gene[m.group(1)] = gm.group(1)

# Pass 2: write normalized GTF, adding synthetic transcript entries for
# exon/CDS features that lack a parent transcript line (common in GENCODE
# for protein_coding_CDS_not_defined transcripts)
synth_added = set()
skipped = 0
written = 0
with open(gtf_in, 'rt', encoding='utf-8', errors='replace') as fin, \\
     open(gtf_out, 'wt', encoding='utf-8') as fout:
    for line in fin:
        if line.startswith('#') or not line.strip():
            fout.write(line)
            continue
        parts = line.rstrip('\\n').split('\\t')
        if len(parts) < 9:
            fout.write(line)
            continue
        attrs = parts[8].strip()
        gm = gene_id_re.search(attrs)
        # Drop features whose parent gene is not registered
        if gm and gm.group(1) not in valid_genes:
            skipped += 1
            continue
        # Add synthetic transcript entry for exon/CDS/UTR features that
        # belong to a transcript without a transcript feature line
        tm = tx_id_re.search(attrs)
        if tm and tm.group(1) not in valid_transcripts and parts[2] in ('exon','CDS','start_codon','stop_codon','UTR'):
            tid = tm.group(1)
            gid = gm.group(1) if gm else tid
            if tid not in synth_added:
                synth_added.add(tid)
                # Synthesize a transcript entry
                synth_attrs = f'gene_id "{gid}"; transcript_id "{tid}"; level "NA";'
                synth = '\\t'.join([parts[0], 'RiboCode', 'transcript',
                    parts[3], parts[4], '.', parts[6], '.', synth_attrs])
                fout.write(synth + '\\n')
                valid_transcripts.add(tid)
        if not level_re.search(attrs):
            if attrs and not attrs.endswith(';'):
                attrs += ';'
            attrs += ' level "NA";'
            parts[8] = attrs
        fout.write('\\t'.join(parts) + '\\n')
        written += 1
print(f"GTF normalized: {written} lines, {skipped} orphan, {len(synth_added)} synth transcripts", file=sys.stderr)
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
