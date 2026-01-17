process DOWNLOAD_ENSEMBL_FILES {
    tag "Ensembl_${ensembl_release}_${genome_assembly}"
    label 'process_low'

    conda "conda-forge::wget=1.21.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/wget:1.21.3' :
        'biocontainers/wget:1.21.3' }"

    input:
    val ensembl_release
    val genome_assembly
    val species  // 'homo_sapiens' or 'mus_musculus'

    output:
    path "Ens${ensembl_release}/", emit: ensembl_dir
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def species_lower = species.toLowerCase()
    def species_caps = species.toLowerCase().split('_').collect { it.capitalize() }.join('_')

    """
    mkdir -p Ens${ensembl_release}
    cd Ens${ensembl_release}

    # Download GTF
    echo "Downloading GTF..."
    wget -q http://ftp.ensembl.org/pub/release-${ensembl_release}/gtf/${species_lower}/${species_caps}.${genome_assembly}.${ensembl_release}.gtf.gz

    # Download proteome
    echo "Downloading proteome..."
    wget -q http://ftp.ensembl.org/pub/release-${ensembl_release}/fasta/${species_lower}/pep/${species_caps}.${genome_assembly}.pep.all.fa.gz

    # Download cDNA
    echo "Downloading cDNA..."
    wget -q http://ftp.ensembl.org/pub/release-${ensembl_release}/fasta/${species_lower}/cdna/${species_caps}.${genome_assembly}.cdna.all.fa.gz

    # Download ncRNA
    echo "Downloading ncRNA..."
    wget -q http://ftp.ensembl.org/pub/release-${ensembl_release}/fasta/${species_lower}/ncrna/${species_caps}.${genome_assembly}.ncrna.fa.gz

    # Decompress all files
    echo "Decompressing files..."
    gunzip -f *.gz

    # Sort GTF file
    echo "Sorting GTF..."
    sort -k1,1 -k4,4n -k5,5n ${species_caps}.${genome_assembly}.${ensembl_release}.gtf > ${species_caps}.${genome_assembly}.sorted.gtf

    # Merge transcriptome sequences (cDNA + ncRNA)
    echo "Merging transcriptome sequences..."
    cat ${species_caps}.${genome_assembly}.cdna.all.fa ${species_caps}.${genome_assembly}.ncrna.fa | \\
        cut -d"." -f1,1 > ${species_caps}.${genome_assembly}.trans.fa

    # Clean up protein IDs (remove version numbers)
    echo "Processing proteome..."
    cut -d"." -f1,1 ${species_caps}.${genome_assembly}.pep.all.fa > tmpfile
    mv tmpfile ${species_caps}.${genome_assembly}.pep.all.fa

    # Download Transcript Support information from Biomart
    echo "Downloading transcript support information..."
    wget -q -O ENST_support.txt 'http://ensembl.org/biomart/martservice?query=<Query virtualSchemaName="default" formatter="TSV" header="1" uniqueRows="0" count="" datasetConfigVersion="0.6"><Dataset name="${species_lower}_gene_ensembl" interface="default"><Attribute name="ensembl_transcript_id"/><Attribute name="transcript_tsl"/><Attribute name="transcript_appris"/></Dataset></Query>' || echo "Warning: Could not download transcript support information"

    cd ..

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ensembl_release: ${ensembl_release}
        genome_assembly: ${genome_assembly}
        species: ${species}
        wget: \$(wget --version 2>&1 | head -n 1 | sed 's/GNU Wget //')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p Ens${ensembl_release}
    touch Ens${ensembl_release}/stub.gtf
    touch Ens${ensembl_release}/stub.fa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ensembl_release: ${ensembl_release}
        genome_assembly: ${genome_assembly}
        species: ${species}
        wget: "1.21.3"
    END_VERSIONS
    """
}

process CALCULATE_PSITE_BED {
    tag "Ensembl_${ensembl_release}"
    label 'process_low'

    conda "conda-forge::python=3.9 conda-forge::biopython=1.81"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.81' :
        'biocontainers/biopython:1.81' }"

    input:
    path ensembl_dir
    val ensembl_release
    val genome_assembly
    val species

    output:
    path "${ensembl_dir}/", emit: ensembl_dir_complete
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def species_lower = species.toLowerCase()
    def species_caps = species.toLowerCase().split('_').collect { it.capitalize() }.join('_')

    """
    # Copy the calculate_frame_bed.py script to working directory
    cp ${moduleDir}/assets/calculate_frame_bed.py .

    # Run the script on the sorted GTF
    cd ${ensembl_dir}
    python3 ../calculate_frame_bed.py ${species_caps}.${genome_assembly}.sorted.gtf
    cd ..

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
        biopython: \$(python3 -c "import Bio; print(Bio.__version__)" 2>/dev/null || echo "1.81")
    END_VERSIONS
    """

    stub:
    """
    mkdir -p ${ensembl_dir}
    touch ${ensembl_dir}/psites.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: "3.9.0"
        biopython: "1.81"
    END_VERSIONS
    """
}
