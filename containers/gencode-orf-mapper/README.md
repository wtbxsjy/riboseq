# gencode-orf-mapper Container

This container provides the [gencode-riboseqORFs](https://github.com/jorruior/gencode-riboseqORFs) tool for unified ORF annotation in the nf-core/riboseq pipeline.

## Contents

- **Python 3** with Biopython 1.81
- **bedtools** 2.30.0
- **gffread** 0.12.7
- **gencode-riboseqORFs** v1.1.0 scripts

## Building

### Docker

```bash
docker build -t nfcore/gencode-orf-mapper:1.1.0 .
```

### Singularity

```bash
singularity build --fakeroot gencode-orf-mapper_1.1.0.sif Singularity.def
```

### Automated Build

```bash
bash build.sh
```

## Testing

### Docker

```bash
# Check versions
docker run --rm nfcore/gencode-orf-mapper:1.1.0 python3 --version
docker run --rm nfcore/gencode-orf-mapper:1.1.0 bedtools --version

# Test main script
docker run --rm nfcore/gencode-orf-mapper:1.1.0 \
    python3 /opt/gencode-riboseqORFs/ORF_mapper_to_GENCODE_v1.1.py --help

# Interactive shell
docker run -it --rm nfcore/gencode-orf-mapper:1.1.0 bash
```

### Singularity

```bash
# Check versions
singularity exec gencode-orf-mapper_1.1.0.sif python3 --version
singularity exec gencode-orf-mapper_1.1.0.sif bedtools --version

# Test main script
singularity exec gencode-orf-mapper_1.1.0.sif \
    python3 /opt/gencode-riboseqORFs/ORF_mapper_to_GENCODE_v1.1.0.py --help

# Interactive shell
singularity shell gencode-orf-mapper_1.1.0.sif
```

## Usage in Nextflow

```groovy
process GENCODE_ORF_MAPPER {
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'file://gencode-orf-mapper_1.1.0.sif' :
        'nfcore/gencode-orf-mapper:1.1.0' }"

    script:
    """
    python3 /opt/gencode-riboseqORFs/ORF_mapper_to_GENCODE_v1.1.py \\
        -d ensembl_dir \\
        -f orfs.fa \\
        -b orfs.bed \\
        -o output_prefix
    """
}
```

## Files Location

- Main script: `/opt/gencode-riboseqORFs/ORF_mapper_to_GENCODE_v1.1.py`
- Helper scripts: `/opt/gencode-riboseqORFs/scripts/`
- Functions: `/opt/gencode-riboseqORFs/functions.py`

## Publishing

### Docker Hub

```bash
# Login
docker login

# Tag
docker tag nfcore/gencode-orf-mapper:1.1.0 nfcore/gencode-orf-mapper:latest

# Push
docker push nfcore/gencode-orf-mapper:1.1.0
docker push nfcore/gencode-orf-mapper:latest
```

### GitHub Container Registry

```bash
# Login
echo $GITHUB_TOKEN | docker login ghcr.io -u USERNAME --password-stdin

# Tag
docker tag nfcore/gencode-orf-mapper:1.1.0 ghcr.io/nf-core/gencode-orf-mapper:1.1.0

# Push
docker push ghcr.io/nf-core/gencode-orf-mapper:1.1.0
```

## Troubleshooting

### Permission Issues (Singularity)

If you encounter permission errors when building:

```bash
# Use sudo
sudo singularity build gencode-orf-mapper_1.1.0.sif Singularity.def

# Or build in a writable sandbox first
singularity build --sandbox gencode-orf-mapper/ Singularity.def
sudo singularity build gencode-orf-mapper_1.1.0.sif gencode-orf-mapper/
```

### Network Issues

If downloads fail, try:

```bash
# Use a proxy
export http_proxy=http://proxy.example.com:8080
export https_proxy=http://proxy.example.com:8080

# Or build with --network=host
docker build --network=host -t nfcore/gencode-orf-mapper:1.1.0 .
```

## References

- gencode-riboseqORFs: https://github.com/jorruior/gencode-riboseqORFs
- Publication: https://doi.org/10.1038/s41587-022-01369-0
- nf-core/riboseq: https://nf-co.re/riboseq
