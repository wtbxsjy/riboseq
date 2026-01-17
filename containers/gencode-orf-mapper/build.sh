#!/bin/bash
set -euo pipefail

echo "=================================================="
echo "Building gencode-orf-mapper containers"
echo "=================================================="

# Check if Docker is available
if command -v docker &> /dev/null; then
    echo ""
    echo "Building Docker container..."
    docker build -t nfcore/gencode-orf-mapper:1.1.0 .

    # Test Docker container
    echo ""
    echo "Testing Docker container..."
    docker run --rm nfcore/gencode-orf-mapper:1.1.0 python3 --version
    docker run --rm nfcore/gencode-orf-mapper:1.1.0 bedtools --version
    docker run --rm nfcore/gencode-orf-mapper:1.1.0 python3 -c "import Bio; print('Biopython OK')"

    echo "✅ Docker container built successfully!"
    echo "   Image: nfcore/gencode-orf-mapper:1.1.0"
    echo ""
    echo "To push to Docker Hub:"
    echo "   docker push nfcore/gencode-orf-mapper:1.1.0"
else
    echo "⚠️  Docker not found, skipping Docker build"
fi

# Check if Singularity/Apptainer is available
if command -v singularity &> /dev/null || command -v apptainer &> /dev/null; then
    SING_CMD=$(command -v singularity || command -v apptainer)
    echo ""
    echo "Building Singularity container with: $SING_CMD"

    # Try with --fakeroot if available, otherwise without
    if $SING_CMD build --help 2>&1 | grep -q "fakeroot"; then
        $SING_CMD build --fakeroot -F gencode-orf-mapper_1.1.0.sif Singularity.def
    else
        $SING_CMD build -F gencode-orf-mapper_1.1.0.sif Singularity.def
    fi

    # Test Singularity container
    echo ""
    echo "Testing Singularity container..."
    $SING_CMD exec gencode-orf-mapper_1.1.0.sif python3 --version
    $SING_CMD exec gencode-orf-mapper_1.1.0.sif bedtools --version
    $SING_CMD exec gencode-orf-mapper_1.1.0.sif python3 -c "import Bio; print('Biopython OK')"

    echo "✅ Singularity container built successfully!"
    echo "   Image: gencode-orf-mapper_1.1.0.sif"
else
    echo "⚠️  Singularity/Apptainer not found, skipping Singularity build"
fi

echo ""
echo "=================================================="
echo "Build completed!"
echo "=================================================="
echo ""
echo "Next steps:"
echo "1. Test the container:"
echo "   docker run --rm nfcore/gencode-orf-mapper:1.1.0 python3 /opt/gencode-riboseqORFs/ORF_mapper_to_GENCODE_v1.1.py --help"
echo ""
echo "2. Push to registry (optional):"
echo "   docker push nfcore/gencode-orf-mapper:1.1.0"
echo ""
echo "3. Use in Nextflow:"
echo "   container 'nfcore/gencode-orf-mapper:1.1.0'"
