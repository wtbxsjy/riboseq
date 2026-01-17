FROM continuumio/miniconda3:23.10.0-1

LABEL \
    authors="nf-core/riboseq" \
    description="Container for gencode-riboseqORFs integration" \
    maintainer="nf-core/riboseq" \
    version="1.1.0"

# Install system dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    wget \
    git \
    procps \
    && rm -rf /var/lib/apt/lists/*

# Install Conda packages
RUN conda install -c bioconda -c conda-forge -y \
    bedtools=2.30.0 \
    gffread=0.12.7 \
    && conda clean -a -y

# Install Python dependencies
RUN pip install --no-cache-dir biopython==1.81

# Clone gencode-riboseqORFs project
WORKDIR /opt
RUN git clone https://github.com/jorruior/gencode-riboseqORFs.git && \
    cd gencode-riboseqORFs && \
    chmod +x scripts/*.sh scripts/*.py ORF_mapper_to_GENCODE_v1.1.py

# Set environment variables
ENV PATH="/opt/gencode-riboseqORFs:/opt/gencode-riboseqORFs/scripts:${PATH}"

# Verify installations
RUN python3 --version && \
    bedtools --version && \
    gffread --version && \
    python3 -c "import Bio; print('Biopython:', Bio.__version__)" && \
    ls -l /opt/gencode-riboseqORFs/ORF_mapper_to_GENCODE_v1.1.py

# Set working directory
WORKDIR /data

# Default command
CMD ["/bin/bash"]
