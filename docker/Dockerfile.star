FROM nfcore/base:1.14

LABEL maintainer="RNA-Seq Pipeline Team" \
      description="STAR alignment container for RNA-Seq Nextflow pipeline"

# Install STAR
RUN apt-get update && apt-get install -y \
    g++ \
    make \
    zlib1g-dev \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

RUN conda install -y star=2.7.9a \
    && conda clean -afy

# Default command
CMD ["STAR", "--version"]
