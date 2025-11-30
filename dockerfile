FROM python:3.10-slim

RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential \
        libhdf5-dev \
        libopenblas-dev \
        liblapack-dev \
        gfortran \
        curl \
        openjdk-21-jre \
        && rm -rf /var/lib/apt/lists/*

# Installing the CPU version of PyTorch makes it faster to build the image and reduces image size. Only works because I have no plans of running scVI with a GPU.
RUN pip install --no-cache-dir \
        torch --extra-index-url https://download.pytorch.org/whl/cpu \ 
        scvi-tools \
        scanpy \
        anndata \
        celltypist \
        scipy \
        pandas

RUN curl -s https://get.nextflow.io | bash && \
    mv nextflow /usr/local/bin/ && chmod +x /usr/local/bin/nextflow

WORKDIR /workspace