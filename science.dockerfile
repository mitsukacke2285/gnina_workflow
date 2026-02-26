FROM nvidia/cuda:12.2.2-runtime-ubuntu22.04

# Install Miniconda
RUN apt-get update && apt-get install -y wget bzip2 && \
    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh && \
    bash Miniconda3-latest-Linux-x86_64.sh -b -p /opt/conda && \
    rm Miniconda3-latest-Linux-x86_64.sh
 RUN apt-get update
 RUN apt-get install -y libxrender1 libxext6 libx11-6

ENV PATH=/opt/conda/bin:$PATH

WORKDIR /app

# Install Python packages
RUN pip install numpy \
    pandas \
    scipy \
    biopython \
    rdkit \
    useful-rdkit-utils \
    openbabel-wheel \
    molscrub \
    openmm \
    pdbfixer \
    mdanalysis

# Copy your workflow
COPY . .

CMD ["bash"]
