FROM continuumio/miniconda3

WORKDIR /app

# Install CUDA
RUN apt-get update
RUN apt-get install -y wget gnupg wget https://developer.download.nvidia.com/compute/cuda/repos/ubuntu2204/x86_64/cuda-keyring_1.1-1_all.deb
RUN dpkg -i cuda-keyring_1.1-1_all.deb
RUN apt-get update
RUN apt-get install -y cuda-runtime-12-2



# Install pip-only packages (note the hyphenated name on PyPI)
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


# Copy your project files
COPY . .

CMD ["bash"]
