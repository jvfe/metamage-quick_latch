FROM 812206152185.dkr.ecr.us-west-2.amazonaws.com/latch-base:6839-main

RUN apt-get install -y curl unzip libz-dev

# BowTie2
RUN curl -L https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.4.4/bowtie2-2.4.4-linux-x86_64.zip/download -o bowtie2-2.4.4.zip &&\
    unzip bowtie2-2.4.4.zip &&\
    mv bowtie2-2.4.4-linux-x86_64 bowtie2

# Samtools
RUN apt-get update -y &&\
    apt-get install -y autoconf samtools

# Get miniconda
RUN curl https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh --output miniconda.sh
ENV CONDA_DIR /opt/conda
RUN bash miniconda.sh -b -p /opt/conda
ENV PATH=$CONDA_DIR/bin:$PATH

# Get Mamba
RUN conda install mamba -n base -c conda-forge

# Get Macrel
RUN mamba install -y -c bioconda macrel

# Get FarGene
RUN mamba create -y -n fargene_env python=2.7
RUN mamba install -y -n fargene_env -c bioconda fargene
ENV PATH=$CONDA_DIR/envs/fargene_env/bin:$PATH

# Get Gecco
RUN pip3 install gecco-tool

# Get MegaHIT and Quast
RUN mamba create -y -n metassembly python=3.6
RUN mamba install -y -n metassembly -c bioconda megahit quast
ENV META_ENV $CONDA_DIR/envs/metassembly/bin

# Create symlink
RUN ln -s $META_ENV/megahit /root/megahit &&\
    ln -s $META_ENV/metaquast.py /root/metaquast.py 

# Install metabat2

RUN mamba install -y -c bioconda metabat2

# STOP HERE:
# The following lines are needed to ensure your build environement works
# correctly with latch.
RUN python3 -m pip install --upgrade latch
COPY wf /root/wf
ARG tag
ENV FLYTE_INTERNAL_IMAGE $tag
WORKDIR /root

