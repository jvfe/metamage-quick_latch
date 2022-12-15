FROM 812206152185.dkr.ecr.us-west-2.amazonaws.com/latch-base:6839-main

RUN apt-get update -y &&\
    apt-get install -y curl unzip libz-dev

# Get BowTie2, Samtools and Prodigal
RUN curl -L https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.4.4/bowtie2-2.4.4-linux-x86_64.zip/download -o bowtie2-2.4.4.zip &&\
    unzip bowtie2-2.4.4.zip &&\
    mv bowtie2-2.4.4-linux-x86_64 bowtie2 &&\
    apt-get update -y &&\
    apt-get install -y autoconf samtools &&\
    curl -L https://github.com/hyattpd/Prodigal/releases/download/v2.6.3/prodigal.linux -o prodigal &&\
    chmod +x prodigal

# Get miniconda
RUN curl https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh --output miniconda.sh
ENV CONDA_DIR /opt/conda
RUN bash miniconda.sh -b -p /opt/conda
ENV PATH=$CONDA_DIR/bin:$PATH

# Get Mamba
RUN conda install mamba -n base -c conda-forge

# Get MegaHIT and Quast
RUN mamba create -y -n metassembly python=3.6 &&\
    mamba install -y -n metassembly -c bioconda megahit quast &&\
    ln -s $CONDA_DIR/envs/metassembly/bin/megahit /root/megahit &&\
    ln -s $CONDA_DIR/envs/metassembly/bin/metaquast.py /root/metaquast.py

# Get metabat2, macrel, Gecco and FarGene

RUN mamba install -y -c bioconda metabat2 &&\
    mamba install -y -c bioconda macrel &&\
    pip3 install gecco-tool &&\
    mamba create -y -n fargene_env python=2.7 &&\
    mamba install -y -n fargene_env -c bioconda fargene

ENV PATH=$CONDA_DIR/envs/fargene_env/bin:$PATH

# Kaiju installation
RUN curl -L \
    https://github.com/bioinformatics-centre/kaiju/releases/download/v1.9.0/kaiju-v1.9.0-linux-x86_64.tar.gz -o kaiju-v1.9.0.tar.gz &&\
    tar -xvf kaiju-v1.9.0.tar.gz

ENV PATH /root/kaiju-v1.9.0-linux-x86_64-static:$PATH

# Krona installation
RUN curl -L \
    https://github.com/marbl/Krona/releases/download/v2.8.1/KronaTools-2.8.1.tar \
    -o KronaTools-2.8.1.tar &&\
    tar -xvf KronaTools-2.8.1.tar --no-same-owner &&\
    cd KronaTools-2.8.1 &&\
    ./install.pl

# STOP HERE:
# The following lines are needed to ensure your build environement works
# correctly with latch.
RUN python3 -m pip install --upgrade latch
COPY wf /root/wf
ARG tag
ENV FLYTE_INTERNAL_IMAGE $tag
WORKDIR /root
