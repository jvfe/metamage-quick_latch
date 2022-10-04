FROM 812206152185.dkr.ecr.us-west-2.amazonaws.com/latch-base:6839-main

RUN apt-get install -y curl unzip libz-dev

# Kaiju installation
RUN curl -L \
    https://github.com/bioinformatics-centre/kaiju/releases/download/v1.9.0/kaiju-v1.9.0-linux-x86_64.tar.gz -o kaiju-v1.9.0.tar.gz &&\
    tar -xvf kaiju-v1.9.0.tar.gz

ENV PATH /root/kaiju-v1.9.0-linux-x86_64-static:$PATH

# Krona installation
RUN curl -L \
    https://github.com/marbl/Krona/releases/download/v2.8.1/KronaTools-2.8.1.tar \
    -o KronaTools-2.8.1.tar &&\
    tar -xvf KronaTools-2.8.1.tar &&\
    cd KronaTools-2.8.1 &&\
    ./install.pl
    