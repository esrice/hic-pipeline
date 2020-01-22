FROM ubuntu:18.04

WORKDIR /app

RUN apt-get -y update
RUN apt-get -y install python2.7 python-pip python3 python3-pip wget \
    build-essential libboost-all-dev git libbz2-dev liblzma-dev zlib1g-dev

# install bwa
RUN git clone https://github.com/lh3/bwa.git
RUN cd bwa && make
RUN ln -s /app/bwa/bwa /usr/bin/bwa

# install samtools
RUN wget https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2
RUN tar xf samtools-1.9.tar.bz2
RUN cd samtools-1.9 && ./configure --without-curses && make && make install

# install bedtools
RUN wget https://github.com/arq5x/bedtools2/releases/download/v2.28.0/bedtools
RUN chmod +x bedtools && mv bedtools /usr/bin/

# install pipeline scripts
RUN pip3 install Cython
RUN pip3 install pysam
RUN git clone https://github.com/esrice/slurm-hic.git
RUN mv slurm-hic/*.py /usr/bin/

# install SALSA and dependencies
RUN pip install numpy networkx==1.1
RUN wget https://github.com/machinegun/SALSA/archive/v2.2.tar.gz
RUN tar xf v2.2.tar.gz && cd SALSA-2.2 && make

ENV SALSA_DIR "/app/SALSA-2.2"
