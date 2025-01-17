#
#

# ARG UBUNTU_VER=18.04
# FROM ubuntu:${UBUNTU_VER}
FROM continuumio/miniconda3:latest
# SHELL ["/bin/bash", "-c"] 

# System packages 
RUN apt-get update && apt-get install -yq wget autoconf automake bzip2 cmake make gcc perl \
    zlib1g-dev libbz2-dev liblzma-dev libcurl4-gnutls-dev libssl-dev libncurses5-dev g++ \
    build-essential libtool git pkg-config openjdk-17-jdk \
    r-base r-base-dev libxml2-dev r-cran-devtools rsync tabix



# samtools 1.21
RUN wget https://github.com/samtools/samtools/releases/download/1.21/samtools-1.21.tar.bz2 -O samtools-1.21.tar.bz2
RUN tar xvjf samtools-1.21.tar.bz2
RUN cd samtools-1.21 && \
    make && \
    make prefix=/usr/local/bin install && \
    ln -s /usr/local/bin/bin/samtools /usr/bin/samtools

# bcftools 1.21
RUN wget https://github.com/samtools/bcftools/releases/download/1.21/bcftools-1.21.tar.bz2 -O bcftools-1.21.tar.bz2
RUN tar -xjvf bcftools-1.21.tar.bz2
RUN cd bcftools-1.21 && \
    make && \
    make prefix=/usr/local/bin install && \
    ln -s /usr/local/bin/bin/bcftools /usr/bin/bcftools

# htslib 1.21
RUN wget https://github.com/samtools/htslib/releases/download/1.21/htslib-1.21.tar.bz2 -O htslib-1.21.tar.bz2
RUN tar -xjvf htslib-1.21.tar.bz2
RUN cd htslib-1.21 && \
    make && \
    make prefix=/usr/local/bin install && \
    ln -s /usr/local/bin/bin/htslib /usr/bin/htslib

# angsd and ngsrelate
RUN git clone --recursive https://github.com/SAMtools/htslib && \
    git clone https://github.com/ANGSD/ngsRelate && \
    cd htslib/;make -j2;cd ../ngsRelate;make HTSSRC=../htslib/ && \
    cd && ln -s /ngsRelate/ngsRelate /usr/bin/ngsRelate
RUN git clone https://github.com/ANGSD/angsd.git && \
    cd htslib;make;cd ../angsd ;make HTSSRC=../htslib && \
    cd && ln -s /angsd/angsd /usr/bin/angsd


# beagle 5.4
RUN wget https://faculty.washington.edu/browning/beagle/beagle.06Aug24.a91.jar -O beagle-5.4.jar 

# glimpse 2.0.1
RUN mkdir glimpse && \
    cd glimpse && \
    wget https://github.com/odelaneau/GLIMPSE/releases/download/v2.0.1/GLIMPSE2_chunk_static && \
    wget https://github.com/odelaneau/GLIMPSE/releases/download/v2.0.1/GLIMPSE2_concordance_static && \
    wget https://github.com/odelaneau/GLIMPSE/releases/download/v2.0.1/GLIMPSE2_ligate_static && \
    wget https://github.com/odelaneau/GLIMPSE/releases/download/v2.0.1/GLIMPSE2_phase_static && \
    wget https://github.com/odelaneau/GLIMPSE/releases/download/v2.0.1/GLIMPSE2_split_reference_static && \
    chmod 777 * && \
    cd && ln -s /glimpse/* /usr/bin/

# stitch 1.7.1
ARG STITCH_version=1.7.1
RUN wget -O STITCH.zip https://github.com/rwdavies/STITCH/archive/refs/tags/${STITCH_version}.zip ## or curl
RUN unzip STITCH.zip && mv STITCH-${STITCH_version} STITCH
RUN cd STITCH && ./scripts/install-dependencies.sh && \
    make install
RUN cd && ln -s /STITCH/STITCH.R /usr/bin/STITCH

# geneimp 1.3.0.9000
RUN Rscript -e 'install.packages("BiocManager")'
COPY scripts/GeneImp_1.3.0.9000.tar.gz .
RUN Rscript -e 'BiocManager::install(c("BiocGenerics", "biganalytics", "bigmemory", "Biostrings", "GenomeInfoDb", "IRanges", "Rsamtools", "S4Vectors", "SummarizedExperiment", "VariantAnnotation"))'
RUN R CMD INSTALL GeneImp_1.3.0.9000.tar.gz

# vcftools 0.1.16
RUN git clone https://github.com/vcftools/vcftools.git && \
    cd vcftools && \
    ./autogen.sh && \
    ./configure && \
    make && \
    make install 

# imputation_accuracy_calculator
RUN git clone https://github.com/TorkamaniLab/imputation_accuracy_calculator.git
RUN pip install pandas numpy==1.26.4 cyvcf2


ENV BCFTOOLS_PLUGINS=/bcftools-1.21/plugins

# QUILT 2.0.1
RUN git clone --recursive https://github.com/rwdavies/QUILT.git && \
    cd QUILT && \
    bash ./scripts/install-dependencies.sh && \
    wget https://github.com/rwdavies/QUILT/releases/download/2.0.1/QUILT_2.0.1.tar.gz && \
    R CMD INSTALL QUILT_2.0.1.tar.gz

RUN apt-get install -yq bc

# R packages
RUN Rscript -e 'install.packages("vcfR")'