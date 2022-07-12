FROM 812206152185.dkr.ecr.us-west-2.amazonaws.com/latch-base:9a7d-main

RUN apt-get install -y curl unzip

# Its easy to build binaries from source that you can later reference as
# subprocesses within your workflow.
# RUN curl -L https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.4.4/bowtie2-2.4.4-linux-x86_64.zip/download -o bowtie2-2.4.4.zip &&\
# unzip bowtie2-2.4.4.zip &&\
#  mv bowtie2-2.4.4-linux-x86_64 bowtie2

# Or use managed library distributions through the container OS's package
# manager.
# RUN apt-get update -y &&\
# apt-get install -y autoconf samtools

#Install conda
ENV PATH="/root/miniconda3/bin:${PATH}"
ARG PATH="/root/miniconda3/bin:${PATH}"
RUN apt-get update
RUN apt-get install -y wget && rm -rf /var/lib/apt/lists/*
RUN wget \
   https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh \
   && mkdir /root/.conda \
   && bash Miniconda3-latest-Linux-x86_64.sh -b \
   && rm -f Miniconda3-latest-Linux-x86_64.sh 
RUN conda config --add channels bioconda --add channels conda-forge
RUN conda install -c bioconda snippy
# install sudo
# RUN apt-get install sudo

# install snippy 
RUN conda update -n base -c defaults conda &&\
   conda config --add channels bioconda
RUN apt-get update -y &&\
   apt-get install -y autoconf samtools
# create snippy env
# RUN . /opt/conda/etc/profile.d/conda.sh &&\
#conda activate base &&\
# conda create -n snippy snippy &&\
#v conda activate snippy
# STOP HERE:
# The following lines are needed to ensure your build environement works
# correctly with latch.
COPY wf /root/wf
ARG tag
ENV FLYTE_INTERNAL_IMAGE $tag
RUN python3 -m pip install --upgrade latch
WORKDIR /root
