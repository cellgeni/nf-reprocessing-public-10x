FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive 

RUN apt-get update && apt-get -y --no-install-recommends -qq install \
    wget gcc build-essential software-properties-common libz-dev \
    git libncurses5-dev libbz2-dev liblzma-dev default-jre bsdmainutils

#Install seqtk
RUN git clone https://github.com/lh3/seqtk.git && \
    mv seqtk /opt && \
    cd /opt/seqtk && \
    make

#Download bamtofastq
RUN wget 'https://github.com/10XGenomics/bamtofastq/releases/download/v1.4.1/bamtofastq_linux' -O bamtofastq && \
    mv bamtofastq /opt/bamtofastq

#Install SRA toolkit
RUN wget 'https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-ubuntu64.tar.gz' -O sratoolkit.tar.gz && \
   tar -xzf sratoolkit.tar.gz -C /opt/

RUN chmod -R 755 /opt

ENV PATH="/opt/seqtk:/opt/sratoolkit.3.0.5-ubuntu64/bin:/opt:${PATH}"
