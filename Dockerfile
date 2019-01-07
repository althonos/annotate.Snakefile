#FROM python:alpine
FROM debian:9-slim

ADD Snakefile Snakefile
ADD COPYING COPYING

ENV LANG=C.UTF-8 BLAST_VERSION=2.8.1 GIT_PYTHON_REFRESH=quiet
RUN apt update \
 && apt install -y upx-ucl wget python3 python3-pip \
 && wget "http://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/${BLAST_VERSION}/ncbi-blast-${BLAST_VERSION}+-x64-linux.tar.gz" -O - \
  | tar xz -C /bin --strip-components 2 -- "ncbi-blast-${BLAST_VERSION}+/bin/makeblastdb" "ncbi-blast-${BLAST_VERSION}+/bin/blastn" \
 && upx /bin/makeblastdb /bin/blastn \
 && mkdir -p /reference/moclo \
 && wget "https://github.com/althonos/moclo/archive/master.tar.gz" -O - \
  | tar xz -C /reference/moclo --strip-components 4 --wildcards '*.gb' \
 && python3 -m pip --no-cache-dir install snakemake biopython \
 && apt remove -y gcc libpython3.5-dev dpkg-dev openssh-client wget upx-ucl \
 && apt autoremove -y

CMD snakemake

