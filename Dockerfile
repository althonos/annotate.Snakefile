FROM python:alpine

ADD Snakefile Snakefile
ADD COPYING COPYING

ENV LANG=C.UTF-8 BLAST_VERSION=2.8.1
RUN apk add --no-cache --update --virtual=.build-dependencies build-base upx tar python3-dev \
 && wget "http://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/${BLAST_VERSION}/ncbi-blast-${BLAST_VERSION}+-x64-linux.tar.gz" -O - \
  | tar xz -C /bin --strip-components 2 -- "ncbi-blast-${BLAST_VERSION}+/bin/makeblastdb" "ncbi-blast-${BLAST_VERSION}+/bin/blastn" \
 && upx /bin/makeblastdb /bin/blastn \
 && mkdir -p /reference/moclo \
 && wget "https://github.com/althonos/moclo/archive/master.tar.gz" -O - \
  | tar xz -C /reference/moclo --strip-components 4 --wildcards '*.gb' \
 && python -m pip --no-cache-dir install snakemake \
 && apk del .build-dependencies

CMD snakemake

