# `annotate.Snakefile`

*A Snakemake pipeline to copy annotations between GenBank files*

[![Source](https://img.shields.io/badge/source-GitHub-303030.svg?maxAge=86400&style=flat-square)](https://github.com/althonos/annotate.Snakefile/)
[![Docker](https://img.shields.io/badge/docker%20build-automatic-blue.svg?maxAge=86400&style=flat-square)](https://hub.docker.com/r/althonos/annotate/)
[![Snakemake](https://img.shields.io/badge/built%20with-snakemake-yellowgreen.svg?maxAge=86400&style=flat-square)](https://snakemake.readthedocs.io/en/stable/)


## Usage

### In-site pipeline

1. Clone the repository
2. Create a `reference` directory and put your reference sequences in it in Genbank format (i.e. `reference/myseqs/seq1.gb`, `reference/myseqs/seq2.gb`, ...)
3. Create an `input` directory and put your sequences in it in Genbank format.
4. Run `snakemake`.
5. Get your annotated sequences in the `output` directory.

### With `docker`

```
docker run -v $(pwd):/input -v $(pwd):/output althonos/annotate
```
will annotate all the sequences in your current directory using the annotated sequences
distributed with the [`moclo`](https://github.com/althonos/moclo) library.


## Dependencies

* `blast+` binaries (`blastn` and `makeblastdb`).
* `snakemake`
* `biopython`


## About

This pipeline was developed by [Martin Larralde](https://github.com/althonos) for the
*data integration* course of the [AMI2B Master's degree](http://www.bibs.u-psud.fr/m2_ami2b.php)
of [Université Paris-Sud](https://www.u-psud.fr).
