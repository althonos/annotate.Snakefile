# `annotate.Snakefile`

*A Snakemake pipeline to copy annotations between GenBank files*

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
