import glob
import os
import shutil

import Bio.SeqIO
from Bio.Blast.Applications import NcbiblastnCommandline

# --- Configuration -------------------------------

REFERENCE = config.get("reference", os.path.basename(next(glob.iglob("reference/*"))))


# --- Rules ---------------------------------------


rule all:
    input:
        [x.replace("input", "output").replace(".gb", ".annotated.gb") for x in glob.iglob("input/*.gb")]

rule clean:
    input:
        filter(os.path.exists, ["blast", "catalog", "db", "output"])
    run:
        for directory in input:
            shutil.rmtree(directory)


rule makeCatalog:
    """Build a catalog of all reference sequences in a single FASTA file.
    """
    input:
        "reference/{sequences}"
    output:
        "catalog/{sequences}.fa"
    run:
        seqs = []
        for f in glob.iglob(os.path.join(input[0], "*.gb")):
            seq = Bio.SeqIO.read(open(f), 'gb')
            seq.id = f
            seqs.append(seq)
        with open(output[0], "w") as out:
            Bio.SeqIO.write(seqs, out, "fasta")

rule makeBlastDb:
    """Build a BLASTn database from a catalog of sequences.
    """
    input:
        "catalog/{sequences}.fa"
    output:
        expand("db/{{sequences}}.{ext}", ext=["nhr", "nin", "nsq"])
    log:
        "db/{sequences}.log"
    wildcard_constraints:
        sequences = ".+"
    shell:
        "makeblastdb -dbtype nucl -out db/{wildcards.sequences} -title 'Annotations' -in {input} 2>&1 >{log}"


rule runBlastN:
    """Run a BLASTn query against a reference BLASTn database.
    """
    input:
        "input/{sequence}.gb", "db/{reference}.nhr"
    output:
        "blast/{sequence}.{reference}.blastxml"
    log:
        "blast/{sequence}.{reference}.log"
    shell:
        "blastn -db db/{wildcards.reference} -query {input[0]}  >{output}  2>{log}"


rule copyAnnotations:
    """Annotate a sequence using the BLASTn result.
    """
    input:
        expand("blast/{{sequence}}.{reference}.blastxml", reference=REFERENCE)
    output:
        "output/{sequence}.annotated.gb"
    run:
        shutil.copy(f"input/{wildcards.sequence}.gb", output[0])

