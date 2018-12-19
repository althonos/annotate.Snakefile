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
        filter(os.path.exists, ("blast", "catalog", "db", "features", "output"))
    run:
        for directory in input:
            shutil.rmtree(directory)


rule makeCatalog:
    input:
        "reference/{sequences}"
    output:
        directory("features/{sequences}"), "catalog/{sequences}.fa"
    run:
        features = {}
        # extract each feature into an individual record
        for sequence in glob.iglob(os.path.join(input[0], "*.gb")):
            record = Bio.SeqIO.read(sequence, 'genbank')
            for idx, feat in enumerate(record.features):
                # NB: do not use `feat.extract` to keep the feature itself
                feat_record = record[feat.location.start:feat.location.end]
                feat_record.features = [f for f in feat_record.features if f.location.start == 0 and f.location.end >= len(feat_record)]
                if feat.strand < 0:
                    feat_record = feat_record.reverse_complement(features=True)
                feat_record.id = feat_record.name = "{}.{}".format(record.id, idx)
                features[str(feat_record.seq)] = feat_record
        # write each feature as an individual genbank file
        os.makedirs(output[0], exist_ok=True)
        for feat_record in features.values():
            with open(os.path.join(output[0], "{}.gb".format(feat_record.id)), 'w') as f:
                Bio.SeqIO.write(feat_record, f, 'genbank')
        # write all features into a single FASTA file
        with open(output[1], 'w') as f:
            Bio.SeqIO.write(features.values(), f, 'fasta')


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

