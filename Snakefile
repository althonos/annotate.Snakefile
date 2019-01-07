import glob
import os
import shutil

import Bio.SeqIO
import Bio.Blast.NCBIXML
from Bio.Blast.Applications import NcbiblastnCommandline

# --- Configuration -------------------------------

REFERENCE = config.get("reference", os.path.basename(next(glob.iglob("reference/*"))))
CLEAR = config.get("clear", False)
IDENTITY_THRESHOLD = float(config.get("identity", 0.98))

assert 0 <= IDENTITY_THRESHOLD <= 1, "`identity` must be between 0 and 1"


# --- Maintenance rules ------------------------------------------------------

rule all:
    """Annotate all sequences in the `input` directory.
    """

    input:
        [x.replace("input", "output").replace(".gb", ".annotated.gb") for x in glob.iglob("input/*.gb")]


rule clean:
    """Remove temporary files.
    """

    input:
        filter(os.path.exists, ("blast", "catalog", "db", "features", "output"))

    run:
        for directory in input:
            shutil.rmtree(directory)


# --- Step 1: Extract features from reference GenBank files ------------------

rule extractFeatures:
    """Build a catalog of features from reference sequences.
    """

    input:
        "reference/{reference}"

    output:
        directory("features/{reference}")

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
                feat_record.description = ""
                features[str(feat_record.seq)] = feat_record
        # write each feature as an individual genbank file
        os.makedirs(output[0], exist_ok=True)
        for feat_record in features.values():
            with open(os.path.join(output[0], "{}.gb".format(feat_record.id)), 'w') as f:
                Bio.SeqIO.write(feat_record, f, 'genbank')


# --- Step 2: Build a FASTA catalog of reference features --------------------

rule makeCatalog:
    """Build a FASTA catalog of sequences from a directory of GenBank files.
    """

    input:
        "features/{reference}"

    output:
        "catalog/{reference}.fa"

    run:
        features = [Bio.SeqIO.read(gb, 'gb') for gb in glob.iglob(os.path.join(input[0], "*.gb"))]
        with open(output[0], 'w') as f:
            Bio.SeqIO.write(features, f, 'fasta')


# --- Step 3: Build a BlastDB of reference features --------------------------

rule makeBlastDb:
    """Build a BLASTn database from a catalog of sequences.
    """

    input:
        "catalog/{reference}.fa"

    output:
        expand("db/{{reference}}.{ext}", ext=["nhr", "nin", "nsq"])

    log:
        "db/{reference}.log"

    shell:
        "makeblastdb -dbtype nucl -out db/{wildcards.reference} -title 'Annotations' -in {input} 2>&1 >{log}"


# --- Step 4: Run each input sequence against the database of features -------

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
        "blastn -task blastn-short -ungapped -outfmt 5 -db db/{wildcards.reference} -query {input[0]} >{output} 2>{log}"


# --- Step 5: Copy features detected with `blastn` to the query --------------

rule copyFeatures:
    """Copy features from reference sequences using BLASTn results as a guide.
    """

    input:
        "input/{sequence}.gb",
        expand("blast/{{sequence}}.{reference}.blastxml", reference=REFERENCE),
        expand("features/{reference}", reference=REFERENCE)

    output:
        "output/{sequence}.annotated.gb"

    run:
        # load the query record input sequence
        record = Bio.SeqIO.read(input[0], "genbank")
        with open(input[1]) as handle:
            result = next(Bio.Blast.NCBIXML.parse(handle))

        # remove existing features if desired
        if CLEAR:
            record.features.clear()

        # FIXME: copy features only if alignment has 100% identity
        for alignment in result.alignments:
            hsp = alignment.hsps[0]
            if (hsp.identities / hsp.align_length) > IDENTITY_THRESHOLD:
                # get the GenBank file corresponding to the feature
                hit_gb = os.path.join(input[2], "{}.gb".format(alignment.hit_def))
                hit_record = Bio.SeqIO.read(hit_gb, "genbank")
                # use the reverse complement if it is on the reverse strand
                if hit_record.seq not in record.seq:
                    hit_record = hit_record.reverse_complement(features=True)
                # position the feature on the query record
                try:
                    hit_feat = hit_record.features[0]
                    hit_feat.location += record.seq.find(hsp.query)
                except IndexError:
                    continue

                # add the feature to the query record
                record.features.append(hit_feat)

        # write the query record
        with open(output[0], 'w') as f:
            Bio.SeqIO.write(record, f, "genbank")
