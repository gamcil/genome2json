"""GFF3 parser."""


import re

from collections import defaultdict

from g2j import fasta
from g2j.classes import Feature, Scaffold, Organism


def parse_attributes(string):
    """Parse keyword:value pairs from the attributes field.

    Uses regular expression to allow for different formatting.

    For example, a GFF3 from JGI is formatted like:
        name "gm1.1_g"; proteinId 412464; exonNumber 3

    From FungiDB:
        ID=AN11191-T-p1-CDS1;Parent=AN11191-T

    From NCBI:
        ID=id1;Parent=rna0;Dbxref=GeneID:36531755
    """
    pattern = re.compile(r' ?(\w+?)[ =]"?(.+?)"?(?:;|$)')
    attributes = {}
    for match in pattern.finditer(string):
        key, value = match.groups()
        attributes.update({key: value})
    return attributes


def parse(gff_handle, fasta_handle, name=None, strain=None):
    """Parse open GFF3 file handle.

    Note that individual rows are treated as separate features. For example, a CDS split
    over 5 rows would result in 5 separate Feature objects.
    """

    scaffolds = defaultdict(list)
    feature = None

    for line in gff_handle:
        if line.startswith("#"):
            continue

        scaffold, *fields = line.strip().split("\t")

        assert len(fields) == 8, "Malformed GFF, field number mismatch"

        source, type, start, end, _, strand, _, attributes = fields

        feature = Feature(type)

        # Location
        feature.location.strand = strand
        feature.location.start.append(int(start))
        feature.location.end.append(int(end))

        # Any extra saved attributes
        feature.qualifiers.update(parse_attributes(attributes))

        # Save to scaffold
        scaffolds[scaffold].append(feature)

    # Get assembly scaffold sequences from FASTA file
    sequences = fasta.parse(fasta_handle)

    # Instantiate and return Organism
    return Organism(
        name,
        strain,
        scaffolds=[
            Scaffold(accession, sequences[accession], features)
            for accession, features in scaffolds.items()
        ],
    )
