"""GFF3 parser."""


import re

from collections import defaultdict

from g2j import fasta
from g2j.classes import Feature, Scaffold, Organism, Interval, Location


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


def set_partiality(feature):
    """

    """
    partial = "partial" in feature.qualifiers
    start_range = "start_range" in feature.qualifiers
    end_range = "end_range" in feature.qualifiers

    if not (partial and (start_range or end_range)):
        return

    if feature.location.strand == "+":
        feature.location.three_is_partial = end_range
        feature.location.five_is_partial = start_range
    else:
        feature.location.three_is_partial = start_range
        feature.location.five_is_partial = end_range


def parse(gff_handle, fasta_handle, name=None, strain=None):
    """Parse open GFF3 file handle.

    Note that individual rows are treated as separate features. For example, a
    CDS split over 5 rows would result in 5 separate Feature objects.
    """

    scaffolds = defaultdict(list)
    feature = None

    for line in gff_handle:
        if line.startswith("#"):
            continue

        scaffold, *fields = line.strip().split("\t")

        assert len(fields) == 8, "Malformed GFF, field number mismatch"

        _, type, start, end, _, strand, phase, attributes = fields

        feature = Feature(type)

        # Location
        if phase.isdigit():
            phase = int(phase)
        interval = Interval(int(start) - 1, int(end), phase=phase)
        feature.location.intervals.append(interval)
        feature.location.strand = strand

        # Any extra saved attributes
        feature.qualifiers.update(parse_attributes(attributes))

        # Partiality, assume non-partial but NCBI files will have partial=true,
        # start/end_range qualifiers
        set_partiality(feature)

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
