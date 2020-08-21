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


def parse(
    gff_handle,
    fasta_handle=None,
    name=None,
    strain=None,
    feature_types=None,
    save_scaffold_sequence=True,
):
    """Parse a GFF3 file.

    Note that individual rows are treated as separate features. For example, a
    CDS split over 5 rows will result in 5 separate Feature objects.
    """

    header = None
    sequence, sequences = None, {}
    fasta_flag = False

    qualifiers = {}
    scaffolds = defaultdict(list)
    feature = None

    for line in gff_handle:
        if line.startswith("###"):
            continue

        # Check for the ##FASTA directive.
        # If present, and save_scaffold_sequence=True, parse scaffold
        # sequences and populate the sequences dictionary.
        if save_scaffold_sequence and line.startswith("##FASTA"):
            fasta_flag = True
            continue
        if fasta_flag:
            if line.startswith(">"):
                if sequence and header:
                    sequences[header] = sequence
                header = line[1:].strip()
                sequence = ""
            else:
                sequence += line.strip()
            continue

        # Check for other official directives, e.g. ##species
        if line.startswith(("##", "#!")):
            key, value = line[2:].strip().split(" ", 1)
            if key in ("sequence-region"):
                continue
            qualifiers[key] = value

        # Ignore comments
        if line.startswith("#"):
            continue

        # Should be at a delimited section now, so parse as normal
        scaffold, *fields = line.strip().split("\t")

        if not len(fields) == 8:
            raise ValueError("Malformed GFF, field number mismatch")

        _, type, start, end, _, strand, phase, attributes = fields

        # Filter out any non-specified feature types
        if feature_types and type not in feature_types:
            continue

        feature = Feature(type)

        # Location
        if phase.isdigit():
            phase = int(phase)
        interval = Interval(int(start) - 1, int(end), phase=phase)
        feature.location.intervals.append(interval)
        feature.location.strand = strand

        # Any extra saved attributes
        attributes = parse_attributes(attributes)
        feature.qualifiers.update(attributes)

        # Partiality, assume non-partial but NCBI files will have partial=true,
        # start/end_range qualifiers
        set_partiality(feature)

        # Save to scaffold
        scaffolds[scaffold].append(feature)

    # Get assembly scaffold sequences from FASTA file
    if save_scaffold_sequence and not sequences:
        if fasta_handle:
            sequences = fasta.parse(fasta_handle)
        else:
            raise ValueError(
                "save_scaffold_sequence=True, but GFF contains no ##FASTA directive"
                " and no FASTA file was provided."
            )

    # Instantiate and return Organism
    return Organism(
        name,
        strain,
        qualifiers=qualifiers,
        scaffolds=[
            Scaffold(
                accession,
                sequences[accession] if save_scaffold_sequence else "",
                features
            )
            for accession, features in scaffolds.items()
        ],
    )
