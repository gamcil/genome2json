"""
Module containing parsers for common file formats.

The goal of this module is to provide a basic, unified interface for
GenBank, GFF3 and FASTA format files (i.e. once parsed, they all end
up looking and acting the same way, regardless of differing origin).

Then, models can be instantiated by calling these parsers; for example,
calling Protein.from_genbank('genome.gbk'), would call the classmethod
GenBank.parse('genome.gbk') and map its attributes to a live Protein
instance ready for insertion.

    Record (Scaffold)
      -> Feature (List of sequence features grouped by size)
           -> Parts (List of locations representing gene model features -
                     genes, mRNA, CDS, etc.)
"""

import re
import warnings

from collections import defaultdict, namedtuple
from itertools import groupby
from operator import itemgetter, attrgetter
from pathlib import Path

from Bio import BiopythonWarning, Seq, SeqIO
from Bio.SeqFeature import FeatureLocation

from helpers import Location

warnings.simplefilter("ignore", BiopythonWarning)


class Parser:
    """Base Parser class to be subclassed for specific file formats."""

    gene_related = {
        "gene",
        "mRNA",
        "rRNA",
        "tRNA",
        "CDS",
        "tRNA",
        "ORF",
        "3'UTR",
        "5'UTR",
        "intron",
        "exon",
        "three_prime_UTR",
        "five_prime_UTR",
    }

    def __init__(self, file):
        if not Path(file).exists():
            raise ValueError("No file found at given path")

        self.file = file
        self.parsed = False
        self.records = defaultdict(dict)

    @classmethod
    def from_file(cls, file):
        """Method for parsing the file."""
        instance = cls(file)
        instance.parse_file()
        return instance

    def get_sequence(self, record):
        return self.records[record].sequence

    def get_features(self, record):
        return self.records[record].features

    def parse_file(self):
        """_parse_file() wrapper to enable setting parsed flag to True."""
        self._parse_file()
        self.parsed = True

    def _parse_file(self):
        """Actual method for parsing open file handles; override this."""
        return NotImplementedError

    @staticmethod
    def features_overlap(one, two):
        """Method for testing feature overlap in grouping.

        Should be overridden in child classes to enable format-specific
        testing of feature overlap. For example, GenBank files parsed
        using SeqIO will return SeqFeatures with SeqLocation objects that
        have start and end attributes, whereas a parsed GFF3 file will
        have a Python list of fields after splitting a feature line.
        """
        raise NotImplementedError

    def group_overlapping_features(self, features):
        """Group overlapping sequence features.

        Given a sorted list of features, and format-specific overrides
        for self.features_overlap(), this generic iterator function will
        group sequence features via location overlap.

        This removes the dependence on shared identifier attributes
        (i.e. locus tags or protein IDs), which are unreliable.
        """
        i, total = 0, len(features)
        while i < total:
            current = features[i]  # grab current
            group = [current]
            for j in range(i + 1, total):  # iterate rest
                future = features[j]
                if self.features_overlap(current, future):
                    group.append(future)  # add if contained
                else:
                    yield group  # else yield group to iterator
                    break
            i += len(group)  # move index ahead of last group
        yield group


class GenBank(Parser):
    """Parse GenBank file."""

    @staticmethod
    def features_overlap(one, two):
        """Check SeqFeature one contains SeqFeature two."""
        return (
            two.location.start >= one.location.start
            and two.location.end <= one.location.end
        )

    def parse_wraparound_features(self, record):
        """Take care of wraparound SeqFeatures on a circular SeqRecord.

        For example, in a circular bacterial record, a gene feature could be split
        across the origin:
        >>> feature.location
        CompoundLocation(
            [
                FeatureLocation(ExactPosition(5226888), ExactPosition(5227080), strand=1),
                FeatureLocation(ExactPosition(0), ExactPosition(1331), strand=1)
            ],
            'join'
        )

        If this special case is not considered, the feature overlap logic will choke
        because the start/end of a SeqFeature reflect the min/max of coordinates, NOT
        the biological start/end of the gene. Thus, the start and end coordinates will
        just be 0 and the total length of the SeqRecord:
        >>> feature.location.start, feature.location.end
        (0, 5227080)

        Then, ALL following features will be considered as overlapping since they all
        lie within this range.

        When parsed by BioPython, wraparound features will always be the first features
        in a SeqRecord following a 'source' feature. This function shifts the origin of
        the SeqRecord to the end of the wraparound gene. This always removes the feature
        directly preceeding the split point, thus wraparound features are saved, have
        their FeatureLocations changed, and are converted to a Feature namedtuple
        separately.
        """
        record_length = len(record)
        wraparounds = []
        for feature in record.features:
            if feature.type not in self.gene_related:
                continue
            if feature.location.start <= record_length <= feature.location.end:
                wraparounds.append(feature)
                continue
            break

        start = min(wraparounds[0].location.parts, key=lambda x: x.start + x.end).end

        for feature in wraparounds:
            feature.location = FeatureLocation(
                record_length - len(feature),
                record_length,
                strand=feature.location.strand,
            )

        shifted = record[start:] + record[:start]
        feature = self.build_feature(wraparounds, shifted.seq) if wraparounds else None
        return feature, shifted

    def _parse_file(self):
        """Parse a BioPython SeqRecord object for features."""
        for record in self.record_iterator:
            if record.annotations["topology"] == "circular":
                wraparound, shifted = self.parse_wraparound_features(record)
                if wraparound:
                    record = shifted

            scaffold, sequence = record.id, record.seq

            # Filter out non-gene related features, sort by start
            features = sorted(
                [f for f in record.features if f.type in self.gene_related],
                key=attrgetter("location.start"),
            )

            # Build Feature objects from filtered SeqFeatures, save Records
            self.records[scaffold] = self.Record(
                scaffold,
                str(sequence),
                [
                    self.build_feature(group, sequence)
                    for group in self.group_overlapping_features(features)
                ],
            )

            if wraparound:
                self.records[scaffold].features.append(wraparound)

    def build_feature(self, group, sequence):
        """Parse a group of overlapping features to build Feature tuples.

        Attempts to retrieve translation directly from feature
        qualifiers if one is stored; otherwise falls back to BioPython.
        """
        identifiers = {"locus_tag": None, "protein_id": None}
        translation = ""
        parts = {}
        notes = {}

        for element in group:
            if element.type not in self.gene_related:
                continue

            parts[element.type] = Location.from_seqlocation(element.location)

            if element.type == "gene":
                notes["pseudo"] = "pseudo" in element.qualifiers

            if element.type == "tRNA":
                notes["tRNA"] = extract(element.qualifiers, "product")

            if element.type == "CDS":
                try:
                    translation = element.qualifiers["translation"][0]
                except KeyError:
                    translation = str(element.extract(sequence).translate())

            if "product" not in notes and "product" in element.qualifiers:
                notes["product"] = element.qualifiers["product"][0]

            # Try and find identifiers if haven't already gotten them
            for key, value in identifiers.items():
                if not value and key in element.qualifiers:
                    identifiers[key] = element.qualifiers[key][0]

        return self.Feature(identifiers, parts, translation, notes)


class GFF3(Parser):
    """Parse features from a GFF3 file.

    Should be used alongside a FASTA parser, such that features in the
    GFF3 file can be mapped to their sequences in the FASTA file:

    Instantiate the FASTA parser
    >>> fasta = FASTA.from_file('genome.fasta')

    Initialise GFF3 directly with the FASTA parser, then parse
    >>> gff3 = GFF3('genome.gff3', fasta_parser=fasta)
    >>> gff3.parse_file()

    Instantiate without, then add, the FASTA parser, and parse
    >>> gff3 = GFF3('genome.gff3')
    >>> gff3.add_fasta(fasta)
    >>> gff3.parse_file()

    Or just use the parse() factory method
    >>> gff3 = GFF3.from_file('genome.gff3', fasta_parser=fasta)

    Note that this class can parse a GFF3 file without a FASTA parser,
    for example:
    >>> gff3 = GFF3.parse('genome.gff3')

    However, the resulting records that it stores will be missing both
    genomic sequence and protein translations.
    """

    def __init__(self, file, fasta_parser=None):
        super().__init__(file)

        if fasta_parser:
            self.set_fasta_parser(fasta_parser)
        else:
            self.fasta = None

    @classmethod
    def from_file(cls, file, fasta_parser=None):
        """Override of from_file() classmethod to allow FASTA parser input."""
        gff3 = cls(file, fasta_parser=fasta_parser)
        gff3.parse_file()
        return gff3

    def set_fasta_parser(self, fasta_parser):
        """Attach a FASTA parser object to this instance."""
        if not isinstance(fasta_parser, FASTA):
            raise ValueError("Expected a FASTA parser object")

        if not fasta_parser.parsed:
            raise ValueError("FASTA parser has not parsed its file")

        self.fasta = fasta_parser

    def _parse_file(self):
        """Parse features from the GFF3 file.

        Have to loop twice: first pass collects raw feature lines on
        each scaffold; second pass groups those features and builds
        consolidated Feature objects.
        """
        raw_features = defaultdict(list)

        # TODO: use Table.from_file() for this
        with open(self.file, "r") as handle:
            for line in handle:
                if line.startswith("#"):
                    continue

                scaffold, *feature_fields = line.strip().split("\t")

                total_fields = len(feature_fields)
                if total_fields != 8:
                    raise ValueError(
                        f"Feature line contains only {total_fields}"
                        " columns. Your GFF3 file may be malformed."
                    )

                raw_features[scaffold].append(feature_fields)

        # Check if this is a JGI file; will have 'jgi' in each line.
        # If so, grab the portal ID and pass to build_feature() to avoid
        # thousands of unecessary regex searches
        first_feature = list(raw_features.values())[0][0][-1]
        if "jgi" in first_feature:
            portal_id = re.search("portal_id=(.+?);", first_feature).group(1)
        else:
            portal_id = None

        for scaffold, feature_fields in raw_features.items():
            # Get scaffold sequence if this object has a FASTA parser
            sequence = self.fasta.get_sequence(scaffold) if self.fasta else None

            # Filter out non-gene related feature types, sort by start
            features = [f for f in feature_fields if f[1] in self.gene_related]

            # Build Feature objects from filtered feature lines
            self.records[scaffold] = self.Record(
                scaffold,
                sequence,
                [
                    self.build_feature(group, scaffold, jgi_portal=portal_id)
                    for group in self.group_overlapping_features(features)
                ],
            )

    @staticmethod
    def features_overlap(one, two):
        """Check feature one contains feature two."""
        return int(two[2]) >= int(one[2]) and int(two[3]) <= int(one[3])

    def build_feature(self, group, scaffold, jgi_portal=None):
        """Build a full Feature object from a group of features.

        If this object has an attached FASTA parser, this function
        will attempt to translate the features CDS region.

        Considerations for different file sources:
        - AspGD/gbrowser derived GFF3s have 'Sequence:' in their
          feature lines, and only a single ID for both gene/protein

        - JGI GFF3 files have 'jgi' in their feature lines, and a
          portal id (e.g. A. alliaceus --> Aspalli1) field. This
          is used to build unique locus tags and protein IDs, since
          JGI stores only generic gene/CDS_0001 strings in the ID
          field, and a plain number in proteinId. A JGI Feature
          object might then resemble:

            Feature(identifiers={'locus_tag': 'Aspalli1_00001',
                                 'protein_id': 'Aspalli1_21251'},
                    parts=[Part(type='gene', ...), ...],
                    translation='M...')

          JGI also stores features in strand orientation; i.e. CDS
          features for a negative strand gene are listed from their
          actual start (e.g. 500-600, 350-450, 100-300, ...).

        - NCBI files actually use 'locus_tag' and 'protein_id' as
          identifier fields. So, if no telltale JGI/gbrowser signs
          are found, this function will eventually try to extract
          those fields specifically.

        Arguments:
            group (list): feature lines from parsed GFF3
            scaffold (str): name of scaffold
            jgi_portal (str): JGI portal ID

        Returns:
            feature (namedtuple): Feature object
        """
        identifiers = {"locus_tag": None, "protein_id": None}
        gene_id, translation = "", ""
        parts = {}
        notes = {}

        # Sort by feature type, ensuring gene features are first
        group.sort(key=lambda f: (f[1] not in {"gene", "ORF"}, f[1]))
        for feature_type, features in groupby(group, itemgetter(1)):
            features.sort(key=itemgetter(2))

            # Get gene ID field; this is used as basis of custom identifiers
            if feature_type in {"gene", "ORF"}:
                gene_id = re.search("ID=(.+?);", features[0][-1]).group(1)
                feature_type = "gene"  # don't want ORF

                if "pseudo" not in notes:  # catch 'pseudogene=TYPE' or 'pseudo=True'
                    notes["pseudo"] = True if "pseudo" in features[0][-1] else False

            if "tRNA" not in notes and feature_type == "tRNA":
                notes["tRNA"] = re.search("product=(.+?)[;$]?").group(1)

            parts[feature_type] = Location.from_gff3_features(features, start_col=2)

            if feature_type == "CDS" and self.fasta:
                translation = Seq.translate(
                    parts["CDS"].extract(self.fasta.get_sequence(scaffold))
                )

            if identifiers["locus_tag"] and identifiers["protein_id"]:
                continue

            if jgi_portal:  # JGI files
                gene_number = gene_id.split("_")[1]  # e.g. gene_1
                protein_id = re.search("proteinId=(.+?);", features[0][-1]).group(1)
                identifiers["locus_tag"] = f"{jgi_portal}_{gene_number}"
                identifiers["protein_id"] = f"{jgi_portal}_{protein_id}"
                continue

            if "Sequence:" in gene_id:  # AspGD/gbrowser files
                gene_id = gene_id.split(":")[1]
                identifiers["locus_tag"] = gene_id
                identifiers["protein_id"] = gene_id  # only single ID
                continue

            for key, value in identifiers.items():  # The remainder
                if not value and key in features[-1]:
                    identifiers[key] = re.search(f"{key}=(.+?);", features[-1]).group(1)

        if not identifiers["locus_tag"] or not identifiers["protein_id"]:
            identifiers["locus_tag"] = gene_id  # if all else fails...
            identifiers["protein_id"] = gene_id

        return self.Feature(identifiers, parts, translation, notes)

    def add_fasta(self, fasta):
        """ Add records in a FASTA parser object to this object."""
        self.set_fasta_parser(fasta)

        if not set(self.records).issubset(self.fasta):
            self.fasta = None
            raise ValueError("Records in this object do not match the FASTA")

        # Using the FASTA records dictionary as a base, add features
        _records = fasta.records
        for scaffold, record in _records.items():
            if scaffold not in self.records:
                continue

            for i, feature in enumerate(self.records[scaffold].features):
                if "CDS" not in feature.parts:
                    continue

                self.records[scaffold].features[i]._replace(
                    translation=Seq.translate(
                        feature.parts["CDS"].extract(record.sequence)
                    )
                )

            # Replace Record on FASTA parser with features from GFF3 parser
            _records[scaffold]._replace(features=self.records[scaffold].features)

        # Switcheroo!
        self.records = _records


class FASTA(Parser):
    """ Parse a FASTA format file."""

    def _parse_file(self):
        """ Parse an open FASTA file handle."""
        header, sequence = "", ""
        with open(self.file) as handle:
            for line in handle:
                line = line.strip()
                if line.startswith(">"):
                    if header and sequence:
                        self.records[header] = self.Record(header, sequence, None)
                        sequence = ""
                    header = line.replace(">", "")
                else:
                    sequence += line

            # So we don't miss the last one
            if header and sequence:
                self.records[header] = self.Record(header, sequence, None)


def find_identifier(attributes, fields):
    """ Try to find an identifier given a list of potential fields.
        List given should be in order of preference; will return first
        non-None value that can be extracted. e.g.:
            Get protein ID from 'protein_id', then 'GenBank', then 'Name'..
    """
    for field in fields:
        try:
            value = attributes[field]
        except KeyError:
            pass
        else:
            return value


def extract(dictionary, field):
    """ Simple function wrapping try/except around dictionary query.
    """
    try:
        value = dictionary[field]
    except KeyError:
        return None

    if not value:
        return None
    if isinstance(value, (list, tuple)):
        return value[0]
    return value
