"""Data structures used by parsers."""


import json

from functools import partialmethod

from g2j import grouping


translation_table = (
    "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    "---M------**--*----M---------------M----------------------------",
    "TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG",
    "TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG",
    "TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG",
)


def get_table():
    """Load NCBI codon translation table."""
    table = {"starts": [], "stops": [], "codons": {}}
    for amino, start, *bases in zip(*translation_table):
        codon = "".join(bases)
        table["codons"][codon] = amino
        if start == "M":
            table["starts"].append(codon)
        elif start == "*":
            table["stops"].append(codon)
    return table


def reverse_complement(sequence):
    table = str.maketrans("acgtACGT", "tgcaTGCA")
    return sequence.translate(table)[::-1]


def translate(sequence, codon_start=0, partial=False, to_stop=False):
    """Translate a given nucleotide sequence."""

    table = get_table()
    length = len(sequence)
    sequence = sequence.upper()

    translation = ""
    for i in range(0 + codon_start, length, 3):
        # If partial=True, finish translation at previous codon
        # If to_stop=True, include * in translation
        # Account for ambiguous bases, translate to 'X'
        codon = sequence[i:i+3]
        if len(codon) != 3:
            if not partial:
                raise ValueError("Reached partial codon and partial=False")
            break
        if codon in table["stops"] and not to_stop:
            break
        if "N" in codon or "X" in codon:
            translation += "X"
        else:
            translation += table["codons"][codon]
    return translation


def merge_qualifiers(a, b):
    """Creates a new dictionary from two merged qualifier dictionaries."""
    merged = a.copy()
    for key, value in b.items():
        if key in merged:
            if merged[key] == value:
                continue
            if isinstance(merged[key], list):
                if isinstance(value, list):
                    merged[key].extend(value)
                else:
                    merged[key].append(value)
            else:
                if isinstance(value, list):
                    merged[key] = [merged[key], *value]
                else:
                    merged[key] = [merged[key], value]
        else:
            merged[key] = value
    return merged


class Serialiser:
    """Mixin for JSON serialisation."""

    def to_dict(self):
        raise NotImplementedError

    @classmethod
    def from_dict(cls, d):
        raise NotImplementedError

    def to_json(self, fp=None, **kwargs):
        if fp:
            json.dump(self.to_dict(), fp, **kwargs)
            return
        return json.dumps(self.to_dict(), **kwargs)

    @classmethod
    def from_json(cls, fp):
        return cls.from_dict(json.load(fp))


class Organism(Serialiser):
    __slots__ = ("name", "strain", "scaffolds")

    def __init__(self, name="", strain="", qualifiers=None, scaffolds=None):
        self.name = name
        self.strain = strain
        self.qualifiers = qualifiers if qualifiers else {}
        self.scaffolds = scaffolds if scaffolds else []

    def __str__(self):
        return f"{self.name} {self.strain} [{len(self.scaffolds)} scaffolds]"

    def _action(self, method):
        for scaffold in self.scaffolds:
            getattr(scaffold, method)()

    sort = partialmethod(_action, "sort")
    group = partialmethod(_action, "group")
    ungroup = partialmethod(_action, "ungroup")
    collapse = partialmethod(_action, "collapse")
    uncollapse = partialmethod(_action, "uncollapse")
    translate = partialmethod(_action, "translate")

    def to_dict(self):
        return {
            "name": self.name,
            "strain": self.strain,
            "qualifiers": self.qualifiers,
            "scaffolds": [scaffold.to_dict() for scaffold in self.scaffolds],
        }

    @classmethod
    def from_dict(cls, d):
        return cls(
            name=d["name"],
            strain=d["strain"],
            qualifiers=d["qualifiers"],
            scaffolds=[Scaffold.from_dict(scaffold) for scaffold in d["scaffolds"]],
        )


class Scaffold(Serialiser):
    __slots__ = ("accession", "sequence", "features")

    def __init__(self, accession, sequence=None, features=None):
        self.accession = accession
        self.sequence = sequence
        self.features = features if features else []
        self._grouped = False

    def __str__(self):
        return f"{self.accession} [{len(self.features)} features]"

    def __len__(self):
        return len(self.sequence)

    def sort(self, reverse=False):
        """Sorts Features on scaffold by start."""
        def func(g):
            if isinstance(g, Feature):
                return g.location.min()
            g.sort(key=lambda f: f.location.min(), reverse=reverse)
            return min(f.location.min() for f in g)
        self.features.sort(key=func, reverse=reverse)

    def group(self):
        """Group sequence features by overlap/shared attributes."""
        if not self._grouped:
            self.features = grouping.group_features(self.features)
            self._grouped = True
        else:
            raise ValueError("Features are already grouped")

    def ungroup(self):
        """Flatten grouped sequence features."""
        if not self._grouped:
            raise ValueError("Features are not grouped")
        features = []
        for feature in self.features:
            if isinstance(feature, list):
                features.extend(feature)
            else:
                features.append(feature)
        self.features = features
        self._grouped = False

    def collapse(self):
        """Collapses same type feature groups (e.g. exons, CDS)."""
        if self._grouped:
            self.features = [
                grouping.collapse_features(group)
                for group in self.features
            ]
        else:
            self.features = grouping.collapse_features(self.features)

    def uncollapse(self):
        """Uncollapses same type feature groups by intervals (e.g. exons, CDS)."""
        self.features = [
            feature.uncollapse()
            if isinstance(feature, Feature)
            else [f.uncollapse() for f in feature]
            for feature in self.features
            if len(feature.location.intervals) > 1
        ]

    def translate(self, partial=True, to_stop=False):
        """Translates all features on this scaffold."""
        for feature in self.features:
            if isinstance(feature, Feature):
                feature.translate(self.sequence, partial=partial, to_stop=to_stop)
            else:
                for f in feature:
                    f.translate(self.sequence, partial=partial, to_stop=to_stop)

    def to_dict(self):
        features = []
        for feature in self.features:
            if isinstance(feature, Feature):
                features.append(feature.to_dict())
            else:
                features.append([f.to_dict() for f in feature])
        return {
            "accession": self.accession,
            "sequence": self.sequence,
            "features": features,
        }

    @classmethod
    def from_dict(cls, d):
        features = []
        for feature in d["features"]:
            if isinstance(feature, dict):
                features.append(Feature.from_dict(feature))
            else:
                features.append([Feature.from_dict(f) for f in feature])
        return cls(d["accession"], d["sequence"], features)


class Feature(Serialiser):
    __slots__ = ("type", "qualifiers", "location")

    def __init__(self, type, qualifiers=None, location=None):
        self.type = type
        self.qualifiers = qualifiers if qualifiers else {}
        self.location = location if location else Location()

    def __str__(self):
        return f"{self.type} {str(self.location)}"

    def __iadd__(self, other):
        self.location += other.location
        self.qualifiers = merge_qualifiers(self.qualifiers, other.qualifiers)
        return self

    def extend(self, others):
        if any(not isinstance(f, Feature) for f in others):
            raise NotImplementedError("Expected Feature object")
        for feature in others:
            self += feature

    def is_collapsed(self):
        return len(self.location.intervals) > 1

    def uncollapse(self):
        """Splits collapsed Features into multiple Features."""
        if not self.is_collapsed():
            raise ValueError("This feature is not collapsed")
        return [
            Feature(
                type=self.type,
                qualifiers=self.qualifiers,
                location=Location(
                    intervals=[interval],
                    strand=self.location.strand,
                    five_is_partial=self.location.five_is_partial,
                    three_is_partial=self.location.three_is_partial,
                )
            )
            for interval in self.location.intervals
        ]

    def translate(self, sequence, partial=True, to_stop=False):
        """Translates parent sequence given this features location."""
        if not self.type == "CDS":
            return
        nts = self.location.slice(sequence)
        if self.location.strand == "-":
            nts = reverse_complement(nts)
        self.qualifiers["translation"] = translate(nts, partial=partial, to_stop=to_stop)

    def is_genic(self):
        return self.type in {
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
            "start_codon",
            "stop_codon",
            "three_prime_UTR",
            "five_prime_UTR",
        }

    def to_dict(self):
        return {
            "type": self.type,
            "qualifiers": self.qualifiers,
            "location": self.location.to_dict(),
        }

    @classmethod
    def from_dict(cls, d):
        return cls(
            d["type"],
            qualifiers=d["qualifiers"],
            location=Location.from_dict(d["location"]),
        )


class Location(Serialiser):
    __slots__ = ("intervals", "strand", "five_is_partial", "three_is_partial")

    def __init__(
        self, intervals=None, strand="+", five_is_partial=False, three_is_partial=False,
    ):
        self.intervals = intervals if intervals else []
        self.strand = strand
        self.five_is_partial = five_is_partial
        self.three_is_partial = three_is_partial

    def __str__(self):
        loc = ",".join(f"{start}..{end}" for start, end in self.iter_intervals())
        return f"{loc}[{self.strand}]"

    def __len__(self):
        return self.max() - self.min()

    def __iadd__(self, other):
        self.intervals.extend(other.intervals)
        self.five_is_partial |= other.five_is_partial
        self.three_is_partial |= other.three_is_partial
        return self

    def iter_intervals(self):
        for interval in sorted(self.intervals, key=lambda i: i.start):
            yield interval.start, interval.end

    def slice(self, sequence):
        return "".join(sequence[start:end] for start, end in self.iter_intervals())

    def starts(self):
        return [i.start for i in self.intervals]

    def ends(self):
        return [i.end for i in self.intervals]

    def min(self):
        return min(self.starts())

    def max(self):
        return max(self.ends())

    def flip(self):
        """Flip strand and five/three end flags."""
        self.strand = "-" if self.strand == "+" else "+"
        self.five_is_partial, self.three_is_partial = (
            self.three_is_partial,
            self.five_is_partial,
        )

    def to_dict(self):
        return {
            "intervals": [i.to_dict() for i in self.intervals],
            "strand": self.strand,
            "five_is_partial": self.five_is_partial,
            "three_is_partial": self.three_is_partial,
        }

    @classmethod
    def from_dict(cls, d):
        return cls(
            intervals=[Interval.from_dict(i) for i in d["intervals"]],
            strand=d["strand"],
            five_is_partial=d["five_is_partial"],
            three_is_partial=d["three_is_partial"],
        )


class Interval(Serialiser):
    def __init__(self, start, end, phase=None):
        self.start = start
        self.end = end
        self.phase = phase

    def __str__(self):
        return f"{self.start}-{self.end} [phase={self.phase}]"

    def to_dict(self):
        return {"start": self.start, "end": self.end, "phase": self.phase}

    @classmethod
    def from_dict(cls, d):
        return cls(start=d["start"], end=d["end"], phase=d["phase"])
