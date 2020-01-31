"""Data structures used by parsers."""


import json

from g2j import grouping


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

    def __init__(self, name="", strain="", scaffolds=None):
        self.name = name
        self.strain = strain
        self.scaffolds = scaffolds if scaffolds else []

    def __str__(self):
        return f"{self.name} {self.strain} [{len(self.scaffolds)} scaffolds]"

    def group(self):
        for scaffold in self.scaffolds:
            scaffold.group()

    def ungroup(self):
        for scaffold in self.scaffolds:
            scaffold.ungroup()

    def to_dict(self):
        return {
            "name": self.name,
            "strain": self.strain,
            "scaffolds": [scaffold.to_dict() for scaffold in self.scaffolds],
        }

    @classmethod
    def from_dict(cls, d):
        return cls(
            name=d["name"],
            strain=d["strain"],
            scaffolds=[Scaffold.from_dict(scaffold) for scaffold in d["scaffolds"]],
        )


class Scaffold(Serialiser):
    __slots__ = ("accession", "sequence", "features")

    def __init__(self, accession, sequence=None, features=None):
        self.accession = accession
        self.sequence = sequence
        self.features = features if features else []

    def __str__(self):
        return f"{self.accession} [{len(self.features)} features]"

    def __len__(self):
        return len(self.sequence)

    def group(self):
        """Group sequence features by overlap/shared attributes."""
        self.features = grouping.group_features(self.features)

    def ungroup(self):
        """Flatten grouped sequence features."""
        features = []
        for feature in self.features:
            if isinstance(feature, list):
                features.extend(feature)
            else:
                features.append(feature)
        self.features = features

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
        self.qualifiers.update(other.qualifiers)
        return self

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
        for interval in self.intervals:
            yield interval.start, interval.end

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
