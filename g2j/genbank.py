"""GenBank parser."""


import re

from collections import defaultdict
from pathlib import Path

from g2j.classes import Organism, Scaffold, Feature, Location, Interval


PATTERNS = {
    "scaffold": re.compile(r"LOCUS\s+?(?P<accession>\b[\w.-]+?)\s.+?\n//", re.DOTALL),
    "organism": re.compile(r"ORGANISM\s+?(?P<organism>\w[\w .-]+?)[\n\r]"),
    "strain": re.compile(r'/strain="(?P<strain>([\w .-]+?))"'),
    "features": re.compile(
        r"FEATURES {13}Location/Qualifiers(.+?)(?:ORIGIN|CONTIG)", re.DOTALL
    ),
    "sequence": re.compile(r"ORIGIN\s+?([^\s].*)", re.DOTALL),
    "identifier": re.compile(r'(protein_id|locus_tag|gene|ID)="([\w.:-]+?)"'),
    "qualifier": re.compile(r'/(\w+?)=(?:"(.+?)"|(\d+?))', re.DOTALL),
}


def parse_location(string, codon_start):
    """Parse GenBank location string."""

    location = Location(
        five_is_partial=True if ">" in string else False,
        three_is_partial=True if "<" in string else False,
    )

    if "complement" in string:
        location.flip()

    for chars in ["complement", "order", "join", "(", ")", "<", ">"]:
        string = string.replace(chars, "")

    for span in string.split(","):
        part = [int(x) for x in span.split("..")]
        start = part[0] - 1
        end = part[0] if len(part) == 1 else part[1]
        location.intervals.append(Interval(start, end))

    set_CDS_interval_phases(location, codon_start)

    return location


def compute_phase(start, end, offset=0):
    leftover = abs(end - start + offset) % 3
    return 0 if leftover == 0 else 2 - leftover


def set_CDS_interval_phases(location, codon_start):
    """Compute phase (a la GFF3 files) of a feature location."""

    # Initial CDS phase is just codon_start - 1 (zero index)
    phases = [codon_start - 1]

    # Compute phases in staggered fashion
    for interval in location.intervals:
        interval.phase = phases[-1]
        phase = compute_phase(interval.start, interval.end, interval.phase)
        phases.append(phase)


def parse_qualifier_block(text):
    """Parse qualifiers from qualifier block.


    Qualifiers are split by newline -> 21 spaces -> slash
    If value, Remove leading/trailing ", leading newline -> 21 spaces
    Otherwise it's boolean, e.g. /pseudo, so set True

    Store qualifiers of same type in lists, otherwise just str
    """
    qualifiers = {}
    gap = "\n" + " " * 21

    for qualifier in text.split(f"{gap}/"):
        key, *value = qualifier.split("=")

        # Determine if boolean or str value
        if value:
            value = value[0].lstrip('"').strip('"\n').replace(gap, "")
        else:
            value = True

        # Save multiple qualifiers of same type in lists
        if key in qualifiers:
            if isinstance(qualifiers[key], list):
                qualifiers[key].append(value)
            else:
                qualifiers[key] = [qualifiers.pop(key), value]
        else:
            qualifiers[key] = value

    return qualifiers


def parse_feature_block(text):
    """Parse features from feature block.

    Block should resemble:
        ^     feature         1..100
        ^                     /qualifier1="value1"
        ^                     /qualifier2="value2"
        ^     feature2        200..500
        ...

    as matched using PATTERNS["features"].

    Separate feature blocks are identified by looking for differences in whitespace;
    new feature lines start with 5 spaces, qualifiers with 21.
    """
    pattern = re.compile(r"^ {5}([\w0-9']+?)\s+([^/]+)", re.M | re.DOTALL)
    features = []
    text_length = len(text)

    # Get positions of each unique sequence feature
    for match in pattern.finditer(text):
        feature = {
            "type": match.group(1),
            "location": match.group(2).replace(" ", "").replace("\n", ""),
            "interval": [match.end(), text_length],
        }
        if features:
            features[-1]["interval"][1] = match.start()

        features.append(feature)

    # Convert above to Feature objects
    for ix, feature in enumerate(features):

        # Get qualifier text corresponding to the current feature
        # The start+1 gets rid of the / in the first qualifier
        start, end = feature["interval"]
        qualifiers = parse_qualifier_block(text[start + 1 : end])

        try:
            codon_start = int(qualifiers["codon_start"])
        except KeyError:
            codon_start = 1

        features[ix] = Feature(
            type=feature["type"],
            location=parse_location(feature["location"], codon_start),
            qualifiers=qualifiers,
        )

    return features


def scaffold_iter(text):
    """Wrapper around scaffold RegEx pattern."""
    yield from PATTERNS["scaffold"].finditer(text)


def get_feature_block(text):
    """Wrapper around scaffold RegEx pattern."""
    match = PATTERNS["features"].search(text)
    if match:
        return match.group(1)


def get_scaffold_sequence(text):
    """Extract scaffold sequence, stripping any digits/whitespace."""
    match = PATTERNS["sequence"].search(text)
    if match:
        return re.sub(r"[^A-Za-z]", "", match.group(1))


def find_pattern(pattern, text):
    """Find a match in a text block for a pattern in PATTERNS.

    If no match is found, as indicated by either IndexError or AttributeError
    being raised, this function will return None.
    """
    if pattern not in PATTERNS:
        raise ValueError("Invalid pattern specified")
    try:
        return PATTERNS[pattern].search(text).groups()[0]
    except (IndexError, AttributeError):
        return None


def parse(handle):
    """Parse a GenBank file."""

    organism = Organism()

    for scaffold in scaffold_iter(handle.read()):

        text = scaffold.group(0)
        accession = scaffold.group("accession")

        if not organism.name:
            organism.name = find_pattern("organism", text)

        if not organism.strain:
            organism.strain = find_pattern("strain", text)

        feature_block = get_feature_block(text)

        if feature_block:
            scaffold = Scaffold(
                accession,
                get_scaffold_sequence(text),
                parse_feature_block(feature_block),
            )
            organism.scaffolds.append(scaffold)

    return organism


def from_path(path):
    with open(path) as handle:
        organism = parse(handle)
    return organism


def get_genbank_paths(folder):
    """Generate a collection of paths to GenBank files in a specified folder."""
    if not Path(folder).is_dir():
        raise ValueError("Expected valid folder")
    valid_extensions = (".gb", ".gbk", ".genbank")
    return [
        file for file in Path(folder).iterdir() if str(file).endswith(valid_extensions)
    ]


if __name__ == "__main__":
    import sys

    with open(sys.argv[1]) as gbk:
        print(f"Parsing: {gbk.name}")
        organism = parse(gbk)

    with open("test.json", "w") as fp:
        organism.to_json(fp=fp, indent=2)
