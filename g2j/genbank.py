"""GenBank parser."""


import re

from pathlib import Path

from g2j.classes import Organism, Scaffold, Feature, Location


PATTERNS = {
    "scaffold": re.compile(r"LOCUS\s+?(?P<accession>\b[\w.-]+?)\s.+?//", re.DOTALL),
    "organism": re.compile(r"ORGANISM\s+?(?P<organism>\w[\w .-]+?)[\n\r]"),
    "strain": re.compile(r'/strain="(?P<strain>([\w .-]+?))"'),
    "features": re.compile(
        r"FEATURES {13}Location/Qualifiers(.+?)(?:ORIGIN|CONTIG)", re.DOTALL
    ),
    "sequence": re.compile(r"ORIGIN\s+?([^\s].*)", re.DOTALL),
    "identifier": re.compile(r'(protein_id|locus_tag|gene|ID)="([\w.:-]+?)"'),
}


def parse_location(string):
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
        location.start.append(part[0] - 1)
        location.end.append(part[0] if len(part) == 1 else part[1])

    return location


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
    # RegEx pattern to extract qualifiers, e.g. /locus_tag="GENE_0001"
    pattern = re.compile(r'/(\w+?)="(.+?)"')

    feature = None
    features = []

    for line in text.split("\n"):
        if not line or line.isspace():
            continue
        if not line.startswith("                     "):
            if feature:
                features.append(feature)
            header, rest = line.strip().split()
            feature = Feature(header, qualifiers=rest)
        else:
            feature.qualifiers += line.strip()
    if feature:
        features.append(feature)

    # Separate location and qualifiers (if any)
    # e.g. 1..100/qualifier1="value"/qualifier2="value"
    for feature in features:
        parts = feature.qualifiers.split("/", 1)
        if len(parts) > 1:
            feature.qualifiers = {
                key: value for key, value in pattern.findall(f"/{parts[1]}")
            }
        feature.location = parse_location(parts[0])

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
