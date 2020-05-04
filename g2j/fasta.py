"""FASTA parser."""


def iter_records(handle):
    """Parse FASTA file."""

    header = ""
    sequence = ""

    for line in handle:
        line = line.strip()
        if line.startswith(">"):
            if header and sequence:
                yield header, sequence
                sequence = ""
            header = line.strip().replace(">", "")
        else:
            sequence += line.strip()
    yield header, sequence


def parse(handle):
    return {header: sequence for header, sequence in iter_records(handle)}
