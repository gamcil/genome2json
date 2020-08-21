import argparse
import logging

from g2j import genbank, gff3, __version__


logging.basicConfig(
    format="[%(asctime)s] %(levelname)s - %(message)s",
    datefmt="%H:%M:%S"
)
LOG = logging.getLogger("g2j")
LOG.setLevel(logging.INFO)


def parse(
    genbank_handle=None,
    gff3_handle=None,
    fasta_handle=None,
    output_handle=None,
    indent=2,
    collapse=False,
    group=False,
    translate=False,
    feature_types=None,
    save_scaffold_sequence=True,
):
    if genbank_handle and not (gff3_handle or fasta_handle):
        LOG.info(f"Parsing GenBank file: {genbank_handle.name}")
        organism = genbank.parse(
            genbank_handle,
            feature_types=feature_types,
            save_scaffold_sequence=save_scaffold_sequence,
        )
    elif gff3_handle:
        LOG.info("Parsing files...")
        LOG.info(f"  GFF: {gff3_handle.name}")
        if fasta_handle:
            LOG.info(f"  FASTA: {fasta_handle.name}")
        organism = gff3.parse(
            gff3_handle,
            fasta_handle=fasta_handle,
            feature_types=feature_types,
            save_scaffold_sequence=save_scaffold_sequence,
        )
    else:
        raise ValueError("Expected GenBank or GFF3+FASTA")

    if collapse:
        LOG.info("Collapsing same type sequence features")
        organism.collapse()

        if translate:
            LOG.info("Translating CDS features")
            organism.translate()

    if group:
        LOG.info("Grouping overlapping sequence features")
        organism.group()

    LOG.info("Sorting sequence features")
    organism.sort()

    if output_handle:
        LOG.info(f"Writing JSON: {output_handle.name}")
        organism.to_json(fp=output_handle, indent=indent)

    return organism


def run():
    args = get_arguments()
    LOG.info("Starting genome2json")
    parse(
        genbank_handle=args.genbank,
        gff3_handle=args.general,
        fasta_handle=args.fasta,
        output_handle=args.output,
        indent=args.indent,
        collapse=args.collapse,
        group=args.group,
        translate=args.translate,
        feature_types=args.feature_types,
    )
    LOG.info("Done!")


def get_arguments():
    parser = argparse.ArgumentParser(
        "Genome2JSON",
        description="Parse genomes in GenBank/GFF3 format, and convert to JSON",
        epilog="Cameron Gilchrist 2020",
    )
    parser.add_argument("--version", action="version", version="%(prog)s " + __version__)

    inputs = parser.add_mutually_exclusive_group(required=True)
    inputs.add_argument("-gbk", "--genbank", help="GenBank file", type=argparse.FileType("r"))
    inputs.add_argument("-gff", "--general", help="GFF3 file", type=argparse.FileType("r"))

    parser.add_argument(
        "-fa",
        "--fasta",
        help="FASTA file, required when parsing GFF",
        type=argparse.FileType("r"),
    )
    parser.add_argument(
        "-o",
        "--output",
        required=True,
        help="Save JSON to file",
        type=argparse.FileType("w"),
    )
    parser.add_argument(
        "-ft",
        "--feature_types",
        help="Feature types to save",
        nargs="+",
    )
    parser.add_argument(
        "-i",
        "--indent",
        help="Number of spaces to indent in JSON",
        type=int,
        default=2,
    )
    parser.add_argument(
        "-c",
        "--collapse",
        help="Collapse same type sequence features (e.g. exons, CDS)",
        action="store_true"
    )
    parser.add_argument(
        "-g",
        "--group",
        help="Group overlapping sequence features",
        action="store_true",
    )
    parser.add_argument(
        "-t",
        "--translate",
        help="Translate CDS sequence features",
        action="store_true",
    )

    args = parser.parse_args()

    return args


if __name__ == "__main__":
    run()
