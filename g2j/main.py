"""Entry point."""

import argparse

from g2j import genbank, gff3, __version__


def parse(
    genbank_handle=None,
    gff3_handle=None,
    fasta_handle=None,
    output_handle=None,
    json_indent=2,
    grouped=False,
):
    if genbank_handle and not (gff3_handle or fasta_handle):
        print("Parsing GenBank file:")
        print(f"  {genbank_handle.name}")
        organism = genbank.parse(genbank_handle)

    elif gff3_handle and fasta_handle:
        print("Parsing GFF3 and FASTA files:")
        print(f"  {gff3_handle.name}")
        print(f"  {fasta_handle.name}")
        organism = gff3.parse(gff3_handle, fasta_handle)

    else:
        raise ValueError("Expected GenBank or GFF3+FASTA")

    if grouped:
        print("\nGrouping overlapping sequence features")
        organism.group()

    if output_handle:
        print("\nWriting JSON: ")
        print(f"  {output_handle.name}")
        organism.to_json(fp=output_handle, indent=json_indent)

    return organism


def run():
    args = get_arguments()
    parse(
        genbank_handle=args.genbank,
        gff3_handle=args.general,
        fasta_handle=args.fasta,
        output_handle=args.output,
        json_indent=args.json_indent,
        grouped=args.grouped,
    )


def get_arguments():
    parser = argparse.ArgumentParser(
        "Genome2JSON",
        description="Parse genomes in GenBank/GFF3 format, and convert to JSON",
        epilog="Cameron Gilchrist 2020",
    )
    parser.add_argument(
        "--version", action="version", version="%(prog)s " + __version__
    )

    inputs = parser.add_mutually_exclusive_group(required=True)
    inputs.add_argument(
        "-gbk", "--genbank", help="GenBank file", type=argparse.FileType("r")
    )
    inputs.add_argument(
        "-gff", "--general", help="GFF3 file", type=argparse.FileType("r")
    )

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
        "--json_indent", help="Number of spaces to indent in JSON", type=int, default=2
    )
    parser.add_argument(
        "--grouped", help="Group overlapping sequence features", action="store_true",
    )

    args = parser.parse_args()

    if args.general and not args.fasta:
        raise ValueError("FASTA file must be provided when parsing GFF3")

    return args


if __name__ == "__main__":
    run()
