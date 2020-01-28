# genome2json
Parse genomes in GenBank/GFF3 format to JSON

Implemented in pure Python using builtin libraries.


```
usage: Genome2JSON [-h] (-gbk GENBANK | -gff GENERAL) [-fa FASTA] -o OUTPUT
                   [--json_indent JSON_INDENT] [--grouped]

Parse genomes in GenBank/GFF3 format, and convert to JSON

optional arguments:
  -h, --help            show this help message and exit
  -gbk GENBANK, --genbank GENBANK
                        GenBank file
  -gff GENERAL, --general GENERAL
                        GFF3 file
  -fa FASTA, --fasta FASTA
                        FASTA file, required when parsing GFF
  -o OUTPUT, --output OUTPUT
                        Save JSON to file
  --json_indent JSON_INDENT
                        Number of spaces to indent in JSON
  --grouped             Group overlapping sequence features

Cameron Gilchrist 2020
```
