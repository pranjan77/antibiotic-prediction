#!/usr/bin/env python3

import sys

from Bio import SeqIO

input_file, output_file = sys.argv[1:]

SeqIO.convert(input_file, "genbank", output_file, "fasta")
