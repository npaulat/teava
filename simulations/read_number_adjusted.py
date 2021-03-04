from Bio import SeqIO
import math
import argparse

parser = argparse.ArgumentParser(description="Arguments for de program")
parser.add_argument("sequence", metavar="sequence", type=argparse.FileType('r'), nargs="?", help="Sequence with Insertions")
parser.add_argument("coverage", metavar="coverage", type=int, nargs='?', help="Coverage, please enter withouth 'X' ")
args = parser.parse_args()

coverage = args.coverage
sequence = SeqIO.read(args.sequence, "fasta")

read_number = len(sequence)*coverage / 250

print(int(read_number))
