'''
Generate semi-synthetic ground truth by shuffling a fasta file and then reinserting specific elements
Essentially, removes biological meaning from everything but desired elements
'''

import argparse
import os

parser = argparse.ArgumentParser(description='Generate semi-synthetic ground truth by shuffling a fasta file and then reinserting specific elements')
parser.add_argument('--fasta', help='Input fasta file')
parser.add_argument('--bed', help='Bed file with elements to reinsert')
parser.add_argument('--output', help='Output fasta file')

#@#@#@#@#@#@#@#@#@#@#@#@#@#@
# Getting argments and validating
#@#@#@#@#@#@#@#@#@#@#@#@#@#@

args = parser.parse_args()
fasta_file = args.fasta
output_file = args.output
bed_file = args.bed

if not os.path.isfile(fasta_file):
    print("Fasta file does not exist")
    exit(1)

if not os.path.isfile(bed_file):
    print("Bed file does not exist")
    exit(1)

#@#@#@#@#@#@#@#@#@#@#@#@#@#@
# Importing modules and setting up
#@#@#@#@#@#@#@#@#@#@#@#@#@#@

from Bio import SeqIO
import random
import re

#@#@#@#@#@#@#@#@#@#@#@#@#@#@
# Classes
#@#@#@#@#@#@#@#@#@#@#@#@#@#@

class Element:
    def __init__(self, chrom, start, end, seq):
        self.chrom = chrom
        self.start = start
        self.end = end
        self.seq = seq

#@#@#@#@#@#@#@#@#@#@#@#@#@#@
# Functions
#@#@#@#@#@#@#@#@#@#@#@#@#@#@
        
def random_nucleotide(match):
    random_value = random.choice(['A', 'T', 'C', 'G'])
    return random_value

#@#@#@#@#@#@#@#@#@#@#@#@#@#@
# Reading elements from bed file and fasta file
#@#@#@#@#@#@#@#@#@#@#@#@#@#@
        
# get sequence from fasta file
seq = str(SeqIO.read(fasta_file, "fasta").seq)

# get elements from bed file
element_list = []
with open(bed_file, "r") as f:
    f.readline() # Skip header
    line = f.readline()
    while line:
        line = line.strip()
        fields = line.split('\t')

        chrom = fields[0]
        start = int(fields[1])
        end = int(fields[2])
        element_seq = seq[start:end]

        element = Element(fields[0], start, end, element_seq)
        element_list.append(element)

        line = f.readline()

#@#@#@#@#@#@#@#@#@#@#@#@#@#@
# Shuffling sequence
#@#@#@#@#@#@#@#@#@#@#@#@#@#@
seq = list(seq)

random.shuffle(seq)
seq = ''.join(seq)

seq = re.sub('N', random_nucleotide, seq)


#@#@#@#@#@#@#@#@#@#@#@#@#@#@
# Insert elements
#@#@#@#@#@#@#@#@#@#@#@#@#@#@
biggest_end = -1
element_list.sort(key=lambda x: x.start)

prev_end = 0
l = []
for i in range(len(element_list)):
    if element_list[i].start > biggest_end:
        l.append(seq[prev_end:element_list[i].start])
        l.append(element_list[i].seq)
        prev_end = element_list[i].end
        biggest_end = element_list[i].end
    
    elif element_list[i].end > biggest_end:
        seq_start = biggest_end - element_list[i].start
        l.append(element_list[i].seq[seq_start:])
        prev_end = element_list[i].end
        biggest_end = element_list[i].end

l.append(seq[prev_end:])
seq = ''.join(l)



#@#@#@#@#@#@#@#@#@#@#@#@#@#@
# Writing output
#@#@#@#@#@#@#@#@#@#@#@#@#@#@

with open(output_file, "w") as f:
    f.write(f">{chrom}\n")
    f.write(seq)
