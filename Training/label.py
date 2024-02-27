'''
Adds labels to consecutive pairs of stretched denoting whether they belong to the same LTR or not.
'''

import argparse
import os
import gc
import numpy as np

parser = argparse.ArgumentParser(description='Label consecutive pairs of stretches')
parser.add_argument('--forward', required=True, help='File containing forward stretches')
parser.add_argument('--backward', required=True, help='File containing backward stretches')
parser.add_argument('--feature', required=True, help='File containing features matrix')
parser.add_argument('--bed', required=True, help='golden standard bed file')
parser.add_argument('--output', required=True, help='Output file')

#@#@#@#@#@#@#@#@#@#@#@#@#@#@
# Get arguments and validate 
#@#@#@#@#@#@#@#@#@#@#@#@#@#@

args = parser.parse_args()
forward_file = args.forward
backward_file = args.backward
feature_file = args.feature
bed_file = args.bed
output_file = args.output

if not os.path.exists(forward_file):
    print("Forward file {} does not exist".format(forward_file))
    exit(1)

if not os.path.exists(backward_file):
    print("Backward file {} does not exist".format(backward_file))
    exit(1)

if not os.path.exists(feature_file):
    print("Feature file {} does not exist".format(feature_file))
    exit(1)

if not os.path.exists(bed_file):
    print("Bed file {} does not exist".format(bed_file))
    exit(1)

#@#@#@#@#@#@#@#@#@#@#@#@#@#@
# Classes
#@#@#@#@#@#@#@#@#@#@#@#@#@#@

class Stretch:
    def __init__(self, start, end, score):
        self.start = start
        self.end = end
        self.score = score

    def is_overlap(self, other, threshold):
        overlap = min(self.end, other.end) - max(self.start, other.start)
        return overlap / (self.end - self.start) > threshold
    

    
class LTR:
    def __init__(self, start, end):
        self.start = start
        self.end = end


#@#@#@#@#@#@#@#@#@#@#@#@#@#@
# Functions
#@#@#@#@#@#@#@#@#@#@#@#@#@#@

def load_stretches(file):
    stretches = []
    with open(file, 'r') as f:
        for line in f:
            stretches.append(Stretch(*[int(x) for x in line.strip().split(',')]))
    return stretches

def convert_to_pairs(stretches):
    pairs = []
    for i in range(len(stretches) - 1):
        pairs.append((stretches[i], stretches[i + 1]))
    return pairs

def load_ltrs(file):
    ltrs = []
    with open(file, 'r') as f:
        f.readline() # skip header
        line = f.readline()
        while line:
            data = line.strip().split()
            l_ltr = LTR(int(data[3]), int(data[4]))
            r_ltr = LTR(int(data[5]), int(data[6]))

            ltrs.append(l_ltr)
            ltrs.append(r_ltr)

            line = f.readline()

    return ltrs

def label_pairs(pairs, ltrs):
    labels = []
    for pair in pairs:
        ltr_overlap = False
        for ltr in ltrs:
            if pair[0].is_overlap(ltr, 0.5) and pair[1].is_overlap(ltr, 0.5):
                ltr_overlap = True
                break
        if ltr_overlap:
            labels.append(1)
        else:
            labels.append(0)
    return labels

def generate_labels(ltrs, stretch_file):
    stretches = load_stretches(stretch_file)
    pairs = convert_to_pairs(stretches)

    print()

    del stretches
    gc.collect()

    labels = label_pairs(pairs, ltrs)
    return labels


#@#@#@#@#@#@#@#@#@#@#@#@#@#@
#  Generate labels
#@#@#@#@#@#@#@#@#@#@#@#@#@#@

ltrs = load_ltrs(bed_file)

labels = generate_labels(ltrs, forward_file)
labels.extend(generate_labels(ltrs, backward_file))

#@#@#@#@#@#@#@#@#@#@#@#@#@#@
#  Update feature matrix with labels
#@#@#@#@#@#@#@#@#@#@#@#@#@#@
    
feature_table = np.loadtxt(feature_file, delimiter=',')

if len(labels) != feature_table.shape[0]:
    print("Number of labels ({}) does not match number of rows in feature table ({})".format(len(labels), feature_table.shape[0]))
    exit(1)

# Add labels to first column
feature_table = np.insert(feature_table, 0, labels, axis=1)

np.savetxt(output_file, feature_table, delimiter=',')

