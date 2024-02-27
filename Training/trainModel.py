'''
a workflow that generates training data then trains a model for Look4LTRs
'''

import argparse
import os
import graphStatistics as gs

parser = argparse.ArgumentParser(description='Train Look4LTRs detector')
parser.add_argument('-fd', '--fastadirs', nargs='+', required=True, help='FASTA file directories containing sequences to be used for training. Each directory should be its own genome.')
parser.add_argument('-bd', '--beddirs', nargs='+', required=True, help='BED file directories containing annotations to be used for training. Each directory should be its own genome. Should correspond to the FASTA file directories.')
parser.add_argument('-o', '--output', required=True, help='Output directory for all generated training data and model. Will be created if it does not exist.')
parser.add_argument('-ge', '--genexe', required=True, help='Path to the generateTrainingData executable')
parser.add_argument('-gg', '--gengraph', required=True, help='Path to the generateGraphData executable')

#@#@#@#@#@#@#@#@#@#@#@#@#@#@
# Get arguments and validate
#@#@#@#@#@#@#@#@#@#@#@#@#@#@

args = parser.parse_args()
fasta_dirs = args.fastadirs
bed_dirs = args.beddirs
output_dir = args.output
genexe_path = args.genexe
gengraph_path = args.gengraph

assert len(fasta_dirs) == len(bed_dirs), "Number of FASTA directories does not match number of BED directories"
assert os.path.exists(genexe_path), "Path to generateTrainingData executable does not exist"
assert os.path.exists(gengraph_path), "Path to generateGraphData executable does not exist"


for fasta_dir, bed_dir in zip(fasta_dirs, bed_dirs):

    assert os.path.exists(fasta_dir), "FASTA directory {} does not exist".format(fasta_dir)

    assert os.path.exists(bed_dir), "BED directory {} does not exist".format(bed_dir)

    assert len(os.listdir(fasta_dir)) == len(os.listdir(bed_dir)), "Number of files in FASTA directory {} does not match number of files in BED directory {}".format(fasta_dir, bed_dir)


assert os.path.exists(os.path.dirname(output_dir)), "Output directory can not be created as the parent directory {} does not exist".format(os.path.dirname(output_dir))

if not os.path.exists(output_dir):
    os.mkdir(output_dir)
else:
    os.system("rm -rf {}/*" .format(output_dir))



#@#@#@#@#@#@#@#@#@#@#@#@#@#@
# Functions
#@#@#@#@#@#@#@#@#@#@#@#@#@#@
def sort_file_key(file_name):
    return os.path.splitext(file_name)[0]

#@#@#@#@#@#@#@#@#@#@#@#@#@#@
# For each genome, generate training data
#@#@#@#@#@#@#@#@#@#@#@#@#@#@

training_matrix_dir = "{}/TrainingMatrices" .format(output_dir)
os.mkdir(training_matrix_dir)

for i, fasta_dir, bed_dir in zip(range(len(fasta_dirs)), fasta_dirs, bed_dirs):
    # Setting up directory

    output_genome_dir = "{}/Genome_{}" .format(output_dir, i)

    os.mkdir(output_genome_dir)


    
    semi_synthetic_dir = "{}/SemiSynthetic".format(output_genome_dir)
    os.mkdir(semi_synthetic_dir)

    training_data_dir = "{}/TrainingData".format(output_genome_dir)
    os.mkdir(training_data_dir)

    # Setting up files

    fasta_list = os.listdir(fasta_dir)
    bed_list = os.listdir(bed_dir)

    fasta_list.sort(key=sort_file_key)
    bed_list.sort(key=sort_file_key)

    for fasta_file, bed_file in zip(fasta_list, bed_list):
        assert os.path.splitext(fasta_file)[0] == os.path.splitext(bed_file)[0], "FASTA file {} does not match BED file {}".format(fasta_file, bed_file)

        fasta_path = "{}/{}" .format(fasta_dir, fasta_file)
        bed_path = "{}/{}" .format(bed_dir, bed_file)
        semi_synthetic_output = "{}/{}.fa".format(semi_synthetic_dir, os.path.splitext(fasta_file)[0])

        # Generate semi-synthetic genomes; call generatesSemiSynthetic.py 
        os.system("python3 generateSemiSynthetic.py --fasta {} --bed {} --output {}" .format(fasta_path, bed_path, semi_synthetic_output))

    # Generate training data; call generateTrainingData
    os.system("{} -s {} -f {} -o {} -b {} -pa 1" .format(genexe_path, semi_synthetic_dir, fasta_dir, training_data_dir, bed_dir))

    for fasta_file, bed_file in zip(fasta_list, bed_list):
        bed_path = "{}/{}" .format(bed_dir, bed_file)

        fasta_name = os.path.splitext(fasta_file)[0]

        feature_file = "{}/{}".format(training_data_dir, fasta_name + ".fmx")
        forward_file = "{}/{}".format(training_data_dir, fasta_name + ".fs")
        backward_file = "{}/{}".format(training_data_dir, fasta_name + ".bs")
        output_file = "{}/{}.fmx".format(training_matrix_dir, fasta_name)

        os.system("python3 label.py --forward {} --backward {} --feature {} --bed {} --output {}" .format(forward_file, backward_file, feature_file, bed_path, output_file))
   
   
#@#@#@#@#@#@#@#@#@#@#@#@#@#@
# Training Detector module
#@#@#@#@#@#@#@#@#@#@#@#@#@#@
     
print("Training model...")
feature_files = [os.path.join(training_matrix_dir, x) for x in os.listdir(training_matrix_dir)]
output_file = "{}/config.txt".format(output_dir)
os.system("python3 detectorTrainer.py --feature_file {} --output {}".format(' '.join(feature_files), output_file))

#@#@#@#@#@#@#@#@#@#@#@#@#@#@
# Training Matcher module
#@#@#@#@#@#@#@#@#@#@#@#@#@#@


# getting interior red scores and finding 95th percentile
training_data_file_list = []
for i, fasta_dir, bed_dir in zip(range(len(fasta_dirs)), fasta_dirs, bed_dirs):
    output_genome_dir = "{}/Genome_{}" .format(output_dir, i)
    training_data_dir = "{}/TrainingData".format(output_genome_dir)

    l = [os.path.join(training_data_dir, x) for x in os.listdir(training_data_dir) if x.endswith('.rsc')]

    training_data_file_list.extend(l)

red_scores = []
for red_score_file in training_data_file_list:
    with open(red_score_file, 'r') as f:
        line = f.readline()
        while line:
            data = line.split()
            scores = [int(x) for x in data[1:]]

            red_scores.append(scores)

            line = f.readline()

# for each list of scores, calculate ratio of number of non-zero scores over number of total scores
ratios = []
for scores in red_scores:
    non_zero = 0
    total = 0
    for score in scores:
        if score != 0:
            non_zero += 1
        total += 1

    ratios.append(non_zero / total)

# sort ratios and find 98th percentile
ratios.sort()
index = int(len(ratios) * 0.02) - 1
red_threshold = ratios[index]


# building graphs
connection_list = []
for i, fasta_dir, bed_dir in zip(range(len(fasta_dirs)), fasta_dirs, bed_dirs):
    output_genome_dir = "{}/Genome_{}" .format(output_dir, i)

    graph_dir = "{}/Graphs".format(output_genome_dir)
    os.mkdir(graph_dir)

    os.system("{} -f {} -c {} -o {} -pa 1" .format(gengraph_path, fasta_dir, output_file, graph_dir))

    bed_list = os.listdir(bed_dir)
    graph_list = os.listdir(graph_dir)

    bed_list.sort(key=sort_file_key)
    graph_list.sort(key=sort_file_key)

    for bed_file, graph_file in zip(bed_list, graph_list):

        bed_path = "{}/{}" .format(bed_dir, bed_file)
        graph_path = "{}/{}" .format(graph_dir, graph_file)

        rt_list = gs.load_gsbed(bed_path)
        graph = gs.load_graph(graph_path)
        ltr_list = []
        for rt in rt_list:
            ltr_list.append(rt.get_left())
            ltr_list.append(rt.get_right())

        # sort graph by direction and start position
        forward_graph = sorted([ele for ele in graph if ele.get_direction() == "+"], key = lambda x : x.get_start())

        backward_graph = sorted([ele for ele in graph if ele.get_direction() == "-"], key = lambda x : x.get_start())

        forward_overlaps = gs.match_overlaps(forward_graph, ltr_list)
        backward_overlaps = gs.match_overlaps(backward_graph, ltr_list)

        for rt in rt_list:
            left = rt.get_left()
            right = rt.get_right()
            if left in forward_overlaps and right in backward_overlaps:
                left_ele = forward_overlaps[left]
                right_ele = backward_overlaps[right]

                left_to_right = right_ele in graph[left_ele]
                right_to_left = left_ele in graph[right_ele]

                if left_to_right and right_to_left and not gs.isZero(graph[left_ele][right_ele]) and not gs.isZero(graph[right_ele][left_ele]):
                        connection_list.append(graph[left_ele][right_ele])
                        connection_list.append(graph[right_ele][left_ele])

# calculate 95th percentile of connection scores
connection_list.sort()
index = int(len(connection_list) * 0.05) - 1  # Adjusted to ensure accurate indexing
connection_threshold = connection_list[index]


#@#@#@#@#@#@#@#@#@#@#@#@#@#@
# Adjusting config file
#@#@#@#@#@#@#@#@#@#@#@#@#@#@

# write threshold to config file
with open(output_file, 'a') as f:
    f.write("red: {}\n".format(red_threshold))

    f.write("connection: {}\n".format(connection_threshold))
    
print("Config file located at {}".format(os.path.normpath(output_file)))



