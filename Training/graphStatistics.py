import sys
import os
import re
import statistics
import gc
from math import isclose


tolerance = 1e-6
def isZero(a):
    return isclose(a, 0.0, abs_tol=tolerance)

class Element:
    def __init__(self, ele_str, ):
        self.ele_str = ele_str
        data = self.ele_str.split(':')
        self.start = int(data[0])
        self.end = int(data[1][:-1])
        self.direction = self.ele_str[-1]

    def get_start(self, min_start = 0):
        return self.start - min_start

    def get_end(self, min_start = 0):
        return self.end - min_start

    def get_direction(self, is_num = False):
        return self.direction if is_num == False else 10000 if self.direction == "+" else 0

    def calc_coverage(self, ltr, respect = True):
        overlap = min(self.end, ltr.get_end()) - max(self.start, ltr.get_start())
        return overlap / (self.end - self.start) if respect else overlap / (ltr.get_size())

    def check_reciprocal(self, ltr, thresh = 0.8):
        return self.calc_coverage(ltr) >= thresh and self.calc_coverage(ltr, False) >= thresh

    def get_node(self):
        return str(self.start) + self.direction

    def __str__(self):
        return self.ele_str

    def __repr__(self):
        return self.__str__()

    def __hash__(self):
        return hash(self.ele_str)

    def __eq__(self, other):
        return self.ele_str == other.ele_str

class LTR:
    def __init__(self, start, end):
        self.start = start
        self.end = end
    
    def get_start(self):
        return self.start

    def get_end(self):
        return self.end

    def get_size(self):
        return self.end - self.start

class RT:
    def __init__(self, line):
        data = line.split()
        self.l = LTR(int(data[3]), int(data[4]))
        self.r = LTR(int(data[5]), int(data[6]))

    def get_left(self):
        return self.l

    def get_right(self):
        return self.r

    def get_lstart(self):
        return self.l.get_start()

    def get_lend(self):
        return self.l.get_end()

    def get_rstart(self):
        return self.r.get_start()

    def get_rend(self):
        return self.r.get_end()

class ConStat:
    def __init__(self, con_list):

        self.median = statistics.median(con_list)
        self.mean = statistics.mean(con_list)
        self.min = min(con_list)
        self.max = max(con_list)
        self.std = statistics.stdev(con_list)

    def get_median(self):
        return self.median

    def get_mean(self):
        return self.mean

    def get_min(self):
        return self.min

    def get_max(self):
        return self.max

    def __str__(self):
        r = []
        r.append(f'Min   : {self.min}')
        r.append(f'Max   : {self.max}')
        r.append(f'Median: {self.median}')
        r.append(f'Mean  : {self.mean}')
        r.append(f'STD  : {self.std}')
        return '\n'.join(r)

    def print(self):
        print(self.__str__())



def to_file(path, con_list):
    with open(path, 'w') as file:
        for tup in con_list:
            file.write('\t'.join(str(x) for x in tup) + "\n")
    

    
def load_gsbed(path):
    rt_list = []
    with open(path, 'r') as file:
        line = file.readline()
        while not line.split()[1].isnumeric():
            line = file.readline()
        while line:
            rt_list.append(RT(line))
            line = file.readline()

    return rt_list

def load_graph(path):
    ele_dict = {}
    with open(path, 'r') as file:
        line = file.readline()
        while line:
            data = line.split('-->')
            left = data[0].split()
            right = data[1].split('}')
            pos = left[1]
            direction = left[-1]
            
            assert direction == '+' or direction == '-'

            match_dict = {}
            for group in right[:-1]:
                group.split('{')[1]
                group.split('{')[1].split()
                group.split('{')[1].split()[0]

                match_pos = group.split('{')[1].split()[0]
                match_direction = '+' if direction == '-' else '-'
                weight = float(group.split()[-1])
                
                match_dict[Element(match_pos + match_direction)] = weight

            ele_dict[Element(pos + direction)] = match_dict

            line = file.readline()

    return ele_dict

def sort_key(file_name):
    data = file_name.split('_')  
    genome_name, chr_name = '_'.join(data[:-1]), data[-1].split('.')[0] # Extract the genome and chromosome names from the file name
    return (genome_name, int(chr_name.replace('chr', '')))  # Sort by genome name, then chromosome number


def match_overlaps(ele_list, ltr_list):
    i = 0
    j = 0
    r = {}

    while i < len(ele_list):
        while j < len(ltr_list):

            if i >= len(ele_list):
                break


            elif ele_list[i].check_reciprocal(ltr_list[j], 0.01):
                r[ltr_list[j]] = ele_list[i] 
            
            if ele_list[i].get_end() < ltr_list[j].get_end():
                i += 1

            elif ele_list[i].get_end() > ltr_list[j].get_end():
                j += 1
            
            elif ele_list[i].get_end() == ltr_list[j].get_end():
                i += 1
                j += 1

            else:
                print("How?")
        break

    return r

            


if __name__ == "__main__":

    if len(sys.argv) != 4 or sys.argv[1] == "-h" or sys.argv[1] == "--help":
        print("Given a gs bed file directory and a graph directory, calculate the statistics of the weights betwen valid LTR nodes")
        print("Arg 1: path to gs bed directory")
        print("Arg 2: path to graph directory")
        print("Arg 3: output directory")
        sys.exit(1)

    gs_dir = sys.argv[1] + "/"
    graph_dir = sys.argv[2] + "/"
    out_dir = sys.argv[3] + "/"
    
    assert os.path.exists(gs_dir)
    assert os.path.exists(graph_dir)
    assert os.path.exists(out_dir)

    gs_list = [gs_dir + x for x in sorted(os.listdir(gs_dir), key = sort_key)]
    graph_list = [graph_dir + x for x in sorted(os.listdir(graph_dir), key = sort_key)]

    assert len(gs_list) == len(graph_list)

    con_list = []
    left_only_list = []
    right_only_list = []
    uncon_list = []
    for gs_path, graph_path in zip(gs_list, graph_list):
        print(f'Loading {gs_path} and {graph_path}')
        rt_list = load_gsbed(gs_path)
        graph = load_graph(graph_path)
        ltr_list = []
        for rt in rt_list:
            ltr_list.append(rt.get_left())
            ltr_list.append(rt.get_right())


        print('Looking for overlaps')
        forward_graph = sorted([ele for ele in graph if ele.get_direction() == "+"], key = lambda x : x.get_start())
        backward_graph = sorted([ele for ele in graph if ele.get_direction() == "-"], key = lambda x : x.get_start())

        forward_overlaps = match_overlaps(forward_graph, ltr_list)
        backward_overlaps = match_overlaps(backward_graph, ltr_list)

        print("Matching connections")
        for rt in rt_list:
            left = rt.get_left()
            right = rt.get_right()
            if left in forward_overlaps and right in backward_overlaps:
                left_ele = forward_overlaps[left]
                right_ele = backward_overlaps[right]

                left_to_right = right_ele in graph[left_ele]
                right_to_left = left_ele in graph[right_ele]

                if left_to_right and right_to_left:
                    if not isZero(graph[left_ele][right_ele]) and not isZero(graph[right_ele][left_ele]):
                        con_list.append((graph[left_ele][right_ele], graph[right_ele][left_ele]))

                elif left_to_right and not right_to_left:
                    if not isZero(graph[left_ele][right_ele]):
                        left_only_list.append((graph[left_ele][right_ele], graph[right_ele]))

                elif not left_to_right and right_to_left:
                    if not isZero(graph[right_ele][left_ele]):
                        right_only_list.append((graph[left_ele], graph[right_ele][left_ele]))

                else:
                    uncon_list.append((graph[left_ele], graph[right_ele]))

    print("Collecting statistics")

    lr_list = [x[0] for x in con_list]
    lr_stat = ConStat(lr_list)
    del lr_list
    gc.collect()

    rl_list = [x[1] for x in con_list]
    rl_stat = ConStat(rl_list)
    del rl_list
    gc.collect()

    avg_list = [(x[0] + x[1]) / 2 for x in con_list]
    avg_stat = ConStat(avg_list)
    del avg_list
    gc.collect()

    rat_list = [min(x[0], x[1]) / max(x[0], x[1]) for x in con_list]
    rat_stat = ConStat(rat_list)
    del rat_list
    gc.collect()


    print("\nLeft to Right")
    lr_stat.print()

    print("\nRight to Left")
    rl_stat.print()

    print("\nAverage of Weights")
    avg_stat.print()

    print("\nRatio of Weights")
    rat_stat.print()

    print(f"\nPercentage of One-way Connections: {(len(left_only_list) + len(right_only_list)) / (len(rt_list))}")

    print(f"\nPercentage of Missed Connections: {len(uncon_list) / len(rt_list)}")

    to_file(f'{out_dir}/TwoWay.xlsl', con_list)
    to_file(f'{out_dir}/LeftOnly.xlsl', left_only_list)
    to_file(f'{out_dir}/RightOnly.xlsl', right_only_list)
    to_file(f'{out_dir}/Unconnected.xlsl', uncon_list)

