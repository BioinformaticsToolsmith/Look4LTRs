import sys
import os
from anytree import Node, RenderTree

if len(sys.argv) != 2:
    print("Usage: findNested.py <path>")
    print("Given a directory or path to an RTR file, find the same-graph nested elements and print them out.")
    print("Argument One: Path to a directory or file")
    sys.exit(1)

path = sys.argv[1]
assert os.path.exists(path), "Path does not exist: %s" % path
file_list = [path] if not os.path.isdir(path) else [os.path.join(path, f) for f in os.listdir(path)]

def find_nested(lines, visited, x, graph_id, data, parent):
    if data[6] != "NA":
        for nested in data[6][1:-1].split(',')[:-1]:
            for y in range(x + 1, len(lines)):
                nest_data = lines[y].split()
                if nest_data[1] == nested and nest_data[14] == graph_id and nest_data[4] != "NA":
                    visited.add(y)
                    node_id = f"{nest_data[0]}_{nest_data[1]}"
                    node = Node(node_id, parent=parent)
                    find_nested(lines, visited, y, graph_id, nest_data, node)

nest_list = []
for file in file_list:
    with open(file, 'r') as f:
        visited = set()
        lines = f.readlines()
        for x in range(1, len(lines)):
            if x in visited:
                continue
            data = lines[x].split()
            if data[4] != "NA":
                graph_id = data[14]
                parent_id = f"{data[0]}_{data[1]}"
                root = Node(parent_id)
                find_nested(lines, visited, x, graph_id, data, root)
                
                if root.children:
                    nest_list.append(root)


def get_max_depth(node):
    """
    Recursive function that calculates the maximum depth of the tree.
    """
    if not node.children:
        return 0
    else:
        max_depth = 0
        for child in node.children:
            child_depth = get_max_depth(child)
            max_depth = max(max_depth, child_depth)
        return max_depth + 1

max_deep = 0
for root in nest_list:
    max_deep = max(max_deep, get_max_depth(root))
    for pre, fill, node in RenderTree(root):
        print("%s%s" % (pre, node.name))

print("Maximum nesting depth: %d" % max_deep)
print("Number of nests: %d" % len(nest_list))
