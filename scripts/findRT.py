import sys
import os

if len(sys.argv) != 2:
    print("Usage: findNested.py <path>")
    print("Given an RTR file, type an ID to get the corresponding RT")
    print("Argument One: Path to a file")
    sys.exit(1)

path = sys.argv[1]
assert os.path.exists(path), "I did not find the file at, "+str(path)

rt_dict = {}
with open(path, 'r') as f:
    line = f.readline()
    line = f.readline()
    while line:
        data = line.split()
        rt_dict[int(data[1])] = line
        line = f.readline()


# Get user input until they type 'exit'
while True:
    id = input("Enter an ID: ")
    if id.isnumeric() == False:
        break
    print(rt_dict[int(id)])

