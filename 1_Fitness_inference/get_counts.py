# Open the parsed file and 
import string
import numpy as np
from scipy import linalg
import sys
import csv
import itertools
import time
import random
import argparse
import os
cwd = os.getcwd()
import psutil
process = psutil.Process(os.getpid())

from argparse import ArgumentParser, SUPPRESS
# Disable default help
parser = ArgumentParser(add_help=False)
required = parser.add_argument_group('required arguments')
optional = parser.add_argument_group('optional arguments')

# Add back help 
optional.add_argument(
    '-h',
    '--help',
    action='help',
    default=SUPPRESS,
    help='show this help message and exit'
)
required.add_argument('-order', help='notation is X1,Y1_X2,Y2 for the index order (col/row)')
required.add_argument('-file', help='filename with corrected barcodes')
args = parser.parse_args()

print(args, file=sys.stderr)

# Now parse the argument
orders = args.order.split("_")

#for i in range(len(orders)):
#	index = orders[i].split(",")

#	col[i] = index[0]
#	row[i] = index[1]

# Now read in the file
print("BC	7	14	28	42	49")
counts = {}
with open(args.file,'r') as readfile:
	for line in readfile:
		line.rstrip()
		fields = line.split("\t")
		if(fields[2] + "," + fields[3] in orders):
			order_index = orders.index(fields[2] + "," + fields[3])
			
			BC = fields[4]

			if BC not in counts:
				counts[BC] = [0] * 5

			counts[BC][order_index] = counts[BC][order_index] + 1
			
for BC in counts:
	print(str(BC),end="\t")
	print(*counts[BC],sep="\t")

