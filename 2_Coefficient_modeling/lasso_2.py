# Open the parsed file and 
import string
import numpy as np
import scipy.special
from sklearn.linear_model import LinearRegression, LassoCV
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
required.add_argument('-i', help="input file")
optional.add_argument('--map', help="Barcode to genotype map", default="/n/home00/nnguyenba/scripts/CRISPR/modeling/BC_geno_map.txt")
optional.add_argument('--order', help="epistasis order to go to", default=1,type=int)
optional.add_argument('--basis', help="Whether to use the orthogonal basis", default=0,type=int)
args = parser.parse_args()

print(args, file=sys.stderr)

map = {}
# Parse the map file
with open(args.map, "r") as readfile:
	for line in readfile:
		line = line.rstrip()

		fields = line.split("\t")
		map[fields[0]] = fields[1]


# Parse the fitness file
genotypes = []
unparsed_genotypes = []
fitnesses = []
errors = []
# Number of 'loci' to consider
num_loci = 0
for i in range(1,args.order+1):
	num_loci = num_loci + int(scipy.special.binom(10,i))

with open(args.i, 'r') as readfile:
	for line in readfile:
		line = line.rstrip()
		fields = line.split("\t")

		if("Lineage" in line or len(fields) <4):
			continue

		geno = map[fields[0]]
		fitness = fields[1]
		SE = fields[2]
		# can filter by standard error here if needed. Skip for now.
		#if '2' not in geno and abs(float(fitness)) < 1:
		if '2' not in geno and abs(float(SE)) < 1:
			fitnesses.append(float(fitness))
			errors.append(float(fields[2]))

			#arr = [0] * num_loci
			arr = []
			for i in range(1,1024):
				binary = format(i,'b').zfill(10)
				count_1 = binary.count("1")
				if(count_1 <= args.order):
					if(args.basis == 0):
						# 0/1 basis
						# int(binary,2) & int(geno,2) this returns decimal representation of binary AND (1 for every matched position, and 0 otherwise).
						#print(str(binary) + "	" + str(geno) + "	" + str(int(binary,2) & int(geno,2)) + "	" + str(int(binary,2)))
						if(int(binary,2) & int(geno,2) == int(binary,2)):					
							arr.append(1)
						else:
							arr.append(0)
					else:
					
						# -1/1 basis
						# So we're going to do this by first bitmasking, then counting the number of 1s in the final mask. 
						mask = bin(int(binary,2) & int(geno,2))[2:]
						mask_1 = mask.count("1")
						#print(str(binary) + "	" + str(geno) + "	" + str(int(binary,2) & int(geno,2)) + "	" + str(int(binary,2)) + "	" + str((count_1 - mask_1) % 2))

						# The value of zeros is count_1 - mask_1
						# If the count is even, then append 1, else append -1
						if((count_1 - mask_1) % 2 == 0):
							# is even
							arr.append(1)
						else:
							arr.append(-1) 					

			genotypes.append(arr)
			unparsed_genotypes.append(geno)


genotypes = np.array(genotypes)

#reg = LinearRegression().fit(genotypes,fitnesses)
reg = LassoCV(cv=5, random_state=0, max_iter=2000).fit(genotypes,fitnesses)
print(reg.score(genotypes,fitnesses))
print(reg.alpha_)

predictions = reg.predict(genotypes)

for i in range(len(predictions)):
	#print(str(genotypes[i]) + "	" + str(fitnesses[i]) + "	" + str(predictions[i]))
	print(str(unparsed_genotypes[i]) + "	" + str(predictions[i]) + "	" + str(fitnesses[i]) + "	" +str(errors[i]) )
print()

iterator = 0
for i in range(1,1024):
	binary = format(i,'b').zfill(10)
	if(binary.count("1") <= args.order):
		print(str(binary.count("1")) + "	" + str(binary) + "	" + str(reg.coef_[iterator]))
		iterator = iterator + 1

exit()