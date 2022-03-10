#! /usr/bin/python3 -u

# ##############################################################################
# Performs parsing of fastq files to identify nextera wgs from the trial run of 
# BBQ.
# Modified by Alex N Nguyen Ba
# nnguyenba@fas.harvard.edu
# 2017
#
# Version 1
#
# LICENCE 
#
# ##############################################################################


import sys
from collections import OrderedDict
sys.path.append('/n/home00/nnguyenba/lib/python')

import regex

complements = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
def rc(seq):
	return "".join(complements.get(base, base) for base in reversed(seq))

def make_corrector(options, num_mismatch=0, rc=0):
	if rc:
		checkers = [regex.compile("("+rc(o)+"){e<=" + str(num_mismatch) + "}") for o in options]
	else:
		checkers = [regex.compile("("+o+"){e<=" + str(num_mismatch) + "}") for o in options]
	def corrector(match):
		current_error_count = 10000
		current_best_i = ""
		for (i,c) in enumerate(checkers):
			m = c.fullmatch(match)
			if m:
				error_count = m.fuzzy_counts[0] + m.fuzzy_counts[1] + m.fuzzy_counts[2]
				if(error_count == 0):
					return i+1
				elif (error_count < current_error_count and error_count <= num_mismatch):
					current_best_i = i+1
					current_error_count = error_count
		if (current_best_i != ""):
			return current_best_i
		else:
			return 0
	return corrector

def OR(xs,rc=0):
	if rc:
		return "(" + "|".join(["(?:"+rc(x)+")" for x in xs]) + ")"
	else:
		return "(" + "|".join(["(?:"+x+")" for x in xs]) + ")"


S5index = [
	"GCGATCTA",
	"ATAGAGAG",
	"AGAGGATA",
	"TCTACTCT",
	"CTCCTTAC",
	"TATGCAGT",
	"TACTCCTT",
	"AGGCTTAG",
	"GAGTAGCC",
	"GTCTGAGG",
	"CGTAAGGA",
	"CCACGCGT",
	"GGAGTTCC",
	"CATGGCCA",
	"AATCTCTC",
	"TAACCGCG",
	"TGGCGGTC",
	"CCATCTTA",
	"ATGTCAAT",
	"AGTTGGCT",
	"ACCTAGTA",
	"AACCGTGA",
	"TCATTACA",
	"CTGACGTG",
	"GAATTCAG",
	"CCGGTACG",
	"CCGTCATC",
	"CGTCTATA",
	"TCAATGAC",
	"AACGATGC",
	"GTCAACCT",
	"CAGTTTCA",
	"TGTGATTG",
	"TTGCATGT",
	"GGCGCGAT",
	"TTAACCGA"]

S5indexDict = {idx:(i+1) for i,idx in enumerate(S5index)}



N7index = [
	"TAAGGCGA",
	"CGTACTAG",
	"AGGCAGAA",
	"TCCTGAGC",
	"GGACTCCT",
	"TAGGCATG",
	"CTCTCTAC",
	"CAGAGAGG",
	"GCTACGCT",
	"CGAGGCTG",
	"AAGAGGCA",
	"GTAGAGGA",
	"ATTGTAAT",
	"GATCATTC",
	"ACCGATCG",
	"CCGTTATT",
	"TTCTTCTA",
	"TACCTGAC",
	"AGGACCGC",
	"GTCCGATT",
	"CACGAGTT",
	"CCACGGCC",
	"ACATGTAA",
	"TGTTAACT"]

N7indexDict = {idx:(i+1) for i,idx in enumerate(N7index)}

CRISPR_col_index = [
	"TTGATCG",
	"GCATTGAT",
	"AAGGCGTCT",
	"CATTCTCGAG",
	"CCTTAGTACTA",
	"GGTATCA",
	"GATGTCCT",
	"CTCATCCAG",
	"GACGGAACTC",
	"ATCGCAGGCAT"]

correct_CRISPR_col = make_corrector(CRISPR_col_index[0:10])

CRISPR_row_index = [
	"TACCTGA",
	"AGAACCAT",
	"ACCGTTGAC",
	"GCAAGCTGTT",
	"TCGGTGGTACG",
	"GTTCAACCGATT"]

correct_CRISPR_row = make_corrector(CRISPR_row_index[0:6])

molecular_barcode = ".{8}"
CRISPR_barcode = "(.{3})(?e)(?:AT){e<=1}(.{5})(?:AT){e<=1}(.{5})(?:AT){e<=1}(.{3})"
P1 = "(?e)(?:AAGGTACGATTCTGACGCA)"

barcode_read_1 = (
	"(" + molecular_barcode + ")" +
	OR(CRISPR_col_index[0:10]) + "{e<=2}" +
	P1 + "{e<=3}" + "(?e)(?:TGCC){e<=1}" +
	CRISPR_barcode +
	"(?e)(?:CGCT){e<=1}")

index_read_2 = (
	"(" + molecular_barcode + ")" +
	OR(CRISPR_row_index[0:6]) + "{e<=2}"
	)


barcode_read_1_re = regex.compile(barcode_read_1)
index_read_2_re = regex.compile(index_read_2)

nreads = 0
ndiscarded = 0

for line in sys.stdin:
	nreads += 1
	row=OrderedDict([
		(x,"") for x in ["illumina_index_1","illumina_index_2","col_index","row_index","barcode","mb1","mb2"]
	])

	splitline = line.split("	")
	row["illumina_index_1"] = N7indexDict[splitline[0]]
	row["illumina_index_2"] = S5indexDict[rc(splitline[1])]
	R1 = splitline[2]
	R2 = splitline[3]

	m1 = barcode_read_1_re.match(R1)
	if m1:	
		m2 = index_read_2_re.match(R2)
		if m2:
			row["mb1"] = m1.groups()[0]
			row["mb2"] = m2.groups()[0]
			row["col_index"] = correct_CRISPR_col(m1.groups()[1])
			row["row_index"] = correct_CRISPR_row(m2.groups()[1])
			row["barcode"] = m1.groups()[2] + m1.groups()[3] + m1.groups()[4] + m1.groups()[5]

			print(*row.values(), sep="	", file=sys.stdout)
		else:
			ndiscarded+=1
			continue
	else:
		ndiscarded+=1
		continue
