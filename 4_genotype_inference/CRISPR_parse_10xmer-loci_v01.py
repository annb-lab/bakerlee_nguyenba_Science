"""
A program to parse single-end BT reads, counting barcodes

Read the file, for each read

    1) check that the read is good:
        - inline index is correct
        - quality of bc region is above some threshold
        - can extract barcode by regex
    2) record the barcode and umi (just first 7 bp), keeping track of counts of each

"""
import gzip
import numpy as np
import csv
import pandas as pd
import os
import subprocess
from collections import defaultdict, Counter
import regex
from milo_tools import FourLineFastq
print('all modules loaded')

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('demult_file', help='demultiplex info file')
parser.add_argument('output_base', help='basic output directory')
parser.add_argument('reads_base', help='directory with read files')
parser.add_argument('demult_id_key', help='key to identify which file to parse (which row of the demult file) - for jabba rays')
parser.add_argument('run_name', help='run name for stats file naming')
parser.add_argument("-quality_cutoff", type=int, default=25, help='quality threshold for bc region')

args = parser.parse_args()

# POSITIONAL ARGS - THESE SHOULD ALWAYS BE THE SAME FOR A RUN
DEMULT_FILE = args.demult_file
OUTPUT_BASE = args.output_base
if OUTPUT_BASE[-1] != '/':
    OUTPUT_BASE += '/'
READS_BASE = args.reads_base
DEMULT_ID_KEY = int(args.demult_id_key)
QUALITY_CUTOFF = int(args.quality_cutoff)
RUN_NAME = args.run_name

def make_dir(path, name=''):
    if not os.path.isdir(path):
        print('Making', name, ':', path)
        subprocess.call(['mkdir', path])
    else:
        print(name, 'directory was already made:', path)

make_dir(OUTPUT_BASE, 'output')
make_dir(OUTPUT_BASE + 'counts', 'count output')
make_dir(OUTPUT_BASE + 'offtarget_indices', 'offtarget_indices output')
make_dir(OUTPUT_BASE + 'umis', 'UMI output')

BC_REGEX_LIST = [regex.compile('(TCTGCC)(\D{22})(CGCTGA)'),
                 regex.compile('(TCTGCC)(\D{20,24})(CGCTGA)'),
                 regex.compile('(TCTGCC){e<=1}(\D{22})(CGCTGA){e<=1}'),
                 regex.compile('(TCTGCC){e<=1}(\D{20,24})(CGCTGA){e<=1}')]

# lists the inline indices we used that the script should look for
INDEX_LIST = ['CTCGTA','CATTGG','AAGGTC','GACTAC','TTCACG','ACACTT','CGAAGG','TCTTGC']

# Set up a table of regex for each locus
# If having issues with read counts, weird alleles, follow up here
regex_table = pd.DataFrame()
regex_table['locus']=['BUL2','FAS1','MKT1','NCS2','PMA1','RHO5','SCH9','WHI2','AKL1','RPI1','HSL7','SPT7','FRS1']
regex_table['left']=['CAACACAA','TACCAGGA','ATAATTGT','ATCCTGAA','TGTTTGTC','GTGTGATA','TCATTTCT','ATTTAATG','CTTGATAT','ATGGAAAG','ATCTGAAT','CAACCAAT','GCCAACCA']
regex_table['right']=['CTAACGTT','TAATTTAG','CATTATGT','CAGAATCT','TAATAGCA','AAAACGTC','GCCATTAA','TTCACCAA','AAGGTAGT','ATTACTAC','GATTTGCC','CATCCCTG','TCAACGGA']
regex_table['bp'] = [20,23,21,20,20,21,21,21,22,20,21,20,20]


class BcCounter:

    def __init__(self, library_name):
        self.library_name = library_name
        # counts for various outcomes for each read
        self.full_count = 0     # any read
        self.quality_count = 0  # quality score in bc region was too low
        self.inline_count = 0   # inline indices were not right
        self.locus_count = 0    # locus couldn't be identified
        self.regex_count = 0    # could not find bc with regex

        # dictionaries with the barcodes as keys like:
        # {bc: [total_counts (int), umis (dict like {umi: count})]
        self.bc_dict = defaultdict(lambda: [0, Counter()])

        # this keeps track of the different inline indices seen, which can provide information about
        # primer cross-contamination
        self.inline_indices_dict = Counter()
        
        # this keeps track of different loci seen
        self.locus_dict = Counter()

        #I've removed the bc_start argument from here, not from below though
    def add_entry(self, r1_seq, r1_qual):
        # input is an entry (4 lines from 1 single-end file)
        # first, identify and store the inline index present
        # Making the search fuzzy now!

        # Set a bc_start_default as 90 (UMI+distance between inline index and BC)
        bc_start_default = 0
        
        # set the inlinechkr parameter to 0, inline datum to None
        inlinechkr = 0
        
        inline = None
        
        # run through the exact matches. Most will only go through this first loop.
        for indi in np.arange(len(INDEX_LIST)):
            if r1_seq[0:5] == INDEX_LIST[indi]:
                inline = INDEX_LIST[indi]
                break
        
        # if no exact matches, try for an inexact match
        if not inline:
            for indi in np.arange(len(INDEX_LIST)):
                # for epochs 0, 2, and 4, allow up to two mismatches
                inline_regex_list = [regex.compile('('+INDEX_LIST[indi]+'){e<=1}')]
                for inreg_search in inline_regex_list:
                    inreg_hit = inreg_search.search(r1_seq[0:5])
                    # if fuzzy match works, use the associated index as inline
                    if inreg_hit:
                        inline = INDEX_LIST[indi]
                        break
        if not inline:
            # don't recognize the inline index
            inlinechkr = 1
            inline = r1_seq[8:20]

        if inlinechkr == 1:
            self.inline_count += 1
            self.inline_indices_dict[inline]
        
        else:
            bc_start = bc_start_default + 6 # 6 is length of inline
            # Check quality of bc region, continue if above threshold
            if np.mean([ord(c)-33 for c in r1_qual[bc_start:bc_start+22]]) < QUALITY_CUTOFF:
                self.quality_count += 1
            else:
                # Identify the locus
                
                thislocus = None
                
                for locus in np.arange(len(regex_table)):
                    if r1_seq.find(regex_table.loc[locus,'left']) > 0:
                        locind = locus
                        thislocus = regex_table.loc[locus,'locus']
                        break
                    elif r1_seq.find(regex_table.loc[locus,'right']) > 0:
                        locind = locus
                        thislocus = regex_table.loc[locus,'locus']
                        break
                
                if not thislocus:
                    # don't recognize this locus
                    self.locus_count +=1
                    self.locus_dict[r1_seq[8:20]]
                
                else:
                    bc = None
                    
                    rel_regex_table = regex_table[(regex_table['locus'] == thislocus)].reset_index(drop=True)
                    left = rel_regex_table['left'][0]
                    right = rel_regex_table['right'][0]
                    bp = rel_regex_table['bp'][0]

                    BC_REGEX_LIST = [regex.compile('('+left+')(\D{'+str(bp)+'})('+right+')'),
                                     regex.compile('('+left+')(\D{'+str(bp-2)+','+str(bp+2)+'})('+right+')'),
                                     regex.compile('('+left+'){e<=1}(\D{'+str(bp)+'})('+right+'){e<=1}'),
                                     regex.compile('('+left+'){e<=1}(\D{'+str(bp-2)+','+str(bp+2)+'})('+right+'){e<=1}')]
                    
                    for reg_search in BC_REGEX_LIST:
                        # Changing BC definition to include locus name and inline index
                        reg_hit = reg_search.search(r1_seq[bc_start:bc_start+106])
                        if reg_hit:
                            bc_alone = reg_hit.group(2)
                            bc = inline+'_'+thislocus+'-'+bc_alone
                            break
                    if not bc:
                        self.regex_count +=1
                    else:
                        # extract UMI (first 8 bp)
                        tmp_UMI = r1_seq[:8]
                        self.bc_dict[bc][0] += 1
                        self.bc_dict[bc][1][tmp_UMI] += 1
                # increment total real count
                self.full_count += 1

    def output_counts(self, fout, fout_UMI_fam, masterout):
        # Variable to keep track of the number of UMIs.  Note that since I just have 8 bp UMIs,
        # I expect a certain amount of UMI repeats just by chance (depending the number of counts) 
        UMI_repeat_count = 0
        # output results like BC, Count, UMI.Count (number of unique UMIs seen with this bc)
        outfile1 = open(fout, 'w')
        writer1 = csv.writer(outfile1)
        writer1.writerow(['BC', 'Reads', 'UMI.Count'])
        # output UMI family sizes (number of reads with the same barcode and UMI)
        outfile2 = open(fout_UMI_fam, 'w')
        writer2 = csv.writer(outfile2)
        writer2.writerow(['BC', 'Reads', 'UMI.fam.size.1', 'UMI.fam.size.2', 'etc'])
        for bc_tmp in self.bc_dict.keys():
            tmp_rec = self.bc_dict[bc_tmp]
            # the number of UMI repeats is the total count minus the number of unique UMIs
            UMI_repeat_count += (tmp_rec[0] - len(tmp_rec[1])) 
            writer1.writerow([bc_tmp, tmp_rec[0], len(tmp_rec[1])]) # Simple output
            # getting UMI family sizes
            tmp_dict = dict() # this dict is like {umi family size: number of fams with that size}
            for umi in tmp_rec[1]:
                UMI_fam_size = tmp_rec[1][umi]
                if UMI_fam_size in tmp_dict:
                    tmp_dict[UMI_fam_size] += 1
                else:
                    tmp_dict[UMI_fam_size] = 1
            max_fam_size = max([i for i in tmp_dict.keys()])
            output_row = [bc_tmp, tmp_rec[0]]
            writer2.writerow(output_row + [tmp_dict.setdefault(i, 0) for i in range(1, max_fam_size+1)])
        
        outfile1.close()
        outfile2.close()
        #Print stats
        print([self.library_name, self.full_count, self.quality_count, self.regex_count,
               UMI_repeat_count])

        # output statistics to a single file (not in a nice order though, will need to sort)
        with open(masterout, 'a') as outfile:
            writer = csv.writer(outfile)
            # there will be many title rows, I will delete all but one when cleaning this up
            writer.writerow(['Library', 'Reads', 'Inline.Index.Excluded', 'Quality.Excluded', 'Regex.Excluded',
                             'UMI.Repeats'])
            writer.writerow([self.library_name, self.full_count, self.inline_count, self.quality_count,
                             self.regex_count, UMI_repeat_count])

    def output_offtarget_inlines(self, inline_out):
        # outputs a csv with the counts for all inline indices seen that did not match the expectation
        with open(inline_out, 'w') as outfile:
            writer = csv.writer(outfile)
            writer.writerow(['R1_inline_index', 'Reads'])
            sorted_indices = sorted([i for i in self.inline_indices_dict], key=lambda x: -1*self.inline_indices_dict[x])
            for index in sorted_indices:
                writer.writerow([index, self.inline_indices_dict[index]])


def count_file(dem_info, sample_id, count_output_base, umi_output_base, stats_output):
    library_name = dem_info[0]
    #inline_index = dem_info[1]
    # the bc_start region depends on the TnF primer offset:
    #bc_start_point = 59
    # initialize bc counting class
    bc_counterator = BcCounter(library_name)
    
    R1 = READS_BASE + library_name

    print('Reading file:', R1)
    with gzip.open(R1, 'rt') as r1in:
        for R1_title, R1_seq, R1_qual in FourLineFastq(r1in):
            bc_counterator.add_entry(R1_seq, R1_qual)

    bc_counterator.output_counts(count_output_base + library_name + '_counts.csv',
                                 umi_output_base + library_name + '_umi_fam_sizes.csv',
                                 stats_output)

    #bc_counterator.output_offtarget_inlines(offind_base + library_name + '_offtarget_indices.csv')


### BEGIN MAIN ###

# Reading demult file, parsing one file based on demult id key
dem_dat = pd.read_csv(DEMULT_FILE)
rc = 1
for row in dem_dat.as_matrix(['Library']):
#    if rc == DEMULT_ID_KEY:
    count_file(row, rc, OUTPUT_BASE + 'counts/', OUTPUT_BASE + 'umis/', 
               OUTPUT_BASE + RUN_NAME + '_library_stats.csv')
rc += 1


#for row in dem_dat.as_matrix(['Library']):
#    if rc == DEMULT_ID_KEY:
#        count_file(row, rc, OUTPUT_BASE + 'counts/', OUTPUT_BASE + 'umis/', 
#                   OUTPUT_BASE + RUN_NAME + '_library_stats.csv')
#    rc += 1 
