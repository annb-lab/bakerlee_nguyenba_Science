#!/bin/bash
#
#SBATCH -p serial_requeue,desai # Partition to submit to (comma separated)
#SBATCH -J 10xmerbatch_v01 # Job name
#SBATCH -n 1 # Number of cores
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -t 0-10:00 # Runtime in D-HH:MM (or use minutes)
#SBATCH --mem 8000 # Memory in MB
#SBATCH -o /dev/null
#SBATCH --error=/n/desai_lab/users/cwbaker/20201102_10xmer-genotyping/"%x.e%j"

CURRDIR=$PWD

module load Anaconda3/5.0.1-fasrc02
source activate DPEbfa

python /n/desai_lab/users/cwbaker/20201102_10xmer-genotyping/CRISPR_parse_10xmer-loci_v01.py /n/desai_lab/users/cwbaker/20201102_10xmer-genotyping/csvcollection/lib${SLURM_ARRAY_TASK_ID}.csv /n/desai_lab/users/cwbaker/20201102_10xmer-genotyping/output_batch /n/desai_lab/users/cwbaker/20201102_10xmer-genotyping/fastq/CRISPR/ 0 CRISPR_10xmer-g