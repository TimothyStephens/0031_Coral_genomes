#!/usr/bin/env bash

# Print all info to log file
exec 1> "${0}.log.$(date +%s)" 2>&1

#### Pre-run setup
source ~/scripts/script_setup.sh
set +eu; conda activate smudgeplot; set -eu
export PYTHONPATH=/home/timothy/miniconda3/envs/smudgeplot/lib/python3.8:/home/timothy/miniconda3/envs/smudgeplot/lib/python3.8/site-packages

export PATH="$PATH:/home/timothy/programs/jellyfish-2.3.0/bin"
export PATH="$PATH:/home/timothy/programs/smudgeplot/bin"
smudgeplot="/home/timothy/programs/smudgeplot/bin"
NCPUS=48
HASHSIZE="32G"

#### Start Script
KMER=21

## 6x is just after minimum between 1st peak and initial error peak
## 136x is 8.5 x 16x (where 16x is the estimated haplid peak by GenomeScope2)
L=6
U=136
N=16
echo "L:$L U:$U # these need to be sane values like 30 800 or so"
run_cmd "jellyfish dump -c -L $L -U $U smudgeplot.jellyfish.$KMER.jf | $smudgeplot/smudgeplot.py hetkmers -o smudgeplot_kmer_pairs_k${KMER}_manual_cutoffs"
run_cmd "$smudgeplot/smudgeplot.py plot smudgeplot_kmer_pairs_k${KMER}_manual_cutoffs_coverages.tsv -o smudgeplot_results_k${KMER}_manual_cutoffs -n $N"


