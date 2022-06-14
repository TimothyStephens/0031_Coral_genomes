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

SRR="DRR194228"
R1_TRIMMED="${SRR}_trimmed_R1.fastq.gz"
R2_TRIMMED="${SRR}_trimmed_R2.fastq.gz"

#### Start Script
KMER=21
run_cmd "jellyfish count -C -m $KMER -s $HASHSIZE -o smudgeplot.jellyfish.$KMER.jf -t $NCPUS -F 2 <(zcat ${R1_TRIMMED}) <(zcat ${R2_TRIMMED})"
run_cmd "jellyfish histo -t $NCPUS smudgeplot.jellyfish.$KMER.jf > smudgeplot.jellyfish.$KMER.histo"

L=$($smudgeplot/smudgeplot.py cutoff smudgeplot.jellyfish.$KMER.histo L)
U=$($smudgeplot/smudgeplot.py cutoff smudgeplot.jellyfish.$KMER.histo U)
echo "L:$L U:$U # these need to be sane values like 30 800 or so"
run_cmd "jellyfish dump -c -L $L -U $U smudgeplot.jellyfish.$KMER.jf | $smudgeplot/smudgeplot.py hetkmers -o smudgeplot_kmer_pairs_k${KMER}"
run_cmd "$smudgeplot/smudgeplot.py plot smudgeplot_kmer_pairs_k${KMER}_coverages.tsv -o smudgeplot_results_k${KMER}"


