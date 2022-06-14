#!/usr/bin/env bash

# Print all info to log file
exec 1> "${0}.log.$(date +%s)" 2>&1

#### Pre-run setup
source ~/scripts/script_setup.sh
set +eu; conda activate cutadapt35; set -eu

SRR="DRR194228"

R1="${SRR}_R1.fastq.gz"
R2="${SRR}_R2.fastq.gz"
R1_TRIMMED="${SRR}_trimmed_R1.fastq.gz"
R2_TRIMMED="${SRR}_trimmed_R2.fastq.gz"
LOG="${SRR}.cutadapt.log"
NCPUS=24

#### Start Script
# Adapter from https://cutadapt.readthedocs.io/en/stable/guide.html#illumina-truseq
run_cmd "cutadapt -q 20 --minimum-length 25 --cores ${NCPUS} -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -o ${R1_TRIMMED} -p ${R2_TRIMMED} ${R1} ${R2} > ${LOG} 2>&1"

# Check raw and trimmed reads for adapters. Set num reads to check to 100 billion (i.e., check all reads)
set +eu; conda activate py27; set -eu
run_cmd "./find_adapters_in_reads.py -n 100000000000 -a adapters.fa -f ${R1} ${R2} ${R1_TRIMMED} ${R2_TRIMMED}"

