#!/usr/bin/env bash

# Print all info to log file
exec 1> "${0}.log.$(date +%s)" 2>&1

#### Pre-run setup
source ~/scripts/script_setup.sh
set +eu; conda activate py27; set -eu

export PATH="$PATH:/home/timothy/programs/FastQC"
export PATH="$PATH:/home/timothy/programs/bbmap"

NCPUS=4

#### Start Script
run_cmd "fastqc --threads ${NCPUS} *.fastq.gz"
run_cmd "statswrapper.sh *.fastq.gz 2>read_stats.tsv"



