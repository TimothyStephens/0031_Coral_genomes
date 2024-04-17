#!/usr/bin/env bash

# Print all info to log file
exec 1> "${0}.log.$(date +%s)" 2>&1

#### Pre-run setup
source ~/scripts/script_setup.sh
set +eu; conda activate py27; set -eu

export PATH="$PATH:/home/timothy/programs/bowtie2-2.4.4-linux-x86_64"
NCPUS=24

PREFIX="Kudoa_iwatai_RSILv1"

#### Start Script
run_cmd "md5sum ${PREFIX}.transcripts.cds.fna | tee ${PREFIX}.bowtie2-build.job_md5sum_list.txt"
run_cmd "bowtie2-build --threads ${NCPUS} ${PREFIX}.transcripts.cds.fna ${PREFIX}.transcripts.cds.fna"


