#!/usr/bin/env bash

# Print all info to log file
exec 1> "${0}.log.$(date +%s)" 2>&1

#### Pre-run setup
source ~/scripts/script_setup.sh
set +eu; conda activate py27; set -eu

export PATH="$PATH:/home/timothy/programs/bbmap"
export PATH="$PATH:/home/timothy/programs/bedtools-2.29.2/bin"

PREFIX="Madracis_auretenra_REEFv1"

#### Start Script
run_cmd "md5sum ${PREFIX}.transcripts.cds.fna | tee ${PREFIX}.CDS_stats.job_md5sum_list.txt"
run_cmd "./CDS_stats ${PREFIX}.transcripts.cds.fna > ${PREFIX}.GeneStats.tsv"


