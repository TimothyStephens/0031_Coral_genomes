#!/usr/bin/env bash

# Print all info to log file
exec 1> "${0}.log.$(date +%s)" 2>&1

#### Pre-run setup
source ~/scripts/script_setup.sh
set +eu; conda activate py27; set -eu

export PATH="$PATH:/home/timothy/programs/ncbi-blast-2.10.1+/bin"

PREFIX="Myxobolus_pendula_ONCAv1"

#### Start Script
run_cmd "md5sum ${PREFIX}.transcripts.cds.fna ${PREFIX}.transcripts.pep.faa | tee ${PREFIX}.makeblastdb.job_md5sum_list.txt"
run_cmd "makeblastdb -dbtype nucl -in ${PREFIX}.transcripts.cds.fna"
run_cmd "makeblastdb -dbtype prot -in ${PREFIX}.transcripts.pep.faa"


