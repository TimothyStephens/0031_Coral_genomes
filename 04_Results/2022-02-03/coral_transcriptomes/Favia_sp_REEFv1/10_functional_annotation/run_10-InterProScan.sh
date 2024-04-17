#!/usr/bin/env bash

# Print all info to log file
exec 1> "${0}.log.$(date +%s)" 2>&1

#### Pre-run setup
source ~/scripts/script_setup.sh
set +eu; conda activate py38; set -eu

INTERPROSCAN="/home/timothy/programs/my_interproscan/interproscan-5.53-87.0/interproscan.sh"
NCPUS=96
PREFIX="Favia_sp_REEFv1"

#### Start Script
run_cmd "md5sum ${PREFIX}.transcripts.pep.faa ../00_databases/${PREFIX}.transcripts.pep.faa | tee ${PREFIX}.transcripts.pep.faa.InterProScan.job_md5sum_list.txt"
run_cmd "${INTERPROSCAN} -dp --goterms --cpu ${NCPUS} --input ${PREFIX}.transcripts.pep.faa --output-file-base ${PREFIX}.transcripts.pep.faa.InterProScan"
run_cmd "rm -fr temp"


