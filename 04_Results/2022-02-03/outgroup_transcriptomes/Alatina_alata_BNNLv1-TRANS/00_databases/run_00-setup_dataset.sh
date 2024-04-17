#!/usr/bin/env bash

# Print all info to log file
exec 1> "${0}.log.$(date +%s)" 2>&1

#### Pre-run setup
source ~/scripts/script_setup.sh
set +eu; conda activate py27; set -eu

export PATH="$PATH:/home/timothy/programs/samtools-1.11/bin"
export PATH="$PATH:/home/timothy/programs/gatk-4.2.0.0"
export PATH="$PATH:/home/timothy/programs/bbmap"

PREFIX="Alatina_alata_BNNLv1-TRANS"

#### Start Script
run_cmd "samtools faidx ${PREFIX}.transcripts.pep.faa"

run_cmd "md5sum ${PREFIX}.transcripts.pep.faa > ${PREFIX}.md5sum_list.txt"


