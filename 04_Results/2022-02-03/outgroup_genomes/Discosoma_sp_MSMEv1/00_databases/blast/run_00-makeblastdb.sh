#!/usr/bin/env bash

# Print all info to log file
exec 1> "${0}.log.$(date +%s)" 2>&1

#### Pre-run setup
source ~/scripts/script_setup.sh
set +eu; conda activate py27; set -eu

export PATH="$PATH:/home/timothy/programs/ncbi-blast-2.10.1+/bin"

PREFIX="Discosoma_sp_MSMEv1"

#### Start Script
run_cmd "md5sum ${PREFIX}.assembly.fasta ${PREFIX}.genes.cds.fna ${PREFIX}.genes.pep.faa | tee ${PREFIX}.makeblastdb.job_md5sum_list.txt"
run_cmd "makeblastdb -dbtype nucl -in ${PREFIX}.assembly.fasta"
run_cmd "makeblastdb -dbtype nucl -in ${PREFIX}.genes.cds.fna"
run_cmd "makeblastdb -dbtype prot -in ${PREFIX}.genes.pep.faa"


