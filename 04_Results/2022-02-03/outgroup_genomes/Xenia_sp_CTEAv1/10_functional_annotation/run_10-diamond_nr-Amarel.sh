#!/usr/bin/env bash

# Print all info to log file
exec 1> "${0}.log.$(date +%s)" 2>&1

#### Pre-run setup
source ~/scripts/script_setup.sh
set +eu; conda activate py27; set -eu
ulimit -n 10000 # Need to set high to avoid 'Too many open files' error

export PATH="$PATH:/home/timothy/programs/DIAMOND_v2.0.15/bin"
DB="/scratch/databases/ncbi/2022_07/nr.dmnd"
NCPUS=24

PREFIX="Xenia_sp_CTEAv1"
QUERY="${PREFIX}.genes.pep.faa"
OUT="${QUERY}.diamond_blastp_nr.outfmt6_short"

if [ ! -f "${OUT}.gz" ];
then
	run_cmd "md5sum ${QUERY} ../00_databases/${PREFIX}.genes.pep.faa | tee ${OUT}.job_md5sum_list.txt"
	run_cmd "diamond blastp --ultra-sensitive --max-target-seqs 0 --evalue 1e-5 --query ${QUERY} --db ${DB} --out ${OUT}.tmp --outfmt 6 qseqid sseqid evalue bitscore staxids --compress 1 --threads ${NCPUS} && mv ${OUT}.tmp.gz ${OUT}.gz"
else
	log "Already finished running this part!"
fi


