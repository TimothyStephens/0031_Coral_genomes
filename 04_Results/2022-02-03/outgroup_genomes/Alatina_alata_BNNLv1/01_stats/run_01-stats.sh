#!/usr/bin/env bash

# Print all info to log file
exec 1> "${0}.log.$(date +%s)" 2>&1

#### Pre-run setup
source ~/scripts/script_setup.sh
set +eu; conda activate py27; set -eu

export PATH="$PATH:/home/timothy/programs/bbmap"
export PATH="$PATH:/home/timothy/programs/bedtools-2.29.2/bin"

PREFIX="Alatina_alata_BNNLv1"

#### Start Script
run_cmd "md5sum ${PREFIX}.assembly.fasta | tee ${PREFIX}.genome_stats.job_md5sum_list.txt"
run_cmd "./genome_stats.tmp ${PREFIX}.assembly.fasta > ${PREFIX}.GeneStats.genome.tsv"
run_cmd "./CDS_stats ${PREFIX}.genes.cds.fna > ${PREFIX}.GeneStats.CDS.tsv"

