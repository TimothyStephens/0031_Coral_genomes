#!/usr/bin/env bash

# Print all info to log file
exec 1> "${0}.log.$(date +%s)" 2>&1

#### Pre-run setup
source ~/scripts/script_setup.sh
set +eu; conda activate eggnog-mapper; set -eu

export PATH="$PATH:/home/timothy/programs/eggnog-mapper-2.1.6"
DATABASES="/scratch/timothy/databases/eggnog-mapper-rel20211209"

NCPUS=96
PREFIX="Acropora_gemmifera_SLJPv1"

#### Start Script
run_cmd "md5sum ${PREFIX}.genes.pep.faa ../00_databases/${PREFIX}.genes.pep.faa | tee ${PREFIX}.genes.pep.faa.emapper.job_md5sum_list.txt"
run_cmd "emapper.py -i ${PREFIX}.genes.pep.faa --output ${PREFIX}.genes.pep.faa --data_dir ${DATABASES} --pfam_realign denovo --report_orthologs --no_file_comments --cpu ${NCPUS} --dbmem"
run_cmd "rm -fr emappertmp_*"


