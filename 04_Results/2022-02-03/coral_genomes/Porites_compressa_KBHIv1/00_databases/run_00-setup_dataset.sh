#!/usr/bin/env bash

# Print all info to log file
exec 1> "${0}.log.$(date +%s)" 2>&1

#### Pre-run setup
source ~/scripts/script_setup.sh
set +eu; conda activate py27; set -eu

export PATH="$PATH:/home/timothy/programs/samtools-1.11/bin"
export PATH="$PATH:/home/timothy/programs/gatk-4.2.0.0"
export PATH="$PATH:/home/timothy/programs/bbmap"

PREFIX="Porites_compressa_KBHIv1"

#### Start Script
run_cmd "samtools faidx ${PREFIX}.assembly.fasta"
run_cmd "samtools faidx ${PREFIX}.genes.cds.fna"
run_cmd "samtools faidx ${PREFIX}.genes.pep.faa"

run_cmd "gatk CreateSequenceDictionary --REFERENCE ${PREFIX}.assembly.fasta"
run_cmd "gatk CreateSequenceDictionary --REFERENCE ${PREFIX}.genes.cds.fna"

run_cmd "stats.sh -Xmx12g ${PREFIX}.assembly.fasta 1>${PREFIX}.assembly.fasta.stats 2>&1"
run_cmd "stats.sh -Xmx12g ${PREFIX}.genes.cds.fna  1>${PREFIX}.genes.cds.fna.stats 2>&1"

run_cmd "md5sum ${PREFIX}.assembly.fasta ${PREFIX}.genes.cds.fna ${PREFIX}.genes.pep.faa ${PREFIX}.genes.gff3 > ${PREFIX}.md5sum_list.txt"


