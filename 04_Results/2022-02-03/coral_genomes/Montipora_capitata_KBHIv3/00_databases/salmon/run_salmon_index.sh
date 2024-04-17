#!/usr/bin/env bash

# Print all info to log file
exec 1> "${0}.log.$(date +%s)" 2>&1

#### Pre-run setup
source ~/scripts/script_setup.sh
set +eu; conda activate py27; set -eu


export PATH="$PATH:/home/timothy/programs/salmon-1.8.0_linux_x86_64/bin"
NCPUS=48

TRANSCRIPTS="Montipora_capitata_KBHIv3.genes.cds.fna"
INDEX="$TRANSCRIPTS.salmon.idx"

#### Start Script
run_cmd "salmon index --transcripts $TRANSCRIPTS --index $INDEX --threads $NCPUS"


