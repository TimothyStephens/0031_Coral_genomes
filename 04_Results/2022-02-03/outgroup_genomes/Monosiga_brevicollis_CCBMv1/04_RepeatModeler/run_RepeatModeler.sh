#!/usr/bin/env bash

# Print all info to log file
exec 1> "${0}.log.$(date +%s)" 2>&1

#### Pre-run setup
source ~/scripts/script_setup.sh
set +eu; conda activate py27; set -eu

export PATH="/home/timothy/programs/RepeatAnalysis/RepeatModeler-2.0.1:$PATH"
export PATH="/home/timothy/miniconda3/envs/py27/bin:$PATH"
export PERL5LIB="/home/timothy/miniconda3/envs/py27/lib/5.26.2/x86_64-linux-thread-multi"

PA=12
PREFIX="Monosiga_brevicollis_CCBMv1"
REF="${PREFIX}.assembly.fasta"

#### Start Script
run_cmd "BuildDatabase -name $REF.RMDB $REF 1> BuildDatabase.log 2>&1"
run_cmd "RepeatModeler -database $REF.RMDB -pa $PA -LTRStruct -genomeSampleSizeMax 27000000 1> RepeatModeler.log 2>&1"


