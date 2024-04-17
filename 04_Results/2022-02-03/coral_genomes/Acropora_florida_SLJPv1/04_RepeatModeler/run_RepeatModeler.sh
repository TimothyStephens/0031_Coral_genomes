#!/usr/bin/env bash

# Print all info to log file
exec 1> "${0}.log.$(date +%s)" 2>&1

#### Pre-run setup
source ~/scripts/script_setup.sh
set +eu; conda activate py27; set -eu

export PATH="$PATH:/scratch/ts942/RepeatAnalysis/RepeatModeler-2.0.1"
PA=5
PREFIX="Acropora_florida_SLJPv1"
REF="${PREFIX}.assembly.fasta"

#### Start Script
#run_cmd "BuildDatabase -name $REF.RMDB $REF 1> BuildDatabase.log 2>&1"
run_cmd "RepeatModeler -database $REF.RMDB -pa $PA -LTRStruct -recoverDir RM_31794.TueMar82337172022 1>> RepeatModeler.log 2>&1"


