#!/usr/bin/env bash

# Print all info to log file
exec 1> "${0}.log.$(date +%s)" 2>&1

#### Pre-run setup
source ~/scripts/script_setup.sh
set +eu; conda activate BUSCO_v5.0.0; set -eu

PREFIX="Pocillopora_damicornis_SIPAv1"
NCPUS=48
TMPDIR="/scratch/timothy/tmp"
DB="/scratch/timothy/databases/busco"

WD="$PWD"
BUSCO_TMPDIR="${TMPDIR}/busco.$(date +%s)"


#### Start Script
run_cmd "md5sum ${PREFIX}.genes.pep.faa | tee ${PREFIX}.busco_protein.job_md5sum_list.txt"

run_cmd "mkdir ${BUSCO_TMPDIR}"
run_cmd "cd ${BUSCO_TMPDIR}"



## BUSCO - protein
F="pep.faa"

run_cmd "ln -s $WD/${PREFIX}.genes.pep.faa $F"

LINEAGE="metazoa_odb10"
run_cmd "busco --cpu ${NCPUS} --offline --lineage_dataset ${DB}/${LINEAGE} --mode protein --in ${F} --out ${F}.busco_${LINEAGE}"

LINEAGE="eukaryota_odb10"
run_cmd "busco --cpu ${NCPUS} --offline --lineage_dataset ${DB}/${LINEAGE} --mode protein --in ${F} --out ${F}.busco_${LINEAGE}"



## Move everything back
run_cmd "mv * ${WD}/"
run_cmd "cd ${WD}"
run_cmd "rm -fr ${BUSCO_TMPDIR}"


