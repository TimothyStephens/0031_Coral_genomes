#!/usr/bin/env bash
  
# Print all info to log file
exec 1> "${0}.log.$(date +%s)" 2>&1

#### Pre-run setup
source ~/scripts/script_setup.sh
set +eu; conda activate RepeatMasker; set -eu

export PATH="$PATH:/home/timothy/programs/RepeatAnalysis/RepeatMasker-4.1.2-p1"
export PATH="$PATH:/home/timothy/programs/bedtools-2.29.2/bin"
export PATH="$PATH:/home/timothy/programs/bbmap"
UTILS="/home/timothy/programs/RepeatAnalysis/RepeatMasker-4.1.2-p1/util"
PA=18
PREFIX="Calvadosia_cruxmelitensis_CRGBv3.2"
LIB="SCOR_consensus_sequences.fa"

#### Start Script
# Get genome with short seq names and short file name (helps prevent RepeatMasker errors)
run_cmd "sed -e 's/${PREFIX}___//' ${PREFIX}.assembly.fasta > genome.fa"

# Run RepeatMasker
run_cmd "RepeatMasker -no_is -nolow -a -x -gff -pa $PA -lib $LIB genome.fa 1> RepeatMasker.log 2>&1"

# Generate repeat landscape plots.
run_cmd "${UTILS}/calcDivergenceFromAlign.pl -s genome.fa.divsum genome.fa.align"

# Convert short names to original long names.
run_cmd "bedtools maskfasta -fi genome.fa -fo genome.fa.softmasked -bed genome.fa.out.gff -soft"
run_cmd "sed -e 's/^>/>${PREFIX}___/' genome.fa.softmasked > ${PREFIX}.assembly.fasta.repmaskSCORs.softmasked"
run_cmd "sed -e 's/^>/>${PREFIX}___/' genome.fa.masked | reformat.sh in=stdin.fa out=stdout.fa fastawrap=0 > ${PREFIX}.assembly.fasta.repmaskSCORs.masked"
run_cmd "awk -F'\t' -vPREFIX="${PREFIX}" '{ if(\$1!~\"^#\") {print PREFIX\"___\"\$0} else {print \$0} }' genome.fa.out.gff > ${PREFIX}.assembly.fasta.repmaskSCORs.out.gff"

run_cmd "cp genome.fa.divsum ${PREFIX}.assembly.fasta.repmaskSCORs.divsum"
run_cmd "cp genome.fa.tbl ${PREFIX}.assembly.fasta.repmaskSCORs.tbl"

# Kimura distance
run_cmd "tail -n 72 ${PREFIX}.assembly.fasta.repmaskSCORs.divsum | sed -e 's/ $//' > ${PREFIX}.assembly.fasta.repmaskSCORs.divsum.Kimura_distance"

# Compress results files
run_cmd "gzip genome.fa.align genome.fa.divsum genome.fa.masked genome.fa.out genome.fa.out.gff genome.fa.softmasked genome.fa.tbl"


