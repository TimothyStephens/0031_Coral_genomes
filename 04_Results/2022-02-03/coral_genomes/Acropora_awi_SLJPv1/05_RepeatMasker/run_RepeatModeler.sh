#!/usr/bin/env bash

# Print all info to log file
exec 1> "${0}.log.$(date +%s)" 2>&1

#### Pre-run setup
source ~/scripts/script_setup.sh
set +eu; conda activate RepeatMasker; set -eu

export PATH="$PATH:/scratch/ts942/RepeatAnalysis/RepeatMasker-4.1.2-p1"
export PATH="$PATH:/home/ts942/PROGRAMS/bedtools2-v2.29.2/bin"
export PATH="$PATH:/home/ts942/PROGRAMS/bbmap"
module load java/14.0.1

UTILS="/scratch/ts942/RepeatAnalysis/RepeatMasker-4.1.2-p1/util"
PA=5
PREFIX="Acropora_awi_SLJPv1"
REF="${PREFIX}.assembly.fasta"
LIB="${PREFIX}.assembly.fasta.RMDB-families.fa"

#### Start Script
# Get genome with short seq names and short file name (helps prevent RepeatMasker errors)
run_cmd "sed -e 's/${PREFIX}___//' ${REF} > genome.fa"

# Combine RepeatModeler repeat library using all repeats from Dfam_3.3 + RepeatMasker basic library.
run_cmd "cat ${LIB} rmlib.fa > custom_repeat_library.fa"

# Run RepeatMasker
run_cmd "RepeatMasker -no_is -a -x -gff -pa $PA -lib custom_repeat_library.fa genome.fa 1> RepeatMasker.log 2>&1"

# Generate repeat landscape plots.
run_cmd "${UTILS}/calcDivergenceFromAlign.pl -s genome.fa.divsum genome.fa.align"
GENOME_SIZE=$(awk -F'\t' 'NR==4{print $2}' ../01_stats/${PREFIX}.GeneStats.tsv)
echo -e "${PREFIX}\t${GENOME_SIZE}"
run_cmd "${UTILS}/createRepeatLandscape.pl -g ${GENOME_SIZE} -div genome.fa.divsum > genome.fa.RepeatLandscape.html"

# Convert short names to original long names.
run_cmd "bedtools maskfasta -fi genome.fa -fo genome.fa.softmasked -bed genome.fa.out.gff -soft"
run_cmd "sed -e \"s/^>/>${PREFIX}___/\" genome.fa.masked | reformat.sh in=stdin.fa out=stdout.fa fastawrap=0 > ${PREFIX}.assembly.fasta.masked"
run_cmd "sed -e \"s/^>/>${PREFIX}___/\" genome.fa.softmasked > ${PREFIX}.assembly.fasta.softmasked"
run_cmd "awk -F'\t' -vPREFIX=\"${PREFIX}\" '{ if(\$1!~\"^#\") {print PREFIX\"___\"\$0} else {print \$0} }' genome.fa.out.gff > ${PREFIX}.assembly.fasta.out.gff"

run_cmd "cp genome.fa.divsum ${PREFIX}.assembly.fasta.divsum"
run_cmd "cp genome.fa.tbl ${PREFIX}.assembly.fasta.tbl"

# Kimura distance
run_cmd "tail -n 72 ${PREFIX}.assembly.fasta.divsum | sed -e 's/ $//' > ${PREFIX}.assembly.fasta.divsum.Kimura_distance"

# Compress results
run_cmd "gzip genome.fa.align genome.fa.divsum genome.fa.masked genome.fa.out genome.fa.out.gff genome.fa.softmasked genome.fa.tbl"


