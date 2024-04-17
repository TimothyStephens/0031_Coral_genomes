


PREFIX="Calvadosia_cruxmelitensis_CRGBv3.2"

cat ${PREFIX}.GeneStats.genome.tsv > ${PREFIX}.GeneStats.tsv

echo -e "
# Gene stats (combined CDS+introns)
Number of genes\t0
Total gene length\t0
Average gene length (bp)\t0
Gene percent GC\t0

# Transcript stats (based on CDS features)
Average transcript length (bp)\t0
Average number of CDS per gene/transcript\t0

# CDS stats" >> ${PREFIX}.GeneStats.tsv

awk 'NR>=5 && NR <= 7' ${PREFIX}.GeneStats.CDS.tsv >> ${PREFIX}.GeneStats.tsv

echo -e "Number of single-CDS transcripts\t0
Percent single-CDS transcripts\t0" >> ${PREFIX}.GeneStats.tsv

awk 'NR==9' ${PREFIX}.GeneStats.CDS.tsv >> ${PREFIX}.GeneStats.tsv

echo -e "
# Intron stats (predicted between CDS features)
Number of introns\t0
Total intron length (bp)\t0
Average intron length (bp)\t0
Intron percent GC\t0

# Intergenic stats (predicted between genes built from CDS features)
Number of intergenic regions\t0
Total intergenic region length (bp)\t0
Average intergenic region length (bp)\t0
Intergenic region percent GC\t0" >> ${PREFIX}.GeneStats.tsv



# Manually change "Number of genes" to be the value of "No. CDS"



