# Download *Montipora* genome data

Download, unpack, cleanup names, and build databases for the avalible *Montipora* genomes.

Raw genome and read data is placed in `01_Data/2021-12-03/${GENOME_NAME}/` and formatted reference data is in `03_Analysis/2021-12-03/${GENOME_NAME}/00_databases`

## SCORs

From: https://doi.org/10.1016/j.ygeno.2017.06.003

```bash
wget https://ars.els-cdn.com/content/image/1-s2.0-S0888754317300447-mmc3.txt
awk '$1!~"^#" && $1!="" { if($1~">"){split($0,a,"#"); sub(">Mcap.","",a[1]); sub("\\.","-",a[1]); print $1"#SCOR/"a[1]} else {print} }' 1-s2.0-S0888754317300447-mmc3.txt > SCOR_consensus_sequences.fa
```

Extract all repeats in Dfam.h5 (version 3.3). 

```bash
../famdb.py -i RepeatMaskerLib.h5 families --format fasta_name --include-class-in-name --ancestors --descendants 'root' > Dfam_3.3_lib.fa
```

Combine the extracted Dfam and basic RepeatMasker repeat libraries together. Will use to supplement the RepeatModeler repeats derived from each genome.

```bash
LIB="/home/timothy/programs/RepeatAnalysis/RepeatMasker-4.1.2-p1/Libraries"
cat "${LIB}/Dfam_3.3_lib.fa" "${LIB}/RepeatMasker.lib" > rmlib.fa
```

There are a total of 280,646 sequences in the `rmlib.fa` file.
