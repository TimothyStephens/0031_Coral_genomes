# Common analysis for all genomes being analyzed across different projects

This project is a repository of analysis run on each coral or our group genome that will/can be used across multiple projects (e.g. `busco` or `RepeatMasker` results). For each genome the analysis was run as described in this document unless otherwise stated.

## Setup analysis directory

Setup bash environment.

```bash
conda activate py27
```

## Download and cleanup genome data

The coral genome assemblies and predicted genes are in `../../01_Data/2021-12-03/coral_genomes/`

The coral out group genome assemblies and predicted genes are in `../../01_Data/2021-12-03/outgroup_genomes/`

The `../../01_Data/2021-12-03/README.md ` file explains where each genome was downloaded from and the steps I took to clean and reformat the data.

## 00_databases

### Setup

Make sequence databases to use for subsequent analysis. 

- `samtools` `*.fai` (genome+CDS+PEP)
- `gatk` `*.dict` (genome+CDS)
- `blast` (genome+CDS+PEP)
- `bowtie2` (genome+CDS)

Set up "00_databases" directory

```bash
while read PREFIX; 
do 
  mkdir -p $PREFIX/00_databases; 
( cd $PREFIX/00_databases; 
  ln -s ../../../../../01_Data/2021-12-03/*_genomes/$PREFIX/genome_assembly/$PREFIX.assembly.fasta; 
  ln -s ../../../../../01_Data/2021-12-03/*_genomes/$PREFIX/genome_assembly/$PREFIX.genes.pep.faa; 
  ln -s ../../../../../01_Data/2021-12-03/*_genomes/$PREFIX/genome_assembly/$PREFIX.genes.cds.fna; 
  ln -s ../../../../../01_Data/2021-12-03/*_genomes/$PREFIX/genome_assembly/$PREFIX.genes.gff3; 
  sed -e 's/PREFIX=""/PREFIX="'$PREFIX'"/' ../../run_00-setup_dataset.sh > run_00-setup_dataset.sh;
  chmod +x run_00-setup_dataset.sh
); done < samples_new.txt
```

```bash
while read PREFIX; do echo "$PREFIX"; ( cd $PREFIX/00_databases; ./run_00-setup_dataset.sh ); done < samples_new.txt
```

### Run bash scripts

Run `bowtie2` on the genome and CDS.

```bash
while read PREFIX; 
do 
  mkdir -p $PREFIX/00_databases/bowtie2; 
( cd $PREFIX/00_databases/bowtie2; 
  ln -s ../$PREFIX.assembly.fasta;  
  ln -s ../$PREFIX.genes.cds.fna; 
  sed -e 's/PREFIX=""/PREFIX="'$PREFIX'"/' ../../../run_00-bowtie2.sh > run_00-bowtie2.sh;
  chmod +x run_00-bowtie2.sh
); done < samples_new.txt
```

```bash
while read PREFIX; do echo "$PREFIX"; ( cd $PREFIX/00_databases/bowtie2; ./run_00-bowtie2.sh ); done < samples_new.txt
```

Run `makeblastdb` on the genome, CDS, and PEP files.

```bash
while read PREFIX; 
do 
  mkdir -p $PREFIX/00_databases/blast; 
( cd $PREFIX/00_databases/blast; 
  ln -s ../$PREFIX.assembly.fasta;
  ln -s ../$PREFIX.genes.pep.faa; 
  ln -s ../$PREFIX.genes.cds.fna; 
  sed -e 's/PREFIX=""/PREFIX="'$PREFIX'"/' ../../../run_00-makeblastdb.sh > run_00-makeblastdb.sh;
  chmod +x run_00-makeblastdb.sh
); done < samples_new.txt
```

```bash
while read PREFIX; do echo "$PREFIX"; ( cd $PREFIX/00_databases/blast; ./run_00-makeblastdb.sh ); done < samples_new.txt
```

### Check for errors in log files

```bash
while read PREFIX;
do 
  echo -e "\n## $PREFIX";
  find ${PREFIX}/ -name "*.sh.log.*" \
    | sort \
    | grep '00_data' \
    | while read LOG;
      do 
        if [ $(grep -i 'error\|warn\|fault' $LOG | grep -v 'default' | wc -l) -eq 0 ];
        then 
          echo "OK: $LOG";
        else 
          echo "FAILED: $LOG";
          grep -i 'error\|warn\|fault' $LOG
        fi; 
      done; 
done < samples.txt
```

### Check md5sum

```bash
while read PREFIX;
do 
  echo -e "\n## $PREFIX";
  find ${PREFIX}/ -name "*md5sum_list.txt*" \
    | sort \
    | grep '00_databases' \
    | while read F;
      do
        echo "## $F"
        ( cd $(dirname $F); md5sum --quiet --strict --check <(grep -v 'split' $(basename $F)) ); 
      done
done < samples.txt
```

### Download results

```bash
WD="timothy@coral.rutgers.edu:/scratch/timothy/projects/0031_Coral_genomes/03_Analysis/2022-02-03"
rsync -zarv --delete --prune-empty-dirs --relative \
 --include="*/" \
 --exclude="*/STAR/*.txt" \
 --include="*.job_md5sum_list.txt" \
 --include="*.sh*" --include="*.log" \
 --exclude="*" \
 "${WD}/./*_genomes/*/00_databases/" \
 . --dry-run
```

---

Build specific databases for certain datasets when required.

### HISAT2

Created for:

- Montipora_capitata_KBHIv3

```bash
ln -s ../Montipora_capitata_KBHIv3.assembly.fasta
/home/timothy/programs/hisat2-2.2.1/hisat2-build Montipora_capitata_KBHIv3.assembly.fasta Montipora_capitata_KBHIv3.assembly.fasta 1>hisat2-build.log 2>&1
```

### Salmon

Created for:

- Montipora_capitata_KBHIv3

```bash
ln -s ../Montipora_capitata_KBHIv3.genes.cds.fna
./run_salmon_index.sh
```





## 01_stats

Generate some basic stats about the genome and gene models using standard tools such as `stats.sh` (from `bbmap`) or `bedtools`.

### Setup

Run STATS on the genome and gene gff file.

```bash
while read PREFIX; 
do 
  mkdir -p $PREFIX/01_stats; 
( cd $PREFIX/01_stats; 
  ln -s ../00_databases/$PREFIX.assembly.fasta;
  ln -s ../00_databases/$PREFIX.genes.gff3;
  cp ../../../../../02_Scripts/genome_stats .
  sed -e 's/PREFIX=""/PREFIX="'$PREFIX'"/' ../../run_01-stats.sh > run_01-stats.sh;
  chmod +x *.sh genome_stats
); done < samples_new.txt
```

### Run bash scripts

```bash
while read PREFIX; do echo "$PREFIX"; ( cd $PREFIX/01_stats; ./run_01-stats.sh ); done < samples_new.txt
```

### Check for errors in log files

```bash
while read PREFIX;
do 
  echo -e "\n## $PREFIX";
  find ${PREFIX}/ -name "*.sh.log.*" \
    | sort \
    | grep '01_stats' \
    | while read LOG;
      do 
        if [ $(grep -i 'error\|warn\|fault' $LOG | wc -l) -eq 0 ];
        then 
          echo "OK: $LOG";
        else 
          echo "FAILED: $LOG";
          grep -i 'error\|warn\|fault' $LOG
        fi; 
      done; 
done < samples.txt
```

The only genome that had errors in its log file was `Montipora_capitata_WTHIv1.1` which was because some of the gene models were from scaffolds that had been removed from the genome during filtering (by the authors of the genome paper). Not a big problem as this only affects 1 gene.

### Check md5sum

```bash
while read PREFIX;
do 
  echo -e "\n## $PREFIX";
  find ${PREFIX}/ -name "*md5sum_list.txt*" \
    | sort \
    | grep '01_stats' \
    | while read F;
      do
        echo "## $F"
        ( cd $(dirname $F); md5sum --quiet --strict --check <(grep -v 'split' $(basename $F)) ); 
      done
done < samples.txt
```

### Download results

```bash
WD="timothy@coral.rutgers.edu:/scratch/timothy/projects/0031_Coral_genomes/03_Analysis/2022-02-03"
rsync -zarv --delete --prune-empty-dirs --relative \
 --include="*/" \
 --include="*.job_md5sum_list.txt" \
 --include="*.sh*" --include="*.log" \
 --include="*.GeneStats.tsv" \
 --exclude="*" \
 "${WD}/./*_genomes/*/01_stats/" \
 . --dry-run
```

## 02_busco

Run `bunco` on the genome and predicted proteins from each isolate. Run using the `metazoa_odb10` (2021-02-24) and `eukaryota_odb10` (2020-09-10) lineages. See https://busco-data.ezlab.org/v5/data/lineages/ for the list of available lineages (version 5).

### Setup

Run BUSCO on the genome and predicted proteins (run each as a separate job/script).

```bash
while read PREFIX; 
do 
  mkdir -p $PREFIX/02_busco; 
( cd $PREFIX/02_busco; 
  ln -s ../00_databases/$PREFIX.assembly.fasta;
  ln -s ../00_databases/$PREFIX.genes.pep.faa;  
  sed -e 's/PREFIX=""/PREFIX="'$PREFIX'"/' ../../run_02-busco_genome.sh > run_02-busco_genome.sh;
  sed -e 's/PREFIX=""/PREFIX="'$PREFIX'"/' ../../run_02-busco_protein.sh > run_02-busco_protein.sh;
  chmod +x *.sh
); done < samples_new.txt
```

### Run bash script

```bash
while read PREFIX; do echo "$PREFIX"; ( cd $PREFIX/02_busco; ./run_02-busco_genome.sh; ./run_02-busco_protein.sh ); done < samples_new.txt
```

Generate BUSCO summary results file for each genome.

```bash
while read PREFIX; 
do 
  for F in $PREFIX/02_busco/pep.faa*_odb10;
  do
    ./../../../02_Scripts/parse_BUSCO_short_summary "${F}"/short_summary.* > "${F}".results.txt
  done
  for F in $PREFIX/02_busco/genome.fa*_odb10;
  do
    ./../../../02_Scripts/parse_BUSCO_short_summary "${F}"/short_summary.* > "${F}".results.txt
  done
done < samples_new.txt
```

### Check for errors in log files

```bash
while read PREFIX;
do 
  echo -e "\n## $PREFIX";
  find ${PREFIX}/ -name "*.sh.log.*" \
    | sort \
    | grep '02_busco' \
    | while read LOG;
      do 
        if [ $(grep -i 'error\|warn\|fault' $LOG | wc -l) -eq 0 ];
        then 
          echo "OK: $LOG";
        else 
          echo "FAILED: $LOG";
          grep -i 'error\|warn\|fault' $LOG
        fi; 
      done; 
done < samples.txt
```

### Cleanup results files

```bash
while read PREFIX;
do
  echo ${PREFIX};
  rm -fr ${PREFIX}/02_busco/busco_downloads/
  (cd ${PREFIX}/02_busco/; 
   D="genome.fa.busco_eukaryota_odb10"; tar -zcf ${D}.tar.gz ${D} && rm -rf ${D};
   D="genome.fa.busco_metazoa_odb10";   tar -zcf ${D}.tar.gz ${D} && rm -rf ${D};
   D="pep.faa.busco_eukaryota_odb10";   tar -zcf ${D}.tar.gz ${D} && rm -rf ${D};
   D="pep.faa.busco_metazoa_odb10";     tar -zcf ${D}.tar.gz ${D} && rm -rf ${D};
  )
done < samples.txt
```

### Check md5sum

```bash
while read PREFIX;
do 
  echo -e "\n## $PREFIX";
  find ${PREFIX}/ -name "*md5sum_list.txt*" \
    | sort \
    | grep '02_busco' \
    | while read F;
      do
        echo "## $F"
        ( cd $(dirname $F); md5sum --quiet --strict --check <(grep -v 'split' $(basename $F)) ); 
      done
done < samples.txt
```

### Download results

```bash
WD="timothy@coral.rutgers.edu:/scratch/timothy/projects/0031_Coral_genomes/03_Analysis/2022-02-03"
rsync -zarv --delete --prune-empty-dirs --relative \
 --include="*/" \
 --exclude="*.headersMap.tsv" --exclude="*_odb10/logs/*" \
 --include="*.job_md5sum_list.txt" \
 --include="*.sh*" --include="*.log" \
 --include="*.txt" --include="*.tsv" \
 --exclude="*" \
 "${WD}/./*_genomes/*/02_busco/" \
 . --dry-run
```

## 03_sequence_contamination

Check to see how contaminated the genomes are with bacteria (using RefSeq Bacteria genomic v211) and Symbiodiniaceae (available genomes; see `/scratch/timothy/databases/Symbiodiniaceae/20220609/`). Filter results, retaining hits with bitscore > 1000 and *e*-value <1e-10.

### Setup

`blastn` of the genomes against RefSeq bacterial genomes and Symbiont genomes (run each as a separate job/script)

```bash
while read PREFIX;
do 
  mkdir -p $PREFIX/03_sequence_contamination; 
( cd $PREFIX/03_sequence_contamination; 
  ln -s ../00_databases/$PREFIX.assembly.fasta;
  ln -s ../00_databases/$PREFIX.assembly.fasta.fai;  
  sed -e 's/PREFIX=""/PREFIX="'$PREFIX'"/' ../../run_03-genome_blastn_bacteria.sh > run_03-genome_blastn_bacteria.sh;
  sed -e 's/PREFIX=""/PREFIX="'$PREFIX'"/' ../../run_03-genome_blastn_Symbiodiniaceae.sh > run_03-genome_blastn_Symbiodiniaceae.sh;
  chmod +x *.sh
); done < samples_new.txt
```

### Run bash scripts

```bash
while read PREFIX; do echo "$PREFIX"; ( cd $PREFIX/03_sequence_contamination; ./run_03-genome_blastn_bacteria.sh ); done < samples_new.txt
```

```bash
while read PREFIX; do echo "$PREFIX"; ( cd $PREFIX/03_sequence_contamination; ./run_03-genome_blastn_Symbiodiniaceae.sh ); done < samples_new.txt
```

Columns in output file (A: Genome scaffold features. B: BLAST hits):

| Column No. | Column Description                                           |
| ---------- | ------------------------------------------------------------ |
| 1          | Scaffold ID                                                  |
| 2          | Start position (always 0 in our case)                        |
| 3          | End position (always length of scaffold in our case)         |
| 4          | The number of features in B that overlapped (by at least one base pair) the A interval. |
| 5          | The number of bases in A that had non-zero coverage from features in B. |
| 6          | The length of the entry in A.                                |
| 7          | The fraction of bases in A that had non-zero coverage from features in B. |

Consider scaffolds with >10% hist coverage as contaminants. Same threshold as (Liu, H., Stephens, T.G., González-Pech, R.A. *et al.* *Symbiodinium* genomes reveal adaptive evolution of functions related to coral-dinoflagellate symbiosis. *Commun Biol* **1,** 95 (2018). https://doi.org/10.1038/s42003-018-0098-3) except they used an *e*-value cutoff of 1e-20 not 1e-10 like here. 

```bash
while read PREFIX;
do
  for F in ${PREFIX}/03_sequence_contamination/*.coverage;
  do
  echo "${PREFIX} --- ${F}"
  awk '
    BEGIN {
      TOTAL_COV=0; 
      TOTAL_COUNT=0; 
      TOTAL_LENGTH=0; 
      FILTERED_COV=0; 
      FILTERED_COUNT=0; 
      FILTERED_LENGTH=0
    } {
      TOTAL_COUNT++; 
      TOTAL_LENGTH=TOTAL_LENGTH+$6;
      TOTAL_COV=TOTAL_COV+$5; 
      if($7>0.1) {
        FILTERED_COUNT++; 
        FILTERED_LENGTH=FILTERED_LENGTH+$6
        FILTERED_COV=FILTERED_COV+$5; 
      }
    } END {
      print "Total number of Scaffolds\tTotal scaffold length (bp)\tTotal contaminant coverage (bp)\tNumber of filtered scaffolds\tLength of filtered scaffolds (bp)\tContaminant coverage of filtered scaffolds (bp)\tPercent number of filtered scaffolds\tPercent length of filtered scaffolds\tPercent contaminant coverage of filtered scaffolds"; 
      print TOTAL_COUNT"\t"TOTAL_LENGTH"\t"TOTAL_COV"\t"FILTERED_COUNT"\t"FILTERED_LENGTH"\t"FILTERED_COV"\t"(FILTERED_COUNT/TOTAL_COUNT)*100"\t"(FILTERED_LENGTH/TOTAL_LENGTH)*100"\t"(FILTERED_COV/TOTAL_LENGTH)*100}' \
    "$F" > "$F.stats_filtered_scaffolds"
  done
done < samples_new.txt
```

### Check for errors in log files

```bash
while read PREFIX;
do 
  echo -e "\n## $PREFIX";
  find ${PREFIX}/ -name "*.sh.log.*" \
    | sort \
    | grep '03_sequence_contamination' \
    | while read LOG;
      do 
        if [ $(grep -i 'error\|warn\|fault' $LOG | wc -l) -eq 0 ];
        then 
          echo "OK: $LOG";
        else 
          echo "FAILED: $LOG";
          grep -i 'error\|warn\|fault' $LOG
        fi; 
      done; 
done < samples.txt
```

### Check md5sum

```bash
while read PREFIX;
do 
  echo -e "\n## $PREFIX";
  find ${PREFIX}/ -name "*md5sum_list.txt*" \
    | sort \
    | grep '03_sequence_contamination' \
    | while read F;
      do
        echo "## $F"
        ( cd $(dirname $F); md5sum --quiet --strict --check <(grep -v 'split' $(basename $F)) ); 
      done
done < samples.txt
```

### Download results

```bash
WD="timothy@coral.rutgers.edu:/scratch/timothy/projects/0031_Coral_genomes/03_Analysis/2022-02-03"
rsync -zarv --delete --prune-empty-dirs --relative \
 --include="*/" \
 --include="*.job_md5sum_list.txt" \
 --include="*.sh*" --include="*.log" \
 --include="*.coverage.stats_filtered_scaffolds" \
 --exclude="*" \
 "${WD}/./*_genomes/*/03_sequence_contamination/" \
 . --dry-run
```

## 04_RepeatModeler

Run `RepeatModeler` to generate a set of *de novo* predicted repeats for each genome.

### Setup

```bash
while read PREFIX; 
do 
  mkdir -p $PREFIX/04_RepeatModeler; 
( cd $PREFIX/04_RepeatModeler; 
  ln -s ../00_databases/$PREFIX.assembly.fasta;
  sed -e 's/PREFIX=""/PREFIX="'$PREFIX'"/' ../../run_04-RepeatModeler.sh > run_RepeatModeler.sh;
  cp ../../run_04-AmarelJob.sh run_AmarelJob.sh;
  chmod +x *.sh
); done < samples_new.txt
```

### Run bash scripts

```bash
## Copy files to Amarel server and 'sbatch' each run_AmarelJob.sh script to submit to queue
rsync -avPL --relative timothy@coral.rutgers.edu:/scratch/timothy/projects/0031_Coral_genomes/03_Analysis/2022-02-03/./*_genomes/*/00_databases/*.assembly.fasta .

rsync -avzP --relative timothy@coral.rutgers.edu:/scratch/timothy/projects/0031_Coral_genomes/03_Analysis/2022-02-03/./*_genomes/*/01_stats/*.GeneStats.tsv .

rsync -avzP --relative timothy@coral.rutgers.edu:/scratch/timothy/projects/0031_Coral_genomes/03_Analysis/2022-02-03/./*_genomes/*/04_RepeatModeler .

sbatch run_AmarelJob.sh
```

### Check for errors in log files

```bash
while read PREFIX;
do 
  echo -e "\n## $PREFIX";
  find ${PREFIX}/ -name "*.sh.log.*" \
    | sort \
    | grep '04_RepeatModeler' \
    | while read LOG;
      do 
        if [ $(grep -i 'error\|warn\|fault' $LOG | grep -v ' Your shell has not been properly configured to use' | wc -l) -eq 0 ];
        then 
          echo "OK: $LOG";
        else 
          echo "FAILED: $LOG";
          grep -i 'error\|warn\|fault' $LOG
        fi; 
      done; 
done < samples.txt

# Target the BuildDatabase.log and RepeatModeler.log log files
while read PREFIX;
do 
  echo -e "\n## $PREFIX";
  find ${PREFIX}/ -name "*log" \
    | sort \
    | grep '04_RepeatModeler' \
    | while read LOG;
      do 
        if [ $(grep -i 'error\|warn\|fault' $LOG | wc -l) -eq 0 ];
        then 
          echo "OK: $LOG";
        else 
          echo "FAILED: $LOG";
          grep -i 'error\|warn\|fault' $LOG
        fi; 
      done; 
done < samples.txt
```

Ignore the "RECON returned a negative offset:" and "Attempt to extract substring outside the range of the sequence" warnings. These appear in a few of the assemblies but RepeatModeler seems to finish just fine.

### Check md5sum

```bash
while read PREFIX;
do
  N1=$(cat ${PREFIX}/04_RepeatModeler/${PREFIX}.assembly.fasta \
         | grep -v '>' | md5sum | awk '{print $1}'); 
  N2=$(/scratch/databases/bin/ncbi-blast-2.13.0+/bin/get_sequences_from_blastDB \
           ${PREFIX}/04_RepeatModeler/${PREFIX}.assembly.fasta.RMDB \
         | ~/programs/bbmap/reformat.sh in=stdin.fa out=stdout.fa fastawrap=0 2>/dev/null \
         | grep -v '>' | md5sum | awk '{print $1}'); 
  if [ "$N1" = "$N2" ];
  then
    echo "OK: ${PREFIX}"
  else
    echo "FAILED: ${PREFIX}"
  fi
done < samples.txt
```

### Download results

```bash
WD="timothy@coral.rutgers.edu:/scratch/timothy/projects/0031_Coral_genomes/03_Analysis/2022-02-03"
rsync -zarv --delete --prune-empty-dirs --relative \
 --include="*/" \
 --include="*.job_md5sum_list.txt" \
 --include="*.sh*" --include="*.log" \
 --include="*-families.fa" \
 --include="AmarelJob*" \
 --exclude="*" \
 "${WD}/./*_genomes/*/04_RepeatModeler/" \
 . --dry-run
```

## 05_RepeatMasker

Run `RepeatMasker` to identify and mask predicted repeats for each genome using the *de novo* predicted repeats and the repeats from the Dfam library. 

### Setup

```bash
while read PREFIX; 
do 
  mkdir -p $PREFIX/05_RepeatMasker; 
( cd $PREFIX/05_RepeatMasker; 
  ln -s ../00_databases/$PREFIX.assembly.fasta;
  ln -s ../05_RepeatMasker/${PREFIX}.assembly.fasta.RMDB-families.fa;
  ln -s ../../../../01_Data/2021-12-03/repeat_library/rmlib.fa;
  sed -e 's/PREFIX=""/PREFIX="'$PREFIX'"/' ../../run_05-RepeatMasker.sh > run_RepeatModeler.sh;
  cp ../../run_05-AmarelJob.sh run_AmarelJob.sh;
  chmod +x *.sh
); done < samples_new.txt
```

### Run bash scripts

```bash
## Copy files to Amarel server and 'sbatch' each run_AmarelJob.sh script to submit to queue
rsync -avzP --relative timothy@coral.rutgers.edu:/scratch/timothy/projects/0031_Coral_genomes/03_Analysis/2022-02-03/./*_genomes/*/05_RepeatMasker .

sbatch run_AmarelJob.sh
```

### Extract Kimura% values

```bash
while read PREFIX; 
do 
( cd $PREFIX/05_RepeatMasker; 
  awk -F'\t' '{ if($1!~"^Jukes" && $1!~"^====" && $1!~"^File" && $1!~"^Weighted" && $1!~"^----" && $0!=""){ if($1~"^Coverage"){exit}else{print} } }' \
      "${PREFIX}.assembly.fasta.divsum" \
    > "${PREFIX}.assembly.fasta.divsum.Weighted_average"
); done < samples_new.txt
```

### Check for errors in log files

```bash
while read PREFIX;
do 
  echo -e "\n## $PREFIX";
  find ${PREFIX}/ -name "*.sh.log.*" \
    | sort \
    | grep '05_RepeatMasker' \
    | while read LOG;
      do 
        if [ $(grep -i 'error\|warn\|fault' $LOG | grep -v 'Your shell has not been properly configured to use' | wc -l) -eq 0 ];
        then 
          echo "OK: $LOG";
        else 
          echo "FAILED: $LOG";
          grep -i 'error\|warn\|fault' $LOG
        fi; 
      done; 
done < samples.txt

# Target the RepeatMasker.log log file
while read PREFIX;
do 
  echo -e "\n## $PREFIX";
  find ${PREFIX}/ -name "*log" \
    | sort \
    | grep '05_RepeatMasker' \
    | while read LOG;
      do 
        if [ $(grep -i 'error\|warn\|fault' $LOG | wc -l) -eq 0 ];
        then 
          echo "OK: $LOG";
        else 
          echo "FAILED: $LOG";
          grep -i 'error\|warn\|fault' $LOG
        fi; 
      done; 
done < samples.txt
```

### Check md5sum

```bash
while read PREFIX;
do
  N1=$(cat ${PREFIX}/05_RepeatMasker/${PREFIX}.assembly.fasta \
         | grep -v '>' | md5sum | awk '{print $1}'); 
  N2=$(cat ${PREFIX}/05_RepeatMasker/genome.fa \
         | grep -v '>' | md5sum | awk '{print $1}');
  if [ "$N1" = "$N2" ];
  then
    echo "OK: ${PREFIX}"
  else
    echo "FAILED: ${PREFIX}"
  fi
done < samples.txt
```

### Download results

```bash
WD="timothy@coral.rutgers.edu:/scratch/timothy/projects/0031_Coral_genomes/03_Analysis/2022-02-03"
rsync -zarv --delete --prune-empty-dirs --relative \
 --include="*/" \
 --include="*.job_md5sum_list.txt" \
 --include="*.sh*" --include="*.log" \
 --include="*.tbl" --include="*.html" \
 --include="*.divsum" \
 --include="*.divsum.Kimura_distance" \
 --include="*.divsum.Weighted_average" \
 --include="AmarelJob*" \
 --exclude="*" \
 "${WD}/./*_genomes/*/05_RepeatMasker/" \
 . --dry-run
```

## 06_RepeatMasker_SCORs

Run `RepeatMasker` to identify and mask predicted repeats for each genome using just the Scleractinia COral-specific Repeat families (SCORs) that were published previously (https://doi.org/10.1016/j.ygeno.2017.06.003).  

### Setup

```bash
while read PREFIX; 
do 
  mkdir -p $PREFIX/06_RepeatMasker_SCORs; 
( cd $PREFIX/06_RepeatMasker_SCORs; 
  ln -s ../00_databases/$PREFIX.assembly.fasta;
  ln -s ../../../../../01_Data/2021-12-03/SCORs/SCOR_consensus_sequences.fa;
  sed -e 's/PREFIX=""/PREFIX="'$PREFIX'"/' ../../run_06-RepeatMasker_SCORs.sh > run_RepeatModeler_SCORs.sh;
  chmod +x *.sh
); done < samples_new.txt
```

### Run bash scripts

```bash
while read PREFIX; do echo "$PREFIX"; ( cd $PREFIX/06_RepeatMasker_SCORs; ./run_RepeatModeler_SCORs.sh ); done < samples_new.txt
```

### Extract Kimura% values

```bash
while read PREFIX; 
do 
( cd $PREFIX/06_RepeatMasker_SCORs; 
  awk -F'\t' '{ if($1!~"^Jukes" && $1!~"^====" && $1!~"^File" && $1!~"^Weighted" && $1!~"^----" && $0!=""){ if($1~"^Coverage"){exit}else{print} } }' \
      "${PREFIX}.assembly.fasta.repmaskSCORs.divsum" \
    > "${PREFIX}.assembly.fasta.repmaskSCORs.divsum.Weighted_average"
); done < samples_new.txt
```

### Check for errors in log files

```bash
while read PREFIX;
do 
  echo -e "\n## $PREFIX";
  find ${PREFIX}/ -name "*.sh.log.*" \
    | sort \
    | grep '06_RepeatMasker_SCORs' \
    | while read LOG;
      do 
        if [ $(grep -i 'error\|warn\|fault' $LOG | wc -l) -eq 0 ];
        then 
          echo "OK: $LOG";
        else 
          echo "FAILED: $LOG";
          grep -i 'error\|warn\|fault' $LOG
        fi; 
      done; 
done < samples.txt

# Target the RepeatMasker.log log file
while read PREFIX;
do 
  echo -e "\n## $PREFIX";
  find ${PREFIX}/ -name "*log" \
    | sort \
    | grep '06_RepeatMasker_SCORs' \
    | while read LOG;
      do 
        if [ $(grep -i 'error\|warn\|fault' $LOG | wc -l) -eq 0 ];
        then 
          echo "OK: $LOG";
        else 
          echo "FAILED: $LOG";
          grep -i 'error\|warn\|fault' $LOG
        fi; 
      done; 
done < samples.txt
```

### Check md5sum

```bash
while read PREFIX;
do
  N1=$(cat ${PREFIX}/06_RepeatMasker_SCORs/${PREFIX}.assembly.fasta \
         | grep -v '>' | md5sum | awk '{print $1}'); 
  N2=$(cat ${PREFIX}/06_RepeatMasker_SCORs/genome.fa \
         | grep -v '>' | md5sum | awk '{print $1}');
  if [ "$N1" = "$N2" ];
  then
    echo "OK: ${PREFIX}"
  else
    echo "FAILED: ${PREFIX}"
  fi
done < samples.txt
```

### Download results

```bash
WD="timothy@coral.rutgers.edu:/scratch/timothy/projects/0031_Coral_genomes/03_Analysis/2022-02-03"
rsync -zarv --delete --prune-empty-dirs --relative \
 --include="*/" \
 --include="*.sh*" --include="*.log" \
 --include="*.tbl" \
 --include="*.divsum" \
 --include="*.divsum.Kimura_distance" \
 --include="*.divsum.Weighted_average" \
 --exclude="*" \
 "${WD}/./*_genomes/*/06_RepeatMasker_SCORs/" \
 . --dry-run
```

## 07_short_read_trimming

Run `cutadapt` on the available Illumina genome sequencing reads. Execute the `run_Cutadapt.sh` script to runthe trimming process and to double check the adapter content using my custom `find_adapters_in_reads.py` script. 

```bash
ln -s ../../../../02_Scripts/adapters.fa
ln -s ../../../../02_Scripts/find_adapters_in_reads.py
```

Run `cutadapt`.

```bash
./run_Cutadapt.sh
```

### Download results

```bash
WD="timothy@coral.rutgers.edu:/scratch/timothy/projects/0031_Coral_genomes/03_Analysis/2022-02-03"
rsync -zarv --delete --prune-empty-dirs --relative \
 --include="*/" \
 --include="*.sh*" --include="*.log" \
 --exclude="*" \
 ${WD}/./*/*/07_short_read_trimming/ \
 . --dry-run
```

## 08_genome_size_estimation

Run `smudgeplot` using the trimmed reads to identify the ploidy of each sample. Execute the `run_smudgeplot.sh` script to run `jellyfish` (k=21) and `smudgeplot` (using automatically chosen cutoffs). After this is done running, download the `smudgeplot.jellyfish.21.histo`file and run the online version of `GenomeScope2` to estimate the genome size, ploidy, and cutoffs. The latter of which we will use for a second round of `smudgeplot` analysis; manually change the cutoffs in `run_smudgeplot_manual_cutoffs.sh`before executing.

After the first script has finished it is safe to remove the `smudgeplot.jellyfish.21.jf` file (only used by `jellyfish` for generation of the `*.histo` file).

```bash
rm */08_genome_size_estimation/smudgeplot.jellyfish.21.jf
```

### Download results

```bash
WD="timothy@coral.rutgers.edu:/scratch/timothy/projects/0031_Coral_genomes/03_Analysis/2022-02-03"
rsync -zarv --prune-empty-dirs --relative \
 --include="*/" \
 --include="*.sh*" --include="*.log" \
 --include="*.histo" \
 --include="*.png" \
 --exclude="*" \
 ${WD}/./*/*/08_genome_size_estimation/ \
 . --dry-run
```

## 09_gene_prediction

























## 10_functional_annotation

Run functional annotation (DIAMOND blastp, eggnog-mapper, interproscan) on the predicted proteins.

### Setup

Split protein file into 100 parts. This is for parallelization of the BLASTP against RefSeq complete analysis (not performed anymore). 

```bash
while read PREFIX; 
do
  mkdir -p $PREFIX/10_functional_annotation; 
( cd $PREFIX/10_functional_annotation; 
  sed -e 's/*/X/g' ../00_databases/$PREFIX.genes.pep.faa > $PREFIX.genes.pep.faa;
  ./../../../../../02_Scripts/fasta-splitter.pl --line-length 0 --n-parts 100 --out-dir $PREFIX.genes.pep.faa.split $PREFIX.genes.pep.faa;
  md5sum ../00_databases/$PREFIX.genes.pep.faa $PREFIX.genes.pep.faa $PREFIX.genes.pep.faa.split/*.faa > md5sum_list.txt
); done < samples_new.txt
```

Setup scripts.

```bash
while read PREFIX; 
do 
  mkdir -p $PREFIX/10_functional_annotation; 
( cd $PREFIX/10_functional_annotation;
  sed -e 's/PREFIX=""/PREFIX="'$PREFIX'"/' ../../run_10-eggnog-mapper.sh > run_10-eggnog-mapper.sh;
  sed -e 's/PREFIX=""/PREFIX="'$PREFIX'"/' ../../run_10-InterProScan.sh > run_10-InterProScan.sh;
  sed -e 's/PREFIX=""/PREFIX="'$PREFIX'"/' ../../run_10-diamond_nr-Amarel.sh > run_10-diamond_nr-Amarel.sh;
  chmod +x *.sh
); done < samples_new.txt
```

##  `eggnog-mapper` and `Interproscan`

### Run bash scripts

Run `eggnog-mapper` and `Interproscan` on the Coral server.

```bash
while read PREFIX; do echo "$PREFIX"; ( cd $PREFIX/10_functional_annotation; ./run_10-eggnog-mapper.sh; ./run_10-InterProScan.sh ); done < samples_new.txt
```

### Check for errors in log files

```bash
while read PREFIX;
do 
  echo -e "\n## $PREFIX";
  find ${PREFIX}/ -name "*.sh.log.*" \
    | sort \
    | grep '10_functional_annotation' \
    | while read LOG;
      do 
        if [ $(grep -i 'error\|warn\|fault' $LOG | grep -v 'this can take a few minutes and load up to 45GB to RAM' | wc -l) -eq 0 ];
        then 
          echo "OK: $LOG";
        else 
          echo "FAILED: $LOG";
          grep -i 'error\|warn\|fault' $LOG
        fi; 
      done; 
done < samples.txt
```

Cleanup large results file that we won't be using but take up 100's of Gigabytes.

```bash
rm */10_functional_annotation/*.InterProScan.{xml,json}
```

## `diamond blastp` against NCBI nr

Copy files to Amarel.

```bash
rsync -avzPL --relative timothy@coral.rutgers.edu:/scratch/timothy/projects/0031_Coral_genomes/03_Analysis/2022-02-03/./*_genomes/*/00_databases/*.faa . --dry-run
rsync -avzPL --relative timothy@coral.rutgers.edu:/scratch/timothy/projects/0031_Coral_genomes/03_Analysis/2022-02-03/./*_genomes/*/10_functional_annotation/run_10-diamond_nr-Amarel.sh . --dry-run
rsync -avzPL --relative timothy@coral.rutgers.edu:/scratch/timothy/projects/0031_Coral_genomes/03_Analysis/2022-02-03/./*_genomes/*/10_functional_annotation/*.faa . --dry-run
```

Run `diamond blastp` against NCBI's nr database on the Amarel server.

```bash
while read PREFIX; do echo "$PREFIX"; ( cd $PREFIX/10_functional_annotation; sbatch run_10-diamond_nr-Amarel.sh ); done < samples_new.txt
```

Copy files back to the Coral server from Amarel.

```bash
rsync -avzP --relative \
  $PWD/./*_genomes/*/10_functional_annotation \
  timothy@coral.rutgers.edu:/scratch/timothy/projects/0031_Coral_genomes/03_Analysis/2022-02-03/ \
  --dry-run
```

### Check for errors in log files

```bash
while read PREFIX;
do 
  echo -e "\n## $PREFIX";
  tail -n 1 $PREFIX/10_functional_annotation/run_10-diamond_nr-Amarel.sh.*
  ls -1 $PREFIX/10_functional_annotation/run_10-diamond_nr-Amarel.sh.* \
    | sort \
    | while read LOG;
      do 
        if [ $(grep -i 'error\|warn\|fault' $LOG | wc -l) -eq 0 ];
        then 
          echo "OK: $LOG";
        else 
          echo "FAILED: $LOG";
          grep -i 'error\|warn\|fault' $LOG
        fi; 
      done; 
done < samples.txt
```

### Chellck md5sum

```bash
while read PREFIX;
do 
  echo -e "\n## $PREFIX";
  find ${PREFIX}/ -name "*md5sum_list.txt*" \
    | sort \
    | grep '10_functional_annotation' \
    | while read F;
      do
        echo "## $F"
        ( cd $(dirname $F); md5sum --quiet --strict --check <(grep -v 'split' $(basename $F)) ); 
      done
done < samples.txt
```





```bash
while read PREFIX;
do 
  ls ${PREFIX}/10_functional_annotation/${PREFIX}*.diamond_blastp_nr.outfmt6_short.gz
done < samples.txt
```

### Download results

```bash
WD="timothy@coral.rutgers.edu:/scratch/timothy/projects/0031_Coral_genomes/03_Analysis/2022-02-03"
rsync -zarv --delete --prune-empty-dirs --relative \
 --include="*/" \
 --include="*.job_md5sum_list.txt" \
 --include="*.sh*" --include="*.log" \
 --exclude="*" \
 "${WD}/./*_genomes/*/10_functional_annotation/" \
 . --dry-run
```













## 11_read_mapping

Map reads to genome assemblies and gene model CDSs.

### RNA-Seq to genome using `HISAT2` and `StringTie2`

Directory: `genome_assembly/RNA-seq_PRJNA694677_HISAT2-StringTie2/`

```bash
ln_loop ../../../00_databases/hisat2/Montipora_capitata_KBHIv3.assembly.fasta*
```

```bash
./run_hisat2.sh
```

### DNA-Seq to genome using `Bowtie2`

Directory: `genome_assembly/DNA-seq_PRJNA509219_Illumina_bowtie2/`

```bash
ln_loop ../../../00_databases/bowtie2/Montipora_capitata_KBHIv3.assembly.fasta*
```

```bash
./run_bowtie2.sh
```

### DNA-Seq to genome using `Minimap2`

Directory: `genome_assembly/DNA-seq_PRJNA509219_PacBioRSII_minimap2/`

```bash
ln -s ../../../00_databases/Montipora_capitata_KBHIv3.assembly.fasta
```

```bash
./run_minimap2.sh
```

### RNA-Seq to gene model CDS using `Salmon`

Directory: `gene_CDS/RNA-seq_*_salmon/`

```bash
ln_loop ../../../00_databases/salmon/Montipora_capitata_KBHIv3.genes.cds.fna*
```

```bash
./run_salmon_quant.sh
```

### 

### 

### 



















