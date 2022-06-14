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

Run blast on the genome, CDS, and PEP files.

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

**Check for errors in log files**

```bash
while read PREFIX;
do 
  echo -e "\n## $PREFIX";
  find ${PREFIX}/ -name "*.sh.log.*" \
    | sort \
    | grep '00_data' \
    | while read LOG;
      do 
        if ! $(grep -qi 'error\|warn\|fault' $LOG);
        then 
          echo "OK: $LOG";
        else 
          echo "FAILED: $LOG";
          grep -i 'error\|warn\|fault' $LOG
        fi; 
      done; 
done < samples_new.txt
```

**Download results**

```bash
WD="timothy@coral.rutgers.edu:/scratch/timothy/projects/0031_Coral_genomes/03_Analysis/2022-02-03"
rsync -zarv --delete --prune-empty-dirs --relative \
 --include="*/" \
 --exclude="*/STAR/*.txt" \
 --include="*.job_md5sum_list.txt" \
 --include="*.sh*" --include="*.log" \
 --exclude="*" \
 ${WD}/./*_genomes/*/00_databases/ \
 . --dry-run
```

## 01_stats

Generate some basic stats about the genome and gene models using standard tools such as `stats.sh` (from `bbmap`) or `bedtools`.

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

```bash
while read PREFIX; do echo "$PREFIX"; ( cd $PREFIX/01_stats; ./run_01-stats.sh ); done < samples_new.txt
```

**Check for errors in log files**

```bash
while read PREFIX;
do 
  echo -e "\n## $PREFIX";
  find ${PREFIX}/ -name "*.sh.log.*" \
    | sort \
    | grep '01_stats' \
    | while read LOG;
      do 
        if ! $(grep -qi 'error\|warn\|fault' $LOG);
        then 
          echo "OK: $LOG";
        else 
          echo "FAILED: $LOG";
          grep -i 'error\|warn\|fault' $LOG
        fi; 
      done; 
done < samples_new.txt
```

The only genome that had errors in its log file was `Montipora_capitata_WTHIv1.1` which was because some of the gene models were from scaffolds that had been removed from the genome during filtering (by the authors of the genome paper). Not a big problem as this only affects 1 gene.

**Download results**

```bash
WD="timothy@coral.rutgers.edu:/scratch/timothy/projects/0031_Coral_genomes/03_Analysis/2022-02-03"
rsync -zarv --delete --prune-empty-dirs --relative \
 --include="*/" \
 --include="*.job_md5sum_list.txt" \
 --include="*.sh*" --include="*.log" \
 --include="*.GeneStats.tsv" \
 --exclude="*" \
 ${WD}/./*_genomes/*/01_stats/ \
 . --dry-run
```

## 02_busco

Run `bunco` on the genome and predicted proteins from each isolate. Run using the `metazoa_odb10` (2021-02-24) and `eukaryota_odb10` (2020-09-10) lineages. See https://busco-data.ezlab.org/v5/data/lineages/ for the list of available lineages (version 5).

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

**Check for errors in log files**

```bash
while read PREFIX;
do 
  echo -e "\n## $PREFIX";
  find ${PREFIX}/ -name "*.sh.log.*" \
    | sort \
    | grep '02_busco' \
    | while read LOG;
      do 
        if ! $(grep -qi 'error\|warn\|fault' $LOG);
        then 
          echo "OK: $LOG";
        else 
          echo "FAILED: $LOG";
          grep -i 'error\|warn\|fault' $LOG
        fi; 
      done; 
done < samples_new.txt
```

**Cleanup results files**

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
done < samples_new.txt
```

**Download results**

```bash
WD="timothy@coral.rutgers.edu:/scratch/timothy/projects/0031_Coral_genomes/03_Analysis/2022-02-03"
rsync -zarv --delete --prune-empty-dirs --relative \
 --include="*/" \
 --exclude="*.headersMap.tsv" --exclude="*_odb10/logs/*" \
 --include="*.job_md5sum_list.txt" \
 --include="*.sh*" --include="*.log" \
 --include="*.txt" --include="*.tsv" \
 --exclude="*" \
 ${WD}/./*_genomes/*/02_busco/ \
 . --dry-run
```



















## 07_short_read_trimming

Run `cutadapt` on the available Illumina genome sequencing reads. Execute the `run_Cutadapt.sh` script to runthe trimming process and to double check the adapter content using my custom `find_adapters_in_reads.py` script. 

```bash
ln -s ../../../../02_Scripts/adapters.fa
ln -s ../../../../02_Scripts/find_adapters_in_reads.py
```

**Download results**

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

**Download results**

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

















## 04_RepeatModeler

Run `RepeatModeler` to generate a set of *de novo* predicted repeats for each genome.

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

```bash
## Copy files to Amarel server and 'sbatch' each run_AmarelJob.sh script to submit to queue
rsync -avPL --relative timothy@coral.rutgers.edu:/scratch/timothy/projects/0031_Coral_genomes/03_Analysis/2022-02-03/./*_genomes/*/00_databases/*.assembly.fasta .

rsync -avzP --relative timothy@coral.rutgers.edu:/scratch/timothy/projects/0031_Coral_genomes/03_Analysis/2022-02-03/./*_genomes/*/01_stats/*.GeneStats.tsv .

rsync -avzP --relative timothy@coral.rutgers.edu:/scratch/timothy/projects/0031_Coral_genomes/03_Analysis/2022-02-03/./*_genomes/./*/04_RepeatModeler .

sbatch run_AmarelJob.sh
```

**Check for errors in log files**

```bash
while read PREFIX;
do 
  echo -e "\n## $PREFIX";
  find ${PREFIX}/ -name "*.sh.log.*" \
    | sort \
    | grep '04_RepeatModeler' \
    | while read LOG;
      do 
        if ! $(grep -qi 'error\|warn\|fault' $LOG);
        then 
          echo "OK: $LOG";
        else 
          echo "FAILED: $LOG";
          grep -i 'error\|warn\|fault' $LOG
        fi; 
      done; 
done < samples_new.txt

# Target the BuildDatabase.log and RepeatModeler.log log files
while read PREFIX;
do 
  echo -e "\n## $PREFIX";
  find ${PREFIX}/ -name "*log" \
    | sort \
    | grep '04_RepeatModeler' \
    | while read LOG;
      do 
        if ! $(grep -qi 'error\|warn\|fault' $LOG);
        then 
          echo "OK: $LOG";
        else 
          echo "FAILED: $LOG";
          grep -i 'error\|warn\|fault' $LOG
        fi; 
      done; 
done < samples_new.txt
```

**Cleanup results files**

```bash
while read PREFIX;
do
  echo ${PREFIX};
  rm -fr ${PREFIX}/04_RepeatModeler/busco_downloads/
  (cd ${PREFIX}/04_RepeatModeler/; 
   D=""; tar -zcf ${D}.tar.gz ${D} && rm -rf ${D};
  )
done < samples_new.txt
```

**Download results**

```bash
WD="timothy@coral.rutgers.edu:/scratch/timothy/projects/0031_Coral_genomes/03_Analysis/2022-02-03"
rsync -zarv --delete --prune-empty-dirs --relative \
 --include="*/" \
 --include="*.job_md5sum_list.txt" \
 --include="*.sh*" --include="*.log" \
 --include="*-families.fa" --include="*-families.stk" \
 --include="AmarelJob*" \
 --exclude="*" \
 ${WD}/./*/04_RepeatModeler/ \
 . --dry-run
```



## 05_RepeatMasker

Run `RepeatMasker` to identify and mask predicted repeats for each genome using the *de novo* predicted repeats and the repeats from the Dfam library. 

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

```bash
## Copy files to Amarel server and 'sbatch' each run_AmarelJob.sh script to submit to queue
rsync -avzP --relative timothy@coral.rutgers.edu:/scratch/timothy/projects/0031_Coral_genomes/03_Analysis/2022-02-03/./*_genomes/*/05_RepeatMasker .

sbatch run_AmarelJob.sh
```

**Check for errors in log files**

```bash
while read PREFIX;
do 
  echo -e "\n## $PREFIX";
  find ${PREFIX}/ -name "*.sh.log.*" \
    | sort \
    | grep '05_RepeatMasker' \
    | while read LOG;
      do 
        if ! $(grep -qi 'error\|warn\|fault' $LOG);
        then 
          echo "OK: $LOG";
        else 
          echo "FAILED: $LOG";
          grep -i 'error\|warn\|fault' $LOG
        fi; 
      done; 
done < samples_new.txt

# Target the RepeatMasker.log log file
while read PREFIX;
do 
  echo -e "\n## $PREFIX";
  find ${PREFIX}/ -name "*log" \
    | sort \
    | grep '05_RepeatMasker' \
    | while read LOG;
      do 
        if ! $(grep -qi 'error\|warn\|fault' $LOG);
        then 
          echo "OK: $LOG";
        else 
          echo "FAILED: $LOG";
          grep -i 'error\|warn\|fault' $LOG
        fi; 
      done; 
done < samples_new.txt
```

**Cleanup results files**

```bash
while read PREFIX;
do
  echo ${PREFIX};
  rm -fr ${PREFIX}/05_RepeatMasker/busco_downloads/
  (cd ${PREFIX}/05_RepeatMasker/; 
   D=""; tar -zcf ${D}.tar.gz ${D} && rm -rf ${D};
  )
done < samples_new.txt
```

**Download results**

```bash
WD="timothy@coral.rutgers.edu:/scratch/timothy/projects/0031_Coral_genomes/03_Analysis/2022-02-03"
rsync -zarv --delete --prune-empty-dirs --relative \
 --include="*/" \
 --include="*.job_md5sum_list.txt" \
 --include="*.sh*" --include="*.log" \
 --include="*.tbl" --include="*.html" \
 --include="*.divsum" \
 --include="*.divsum.Kimura_distance" \
 --include="AmarelJob*" \
 --exclude="*" \
 ${WD}/./*/05_RepeatMasker/ \
 . --dry-run
```

























## 10_functional_annotation

Run functional annotation (BLASTP RefSeq complete, BLASTP UniProt, eggnog-mapper, interproscan) on the predicted proteins.

Split protein file into 100 parts. This is for parallelization of the BLASTP against RefSeq complete analysis. 

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
  sed -e 's/PREFIX=""/PREFIX="'$PREFIX'"/' ../../run_10-blastp_UniProt.sh > run_10-blastp_UniProt.sh;
  sed -e 's/PREFIX=""/PREFIX="'$PREFIX'"/' ../../run_10-blastp_RefSeqComplete.sh > run_10-blastp_RefSeqComplete.sh;
  chmod +x *.sh
); done < samples_new.txt
```

```bash
while read PREFIX; do echo "$PREFIX"; ( cd $PREFIX/10_functional_annotation; ./run_10-eggnog-mapper.sh; ./run_10-InterProScan.sh ); done < samples_new.txt
```

Copy files to Amarel so we can run the BLASTP against RefSeq complete in a highly parallel fashion. 

```bash
rsync -avzPL --relative timothy@coral.rutgers.edu:/scratch/timothy/projects/0031_Coral_genomes/03_Analysis/2022-02-03/./*/*/00_databases/*.faa . --dry-run
rsync -avzPL --relative timothy@coral.rutgers.edu:/scratch/timothy/projects/0031_Coral_genomes/03_Analysis/2022-02-03/./*/*/10_functional_annotation . --dry-run

find . -name "md5sum_list.txt" | while read F; do echo $F; ( cd $(dirname $F); md5sum --quiet --strict --check $(basename $F) ); done
```

**Check for errors in log files**

```bash
while read PREFIX;
do 
  echo -e "\n## $PREFIX";
  find ${PREFIX}/ -name "*.sh.log.*" \
    | sort \
    | grep '10_functional_annotation' \
    | while read LOG;
      do 
        if ! $(grep -qi 'error\|warn\|fault' $LOG);
        then 
          echo "OK: $LOG";
        else 
          echo "FAILED: $LOG";
          grep -i 'error\|warn\|fault' $LOG
        fi; 
      done; 
done < samples_new.txt
```









## 03_sequence_contamination

Check to see how contaminated the genomes are wit bacteria (using RefSeq Bacteria genomic v204) and Symbiodiniaceae (available genomes; see `/scratch/timothy/databases/Symbiodiniaceae/genomes/`). Filter results, retaining hits with bitscore > 1000 and *e*-value <1e-10.

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
for F in */03_sequence_contamination/*/*.coverage;
do
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
```

**Download results**

```bash
WD="timothy@coral.rutgers.edu:/scratch/timothy/projects/0031_Coral_genomes/03_Analysis/2022-02-03"
rsync -zarv --delete --prune-empty-dirs --relative \
 --include="*/" \
 --include="*.job_md5sum_list.txt" \
 --include="*.sh*" --include="*.log" \
 --include="*.coverage.stats_filtered_scaffolds" \
 --exclude="*" \
 ${WD}/./*/03_sequence_contamination/ \
 . --dry-run
```

## 06_RepeatMasker_SCORs

Run `RepeatMasker` to identify and mask predicted repeats for each genome using just the Scleractinia COral-specific Repeat families (SCORs) that were published previously (https://doi.org/10.1016/j.ygeno.2017.06.003).  

```bash
./run-06_RepeatMasker_SCORs.sh samples.txt
```

**Download results**

```bash
WD="timothy@coral.rutgers.edu:/scratch/timothy/projects/0031_Coral_genomes/03_Analysis/2022-02-03"
rsync -zarv --delete --prune-empty-dirs --relative \
 --include="*/" \
 --include="*.job_md5sum_list.txt" \
 --include="*.sh*" --include="*.log" \
 --include="*.tbl" \
 --include="*.divsum" \
 --include="*.divsum.Kimura_distance" \
 --exclude="*" \
 ${WD}/./*/06_RepeatMasker_SCORs/ \
 . --dry-run
```

## 09_gene_prediction









