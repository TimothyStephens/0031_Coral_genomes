# Common analysis for all transcriptomes being analyzed across different projects

This project is a repository of analysis run on each coral or our group transcriptome that will/can be used across multiple projects (e.g. `busco` or `eggnog` results). For each genome the analysis was run as described in this document unless otherwise stated.

## Setup analysis directory

Setup bash environment.

```bash
conda activate py27
```

## Download and cleanup transcriptome data

The coral transcriptome assemblies and predicted proteins are in `../../01_Data/2021-12-03/coral_transcriptomes/`

The coral out group transcriptome assemblies and predicted proteins are in `../../01_Data/2021-12-03/outgroup_transcriptomes/`

The `../../01_Data/2021-12-03/README.md ` file explains where each genome was downloaded from and the steps I took to clean and reformat the data.

## 00_databases

Make sequence databases to use for subsequent analysis. 

- `samtools` `*.fai` (CDS+PEP)
- `gatk` `*.dict` (CDS)
- `blast` (CDS+PEP)
- `bowtie2` (CDS)

Set up "00_databases" directory

```bash
while read PREFIX; 
do 
  mkdir -p $PREFIX/00_databases; 
( cd $PREFIX/00_databases; 
  ln -s ../../../../../01_Data/2021-12-03/*_genomes/$PREFIX/transcriptome_assembly/$PREFIX.transcripts.pep.faa; 
  ln -s ../../../../../01_Data/2021-12-03/*_genomes/$PREFIX/transcriptome_assembly/$PREFIX.transcripts.cds.fna;
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
rsync -avzPL --relative timothy@coral.rutgers.edu:/scratch/timothy/projects/0031_Coral_genomes/03_Analysis/2022-02-03/coral_genomes/./*/00_databases/*.faa . --dry-run
rsync -avzPL --relative timothy@coral.rutgers.edu:/scratch/timothy/projects/0031_Coral_genomes/03_Analysis/2022-02-03/coral_genomes/./*/10_functional_annotation . --dry-run
rsync -avzPL --relative timothy@coral.rutgers.edu:/scratch/timothy/projects/0031_Coral_genomes/03_Analysis/2022-02-03/coral_genomes/./*/10_functional_annotation/run_10-blastp_RefSeqComplete.sh . --dry-run

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






