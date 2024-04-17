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
  ln -s ../../../../../01_Data/2021-12-03/*_transcriptomes/$PREFIX/transcriptome_assembly/$PREFIX.transcripts.pep.faa; 
  ln -s ../../../../../01_Data/2021-12-03/*_transcriptomes/$PREFIX/transcriptome_assembly/$PREFIX.transcripts.cds.fna;
  sed -e 's/PREFIX=""/PREFIX="'$PREFIX'"/' ../../run_00-setup_dataset.sh > run_00-setup_dataset.sh;
  chmod +x run_00-setup_dataset.sh
); done < samples_new.txt
```

```bash
while read PREFIX; do echo "$PREFIX"; ( cd $PREFIX/00_databases; ./run_00-setup_dataset.sh ); done < samples_new.txt
```

Run `bowtie2` on the transcript CDS.

```bash
while read PREFIX; 
do 
  mkdir -p $PREFIX/00_databases/bowtie2; 
( cd $PREFIX/00_databases/bowtie2; 
  ln -s ../$PREFIX.transcripts.cds.fna;
  sed -e 's/PREFIX=""/PREFIX="'$PREFIX'"/' ../../../run_00-bowtie2.sh > run_00-bowtie2.sh;
  chmod +x run_00-bowtie2.sh
); done < samples_new.txt
```

```bash
while read PREFIX; do echo "$PREFIX"; ( cd $PREFIX/00_databases/bowtie2; ./run_00-bowtie2.sh ); done < samples_new.txt
```

Run `blast` on the transcript CDS, and PEP files.

```bash
while read PREFIX; 
do 
  mkdir -p $PREFIX/00_databases/blast; 
( cd $PREFIX/00_databases/blast; 
  ln -s ../$PREFIX.transcripts.pep.faa;
  ln -s ../$PREFIX.transcripts.cds.fna;
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
    | grep '00_data' \
    | while read F;
      do
        echo "## $F"
        ( cd $(dirname $F); md5sum --quiet --strict --check $(basename $F) ); 
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
 "${WD}/./*_transcriptomes/*/00_databases/" \
 . --dry-run
```

## 01_stats

Generate some basic stats about the genome and gene models using standard tools such as `stats.sh` (from `bbmap`).

Run STATS on the genome and gene gff file.

```bash
while read PREFIX; 
do 
  mkdir -p $PREFIX/01_stats; 
( cd $PREFIX/01_stats; 
  ln -s ../00_databases/$PREFIX.transcripts.cds.fna;
  cp ../../../../../02_Scripts/CDS_stats .
  sed -e 's/PREFIX=""/PREFIX="'$PREFIX'"/' ../../run_01-stats.sh > run_01-stats.sh;
  chmod +x *.sh CDS_stats
); done < samples_new.txt
```

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
        ( cd $(dirname $F); md5sum --quiet --strict --check $(basename $F) ); 
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
 "${WD}/./*_transcriptomes/*/01_stats/" \
 . --dry-run
```

## 02_busco

Run `bunco` on the genome and predicted proteins from each isolate. Run using the `metazoa_odb10` (2021-02-24) and `eukaryota_odb10` (2020-09-10) lineages. See https://busco-data.ezlab.org/v5/data/lineages/ for the list of available lineages (version 5).

Run BUSCO on the predicted proteins.

```bash
while read PREFIX; 
do 
  mkdir -p $PREFIX/02_busco; 
( cd $PREFIX/02_busco; 
  ln -s ../00_databases/$PREFIX.transcripts.pep.faa;  
  sed -e 's/PREFIX=""/PREFIX="'$PREFIX'"/' ../../run_02-busco_protein.sh > run_02-busco_protein.sh;
  chmod +x *.sh
); done < samples_new.txt
```

```bash
while read PREFIX; do echo "$PREFIX"; ( cd $PREFIX/02_busco; ./run_02-busco_protein.sh ); done < samples_new.txt
```

Generate BUSCO summary results file for each transcriptome.

```bash
while read PREFIX; 
do 
  for F in $PREFIX/02_busco/pep.faa*_odb10;
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
   D="pep.faa.busco_eukaryota_odb10";   tar -zcf ${D}.tar.gz ${D} && rm -rf ${D};
   D="pep.faa.busco_metazoa_odb10";     tar -zcf ${D}.tar.gz ${D} && rm -rf ${D};
  )
done < samples_new.txt
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
        ( cd $(dirname $F); md5sum --quiet --strict --check $(basename $F) ); 
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
 "${WD}/./*_transcriptomes/*/02_busco/" \
 . --dry-run
```

## 03_sequence_contamination

Not applicable.

## 04_RepeatModeler

Not applicable.

## 05_RepeatMasker

Not applicable.

## 06_RepeatMasker_SCORs

Not applicable.

## 07_short_read_trimming

Not applicable.

## 08_genome_size_estimation

Not applicable.

## 09_gene_prediction

Not applicable.

## 10_functional_annotation

Run functional annotation (BLASTP RefSeq complete, BLASTP UniProt, eggnog-mapper, interproscan) on the predicted proteins.

Split protein file into 100 parts. This is for parallelization of the BLASTP against RefSeq complete analysis. 

```bash
while read PREFIX; 
do
  mkdir -p $PREFIX/10_functional_annotation; 
( cd $PREFIX/10_functional_annotation; 
  sed -e 's/*/X/g' ../00_databases/$PREFIX.transcripts.pep.faa > $PREFIX.transcripts.pep.faa;
  ./../../../../../02_Scripts/fasta-splitter.pl --line-length 0 --n-parts 100 --out-dir $PREFIX.transcripts.pep.faa.split $PREFIX.transcripts.pep.faa;
  md5sum ../00_databases/$PREFIX.transcripts.pep.faa $PREFIX.transcripts.pep.faa $PREFIX.transcripts.pep.faa.split/*.faa > md5sum_list.txt
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

---

## `eggnog-mapper` and `Interproscan`

### Run bash scripts

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
        if [ $(grep -i 'error\|warn\|fault' $LOG | grep -v 'this can take a few minutes and load up to' | wc -l) -eq 0 ];
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
rsync -avzPL --relative timothy@coral.rutgers.edu:/scratch/timothy/projects/0031_Coral_genomes/03_Analysis/2022-02-03/./*_transcriptomes/*/00_databases/*.faa . --dry-run
rsync -avzPL --relative timothy@coral.rutgers.edu:/scratch/timothy/projects/0031_Coral_genomes/03_Analysis/2022-02-03/./*_transcriptomes/*/10_functional_annotation/run_10-diamond_nr-Amarel.sh . --dry-run
rsync -avzPL --relative timothy@coral.rutgers.edu:/scratch/timothy/projects/0031_Coral_genomes/03_Analysis/2022-02-03/./*_transcriptomes/*/10_functional_annotation/*.faa . --dry-run
```

Run `diamond blastp` against NCBI's nr database on the Amarel server.

```bash
while read PREFIX; do echo "$PREFIX"; ( cd $PREFIX/10_functional_annotation; sbatch run_10-diamond_nr-Amarel.sh ); done < samples_new.txt
```

Copy files back to the Coral server from Amarel.

```bash
rsync -avzP --relative \
  $PWD/./*_transcriptomes/*/10_functional_annotation \
  timothy@coral.rutgers.edu:/scratch/timothy/projects/0031_Coral_genomes/03_Analysis/2022-02-03/ \
  --dry-run
```



### Check for errors in log files

```bash
while read PREFIX;
do 
  echo -e "\n## $PREFIX";
  cat $PREFIX/10_functional_annotation/run_10-diamond_nr-Amarel.sh.* | grep 'queries aligned'
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



### Download results

```bash
WD="timothy@coral.rutgers.edu:/scratch/timothy/projects/0031_Coral_genomes/03_Analysis/2022-02-03"
rsync -zarv --delete --prune-empty-dirs --relative \
 --include="*/" \
 --include="*.job_md5sum_list.txt" \
 --include="*.sh*" --include="*.log" \
 --exclude="*" \
 "${WD}/./*_transcriptomes/*/10_functional_annotation/" \
 . --dry-run
```





