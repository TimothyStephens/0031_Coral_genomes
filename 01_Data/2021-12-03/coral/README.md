# Download and reformat coral genome data

Download coral genome assembly + predicted genes (`gff` file + `CDS` and/or `PEP` sequences if available).

Rename the files, and sequences in the files, so they both use the same informative prefix. Reformat genome+CDS+PEP `fasta` files using unique prefix; also reformat sequence by removing description, non-ATGCN characters (genome+CDS only), making isingle line `fasta`, and all uppercase characters.

```bash
conda activate py27
export PATH="$PATH:/home/timothy/programs/samtools-1.11/bin"
export PATH="$PATH:/home/timothy/programs/seqkit_v0.15.0"
export PATH="$PATH:/home/timothy/programs/bbmap"
export PATH="$PATH:/home/timothy/programs/sratoolkit.2.9.6-1-centos_linux64/bin"
export PATH="$PATH:/home/timothy/programs/gffread-0.11.6.Linux_x86_64"
export PATH="$PATH:/home/timothy/programs/TransDecoder-v5.5.0"
```

## **NOTES**

- *Montipora capitata* WTHIv1.1:
    - xfSc0014868 appears to be in gff file but is missing from the genome assembly and when you use the gff to translate some of the genes into proteins ~100 of them have in-frame stop codons. This suggests that the version of the genome available online is different (only slightly) to the one used for the gene prediction.
- *Stylophora pistillata* GAJOv1:
    - There are 39453 genes reported in the paper and in the `gff` file but there are 39511 in the supplied protein and CDS `fasta` files. There are 58 proteins (named "Prus_finalg...") which are in the supplied `fasta` files but not the `gff` file.
- *Exaiptasia diaphana* (*Aiptasia* strain CC7)
    - There are a number of gene models with duplicate names in the input `gff` file. Remove the duplicates which have only 'exon' features.
    - There looks to be something weird with the phase of some of the genes in the `gff` that causes a large number of inflame stop codons. Use the original protein file available through the project website - it should mirror up with the CDS file we just generated well enough for any downstream analysis.



**Internal stop codons**

| Isolate                         | No. internal stop codons | Notes        |
| ------------------------------- | ------------------------ | ------------ |
| Acropora_acuminata_SLJPv1       | 136                      |              |
| Acropora_awi_SLJPv1             | 197                      |              |
| Acropora_cytherea_SLJPv1        | 156                      |              |
| Acropora_digitifera_SLJPv2      | 167                      |              |
| Acropora_echinata_SLJPv1        | 203                      |              |
| Acropora_florida_SLJPv1         | 164                      |              |
| Acropora_gemmifera_SLJPv1       | 161                      |              |
| Acropora_hyacinthus_SLJPv1      | 191                      |              |
| Acropora_intermedia_SLJPv1      | 167                      |              |
| Acropora_microphthalma_SLJPv1   | 171                      |              |
| Acropora_muricata_SLJPv1        | 168                      |              |
| Acropora_nasuta_SLJPv1          | 153                      |              |
| Acropora_selago_SLJPv1          | 169                      |              |
| Acropora_yongei_SLJPv1          | 199                      |              |
| Astreopora_myriophthalma_SLJPv1 | 243                      |              |
| Montipora_cactus_SLJPv2         | 1504                     |              |
| Montipora_efflorescens_SLJPv2   | 1443                     |              |
| Pocillopora_verrucosa_RSSAv1    | 248                      | PEP provided |
| Porites_australiensis_SLJPv1    | 276                      |              |
| Stylophora_pistillata_GAJOv1    | 768                      | PEP provided |

# Coral Genomes

## *Acropora acuminata* SLJPv1

From Sekisei Lagoon, Okinawa, Japan

Version: 1.0

Paper: https://doi.org/10.1093/molbev/msaa216

Website: https://marinegenomics.oist.jp/aacu/viewer/download?project_id=84

```bash
wget https://marinegenomics.oist.jp/aacu/download/aacu.gff.gz
wget https://marinegenomics.oist.jp/aacu/download/aacu.fasta.gz
```

```bash
PREFIX="Acropora_acuminata_SLJPv1"
```

**Genome**

```bash
zcat aacu.fasta.gz \
  | sed -e 's/_pilon$//' \
  | sed -e "s/>/>${PREFIX}___/" \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.assembly.fasta
```

**GFF**

```bash
zcat aacu.gff.gz \
  | sed -e "s/^/${PREFIX}___/" -e 's/_pilon//' \
  | sed -e "s/ID=aacu_/ID=${PREFIX}___/" \
  | sed -e "s/Parent=aacu_/Parent=${PREFIX}___/" \
  | sed -e "s/Name=aacu_/Name=${PREFIX}___/" \
  | awk -F'\t' '$3=="mRNA" || $3=="CDS" || $3=="exon"' \
  | gffread \
  | awk '$1!~"^#"' \
  | sed -e 's/;geneID=.*//' \
  > ${PREFIX}.genes.all.gff3

cat ${PREFIX}.genes.all.gff3 \
  | awk -F'\t' '$3=="CDS" {
      split($9,a,"="); 
      CDS[a[2]]=CDS[a[2]]+(($5-$4)+1)
    } END {
      for (i in CDS){
        split(i,ii,"\\.t")
        print ii[1]"\t"i"\t"CDS[i]
    } }' \
  > ${PREFIX}.genes.all.gff3.CDS_lengths.tsv

cat ${PREFIX}.genes.all.gff3.CDS_lengths.tsv \
  | sort -k1,1 -k3,3nr \
  | awk '!seen[$1]++{print $2}' \
  > ${PREFIX}.genes.all.gff3.longest.txt

cat ${PREFIX}.genes.all.gff3 \
  | awk -F'\t' '{split($9,a,"="); print $0"\t"a[2]}' \
  | ~/scripts/grepf_column.py -c 10 -f ${PREFIX}.genes.all.gff3.longest.txt \
  | cut -f 1-9 \
  > ${PREFIX}.genes.gff3

gffread -S \
  -g ${PREFIX}.assembly.fasta \
  -x ${PREFIX}.genes.cds.fna.tmp \
  -y ${PREFIX}.genes.pep.faa.tmp \
  ${PREFIX}.genes.gff3
```

**CDS**

```bash
cat ${PREFIX}.genes.cds.fna.tmp \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.genes.cds.fna
rm ${PREFIX}.genes.cds.fna.tmp
```

**PEP**

```bash
cat ${PREFIX}.genes.pep.faa.tmp \
  | reformat.sh fastawrap=0 iupacToN=f ignorejunk=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.genes.pep.faa
rm ${PREFIX}.genes.pep.faa.tmp
```

## *Acropora awi* SLJPv1

From Sekisei Lagoon, Okinawa, Japan

Version: 1.0

Paper: https://doi.org/10.1093/molbev/msaa216

Website: https://marinegenomics.oist.jp/aawi/viewer/download?project_id=85

```bash
wget https://marinegenomics.oist.jp/aawi/download/aawi.gff.gz
wget https://marinegenomics.oist.jp/aawi/download/aawi.fasta.gz
```

```bash
PREFIX="Acropora_awi_SLJPv1"
```

**Genome**

```bash
zcat aawi.fasta.gz \
  | sed -e 's/_pilon$//' \
  | sed -e "s/>/>${PREFIX}___/" \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.assembly.fasta
```

**GFF**

```bash
zcat aawi.gff.gz \
  | sed -e "s/^/${PREFIX}___/" -e 's/_pilon//' \
  | sed -e "s/ID=aawi_/ID=${PREFIX}___/" \
  | sed -e "s/Parent=aawi_/Parent=${PREFIX}___/" \
  | sed -e "s/Name=aawi_/Name=${PREFIX}___/" \
  | awk -F'\t' '$3=="mRNA" || $3=="CDS" || $3=="exon"' \
  | gffread \
  | awk '$1!~"^#"' \
  | sed -e 's/;geneID=.*//' \
  > ${PREFIX}.genes.all.gff3

cat ${PREFIX}.genes.all.gff3 \
  | awk -F'\t' '$3=="CDS" {
      split($9,a,"="); 
      CDS[a[2]]=CDS[a[2]]+(($5-$4)+1)
    } END {
      for (i in CDS){
        split(i,ii,"\\.t")
        print ii[1]"\t"i"\t"CDS[i]
    } }' \
  > ${PREFIX}.genes.all.gff3.CDS_lengths.tsv

cat ${PREFIX}.genes.all.gff3.CDS_lengths.tsv \
  | sort -k1,1 -k3,3nr \
  | awk '!seen[$1]++{print $2}' \
  > ${PREFIX}.genes.all.gff3.longest.txt

cat ${PREFIX}.genes.all.gff3 \
  | awk -F'\t' '{split($9,a,"="); print $0"\t"a[2]}' \
  | ~/scripts/grepf_column.py -c 10 -f ${PREFIX}.genes.all.gff3.longest.txt \
  | cut -f 1-9 \
  > ${PREFIX}.genes.gff3

gffread -S \
  -g ${PREFIX}.assembly.fasta \
  -x ${PREFIX}.genes.cds.fna.tmp \
  -y ${PREFIX}.genes.pep.faa.tmp \
  ${PREFIX}.genes.gff3
```

**CDS**

```bash
cat ${PREFIX}.genes.cds.fna.tmp \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.genes.cds.fna
rm ${PREFIX}.genes.cds.fna.tmp
```

**PEP**

```bash
cat ${PREFIX}.genes.pep.faa.tmp \
  | reformat.sh fastawrap=0 iupacToN=f ignorejunk=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.genes.pep.faa
rm ${PREFIX}.genes.pep.faa.tmp
```

## *Acropora cytherea* SLJPv1

From Sekisei Lagoon, Okinawa, Japan

Version: 1.0

Paper: https://doi.org/10.1093/molbev/msaa216

Website: https://marinegenomics.oist.jp/acyt/viewer/download?project_id=86

```bash
wget https://marinegenomics.oist.jp/acyt/download/acyt.fasta.gz
wget https://marinegenomics.oist.jp/acyt/download/acyt.gff.gz
```

```bash
PREFIX="Acropora_cytherea_SLJPv1"
```

**Genome**

```bash
zcat acyt.fasta.gz \
  | sed -e 's/_pilon$//' \
  | sed -e "s/>/>${PREFIX}___/" \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.assembly.fasta
```

**GFF**

```bash
zcat acyt.gff.gz \
  | sed -e "s/^/${PREFIX}___/" -e 's/_pilon//' \
  | sed -e "s/ID=acyt_/ID=${PREFIX}___/" \
  | sed -e "s/Parent=acyt_/Parent=${PREFIX}___/" \
  | sed -e "s/Name=acyt_/Name=${PREFIX}___/" \
  | awk -F'\t' '$3=="mRNA" || $3=="CDS" || $3=="exon"' \
  | gffread \
  | awk '$1!~"^#"' \
  | sed -e 's/;geneID=.*//' \
  > ${PREFIX}.genes.all.gff3

cat ${PREFIX}.genes.all.gff3 \
  | awk -F'\t' '$3=="CDS" {
      split($9,a,"="); 
      CDS[a[2]]=CDS[a[2]]+(($5-$4)+1)
    } END {
      for (i in CDS){
        split(i,ii,"\\.t")
        print ii[1]"\t"i"\t"CDS[i]
    } }' \
  > ${PREFIX}.genes.all.gff3.CDS_lengths.tsv

cat ${PREFIX}.genes.all.gff3.CDS_lengths.tsv \
  | sort -k1,1 -k3,3nr \
  | awk '!seen[$1]++{print $2}' \
  > ${PREFIX}.genes.all.gff3.longest.txt

cat ${PREFIX}.genes.all.gff3 \
  | awk -F'\t' '{split($9,a,"="); print $0"\t"a[2]}' \
  | ~/scripts/grepf_column.py -c 10 -f ${PREFIX}.genes.all.gff3.longest.txt \
  | cut -f 1-9 \
  > ${PREFIX}.genes.gff3

gffread -S \
  -g ${PREFIX}.assembly.fasta \
  -x ${PREFIX}.genes.cds.fna.tmp \
  -y ${PREFIX}.genes.pep.faa.tmp \
  ${PREFIX}.genes.gff3
```

**CDS**

```bash
cat ${PREFIX}.genes.cds.fna.tmp \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.genes.cds.fna
rm ${PREFIX}.genes.cds.fna.tmp
```

**PEP**

```bash
cat ${PREFIX}.genes.pep.faa.tmp \
  | reformat.sh fastawrap=0 iupacToN=f ignorejunk=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.genes.pep.faa
rm ${PREFIX}.genes.pep.faa.tmp
```

## *Acropora digitifera* SLJPv2

From Sekisei Lagoon, Okinawa, Japan

Version: 2.0

Paper: https://doi.org/10.1093/molbev/msaa216

Website: https://marinegenomics.oist.jp/adig/viewer/download?project_id=87

```bash
wget https://marinegenomics.oist.jp/adig/download/adig.fasta.gz
wget https://marinegenomics.oist.jp/adig/download/adig.gff.gz
```

```bash
PREFIX="Acropora_digitifera_SLJPv2"
```

**Genome**

```bash
zcat adig.fasta.gz \
  | sed -e 's/_arrow_pilon$//' \
  | sed -e "s/>/>${PREFIX}___/" \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.assembly.fasta
```

**GFF**

```bash
zcat adig.gff.gz \
  | sed -e "s/^/${PREFIX}___/" -e 's/_arrow_pilon//' \
  | sed -e "s/ID=adig_/ID=${PREFIX}___/" \
  | sed -e "s/Parent=adig_/Parent=${PREFIX}___/" \
  | sed -e "s/Name=adig_/Name=${PREFIX}___/" \
  | awk -F'\t' '$3=="mRNA" || $3=="CDS" || $3=="exon"' \
  | gffread \
  | awk '$1!~"^#"' \
  | sed -e 's/;geneID=.*//' \
  > ${PREFIX}.genes.all.gff3

cat ${PREFIX}.genes.all.gff3 \
  | awk -F'\t' '$3=="CDS" {
      split($9,a,"="); 
      CDS[a[2]]=CDS[a[2]]+(($5-$4)+1)
    } END {
      for (i in CDS){
        split(i,ii,"\\.t")
        print ii[1]"\t"i"\t"CDS[i]
    } }' \
  > ${PREFIX}.genes.all.gff3.CDS_lengths.tsv

cat ${PREFIX}.genes.all.gff3.CDS_lengths.tsv \
  | sort -k1,1 -k3,3nr \
  | awk '!seen[$1]++{print $2}' \
  > ${PREFIX}.genes.all.gff3.longest.txt

cat ${PREFIX}.genes.all.gff3 \
  | awk -F'\t' '{split($9,a,"="); print $0"\t"a[2]}' \
  | ~/scripts/grepf_column.py -c 10 -f ${PREFIX}.genes.all.gff3.longest.txt \
  | cut -f 1-9 \
  > ${PREFIX}.genes.gff3

gffread -S \
  -g ${PREFIX}.assembly.fasta \
  -x ${PREFIX}.genes.cds.fna.tmp \
  -y ${PREFIX}.genes.pep.faa.tmp \
  ${PREFIX}.genes.gff3
```

**CDS**

```bash
cat ${PREFIX}.genes.cds.fna.tmp \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.genes.cds.fna
rm ${PREFIX}.genes.cds.fna.tmp
```

**PEP**

```bash
cat ${PREFIX}.genes.pep.faa.tmp \
  | reformat.sh fastawrap=0 iupacToN=f ignorejunk=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.genes.pep.faa
rm ${PREFIX}.genes.pep.faa.tmp
```

## *Acropora echinata* SLJPv1

From Sekisei Lagoon, Okinawa, Japan

Version: 1.0

Paper: https://doi.org/10.1093/molbev/msaa216

Website: https://marinegenomics.oist.jp/aech/viewer/download?project_id=88

```bash
wget https://marinegenomics.oist.jp/aech/download/aech.fasta.gz
wget https://marinegenomics.oist.jp/aech/download/aech.gff.gz
```

```bash
PREFIX="Acropora_echinata_SLJPv1"
```

**Genome**

```bash
zcat aech.fasta.gz \
  | sed -e 's/_pilon$//' \
  | sed -e "s/>/>${PREFIX}___/" \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.assembly.fasta
```

**GFF**

```bash
zcat aech.gff.gz \
  | sed -e "s/^/${PREFIX}___/" -e 's/_pilon//' \
  | sed -e "s/ID=aech_/ID=${PREFIX}___/" \
  | sed -e "s/Parent=aech_/Parent=${PREFIX}___/" \
  | sed -e "s/Name=aech_/Name=${PREFIX}___/" \
  | awk -F'\t' '$3=="mRNA" || $3=="CDS" || $3=="exon"' \
  | gffread \
  | awk '$1!~"^#"' \
  | sed -e 's/;geneID=.*//' \
  > ${PREFIX}.genes.all.gff3

cat ${PREFIX}.genes.all.gff3 \
  | awk -F'\t' '$3=="CDS" {
      split($9,a,"="); 
      CDS[a[2]]=CDS[a[2]]+(($5-$4)+1)
    } END {
      for (i in CDS){
        split(i,ii,"\\.t")
        print ii[1]"\t"i"\t"CDS[i]
    } }' \
  > ${PREFIX}.genes.all.gff3.CDS_lengths.tsv

cat ${PREFIX}.genes.all.gff3.CDS_lengths.tsv \
  | sort -k1,1 -k3,3nr \
  | awk '!seen[$1]++{print $2}' \
  > ${PREFIX}.genes.all.gff3.longest.txt

cat ${PREFIX}.genes.all.gff3 \
  | awk -F'\t' '{split($9,a,"="); print $0"\t"a[2]}' \
  | ~/scripts/grepf_column.py -c 10 -f ${PREFIX}.genes.all.gff3.longest.txt \
  | cut -f 1-9 \
  > ${PREFIX}.genes.gff3

gffread -S \
  -g ${PREFIX}.assembly.fasta \
  -x ${PREFIX}.genes.cds.fna.tmp \
  -y ${PREFIX}.genes.pep.faa.tmp \
  ${PREFIX}.genes.gff3
```

**CDS**

```bash
cat ${PREFIX}.genes.cds.fna.tmp \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.genes.cds.fna
rm ${PREFIX}.genes.cds.fna.tmp
```

**PEP**

```bash
cat ${PREFIX}.genes.pep.faa.tmp \
  | reformat.sh fastawrap=0 iupacToN=f ignorejunk=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.genes.pep.faa
rm ${PREFIX}.genes.pep.faa.tmp
```

## *Acropora florida* SLJPv1

From Sekisei Lagoon, Okinawa, Japan

Version: 1.0

Paper: https://doi.org/10.1093/molbev/msaa216

Website: https://marinegenomics.oist.jp/aflo/viewer/download?project_id=89

```bash
wget https://marinegenomics.oist.jp/aflo/download/aflo.fasta.gz
wget https://marinegenomics.oist.jp/aflo/download/aflo.gff.gz
```

```bash
PREFIX="Acropora_florida_SLJPv1"
```

**Genome**

```bash
zcat aflo.fasta.gz \
  | sed -e 's/_pilon$//' \
  | sed -e "s/>/>${PREFIX}___/" \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.assembly.fasta
```

**GFF**

```bash
zcat aflo.gff.gz \
  | sed -e "s/^/${PREFIX}___/" -e 's/_pilon//' \
  | sed -e "s/ID=aflo_/ID=${PREFIX}___/" \
  | sed -e "s/Parent=aflo_/Parent=${PREFIX}___/" \
  | sed -e "s/Name=aflo_/Name=${PREFIX}___/" \
  | awk -F'\t' '$3=="mRNA" || $3=="CDS" || $3=="exon"' \
  | gffread \
  | awk '$1!~"^#"' \
  | sed -e 's/;geneID=.*//' \
  > ${PREFIX}.genes.all.gff3

cat ${PREFIX}.genes.all.gff3 \
  | awk -F'\t' '$3=="CDS" {
      split($9,a,"="); 
      CDS[a[2]]=CDS[a[2]]+(($5-$4)+1)
    } END {
      for (i in CDS){
        split(i,ii,"\\.t")
        print ii[1]"\t"i"\t"CDS[i]
    } }' \
  > ${PREFIX}.genes.all.gff3.CDS_lengths.tsv

cat ${PREFIX}.genes.all.gff3.CDS_lengths.tsv \
  | sort -k1,1 -k3,3nr \
  | awk '!seen[$1]++{print $2}' \
  > ${PREFIX}.genes.all.gff3.longest.txt

cat ${PREFIX}.genes.all.gff3 \
  | awk -F'\t' '{split($9,a,"="); print $0"\t"a[2]}' \
  | ~/scripts/grepf_column.py -c 10 -f ${PREFIX}.genes.all.gff3.longest.txt \
  | cut -f 1-9 \
  > ${PREFIX}.genes.gff3

gffread -S \
  -g ${PREFIX}.assembly.fasta \
  -x ${PREFIX}.genes.cds.fna.tmp \
  -y ${PREFIX}.genes.pep.faa.tmp \
  ${PREFIX}.genes.gff3
```

**CDS**

```bash
cat ${PREFIX}.genes.cds.fna.tmp \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.genes.cds.fna
rm ${PREFIX}.genes.cds.fna.tmp
```

**PEP**

```bash
cat ${PREFIX}.genes.pep.faa.tmp \
  | reformat.sh fastawrap=0 iupacToN=f ignorejunk=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.genes.pep.faa
rm ${PREFIX}.genes.pep.faa.tmp
```

## *Acropora gemmifera* SLJPv1

From Sekisei Lagoon, Okinawa, Japan

Version: 1.0

Paper: https://doi.org/10.1093/molbev/msaa216

Website: https://marinegenomics.oist.jp/agem/viewer/download?project_id=90

```bash
wget https://marinegenomics.oist.jp/agem/download/agem.fasta.gz
wget https://marinegenomics.oist.jp/agem/download/agem.gff.gz
```

```bash
PREFIX="Acropora_gemmifera_SLJPv1"
```

**Genome**

```bash
zcat agem.fasta.gz \
  | sed -e 's/_pilon$//' \
  | sed -e "s/>/>${PREFIX}___/" \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.assembly.fasta
```

**GFF**

```bash
zcat agem.gff.gz \
  | sed -e "s/^/${PREFIX}___/" -e 's/_pilon//' \
  | sed -e "s/ID=agem_/ID=${PREFIX}___/" \
  | sed -e "s/Parent=agem_/Parent=${PREFIX}___/" \
  | sed -e "s/Name=agem_/Name=${PREFIX}___/" \
  | awk -F'\t' '$3=="mRNA" || $3=="CDS" || $3=="exon"' \
  | gffread \
  | awk '$1!~"^#"' \
  | sed -e 's/;geneID=.*//' \
  > ${PREFIX}.genes.all.gff3

cat ${PREFIX}.genes.all.gff3 \
  | awk -F'\t' '$3=="CDS" {
      split($9,a,"="); 
      CDS[a[2]]=CDS[a[2]]+(($5-$4)+1)
    } END {
      for (i in CDS){
        split(i,ii,"\\.t")
        print ii[1]"\t"i"\t"CDS[i]
    } }' \
  > ${PREFIX}.genes.all.gff3.CDS_lengths.tsv

cat ${PREFIX}.genes.all.gff3.CDS_lengths.tsv \
  | sort -k1,1 -k3,3nr \
  | awk '!seen[$1]++{print $2}' \
  > ${PREFIX}.genes.all.gff3.longest.txt

cat ${PREFIX}.genes.all.gff3 \
  | awk -F'\t' '{split($9,a,"="); print $0"\t"a[2]}' \
  | ~/scripts/grepf_column.py -c 10 -f ${PREFIX}.genes.all.gff3.longest.txt \
  | cut -f 1-9 \
  > ${PREFIX}.genes.gff3

gffread -S \
  -g ${PREFIX}.assembly.fasta \
  -x ${PREFIX}.genes.cds.fna.tmp \
  -y ${PREFIX}.genes.pep.faa.tmp \
  ${PREFIX}.genes.gff3
```

**CDS**

```bash
cat ${PREFIX}.genes.cds.fna.tmp \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.genes.cds.fna
rm ${PREFIX}.genes.cds.fna.tmp
```

**PEP**

```bash
cat ${PREFIX}.genes.pep.faa.tmp \
  | reformat.sh fastawrap=0 iupacToN=f ignorejunk=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.genes.pep.faa
rm ${PREFIX}.genes.pep.faa.tmp
```

## *Acropora hyacinthus* SLJPv1

From Sekisei Lagoon, Okinawa, Japan

Version: 1.0

Paper: https://doi.org/10.1093/molbev/msaa216

Website: https://marinegenomics.oist.jp/ahya/viewer/download?project_id=91

```bash
wget https://marinegenomics.oist.jp/ahya/download/ahya.fasta.gz
wget https://marinegenomics.oist.jp/ahya/download/ahya.gff.gz
```

```bash
PREFIX="Acropora_hyacinthus_SLJPv1"
```

**Genome**

```bash
zcat ahya.fasta.gz \
  | sed -e 's/_pilon$//' \
  | sed -e "s/>/>${PREFIX}___/" \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.assembly.fasta
```

**GFF**

```bash
zcat ahya.gff.gz \
  | sed -e "s/^/${PREFIX}___/" -e 's/_pilon//' \
  | sed -e "s/ID=ahya_/ID=${PREFIX}___/" \
  | sed -e "s/Parent=ahya_/Parent=${PREFIX}___/" \
  | sed -e "s/Name=ahya_/Name=${PREFIX}___/" \
  | awk -F'\t' '$3=="mRNA" || $3=="CDS" || $3=="exon"' \
  | gffread \
  | awk '$1!~"^#"' \
  | sed -e 's/;geneID=.*//' \
  > ${PREFIX}.genes.all.gff3

cat ${PREFIX}.genes.all.gff3 \
  | awk -F'\t' '$3=="CDS" {
      split($9,a,"="); 
      CDS[a[2]]=CDS[a[2]]+(($5-$4)+1)
    } END {
      for (i in CDS){
        split(i,ii,"\\.t")
        print ii[1]"\t"i"\t"CDS[i]
    } }' \
  > ${PREFIX}.genes.all.gff3.CDS_lengths.tsv

cat ${PREFIX}.genes.all.gff3.CDS_lengths.tsv \
  | sort -k1,1 -k3,3nr \
  | awk '!seen[$1]++{print $2}' \
  > ${PREFIX}.genes.all.gff3.longest.txt

cat ${PREFIX}.genes.all.gff3 \
  | awk -F'\t' '{split($9,a,"="); print $0"\t"a[2]}' \
  | ~/scripts/grepf_column.py -c 10 -f ${PREFIX}.genes.all.gff3.longest.txt \
  | cut -f 1-9 \
  > ${PREFIX}.genes.gff3

gffread -S \
  -g ${PREFIX}.assembly.fasta \
  -x ${PREFIX}.genes.cds.fna.tmp \
  -y ${PREFIX}.genes.pep.faa.tmp \
  ${PREFIX}.genes.gff3
```

**CDS**

```bash
cat ${PREFIX}.genes.cds.fna.tmp \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.genes.cds.fna
rm ${PREFIX}.genes.cds.fna.tmp
```

**PEP**

```bash
cat ${PREFIX}.genes.pep.faa.tmp \
  | reformat.sh fastawrap=0 iupacToN=f ignorejunk=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.genes.pep.faa
rm ${PREFIX}.genes.pep.faa.tmp
```

## *Acropora intermedia* SLJPv1

From Sekisei Lagoon, Okinawa, Japan

Version: 1.0

Paper: https://doi.org/10.1093/molbev/msaa216

Website: https://marinegenomics.oist.jp/aint/viewer/download?project_id=92

```bash
wget https://marinegenomics.oist.jp/aint/download/aint.fasta.gz
wget https://marinegenomics.oist.jp/aint/download/aint.gff.gz
```

```bash
PREFIX="Acropora_intermedia_SLJPv1"
```

**Genome**

```bash
zcat aint.fasta.gz \
  | sed -e 's/_pilon$//' \
  | sed -e "s/>/>${PREFIX}___/" \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.assembly.fasta
```

**GFF**

```bash
zcat aint.gff.gz \
  | sed -e "s/^/${PREFIX}___/" -e 's/_pilon//' \
  | sed -e "s/ID=aint_/ID=${PREFIX}___/" \
  | sed -e "s/Parent=aint_/Parent=${PREFIX}___/" \
  | sed -e "s/Name=aint_/Name=${PREFIX}___/" \
  | awk -F'\t' '$3=="mRNA" || $3=="CDS" || $3=="exon"' \
  | gffread \
  | awk '$1!~"^#"' \
  | sed -e 's/;geneID=.*//' \
  > ${PREFIX}.genes.all.gff3

cat ${PREFIX}.genes.all.gff3 \
  | awk -F'\t' '$3=="CDS" {
      split($9,a,"="); 
      CDS[a[2]]=CDS[a[2]]+(($5-$4)+1)
    } END {
      for (i in CDS){
        split(i,ii,"\\.t")
        print ii[1]"\t"i"\t"CDS[i]
    } }' \
  > ${PREFIX}.genes.all.gff3.CDS_lengths.tsv

cat ${PREFIX}.genes.all.gff3.CDS_lengths.tsv \
  | sort -k1,1 -k3,3nr \
  | awk '!seen[$1]++{print $2}' \
  > ${PREFIX}.genes.all.gff3.longest.txt

cat ${PREFIX}.genes.all.gff3 \
  | awk -F'\t' '{split($9,a,"="); print $0"\t"a[2]}' \
  | ~/scripts/grepf_column.py -c 10 -f ${PREFIX}.genes.all.gff3.longest.txt \
  | cut -f 1-9 \
  > ${PREFIX}.genes.gff3

gffread -S \
  -g ${PREFIX}.assembly.fasta \
  -x ${PREFIX}.genes.cds.fna.tmp \
  -y ${PREFIX}.genes.pep.faa.tmp \
  ${PREFIX}.genes.gff3
```

**CDS**

```bash
cat ${PREFIX}.genes.cds.fna.tmp \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.genes.cds.fna
rm ${PREFIX}.genes.cds.fna.tmp
```

**PEP**

```bash
cat ${PREFIX}.genes.pep.faa.tmp \
  | reformat.sh fastawrap=0 iupacToN=f ignorejunk=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.genes.pep.faa
rm ${PREFIX}.genes.pep.faa.tmp
```

## *Acropora microphthalma* SLJPv1

From Sekisei Lagoon, Okinawa, Japan

Version: 1.0

Paper: https://doi.org/10.1093/molbev/msaa216

Website: https://marinegenomics.oist.jp/amic/viewer/download?project_id=93

```bash
wget https://marinegenomics.oist.jp/amic/download/amic.fasta.gz
wget https://marinegenomics.oist.jp/amic/download/amic.gff.gz
```

```bash
PREFIX="Acropora_microphthalma_SLJPv1"
```

**Genome**

```bash
zcat amic.fasta.gz \
  | sed -e 's/_pilon$//' \
  | sed -e "s/>/>${PREFIX}___/" \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.assembly.fasta
```

**GFF**

```bash
zcat amic.gff.gz \
  | sed -e "s/^/${PREFIX}___/" -e 's/_pilon//' \
  | sed -e "s/ID=amic_/ID=${PREFIX}___/" \
  | sed -e "s/Parent=amic_/Parent=${PREFIX}___/" \
  | sed -e "s/Name=amic_/Name=${PREFIX}___/" \
  | awk -F'\t' '$3=="mRNA" || $3=="CDS" || $3=="exon"' \
  | gffread \
  | awk '$1!~"^#"' \
  | sed -e 's/;geneID=.*//' \
  > ${PREFIX}.genes.all.gff3

cat ${PREFIX}.genes.all.gff3 \
  | awk -F'\t' '$3=="CDS" {
      split($9,a,"="); 
      CDS[a[2]]=CDS[a[2]]+(($5-$4)+1)
    } END {
      for (i in CDS){
        split(i,ii,"\\.t")
        print ii[1]"\t"i"\t"CDS[i]
    } }' \
  > ${PREFIX}.genes.all.gff3.CDS_lengths.tsv

cat ${PREFIX}.genes.all.gff3.CDS_lengths.tsv \
  | sort -k1,1 -k3,3nr \
  | awk '!seen[$1]++{print $2}' \
  > ${PREFIX}.genes.all.gff3.longest.txt

cat ${PREFIX}.genes.all.gff3 \
  | awk -F'\t' '{split($9,a,"="); print $0"\t"a[2]}' \
  | ~/scripts/grepf_column.py -c 10 -f ${PREFIX}.genes.all.gff3.longest.txt \
  | cut -f 1-9 \
  > ${PREFIX}.genes.gff3

gffread -S \
  -g ${PREFIX}.assembly.fasta \
  -x ${PREFIX}.genes.cds.fna.tmp \
  -y ${PREFIX}.genes.pep.faa.tmp \
  ${PREFIX}.genes.gff3
```

**CDS**

```bash
cat ${PREFIX}.genes.cds.fna.tmp \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.genes.cds.fna
rm ${PREFIX}.genes.cds.fna.tmp
```

**PEP**

```bash
cat ${PREFIX}.genes.pep.faa.tmp \
  | reformat.sh fastawrap=0 iupacToN=f ignorejunk=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.genes.pep.faa
rm ${PREFIX}.genes.pep.faa.tmp
```

## *Acropora millepora* JSIDv2.1

From Indonesia (isolate JS-1)

Version: 2.1

Paper: https://doi.org/10.1126/science.aba4674

Website: https://www.ncbi.nlm.nih.gov/assembly/GCF_013753865.1

```bash
# Download genome+gff3+CDS+PEP from NCBI
tar -xvf genome_assemblies_cds_fasta.tar
tar -xvf genome_assemblies_prot_fasta.tar
tar -xvf genome_assemblies_genome_gff.tar
tar -xvf genome_assemblies_genome_fasta.tar
```

```bash
PREFIX="Acropora_millepora_JSIDv2.1"
```

**Genome**

```bash
zcat ncbi-genomes-2022-02-07/GCF_013753865.1_Amil_v2.1_genomic.fna.gz \
  | sed -e "s/>/>${PREFIX}___/" \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.assembly.fasta
```

**Get longest ORF per gene**

```bash
zcat ncbi-genomes-2022-02-07/GCF_013753865.1_Amil_v2.1_genomic.gff.gz \
  | awk -F'\t' '$3=="CDS"' \
  | sed -e "s/[^\t]*Dbxref=GeneID:\([^,]*\).*Name=\([^;]*\).*/${PREFIX}___\2\t\1/" \
  > ncbi-genomes-2022-02-07.CDS_features.tsv
  
cat ncbi-genomes-2022-02-07.CDS_features.tsv \
  | awk -F'\t' '$3=="CDS" {
      CDS[$9]=CDS[$9]+(($5-$4)+1)
      GENE[$9]=$10
    } END {
      for (i in CDS){
        print GENE[i]"\t"i"\t"CDS[i]
    } }' \
  > ncbi-genomes-2022-02-07.CDS_lengths.tsv

cat ncbi-genomes-2022-02-07.CDS_lengths.tsv \
  | sort -k1,1 -k3,3nr \
  | awk '!seen[$1]++{print $2}' \
  > ncbi-genomes-2022-02-07.longest.txt
```

**CDS**

```bash
zcat ncbi-genomes-2022-02-07/GCF_013753865.1_Amil_v2.1_cds_from_genomic.fna.gz \
  | sed -e "s/.*protein_id=\([^]]*\)\].*/>${PREFIX}___\1/" \
  | ~/scripts/grepf_fasta.py -f ncbi-genomes-2022-02-07.longest.txt \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.genes.cds.fna
```

**PEP**

```bash
zcat ncbi-genomes-2022-02-07/GCF_013753865.1_Amil_v2.1_protein.faa.gz \
  | sed -e "s/>/>${PREFIX}___/" \
  | ~/scripts/grepf_fasta.py -f ncbi-genomes-2022-02-07.longest.txt \
  | reformat.sh fastawrap=0 iupacToN=f ignorejunk=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.genes.pep.faa
```

**GFF**

```bash
zcat ncbi-genomes-2022-02-07/GCF_013753865.1_Amil_v2.1_genomic.gff.gz \
  | awk -F'\t' '$3=="CDS"' \
  | sed -e "s/[^\t]*Name=\([^;]*\).*/Parent=${PREFIX}___\1\t${PREFIX}___\1/" \
  | sed -e "s/^/${PREFIX}___/" \
  | ~/scripts/grepf_column.py -c 10 -f ncbi-genomes-2022-02-07.longest.txt \
  | cut -f 1-9 \
  | gffread \
  | awk '$1!~"^#"' \
  | awk 'BEGIN{OFS=FS="\t"} { if($3!="CDS"){print} else {print; $3="exon"; print} }' \
  > ${PREFIX}.genes.gff3
```

## *Acropora muricata* SLJPv1

From Sekisei Lagoon, Okinawa, Japan

Version: 1.0

Paper: https://doi.org/10.1093/molbev/msaa216

Website: https://marinegenomics.oist.jp/amur/viewer/download?project_id=94

```bash
wget https://marinegenomics.oist.jp/amur/download/amur.fasta.gz
wget https://marinegenomics.oist.jp/amur/download/amur.gff.gz
```

```bash
PREFIX="Acropora_muricata_SLJPv1"
```

**Genome**

```bash
zcat amur.fasta.gz \
  | sed -e 's/_pilon$//' \
  | sed -e "s/>/>${PREFIX}___/" \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.assembly.fasta
```

**GFF**

```bash
zcat amur.gff.gz \
  | sed -e "s/^/${PREFIX}___/" -e 's/_pilon//' \
  | sed -e "s/ID=amur_/ID=${PREFIX}___/" \
  | sed -e "s/Parent=amur_/Parent=${PREFIX}___/" \
  | sed -e "s/Name=amur_/Name=${PREFIX}___/" \
  | awk -F'\t' '$3=="mRNA" || $3=="CDS" || $3=="exon"' \
  | gffread \
  | awk '$1!~"^#"' \
  | sed -e 's/;geneID=.*//' \
  > ${PREFIX}.genes.all.gff3

cat ${PREFIX}.genes.all.gff3 \
  | awk -F'\t' '$3=="CDS" {
      split($9,a,"="); 
      CDS[a[2]]=CDS[a[2]]+(($5-$4)+1)
    } END {
      for (i in CDS){
        split(i,ii,"\\.t")
        print ii[1]"\t"i"\t"CDS[i]
    } }' \
  > ${PREFIX}.genes.all.gff3.CDS_lengths.tsv

cat ${PREFIX}.genes.all.gff3.CDS_lengths.tsv \
  | sort -k1,1 -k3,3nr \
  | awk '!seen[$1]++{print $2}' \
  > ${PREFIX}.genes.all.gff3.longest.txt

cat ${PREFIX}.genes.all.gff3 \
  | awk -F'\t' '{split($9,a,"="); print $0"\t"a[2]}' \
  | ~/scripts/grepf_column.py -c 10 -f ${PREFIX}.genes.all.gff3.longest.txt \
  | cut -f 1-9 \
  > ${PREFIX}.genes.gff3

gffread -S \
  -g ${PREFIX}.assembly.fasta \
  -x ${PREFIX}.genes.cds.fna.tmp \
  -y ${PREFIX}.genes.pep.faa.tmp \
  ${PREFIX}.genes.gff3
```

**CDS**

```bash
cat ${PREFIX}.genes.cds.fna.tmp \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.genes.cds.fna
rm ${PREFIX}.genes.cds.fna.tmp
```

**PEP**

```bash
cat ${PREFIX}.genes.pep.faa.tmp \
  | reformat.sh fastawrap=0 iupacToN=f ignorejunk=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.genes.pep.faa
rm ${PREFIX}.genes.pep.faa.tmp
```

## *Acropora nasuta* SLJPv1

From Sekisei Lagoon, Okinawa, Japan

Version: 1.0

Paper: https://doi.org/10.1093/molbev/msaa216

Website: https://marinegenomics.oist.jp/anas/viewer/download?project_id=95

```bash
wget https://marinegenomics.oist.jp/anas/download/anas.fasta.gz
wget https://marinegenomics.oist.jp/anas/download/anas.gff.gz
```

```bash
PREFIX="Acropora_nasuta_SLJPv1"
```

**Genome**

```bash
zcat anas.fasta.gz \
  | sed -e 's/_pilon$//' \
  | sed -e "s/>/>${PREFIX}___/" \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.assembly.fasta
```

**GFF**

```bash
zcat anas.gff.gz \
  | sed -e "s/^/${PREFIX}___/" -e 's/_pilon//' \
  | sed -e "s/ID=anas_/ID=${PREFIX}___/" \
  | sed -e "s/Parent=anas_/Parent=${PREFIX}___/" \
  | sed -e "s/Name=anas_/Name=${PREFIX}___/" \
  | awk -F'\t' '$3=="mRNA" || $3=="CDS" || $3=="exon"' \
  | gffread \
  | awk '$1!~"^#"' \
  | sed -e 's/;geneID=.*//' \
  > ${PREFIX}.genes.all.gff3

cat ${PREFIX}.genes.all.gff3 \
  | awk -F'\t' '$3=="CDS" {
      split($9,a,"="); 
      CDS[a[2]]=CDS[a[2]]+(($5-$4)+1)
    } END {
      for (i in CDS){
        split(i,ii,"\\.t")
        print ii[1]"\t"i"\t"CDS[i]
    } }' \
  > ${PREFIX}.genes.all.gff3.CDS_lengths.tsv

cat ${PREFIX}.genes.all.gff3.CDS_lengths.tsv \
  | sort -k1,1 -k3,3nr \
  | awk '!seen[$1]++{print $2}' \
  > ${PREFIX}.genes.all.gff3.longest.txt

cat ${PREFIX}.genes.all.gff3 \
  | awk -F'\t' '{split($9,a,"="); print $0"\t"a[2]}' \
  | ~/scripts/grepf_column.py -c 10 -f ${PREFIX}.genes.all.gff3.longest.txt \
  | cut -f 1-9 \
  > ${PREFIX}.genes.gff3

gffread -S \
  -g ${PREFIX}.assembly.fasta \
  -x ${PREFIX}.genes.cds.fna.tmp \
  -y ${PREFIX}.genes.pep.faa.tmp \
  ${PREFIX}.genes.gff3
```

**CDS**

```bash
cat ${PREFIX}.genes.cds.fna.tmp \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.genes.cds.fna
rm ${PREFIX}.genes.cds.fna.tmp
```

**PEP**

```bash
cat ${PREFIX}.genes.pep.faa.tmp \
  | reformat.sh fastawrap=0 iupacToN=f ignorejunk=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.genes.pep.faa
rm ${PREFIX}.genes.pep.faa.tmp
```

## *Acropora selago* SLJPv1

From Sekisei Lagoon, Okinawa, Japan

Version: 1.0

Paper: https://doi.org/10.1093/molbev/msaa216

Website: https://marinegenomics.oist.jp/asel/viewer/info?project_id=96

```bash
wget https://marinegenomics.oist.jp/asel/download/asel.fasta.gz
wget https://marinegenomics.oist.jp/asel/download/asel.gff.gz
```

```bash
PREFIX="Acropora_selago_SLJPv1"
```

**Genome**

```bash
zcat asel.fasta.gz \
  | sed -e 's/_pilon$//' \
  | sed -e "s/>/>${PREFIX}___/" \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.assembly.fasta
```

**GFF**

```bash
zcat asel.gff.gz \
  | sed -e "s/^/${PREFIX}___/" -e 's/_pilon//' \
  | sed -e "s/ID=asel_/ID=${PREFIX}___/" \
  | sed -e "s/Parent=asel_/Parent=${PREFIX}___/" \
  | sed -e "s/Name=asel_/Name=${PREFIX}___/" \
  | awk -F'\t' '$3=="mRNA" || $3=="CDS" || $3=="exon"' \
  | gffread \
  | awk '$1!~"^#"' \
  | sed -e 's/;geneID=.*//' \
  > ${PREFIX}.genes.all.gff3

cat ${PREFIX}.genes.all.gff3 \
  | awk -F'\t' '$3=="CDS" {
      split($9,a,"="); 
      CDS[a[2]]=CDS[a[2]]+(($5-$4)+1)
    } END {
      for (i in CDS){
        split(i,ii,"\\.t")
        print ii[1]"\t"i"\t"CDS[i]
    } }' \
  > ${PREFIX}.genes.all.gff3.CDS_lengths.tsv

cat ${PREFIX}.genes.all.gff3.CDS_lengths.tsv \
  | sort -k1,1 -k3,3nr \
  | awk '!seen[$1]++{print $2}' \
  > ${PREFIX}.genes.all.gff3.longest.txt

cat ${PREFIX}.genes.all.gff3 \
  | awk -F'\t' '{split($9,a,"="); print $0"\t"a[2]}' \
  | ~/scripts/grepf_column.py -c 10 -f ${PREFIX}.genes.all.gff3.longest.txt \
  | cut -f 1-9 \
  > ${PREFIX}.genes.gff3

gffread -S \
  -g ${PREFIX}.assembly.fasta \
  -x ${PREFIX}.genes.cds.fna.tmp \
  -y ${PREFIX}.genes.pep.faa.tmp \
  ${PREFIX}.genes.gff3
```

**CDS**

```bash
cat ${PREFIX}.genes.cds.fna.tmp \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.genes.cds.fna
rm ${PREFIX}.genes.cds.fna.tmp
```

**PEP**

```bash
cat ${PREFIX}.genes.pep.faa.tmp \
  | reformat.sh fastawrap=0 iupacToN=f ignorejunk=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.genes.pep.faa
rm ${PREFIX}.genes.pep.faa.tmp
```

## *Acropora tenuis* SLJPv1

From Sekisei Lagoon, Okinawa, Japan

Version: 1.0

Paper: https://doi.org/10.1093/molbev/msaa216

Website: https://marinegenomics.oist.jp/aten/viewer/info?project_id=97

```bash
wget https://marinegenomics.oist.jp/aten/download/aten.fasta.gz
wget https://marinegenomics.oist.jp/aten/download/aten.gff.gz
wget https://marinegenomics.oist.jp/aten/download/aten.mrna.gz
wget https://marinegenomics.oist.jp/aten/download/aten.prot.t1.fa.gz
```

```bash
PREFIX="Acropora_tenuis_SLJPv1"
```

**Genome**

```bash
zcat aten.fasta.gz \
  | sed -e 's/_pilon$//' \
  | sed -e "s/>/>${PREFIX}___/" \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.assembly.fasta
```

**GFF**

```bash
zcat aten.gff.gz \
  | sed -e "s/^/${PREFIX}___/" -e 's/_pilon//' \
  | sed -e "s/ID=aten_/ID=${PREFIX}___/" \
  | sed -e "s/Parent=aten_/Parent=${PREFIX}___/" \
  | sed -e "s/Name=aten_/Name=${PREFIX}___/" \
  | awk -F'\t' '$3=="mRNA" || $3=="CDS" || $3=="exon"' \
  | gffread \
  | awk '$1!~"^#"' \
  | sed -e 's/;geneID=.*//' \
    > ${PREFIX}.genes.all.gff3

cat ${PREFIX}.genes.all.gff3 \
  | awk -F'\t' '$3=="CDS" {
      split($9,a,"="); 
      CDS[a[2]]=CDS[a[2]]+(($5-$4)+1)
    } END {
      for (i in CDS){
        split(i,ii,"\\.t")
        print ii[1]"\t"i"\t"CDS[i]
    } }' \
  > ${PREFIX}.genes.all.gff3.CDS_lengths.tsv

cat ${PREFIX}.genes.all.gff3.CDS_lengths.tsv \
  | sort -k1,1 -k3,3nr \
  | awk '!seen[$1]++{print $2}' \
  > ${PREFIX}.genes.all.gff3.longest.txt

cat ${PREFIX}.genes.all.gff3 \
  | awk -F'\t' '{split($9,a,"="); print $0"\t"a[2]}' \
  | ~/scripts/grepf_column.py -c 10 -f ${PREFIX}.genes.all.gff3.longest.txt \
  | cut -f 1-9 \
  > ${PREFIX}.genes.gff3

gffread -S \
  -g ${PREFIX}.assembly.fasta \
  -x ${PREFIX}.genes.cds.fna.tmp \
  -y ${PREFIX}.genes.pep.faa.tmp \
  ${PREFIX}.genes.gff3
```

**CDS**

```bash
cat ${PREFIX}.genes.cds.fna.tmp \
  | sed -e "s/>aten_/>${PREFIX}___/" \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.genes.cds.fna
rm ${PREFIX}.genes.cds.fna.tmp
```

**PEP**

```bash
cat ${PREFIX}.genes.pep.faa.tmp \
  | sed -e "s/>aten_/>${PREFIX}___/" \
  | reformat.sh fastawrap=0 iupacToN=f ignorejunk=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.genes.pep.faa
rm ${PREFIX}.genes.pep.faa.tmp
```

## *Acropora yongei* SLJPv1

From Sekisei Lagoon, Okinawa, Japan

Version: 1.0

Paper: https://doi.org/10.1093/molbev/msaa216

Website: https://marinegenomics.oist.jp/ayon/viewer/info?project_id=98

```bash
wget https://marinegenomics.oist.jp/ayon/download/ayon.fasta.gz
wget https://marinegenomics.oist.jp/ayon/download/ayon.gff.gz
```

```bash
PREFIX="Acropora_yongei_SLJPv1"
```

**Genome**

```bash
zcat ayon.fasta.gz \
  | sed -e 's/_pilon$//' \
  | sed -e "s/>/>${PREFIX}___/" \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.assembly.fasta
```

**GFF**

```bash
zcat ayon.gff.gz \
  | sed -e "s/^/${PREFIX}___/" -e 's/_pilon//' \
  | sed -e "s/ID=ayon_/ID=${PREFIX}___/" \
  | sed -e "s/Parent=ayon_/Parent=${PREFIX}___/" \
  | sed -e "s/Name=ayon_/Name=${PREFIX}___/" \
  | awk -F'\t' '$3=="mRNA" || $3=="CDS" || $3=="exon"' \
  | gffread \
  | awk '$1!~"^#"' \
  | sed -e 's/;geneID=.*//' \
  > ${PREFIX}.genes.all.gff3

cat ${PREFIX}.genes.all.gff3 \
  | awk -F'\t' '$3=="CDS" {
      split($9,a,"="); 
      CDS[a[2]]=CDS[a[2]]+(($5-$4)+1)
    } END {
      for (i in CDS){
        split(i,ii,"\\.t")
        print ii[1]"\t"i"\t"CDS[i]
    } }' \
  > ${PREFIX}.genes.all.gff3.CDS_lengths.tsv

cat ${PREFIX}.genes.all.gff3.CDS_lengths.tsv \
  | sort -k1,1 -k3,3nr \
  | awk '!seen[$1]++{print $2}' \
  > ${PREFIX}.genes.all.gff3.longest.txt

cat ${PREFIX}.genes.all.gff3 \
  | awk -F'\t' '{split($9,a,"="); print $0"\t"a[2]}' \
  | ~/scripts/grepf_column.py -c 10 -f ${PREFIX}.genes.all.gff3.longest.txt \
  | cut -f 1-9 \
  > ${PREFIX}.genes.gff3

gffread -S \
  -g ${PREFIX}.assembly.fasta \
  -x ${PREFIX}.genes.cds.fna.tmp \
  -y ${PREFIX}.genes.pep.faa.tmp \
  ${PREFIX}.genes.gff3
```

**CDS**

```bash
cat ${PREFIX}.genes.cds.fna.tmp \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.genes.cds.fna
rm ${PREFIX}.genes.cds.fna.tmp
```

**PEP**

```bash
cat ${PREFIX}.genes.pep.faa.tmp \
  | reformat.sh fastawrap=0 iupacToN=f ignorejunk=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.genes.pep.faa
rm ${PREFIX}.genes.pep.faa.tmp
```

## *Astreopora myriophthalma* SLJPv2

From Sekisei Lagoon, Okinawa, Japan

Version: 2.0

Paper (v1; originally published): https://doi.org/10.1093/molbev/msaa216

Website: https://marinegenomics.oist.jp/astreopora/viewer/info?project_id=99

Paper (v2; new gene models): https://doi.org/10.21203/rs.3.rs-944849/v1

Website: Supp Data associated with paper

```bash
wget https://marinegenomics.oist.jp/astreopora/download/asteropora.fasta.gz
wget https://marinegenomics.oist.jp/astreopora/download/asteropora.gff.gz
# Download Supplementary Data 5 from paper
# Supplementary Data S5. Gene models for Astreopora myriophthalma in GTF format.
wget https://assets.researchsquare.com/files/rs-944849/v1/106ef38c9bf204c05f5342eb.gz
```

```bash
PREFIX="Astreopora_myriophthalma_SLJPv2"
```

**Genome**

```bash
zcat  asteropora.fasta.gz \
  | sed -e 's/_pilon$//' \
  | sed -e "s/>/>${PREFIX}___/" \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.assembly.fasta
```

**GFF**

```bash
zcat SupplementaryData5.gz \
  | awk -F'\t' '$3=="transcript" || $3=="CDS" || $3=="exon"' \
  | sed -e "s/^s/${PREFIX}___sc000/" -e "s/astr_s/${PREFIX}___sc000/g" \
  | gffread \
  | awk '$1!~"^#"' \
  | sed -e 's/;geneID=.*//' \
  > ${PREFIX}.genes.gff3

gffread -S \
  -g ${PREFIX}.assembly.fasta \
  -x ${PREFIX}.genes.cds.fna.tmp \
  -y ${PREFIX}.genes.pep.faa.tmp \
  ${PREFIX}.genes.gff3
```

**CDS**

```bash
cat ${PREFIX}.genes.cds.fna.tmp \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.genes.cds.fna
rm ${PREFIX}.genes.cds.fna.tmp
```

**PEP**

```bash
cat ${PREFIX}.genes.pep.faa.tmp \
  | reformat.sh fastawrap=0 iupacToN=f ignorejunk=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.genes.pep.faa
rm ${PREFIX}.genes.pep.faa.tmp
```

## *Fungia* spp. REEFv1

From: Unknown

Version: 1.0

Paper: https://doi.org/10.3389/fmars.2015.00068

Website: http://ffun.reefgenomics.org

```bash
wget http://ffun.reefgenomics.org/download/ffun_final_1.0.fasta.gz
wget http://ffun.reefgenomics.org/download/ffun_1.0.transcripts.fasta.gz
wget http://ffun.reefgenomics.org/download/ffun_1.0.proteins.fasta.gz
wget http://ffun.reefgenomics.org/download/ffun_1.0.genes.gff3.gz
```

```bash
PREFIX="Fungia_spp_REEFv1"
```

**Genome**

```bash
zcat ffun_final_1.0.fasta.gz \
  | sed -e "s/>/>${PREFIX}___/" \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.assembly.fasta
```

**CDS**

```bash
zcat ffun_1.0.transcripts.fasta.gz \
  | sed -e "s/>ffun1./>${PREFIX}___/" \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.genes.cds.fna
```

**PEP**

```bash
zcat ffun_1.0.proteins.fasta.gz \
  | sed -e "s/>ffun1./>${PREFIX}___/" \
  | reformat.sh fastawrap=0 iupacToN=f ignorejunk=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.genes.pep.faa
```

**GFF**

```bash
zcat ffun_1.0.genes.gff3.gz \
  | awk '$1!~"^#"' \
  | sed -e "s/^/${PREFIX}___/" \
  | sed -e "s/ID=ffun1./ID=${PREFIX}___/" \
  | sed -e "s/Parent=ffun1./Parent=${PREFIX}___/" \
  | sed -e "s/Name=ffun1./Name=${PREFIX}___/" \
  | awk -F'\t' '$3=="mRNA" || $3=="CDS" || $3=="exon"' \
  | gffread \
  | awk '$1!~"^#"' \
  | awk 'BEGIN{OFS=FS="\t"} { if($3=="mRNA") {$3="transcript"; print} else {print} }' \
  > ${PREFIX}.genes.gff3
```

## *Galaxea fascicularis* REEFv1

From: Unknown

Version: 1.0

Paper: https://doi.org/10.3389/fmars.2015.00068

Website: http://gfas.reefgenomics.org

```bash
wget http://gfas.reefgenomics.org/download/gfas_final_1.0.fasta.gz
wget http://gfas.reefgenomics.org/download/gfas_1.0.transcripts.fasta.gz
wget http://gfas.reefgenomics.org/download/gfas_1.0.proteins.fasta.gz
wget http://gfas.reefgenomics.org/download/gfas_1.0.genes.gff3.gz
```

```bash
PREFIX="Galaxea_fascicularis_REEFv1"
```

**Genome**

```bash
zcat gfas_final_1.0.fasta.gz \
  | sed -e "s/>/>${PREFIX}___/" \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.assembly.fasta
```

**CDS**

```bash
zcat gfas_1.0.transcripts.fasta.gz \
  | sed -e "s/>gfas1./>${PREFIX}___/" \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.genes.cds.fna
```

**PEP**

```bash
zcat gfas_1.0.proteins.fasta.gz \
  | sed -e "s/>gfas1./>${PREFIX}___/" \
  | reformat.sh fastawrap=0 iupacToN=f ignorejunk=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.genes.pep.faa
```

**GFF**

```bash
zcat gfas_1.0.genes.gff3.gz \
  | awk '$1!~"^#"' \
  | sed -e "s/^/${PREFIX}___/" \
  | sed -e "s/ID=gfas1./ID=${PREFIX}___/" \
  | sed -e "s/Parent=gfas1./Parent=${PREFIX}___/" \
  | sed -e "s/Name=gfas1./Name=${PREFIX}___/" \
  | awk -F'\t' '$3=="mRNA" || $3=="CDS" || $3=="exon"' \
  | gffread \
  | awk '$1!~"^#"' \
  | awk 'BEGIN{OFS=FS="\t"} { if($3=="mRNA") {$3="transcript"; print} else {print} }' \
  > ${PREFIX}.genes.gff3
```

## *Goniastrea aspera* REEFv1

From: Unknown

Version: 1.0

Paper: https://doi.org/10.3389/fmars.2015.00068

Website: http://gasp.reefgenomics.org

```bash
wget http://gasp.reefgenomics.org/download/gasp_final_1.0.fasta.gz
wget http://gasp.reefgenomics.org/download/gasp_1.0.transcripts.fasta.gz
wget http://gasp.reefgenomics.org/download/gasp_1.0.proteins.fasta.gz
wget http://gasp.reefgenomics.org/download/gasp_1.0.genes.gff3.gz
```

```bash
PREFIX="Goniastrea_aspera_REEFv1"
```

**Genome**

```bash
zcat gasp_final_1.0.fasta.gz \
  | sed -e "s/>/>${PREFIX}___/" \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.assembly.fasta
```

**CDS**

```bash
zcat gasp_1.0.transcripts.fasta.gz \
  | sed -e "s/>gasp1./>${PREFIX}___/" \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.genes.cds.fna
```

**PEP**

```bash
zcat gasp_1.0.proteins.fasta.gz \
  | sed -e "s/>gasp1./>${PREFIX}___/" \
  | reformat.sh fastawrap=0 iupacToN=f ignorejunk=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.genes.pep.faa
```

**GFF**

```bash
zcat gasp_1.0.genes.gff3.gz \
  | awk '$1!~"^#"' \
  | sed -e "s/^/${PREFIX}___/" \
  | sed -e "s/ID=gasp1./ID=${PREFIX}___/" \
  | sed -e "s/Parent=gasp1./Parent=${PREFIX}___/" \
  | sed -e "s/Name=gasp1./Name=${PREFIX}___/" \
  | awk -F'\t' '$3=="mRNA" || $3=="CDS" || $3=="exon"' \
  | gffread \
  | awk '$1!~"^#"' \
  | awk 'BEGIN{OFS=FS="\t"} { if($3=="mRNA") {$3="transcript"; print} else {print} }' \
  > ${PREFIX}.genes.gff3
```

## *Montastrea cavernosa* FGUSv2

From: Flower Garden Banks, United States

Version: 2.0

Paper: https://matzlab.weebly.com/data--code.html

Website: https://www.dropbox.com/s/yfqefzntt896xfz/Mcavernosa_genome.tgz

```bash
# Download tar-zip file from Dropbox: https://www.dropbox.com/s/yfqefzntt896xfz/Mcavernosa_genome.tgz
```

```bash
PREFIX="Montastrea_cavernosa_FGUSv2"
```

**Genome**

```bash
cat Mcav_genome/Mcavernosa_July2018.fasta \
  | sed -e "s/>/>${PREFIX}___/" -e 's/ .*//' \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.assembly.fasta
```

**GFF**

Need to fix the gff records for one gene ("Mcavernosa33090-RA") which appears to be missing the scaffold ID for some of its features.

**FROM**

```bash
maker	mRNA	1170	1	1	.	ID=Mcavernosa33090-RA;Parent=Mcavernosa33090;Name=Mcavernosa33090-RA;Alias=maker-NewChr39-augustus-gene-95.19-mRNA-1;_AED=0.25;_QI=0|0|0|1|0.5|0.33|3|0|259;_eAED=0.65;
xfSc0003276	maker	exon	873	1171	.	-	.	ID=Mcavernosa33090-RA:exon:1901;Parent=Mcavernosa33090-RA;
xfSc0003276	maker	exon	459	559	.	-	.	ID=Mcavernosa33090-RA:exon:1900;Parent=Mcavernosa33090-RA;
maker	exon	210	1	1	.	ID=Mcavernosa33090-RA:exon:1899;Parent=Mcavernosa33090-RA;
xfSc0003276	maker	CDS	873	1171	.	-	0	ID=Mcavernosa33090-RA:cds;Parent=Mcavernosa33090-RA;
xfSc0003276	maker	CDS	459	559	.	-	1	ID=Mcavernosa33090-RA:cds;Parent=Mcavernosa33090-RA;
maker	CDS	210	1	1	2	ID=Mcavernosa33090-RA:cds;Parent=Mcavernosa33090-RA;
```

**TO**

```bash
xfSc0003276     maker   gene    1       1170    .       -       .       ID=Mcavernosa33090;Name=Mcavernosa33090;Alias=maker-NewChr39-augustus-gene-95.19;
xfSc0003276     maker   mRNA    1       1170    .       -       .       ID=Mcavernosa33090-RA;Parent=Mcavernosa33090;Name=Mcavernosa33090-RA;Alias=maker-NewChr39-augustus-gene-95.19-mRNA-1;_AED=0.25;_QI=0|0|0|1|0.5|0.33|3|0|259;_eAED=0.65;
xfSc0003276     maker   exon    873     1171    .       -       .       ID=Mcavernosa33090-RA:exon:1901;Parent=Mcavernosa33090-RA;
xfSc0003276     maker   exon    459     559     .       -       .       ID=Mcavernosa33090-RA:exon:1900;Parent=Mcavernosa33090-RA;
xfSc0003276     maker   exon    1       210     .       -       .       ID=Mcavernosa33090-RA:exon:1899;Parent=Mcavernosa33090-RA;
xfSc0003276     maker   CDS     873     1171    .       -       0       ID=Mcavernosa33090-RA:cds;Parent=Mcavernosa33090-RA;
xfSc0003276     maker   CDS     459     559     .       -       1       ID=Mcavernosa33090-RA:cds;Parent=Mcavernosa33090-RA;
xfSc0003276     maker   CDS     1       210     .       -       2       ID=Mcavernosa33090-RA:cds;Parent=Mcavernosa33090-RA;
```

```bash
#cp Mcav_genome/Mcavernosa_annotation/Mcavernosa.maker.coding.gff3 Mcavernosa.maker.coding.fixed.gff3

cat Mcavernosa.maker.coding.fixed.gff3 \
  | awk '$1!~"^#"' \
  | sed -e "s/^/${PREFIX}___/" \
  | sed -e "s/ID=Mcavernosa/ID=${PREFIX}___/" \
  | sed -e "s/Parent=Mcavernosa/Parent=${PREFIX}___/" \
  | sed -e "s/Name=Mcavernosa/Name=${PREFIX}___/" \
  | awk -F'\t' '$3=="mRNA" || $3=="CDS" || $3=="exon"' \
  | gffread \
  | awk '$1!~"^#"' \
  | awk 'BEGIN{OFS=FS="\t"} { if($3=="mRNA") {$3="transcript"; print} else {print} }' \
  > ${PREFIX}.genes.gff3
```

**CDS**

```bash
cat Mcav_genome/Mcavernosa_annotation/Mcavernosa.maker.transcripts.fasta \
  | sed -e "s/>Mcavernosa/>${PREFIX}___/" -e 's/ .*//' \
  | ~/scripts/grepf_fasta.py -f <(cat ${PREFIX}.genes.gff3 | awk -F'\t' '$3=="transcript" {print $9}' | sed -e 's/ID=//') \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.genes.cds.fna
```

**PEP**

```bash
cat Mcav_genome/Mcavernosa_annotation/Mcavernosa.maker.proteins.fasta \
  | sed -e "s/>Mcavernosa/>${PREFIX}___/" -e 's/ .*//' \
  | ~/scripts/grepf_fasta.py -f <(cat ${PREFIX}.genes.gff3 | awk -F'\t' '$3=="transcript" {print $9}' | sed -e 's/ID=//') \
  | reformat.sh fastawrap=0 iupacToN=f ignorejunk=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.genes.pep.faa
```

## *Montipora cactus* SLJPv2

From Sekisei Lagoon, Okinawa, Japan

Version: 2.0

Paper (v1; originally published): https://doi.org/10.1093/molbev/msaa216

Website: https://marinegenomics.oist.jp/mcac/viewer/info?project_id=100

Paper (v2; haplotig reduced assembly + new gene models): https://doi.org/10.21203/rs.3.rs-944849/v1

Website: Supp Data associated with paper

```bash
wget https://marinegenomics.oist.jp/mcac/download/mcac.fasta.gz
wget https://marinegenomics.oist.jp/mcac/download/mcac.gff.gz
# Download Supplementary Data 2 & 3 from paper
# Supplementary Data S2. Summary of retained scaffolds after genome assembly curation for Montipora cactus and M. efflorescens.
# Supplementary Data S3. Gene models for M. cactus in GTF format.
wget https://assets.researchsquare.com/files/rs-944849/v1/ba942c9d56b998ca59001f33.xlsx
wget https://assets.researchsquare.com/files/rs-944849/v1/29fad061ca05feb6f2368592.gz
```

```bash
PREFIX="Montipora_cactus_SLJPv2"
```

**Genome**

Extract just the scaffolds that survived haplotig purging (from `SupplementaryData2-Mcactus.txt`)

```bash
gunzip -c  mcac.fasta.gz \
  | sed -e 's/_pilon$//' \
  | ./../../../../02_Scripts/grepf_fasta.py -f <(awk -F'\t' 'NR>1 {print $1}' \
       SupplementaryData2-Mcactus.txt) \
  | sed -e "s/>/>${PREFIX}___/" \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.assembly.fasta
```

**GFF**

```bash
zcat SupplementaryData3.gz \
  | awk -F'\t' '$3=="transcript" || $3=="CDS" || $3=="exon"' \
  | sed -e "s/^s/${PREFIX}___sc000/" -e "s/mcac_s/${PREFIX}___sc000/g" \
  | gffread \
  | awk '$1!~"^#"' \
  | sed -e 's/;geneID=.*//' \
  > ${PREFIX}.genes.gff3

gffread -S \
  -g ${PREFIX}.assembly.fasta \
  -x ${PREFIX}.genes.cds.fna.tmp \
  -y ${PREFIX}.genes.pep.faa.tmp \
  ${PREFIX}.genes.gff3
```

**CDS**

```bash
cat ${PREFIX}.genes.cds.fna.tmp \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.genes.cds.fna
rm ${PREFIX}.genes.cds.fna.tmp
```

**PEP**

```bash
cat ${PREFIX}.genes.pep.faa.tmp \
  | reformat.sh fastawrap=0 iupacToN=f ignorejunk=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.genes.pep.faa
rm ${PREFIX}.genes.pep.faa.tmp
```

**Reads**

Download data from BioProject PRJDB8519 BioSample SAMD00186143. Five Runs are associated with the BioSample but four of them are mate-pair libraries. Download one Run (DRR194228) that is Illumina paired-end data (generated on a HiSeq 2500 platform). 

```bash
SRR="DRR194228"
fasterq-dump --threads 6 --split-3 --skip-technical --progress --outdir . ${SRR}
reformat.sh deleteinput=t verifypaired=t trimreaddescription=t addslash=t spaceslash=f \
  in=${SRR}_1.fastq in2=${SRR}_2.fastq out=${SRR}_R1.fastq.gz out2=${SRR}_R2.fastq.gz
./statswrapper.sh ${SRR}_R1.fastq.gz ${SRR}_R2.fastq.gz 2>read_stats.txt
```

## *Montipora capitata* KBHIv3

From Kāne‘ohe Bay, Oʻahu, Hawaii

Version: 3.0

Genome files from JunMo (2021/12/06)

```bash
PREFIX="Montipora_capitata_KBHIv3"
```

**Genome**

```bash
cat Montipora_capitata_HIv1_assembly.fasta \
  | sed -e "s/Montipora_capitata_HIv1___/${PREFIX}___/" -e 's/___length.*//' \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.assembly.fasta
samtools faidx ${PREFIX}.assembly.fasta
```

**Separate chromosome**

Separate scaffolds into putative chromosomal and extrachromosomal sequences.

```bash
# Chromosomes
seqkit faidx ${PREFIX}.assembly.fasta --line-width 0 \
  --infile-list <(awk 'NR<=14{print $1}' ${PREFIX}.assembly.fasta.fai) \
  > ${PREFIX}.chromosomes.assembly.fasta
# Extrachromosomal
seqkit faidx ${PREFIX}.assembly.fasta --line-width 0 \
  --infile-list <(awk 'NR>14{print $1}' ${PREFIX}.assembly.fasta.fai) \
  > ${PREFIX}.extrachromosomal.assembly.fasta
```

**CDS**

```bash
cat Montipora_capitata_CDS.fasta \
  | sed -e "s/Montipora_capitata../${PREFIX}___/" \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.genes.cds.fna
```

**PEP**

```bash
cat Montipora_capitata_proteins.fasta \
  | sed -e "s/Montipora_capitata../${PREFIX}___/" \
  | reformat.sh fastawrap=0 iupacToN=f ignorejunk=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.genes.pep.faa
```

**GFF**

```bash
cat Montipora_capitata_genemodels.gff \
  | sed -e "s/Montipora_capitata_HIv1___/${PREFIX}___/" -e 's/___length[^\t]*\t/\t/' \
  | sed -e "s/ID=/ID=${PREFIX}___/" \
  | sed -e "s/Parent=/Parent=${PREFIX}___/" \
  | awk -F'\t' '$3=="transcript" || $3=="CDS" || $3=="exon"' \
  | gffread \
  | awk '$1!~"^#"' \
  | awk 'BEGIN{OFS=FS="\t"} { if($3!="CDS"){print} else {print; $3="exon"; print} }' \
  > ${PREFIX}.genes.gff3
```

**Reads**

Download data from BioProject PRJNA509219. Data is short read Illumina data from a HiSeq 2500 platform. 

```bash
SRR="SRR8497577"
fasterq-dump --threads 6 --split-3 --skip-technical --progress --outdir . ${SRR}
reformat.sh deleteinput=t verifypaired=t trimreaddescription=t addslash=t spaceslash=f \
  in=${SRR}_1.fastq in2=${SRR}_2.fastq out=${SRR}_R1.fastq.gz out2=${SRR}_R2.fastq.gz
statswrapper.sh ${SRR}_R1.fastq.gz ${SRR}_R2.fastq.gz 1>read_stats.txt
```

## *Montipora capitata* WTHIv1.1

From Waiopae tide pools, Hawaii Island, Hawaii

Version: 1.1

Paper: https://doi.org/10.1093/gbe/evz135

Website: Genome: (NCBI: RDEB00000000); Gene models (USDA Ag Data Commons, doi: 10.15482/USDA.ADC/1503958)

```bash
# Genome assembly from NCBI
wget https://data.nal.usda.gov/system/files/Mcap_genemodels_1.1_aa.fas
wget https://data.nal.usda.gov/system/files/Mcap_genemodels_1.1.gff
wget https://data.nal.usda.gov/system/files/Mcap_genemodels_1.1_cds.fas
```

```bash
PREFIX="Montipora_capitata_WTHIv1.1"
```

**Genome**

```bash
gunzip -c GCA_006542545.1_Mcap_UHH_1.1_genomic.fna.gz \
  | awk '{ if($1~"^>") {gsub(",","",$7); print ">"$7} else {print} }' \
  | sed -e "s/>/>${PREFIX}___/" \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.assembly.fasta
```

**GFF**

Only take "CDS" features and let `gffread` add "transcript" features. Have to do this as "transcript" features arent formated correctly and `gffread` ignores them anyway. 

```bash
cat Mcap_genemodels_1.1.gff \
  | awk -F'\t' '$1!~"^#" && $2=="AUGUSTUS"' \
  | awk -F'\t' '$3=="CDS"' \
  | gffread \
  | awk '$1!~"^#"' \
  | sed -e "s/^/${PREFIX}___/" -e "s/=/=${PREFIX}___/g" \
  | awk 'BEGIN{OFS=FS="\t"} { if($3!="CDS"){print} else {print; $3="exon"; print} }' \
  | sed -e 's/;geneID=.*//' \
  > ${PREFIX}.genes.gff3
```

NOTE: xfSc0014868 appears to be in gff file but is missing from the genome assembly and when you use the gff to translate some of the genes into proteins ~100 of them have in-frame stop codons. This suggests that the version of the genome available online is different (only slightly) to the one used for the gene prediction.

**CDS**

```bash
cat Mcap_genemodels_1.1_cds.fas \
  | sed -e "s/>[^\.]*./>${PREFIX}___/" \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.genes.cds.fna
```

**PEP**

```bash
cat Mcap_genemodels_1.1_aa.fas \
  | sed -e "s/>/>${PREFIX}___/" \
  | reformat.sh fastawrap=0 iupacToN=f ignorejunk=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.genes.pep.faa
```

**Reads**

Download data from BioProject PRJNA495325. Data is short read Illumina data but generated from a Chromium 10x machine. 

```bash
SRR="SRR8068582"
fasterq-dump --threads 6 --split-3 --skip-technical --progress --outdir . ${SRR}
reformat.sh deleteinput=t verifypaired=t trimreaddescription=t addslash=t spaceslash=f \
  in=${SRR}_1.fastq in2=${SRR}_2.fastq out=${SRR}_R1.fastq.gz out2=${SRR}_R2.fastq.gz
statswrapper.sh ${SRR}_R1.fastq.gz ${SRR}_R2.fastq.gz 2>read_stats.txt
```

## *Montipora efflorescens* SLJPv2

From Sekisei Lagoon, Okinawa, Japan

Version: 2.0

Paper (v1; originally published): https://doi.org/10.1093/molbev/msaa216

Website: https://marinegenomics.oist.jp/meff/viewer/info?project_id=101

Paper (v2; haplotig reduced assembly + new gene models): https://doi.org/10.21203/rs.3.rs-944849/v1

Website: Supp Data associated with paper

```bash
wget https://marinegenomics.oist.jp/meff/download/meff.fasta.gz
wget https://marinegenomics.oist.jp/meff/download/meff.gff.gz
# Download Supplementary Data 2 & 4 from paper
# Supplementary Data S2. Summary of retained scaffolds after genome assembly curation for Montipora cactus and M. efflorescens.
# Supplementary Data S4. Gene models for M. efflorescens in GTF format.
wget https://assets.researchsquare.com/files/rs-944849/v1/ba942c9d56b998ca59001f33.xlsx
wget https://assets.researchsquare.com/files/rs-944849/v1/7e3d9454dd6bd4a52b5ac0d1.gz
```

```bash
PREFIX="Montipora_efflorescens_SLJPv2"
```

**Genome**

```bash
gunzip -c meff.fasta.gz \
  | sed -e 's/_pilon$//' \
  | ./../../../../02_Scripts/grepf_fasta.py -f <(awk -F'\t' 'NR>1 {print $1}' SupplementaryData2-Mefflorescens.txt) \
  | sed -e "s/>/>${PREFIX}___/" \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.assembly.fasta
```

**GFF**

```bash
zcat SupplementaryData4.gz \
  | sed -e "s/^s/${PREFIX}___sc000/" -e "s/meff_s/${PREFIX}___sc000/g" \
  | gffread \
  | awk '$1!~"^#"' \
  | sed -e 's/;geneID=.*//' \
  > ${PREFIX}.genes.gff3

gffread -S \
  -g ${PREFIX}.assembly.fasta \
  -x ${PREFIX}.genes.cds.fna.tmp \
  -y ${PREFIX}.genes.pep.faa.tmp \
  ${PREFIX}.genes.gff3
```

**CDS**

```bash
cat ${PREFIX}.genes.cds.fna.tmp \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.genes.cds.fna
rm ${PREFIX}.genes.cds.fna.tmp
```

**PEP**

```bash
cat ${PREFIX}.genes.pep.faa.tmp \
  | reformat.sh fastawrap=0 iupacToN=f ignorejunk=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.genes.pep.faa
rm ${PREFIX}.genes.pep.faa.tmp
```

**Reads**

Download data from BioProject PRJDB8519 BioSample SAMD00186142. Five Runs are associated with the BioSample but four of them are mate-pair libraries. Download one Run (DRR194233) that is Illumina paired-end data (generated on a HiSeq 2500 platform). 

```bash
SRR="DRR194233"
fasterq-dump --threads 6 --split-3 --skip-technical --progress --outdir . ${SRR}
reformat.sh deleteinput=t verifypaired=t trimreaddescription=t addslash=t spaceslash=f \
  in=${SRR}_1.fastq in2=${SRR}_2.fastq out=${SRR}_R1.fastq.gz out2=${SRR}_R2.fastq.gz
statswrapper.sh ${SRR}_R1.fastq.gz ${SRR}_R2.fastq.gz 2>read_stats.txt
```

## *Montipora* sp 1 aff. capitata ULFMv1

From Ulithi, Yap State, Federated States of Micronesia

Version: 1.0

Unpublished data from collaborators

```bash
# Genome fasta file from collaborators.
```

```bash
PREFIX="Montipora_sp1_aff_capitata_ULFMv1"
```

**Genome**

```bash
cat Montipora_sp_1_aff_capitata.fasta \
  | sed -e "s/>/>${PREFIX}___/" \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.assembly.fasta
```

**Rename BRAKER results**

```bash
sed -e 's/file_[0-9]*_file_[0-9]*_/Montipora_sp1_aff_capitata_ULFMv1___/g' \
  braker.gtf  > braker.renamed.gtf
sed -e 's/file_[0-9]*_file_[0-9]*_/Montipora_sp1_aff_capitata_ULFMv1___/g' \
  braker.gff3 > braker.renamed.gff3

conda activate BRAKER
~/programs/BRAKER/Augustus/scripts/getAnnoFastaFromJoingenes.py \
  --genome "${PREFIX}.assembly.fasta" \
  --gtf braker.renamed.gtf \
  --out braker.renamed
```

**GFF**

```bash
cat braker.renamed.gff3 \
  | awk -F'\t' '$3=="CDS"' \
  | gffread \
  | awk '$1!~"^#"' \
  | awk 'BEGIN{OFS=FS="\t"} { if($3!="CDS"){print} else {print; $3="exon"; print} }' \
  > ${PREFIX}.genes.gff3
```

**CDS**

```bash
cat braker.renamed.codingseq \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.genes.cds.fna
```

**PEP**

```bash
cat braker.renamed.aa \
  | sed -e 's/*$//' \
  | reformat.sh fastawrap=0 iupacToN=f ignorejunk=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.genes.pep.faa
```

## *Orbicella faveolata* FLUSv1

From Florida, United States (24.812697° N; 80.66925° W)

Version: 1.0

Paper: http://dx.doi.org/10.1016/j.cub.2016.09.039

Website: https://www.ncbi.nlm.nih.gov/assembly/GCF_002042975.1

```bash
# Download genome+gff3+CDS+PEP from NCBI
tar -xvf genome_assemblies_cds_fasta.tar
tar -xvf genome_assemblies_prot_fasta.tar
tar -xvf genome_assemblies_genome_gff.tar
tar -xvf genome_assemblies_genome_fasta.tar
```

```bash
PREFIX="Orbicella_faveolata_FLUSv1"
```

**Genome**

```bash
zcat ncbi-genomes-2022-02-08/GCF_002042975.1_ofav_dov_v1_genomic.fna.gz \
  | sed -e "s/>/>${PREFIX}___/" \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.assembly.fasta
```

**Get longest ORF per gene**

```bash
zcat ncbi-genomes-2022-02-08/GCF_002042975.1_ofav_dov_v1_genomic.gff.gz \
  | awk -F'\t' '$3=="CDS"' \
  | sed -e "s/[^\t]*GeneID:\([^,;]*\).*Name=\([^;]*\).*/${PREFIX}___\2\t\1/" \
  > ncbi-genomes-2022-02-08.CDS_features.tsv
  
cat ncbi-genomes-2022-02-08.CDS_features.tsv \
  | awk -F'\t' '$3=="CDS" {
      CDS[$9]=CDS[$9]+(($5-$4)+1)
      GENE[$9]=$10
    } END {
      for (i in CDS){
        print GENE[i]"\t"i"\t"CDS[i]
    } }' \
  > ncbi-genomes-2022-02-08.CDS_lengths.tsv

cat ncbi-genomes-2022-02-08.CDS_lengths.tsv \
  | sort -k1,1 -k3,3nr \
  | awk '!seen[$1]++{print $2}' \
  > ncbi-genomes-2022-02-08.longest.txt
```

**CDS**

```bash
zcat ncbi-genomes-2022-02-08/GCF_002042975.1_ofav_dov_v1_cds_from_genomic.fna.gz \
  | sed -e "s/.*protein_id=\([^]]*\)\].*/>${PREFIX}___\1/" \
  | ~/scripts/grepf_fasta.py -f ncbi-genomes-2022-02-08.longest.txt \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.genes.cds.fna
```

**PEP**

```bash
zcat ncbi-genomes-2022-02-08/GCF_002042975.1_ofav_dov_v1_protein.faa.gz \
  | sed -e "s/>/>${PREFIX}___/" \
  | ~/scripts/grepf_fasta.py -f ncbi-genomes-2022-02-08.longest.txt \
  | reformat.sh fastawrap=0 iupacToN=f ignorejunk=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.genes.pep.faa
```

**GFF**

```bash
zcat ncbi-genomes-2022-02-08/GCF_002042975.1_ofav_dov_v1_genomic.gff.gz \
  | awk -F'\t' '$3=="CDS"' \
  | sed -e "s/[^\t]*Name=\([^;]*\).*/Parent=${PREFIX}___\1\t${PREFIX}___\1/" \
  | sed -e "s/^/${PREFIX}___/" \
  | ~/scripts/grepf_column.py -c 10 -f ncbi-genomes-2022-02-08.longest.txt \
  | cut -f1-9 \
  | gffread \
  | awk '$1!~"^#"' \
  | awk 'BEGIN{OFS=FS="\t"} { if($3!="CDS"){print} else {print; $3="exon"; print} }' \
  > ${PREFIX}.genes.gff3
```

## *Pachyseris speciosa* OIAUv0.12

From: Orpheus Island Research Station, Australia

Version: 0.12

Paper: https://doi.org/10.1016/j.cub.2021.03.028

Website: http://pspe.reefgenomics.org

```bash
wget http://pspe.reefgenomics.org/download/pspe_final_0.12.fasta.gz
wget http://pspe.reefgenomics.org/download/pspe_0.12.maker_post_002.transcripts.fasta.gz
wget http://pspe.reefgenomics.org/download/pspe_0.12.maker_post_002.proteins.fasta.gz
wget http://pspe.reefgenomics.org/download/pspe_0.12.maker_post_002.genes.gff3.gz
```

```bash
PREFIX="Pachyseris_speciosa_OIAUv0.12"
```

**Genome**

```bash
zcat pspe_final_0.12.fasta.gz \
  | sed -e "s/>/>${PREFIX}___/" \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.assembly.fasta
```

**CDS**

```bash
zcat pspe_0.12.maker_post_002.transcripts.fasta.gz \
  | sed -e "s/>pspe_/>${PREFIX}___/" \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.genes.cds.fna
```

**PEP**

```bash
zcat pspe_0.12.maker_post_002.proteins.fasta.gz \
  | sed -e "s/>pspe_/>${PREFIX}___/" \
  | reformat.sh fastawrap=0 iupacToN=f ignorejunk=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.genes.pep.faa
```

**GFF**

```bash
zcat pspe_0.12.maker_post_002.genes.gff3.gz \
  | awk '$1!~"^#"' \
  | sed -e "s/^/${PREFIX}___/" \
  | sed -e "s/ID=pspe_/ID=${PREFIX}___/" \
  | sed -e "s/Parent=pspe_/Parent=${PREFIX}___/" \
  | sed -e "s/Name=pspe_/Name=${PREFIX}___/" \
  | awk -F'\t' '$3=="mRNA" || $3=="CDS" || $3=="exon"' \
  | gffread \
  | awk '$1!~"^#"' \
  | awk 'BEGIN{OFS=FS="\t"} { if($3=="mRNA") {$3="transcript"; print} else {print} }' \
  > ${PREFIX}.genes.gff3
```

## *Platygyra sinensis* SDTHv1

From: inshore reef in Sattahip district (Samaesarn subdistrict) located in the Gulf of Thailand (12°35′556″N, 100°57′508″E)

Version: 1.0

Paper: https://doi.org/10.3389/fmars.2021.732650

Website: Genome: (NCBI: JAHPZR000000000); Genes (https://www.nstda.or.th/noc/research-development/80-about-us/102-our-genome-assemblies.html)

```bash
# Genome assembly from NCBI
wget https://www.nstda.or.th/noc/images/GenomeData/coral/coral.gff3
wget https://www.nstda.or.th/noc/images/GenomeData/coral/coral_protein.fasta
wget https://www.nstda.or.th/noc/images/GenomeData/coral/coral_CDS.fasta
```

```bash
PREFIX="Platygyra_sinensis_SDTHv1"
```

**Genome**

```bash
gunzip -c ncbi-genomes-2022-04-17/GCA_019787425.1_ASM1978742v1_genomic.fna.gz \
  | awk '{ if($1~"^>") {gsub(",","",$6); print ">"$6} else {print} }' \
  | sed -e "s/>/>${PREFIX}___/" \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.assembly.fasta
```

**GFF**

Only take "CDS" features and let `gffread` add "transcript" features. Have to do this as "transcript" features arent formated correctly and `gffread` ignores them anyway. 

```bash
cat coral.gff3 \
  | awk -F'\t' '$1!~"^#"' \
  | awk -F'\t' '$3=="CDS"' \
  | gffread \
  | awk '$1!~"^#"' \
  | sed -e "s/^/${PREFIX}___/" -e "s/=/=${PREFIX}___/g" \
  | awk 'BEGIN{OFS=FS="\t"} { if($3!="CDS"){print} else {print; $3="exon"; print} }' \
  | sed -e 's/;geneID=.*//' \
  > ${PREFIX}.genes.gff3
```

**CDS**

```bash
cat coral_CDS.fasta \
  | sed -e "s/>/>${PREFIX}___/" \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.genes.cds.fna
```

**PEP**

```bash
cat coral_protein.fasta \
  | sed -e "s/>/>${PREFIX}___/" -e 's/*$//' \
  | reformat.sh fastawrap=0 iupacToN=f ignorejunk=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.genes.pep.faa
```

## *Pocillopora acuta* KBHIv2

From Kāne‘ohe Bay, Oʻahu, Hawaii

Version: 2.0

Genome files from JunMo (2021/12/06)

```bash
PREFIX="Pocillopora_acuta_KBHIv2"
```

**Genome**

```bash
cat Pocillopora_acuta_CANU_HM2x3_ref_assembly.fasta \
  | sed -e "s/>/>${PREFIX}___/" \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.assembly.fasta
```

**CDS**

```bash
cat Pocillopora_acuta_CDS.fasta \
  | sed -e "s/Pocillopora_acuta../${PREFIX}___/" \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.genes.cds.fna
```

**PEP**

```bash
cat Pocillopora_acuta_proteins.fasta \
  | sed -e "s/Pocillopora_acuta../${PREFIX}___/" \
  | reformat.sh fastawrap=0 iupacToN=f ignorejunk=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.genes.pep.faa
```

**GFF**

```bash
cat Pocillopora_acuta_genemodels.gff \
  | sed -e "s/^/${PREFIX}___/" \
  | sed -e "s/ID=/ID=${PREFIX}___/" \
  | sed -e "s/Parent=/Parent=${PREFIX}___/" \
  | awk -F'\t' '$3=="gene" || $3=="CDS"' \
  | gffread \
  | awk '$1!~"^#"' \
  | awk 'BEGIN{OFS=FS="\t"} { if($3=="gene"){$3="transcript"; print} else if ($3=="CDS") {print; $3="exon"; print} else {print} }' \
  | sed -e 's/;geneID=.*//' \
  > ${PREFIX}.genes.gff3
```

**Reads**

Download data from BioProject PRJNA761443. Data is short read Illumina data from a NovaSeq 6000 platform. 

```bash
SRR="SRR16077714"
fasterq-dump --threads 6 --split-3 --skip-technical --progress --outdir . ${SRR}
reformat.sh deleteinput=t verifypaired=t trimreaddescription=t addslash=t spaceslash=f \
  in=${SRR}_1.fastq in2=${SRR}_2.fastq out=${SRR}_R1.fastq.gz out2=${SRR}_R2.fastq.gz
statswrapper.sh ${SRR}_R1.fastq.gz ${SRR}_R2.fastq.gz 1>read_stats.txt
```

## *Pocillopora acuta* LBIDv1

From Lombok, Indonesia

Version: 1.0

Paper: https://doi.org/10.1101/698688

Website: http://ihpe.univ-perp.fr/acces-aux-donnees/

```bash
wget http://ihpe.univ-perp.fr/telechargement/Data_to_downoload.rar
# Unzip on mac using unarchiver
```

```bash
PREFIX="Pocillopora_acuta_LBIDv1"
```

**Genome**

```bash
cat Data_to_downoload/Pocillopora_acuta_genome_v1.fasta \
  | sed -e "s/>/>${PREFIX}___/" -e 's/_cov.*//' \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.assembly.fasta
samtools faidx ${PREFIX}.assembly.fasta
```

**Parse "Experimental" `gff`** - Retaining only the longest transcript per locus.

Reformat `gff` -**NOTE**: Set genes without strand information (looks to be only single exon genes) as being on the forward/plus strand. TransDecoder will correct the strand if needed once it finds a good candidate ORF.

```bash
cat Data_to_downoload/Structural_annotation_experimental.gff \
  | awk -F'\t' '{print "Pocillopora_acuta_LBIDv1___"$0}' \
  | gffread - -T \
  | sed -e 's/transcript_id "/transcript_id "Pocillopora_acuta_LBIDv1___/' \
        -e 's/gene_id "/gene_id "Pocillopora_acuta_LBIDv1___/' \
        -e 's/_cov[^\t]*//' \
  | awk 'BEGIN{OFS=FS="\t"} { if($7=="."){$7="+"};print }' \
  > ${PREFIX}.genes_experimental.gtf

## GTF to "alignment" gff3
~/programs/TransDecoder-v5.5.0/util/gtf_to_alignment_gff3.pl ${PREFIX}.genes_experimental.gtf \
  > ${PREFIX}.genes_experimental.alignment.gff3
```

Get transcripts

```bash
gffread ${PREFIX}.genes_experimental.gtf \
	-g ${PREFIX}.assembly.fasta \
	-w ${PREFIX}.genes_experimental.transcripts.fna
samtools faidx ${PREFIX}.genes_experimental.transcripts.fna
```

Get transcript <-> loco <-> transcript_length info in same file

```bash
awk '$3=="transcript" {print $10"\t"$12}' ${PREFIX}.genes_experimental.gtf \
  | sed -e 's/"//g' -e 's/;//g' \
  | ./add_value_to_table_SQLite3.py -a <(cut -f1,2 ${PREFIX}.genes_experimental.transcripts.fna.fai) \
  > ${PREFIX}.genes_experimental.locus
```

Get just the longest transcript at each locus

```bash
sort -k2,2 -k3,3nr ${PREFIX}.genes_experimental.locus \
  | awk -F'\t' '!seen[$2]++' \
  > ${PREFIX}.genes_experimental.locus.longest
```

Extract `gtf` features for the longest transcripts

```bash
awk '{print $0"\t"$10}' ${PREFIX}.genes_experimental.gtf \
  | awk -F'\t' 'BEGIN{OFS=FS="\t"} {gsub("\"","",$10); gsub(";","",$10); print}' \
  | ./grepf_column.py -c 10 -f <(cut -f1 ${PREFIX}.genes_experimental.locus.longest) \
  | cut -f1-9 \
  > ${PREFIX}.genes_experimental.longest.gtf
```

Filter transcripts.

```bash
./grepf_fasta.py \
  -i ${PREFIX}.genes_experimental.transcripts.fna \
  -o ${PREFIX}.genes_experimental.longest.transcripts.fna \
  -f <(cut -f1 ${PREFIX}.genes_experimental.locus.longest)
```

Execute the `run_Transdecoder.sh` script to run TransDecoder on the longest transcripts.

Allow for multiple ORFs per transcripts as some of the RNA-seq based gene models might be chimeric/fused. 

**GFF**

```bash
cat ${PREFIX}.genes_experimental.longest.transcripts.fna.transdecoder.genome.gff3 \
  | awk -F'\t' '$3=="mRNA" || $3=="CDS"' \
  | gffread \
  | awk '$1!~"^#"' \
  | awk 'BEGIN{OFS=FS="\t"} { if($3=="mRNA"){$3="transcript"; print} else if ($3=="CDS") {print; $3="exon"; $8="."; print} else {print} }' \
  | sed -e 's/;geneID=.*//' \
  > ${PREFIX}.genes.gff3
```

Extract the IDs of the genes with models that could be propagated back to the genome (i.e., ORFs on the correct strand of multi-exon genes).

```bash
awk -F'\t' '$3=="mRNA"' ${PREFIX}.genes_experimental.longest.transcripts.fna.transdecoder.genome.gff3 \
  | sed -e 's/.*ID=\([^;]*\).*/\1/' \
  > ${PREFIX}.genes.names.txt
```

**CDS**

```bash
cat ${PREFIX}.genes_experimental.longest.transcripts.fna.transdecoder.cds \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  | ./grepf_fasta.py -f ${PREFIX}.genes.names.txt \
  > ${PREFIX}.genes.cds.fna
```

**PEP**

```bash
cat ${PREFIX}.genes_experimental.longest.transcripts.fna.transdecoder.pep \
  | reformat.sh fastawrap=0 iupacToN=f ignorejunk=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  | ./grepf_fasta.py -f ${PREFIX}.genes.names.txt \
  | sed -e 's/*$//' \
  > ${PREFIX}.genes.pep.faa
```

Correct gene models on wrong strand.

It appears that TransDecoder v5.5.0, when allowed to predict multiple genes per transcript, can produce gene models with the wrong strand information when the multiple gene models are predicted on opposite strands on a single exon gene. Without introns to hard enforce a strand (i.e., to make TransDecoder discard genes on the opposite strand) the program gets confused and makes all ORFs on one of the strands. Leading to gene models that when translated are in the wrong direction and are full if stop codons. 

```bash
gffread -S \
  -g ${PREFIX}.assembly.fasta \
  -y /dev/stdout \
  ${PREFIX}.genes.gff3 \
  | reformat.sh fastawrap=0 iupacToN=f ignorejunk=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.genes.pep.faa.tmp


cat ${PREFIX}.genes.gff3 \
  | sed -e 's/.*=\(.*\)/\0\t\1/' \
  | ./add_value_to_table_SQLite3.py -c 10 -a <(grep -B 1 '*' ${PREFIX}.genes.pep.faa.tmp | grep '>' \
  | sed -e 's/>//' -e 's/$/\tIntStop/') \
  | awk 'BEGIN{OFS=FS="\t"} { if($11=="IntStop"){ if($7=="+"){$7="-"}else{$7="+"};print }else{print} }' \
  | cut -f1-9 \
  > Pocillopora_acuta_LBIDv1.genes.tmp.gff3 


gffread -S \
  -g ${PREFIX}.assembly.fasta \
  -y /dev/stdout \
  ${PREFIX}.genes.tmp.gff3 \
  | reformat.sh fastawrap=0 iupacToN=f ignorejunk=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.genes.pep.faa.tmp2
```

Also had to edit 3 genes manually that didn't have internal stop codons but were still placed on the incorrect strand by TransDecoder.

```bash
Pocillopora_acuta_LBIDv1___TCONS_00000959.p3
Pocillopora_acuta_LBIDv1___TCONS_00010495.p2
Pocillopora_acuta_LBIDv1___TCONS_00047697.p2
```

Check the difference between proteins extracted using `gffread` and the original TransDecoder proteins.

```bash
diff <(seqkit fx2tab Pocillopora_acuta_LBIDv1.genes.pep.faa | sort -k1,1) <(seqkit fx2tab Pocillopora_acuta_LBIDv1.genes.pep.faa.tmp2 | sort -k1,1) | less
```

Overwrite the original `gff` file now that we have corrected the strand for some of the problem gene models. 

```bash
mv Pocillopora_acuta_LBIDv1.genes.tmp.gff3 Pocillopora_acuta_LBIDv1.genes.gff3
rm Pocillopora_acuta_LBIDv1.genes.pep.faa.tmp Pocillopora_acuta_LBIDv1.genes.pep.faa.tmp2
```

**Reads**

Download data from BioProject PRJNA342887. Data is short read Illumina data from a HiSeq 2000 platform. 

```bash
SRR="SRR4254617"
fasterq-dump --threads 6 --split-3 --skip-technical --progress --outdir . ${SRR}
reformat.sh deleteinput=t verifypaired=t trimreaddescription=t addslash=t spaceslash=f \
  in=${SRR}_1.fastq in2=${SRR}_2.fastq out=${SRR}_R1.fastq.gz out2=${SRR}_R2.fastq.gz
statswrapper.sh ${SRR}_R1.fastq.gz ${SRR}_R2.fastq.gz 1>read_stats.txt
```

## *Pocillopora damicornis* SIPAv1

From: Saboga Is., Panama

Version: 1.0

Paper: https://doi.org/10.1038/s41598-018-34459-8

Website: http://pdam.reefgenomics.org/

```bash
wget http://pdam.reefgenomics.org/download/pdam_scaffolds.fasta.gz
wget http://pdam.reefgenomics.org/download/pdam_transcripts.fasta.gz
wget http://pdam.reefgenomics.org/download/pdam_proteins.fasta.gz
wget http://pdam.reefgenomics.org/download/pdam_annotation.gff3.gz
```

```bash
PREFIX="Pocillopora_damicornis_SIPAv1"
```

**Genome**

```bash
zcat pdam_scaffolds.fasta.gz \
  | sed -e "s/>/>${PREFIX}___/" \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.assembly.fasta
```

**CDS**

```bash
zcat pdam_transcripts.fasta.gz \
  | sed -e "s/>pdam_/>${PREFIX}___/" \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.genes.cds.fna
```

**PEP**

```bash
zcat pdam_proteins.fasta.gz \
  | sed -e "s/>pdam_/>${PREFIX}___/" \
  | reformat.sh fastawrap=0 iupacToN=f ignorejunk=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.genes.pep.faa
```

**GFF**

```bash
zcat pdam_annotation.gff3.gz \
  | sed -e "s/^/${PREFIX}___/" \
  | sed -e "s/ID=pdam_/ID=${PREFIX}___/" \
  | sed -e "s/Parent=pdam_/Parent=${PREFIX}___/" \
  | awk -F'\t' '$3=="mRNA" || $3=="CDS" || $3=="exon"' \
  | gffread \
  | awk '$1!~"^#"' \
  | awk 'BEGIN{OFS=FS="\t"} { if($3=="mRNA") {$3="transcript"; print} else {print} }' \
  > ${PREFIX}.genes.gff3
```

## *Pocillopora meandrina* KBHIv1

From Kāne‘ohe Bay, Oʻahu, Hawaii

Version: 1.0

Genome files from JunMo (2021/12/06)

```bash
PREFIX="Pocillopora_meandrina_KBHIv1"
```

**Genome**

```bash
cat Pocillopora_meandrina_CANU_HM2_ref_assembly.fasta \
  | sed -e "s/>/>${PREFIX}___/" \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.assembly.fasta
```

**CDS**

```bash
cat Pocillopora_meandrina_CDS.fasta \
  | sed -e "s/Pocillopora_meandrina../${PREFIX}___/" \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.genes.cds.fna
```

**PEP**

```bash
cat Pocillopora_meandrina_proteins.fasta \
  | sed -e "s/Pocillopora_meandrina../${PREFIX}___/" \
  | reformat.sh fastawrap=0 iupacToN=f ignorejunk=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.genes.pep.faa
```

**GFF**

```bash
cat Pocillopora_meandrina_genemodels.gff \
  | sed -e "s/^/${PREFIX}___/" \
  | sed -e "s/ID=/ID=${PREFIX}___/" \
  | sed -e "s/Parent=/Parent=${PREFIX}___/" \
  | awk -F'\t' '$3=="gene" || $3=="CDS"' \
  | gffread \
  | awk '$1!~"^#"' \
  | awk 'BEGIN{OFS=FS="\t"} { if($3=="gene"){$3="transcript"; print} else if ($3=="CDS") {print; $3="exon"; print} else {print} }' \
  | sed -e 's/;geneID=.*//' \
  > ${PREFIX}.genes.gff3
```

**Reads**

Download data from BioProject PRJNA761443. Data is short read Illumina data from a NovaSeq 6000 platform. 

```bash
SRR="SRR16077712"
fasterq-dump --threads 6 --split-3 --skip-technical --progress --outdir . ${SRR}
reformat.sh deleteinput=t verifypaired=t trimreaddescription=t addslash=t spaceslash=f \
  in=${SRR}_1.fastq in2=${SRR}_2.fastq out=${SRR}_R1.fastq.gz out2=${SRR}_R2.fastq.gz
statswrapper.sh ${SRR}_R1.fastq.gz ${SRR}_R2.fastq.gz 1>read_stats.txt
```

## *Pocillopora verrucosa* RSSAv1

From Red Sea Al Fahal, Saudi Arabia

Version: 1.0

Paper: https://doi.org/10.1093/gbe/evaa184

Website: http://pver.reefgenomics.org/

```bash
wget http://pver.reefgenomics.org/download/Pver_genome_assembly_v1.0.fasta.gz
wget http://pver.reefgenomics.org/download/Pver_genes_names_v1.0.fna.gz
wget http://pver.reefgenomics.org/download/Pver_proteins_names_v1.0.faa.gz
wget http://pver.reefgenomics.org/download/Pver_genome_assembly_v1.0.gff3.gz
```

```bash
PREFIX="Pocillopora_verrucosa_RSSAv1"
```

**Genome**

```bash
zcat Pver_genome_assembly_v1.0.fasta.gz \
  | sed -e "s/>Pver_/>${PREFIX}___/" -e 's/_size.*//' \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.assembly.fasta
```

**CDS**

```bash
zcat Pver_genes_names_v1.0.fna.gz \
  | sed -e "s/>Pver_/>${PREFIX}___/" \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.genes.cds.fna
```

**PEP**

```bash
zcat Pver_proteins_names_v1.0.faa.gz \
  | sed -e "s/>Pver_/>${PREFIX}___/" \
  | reformat.sh fastawrap=0 iupacToN=f ignorejunk=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  | sed -e 's/*$//' \
  > ${PREFIX}.genes.pep.faa
```

**GFF**

```bash
zcat Pver_genome_assembly_v1.0.gff3.gz \
  | cut -f1-9 \
  | sed -e "s/^Pver_/${PREFIX}___/" -e 's/_size[^\t]*//' \
  | sed -e "s/ID=Pver_/ID=${PREFIX}___/" \
  | sed -e "s/Parent=Pver_/Parent=${PREFIX}___/" \
  | awk -F'\t' '$3=="CDS" || $3=="exon"' \
  | gffread \
  | awk '$1!~"^#"' \
  > ${PREFIX}.genes.gff3
```

**Reads**

Download data from BioProject PRJNA551401 BioSample SAMN12147516. Two Runs are associated with the BioProject (one WGS and one RNA-seq), download the one WGS Run (SRR11880677) that is Illumina paired-end data (generated on a HiSeq 2500 platform). 

```bash
SRR="SRR11880677"
fasterq-dump --threads 6 --split-3 --skip-technical --progress --outdir . ${SRR}
reformat.sh deleteinput=t verifypaired=t trimreaddescription=t addslash=t spaceslash=f \
  in=${SRR}_1.fastq in2=${SRR}_2.fastq out=${SRR}_R1.fastq.gz out2=${SRR}_R2.fastq.gz
statswrapper.sh ${SRR}_R1.fastq.gz ${SRR}_R2.fastq.gz 2>read_stats.txt
```

## *Porites australiensis* SLJPv1

From Sesoko Island, Okinawa, Japan

Version: 1.0

Paper: https://dx.doi.org/10.1093%2Fgbe%2Fevab270

Website: https://drive.google.com/drive/folders/1xZCIAHLmuCcdIQlvDHbN7GGJzxp6VIs0

```bash
# Download mRNA+gff+assembly manually from Google Drive
```

```bash
PREFIX="Porites_australiensis_SLJPv1"
```

**Genome**

```bash
zcat paus_genome-assembly.fasta.gz \
  | sed -e "s/>/>${PREFIX}___/" \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.assembly.fasta
```

**GFF**

```bash
zcat paus_mRNA.gff.gz \
  | sed -e "s/^/${PREFIX}___/" -e "s/paus_/${PREFIX}___/g" \
  | awk -F'\t' '$3=="transcript" || $3=="CDS" || $3=="exon"' \
  | gffread \
  | awk '$1!~"^#"' \
  > ${PREFIX}.genes.all.gff3

cat ${PREFIX}.genes.all.gff3 \
  | awk -F'\t' '$3=="CDS" {
      split($9,a,"="); 
      CDS[a[2]]=CDS[a[2]]+(($5-$4)+1)
    } END {
      for (i in CDS){
        split(i,ii,"\\.t")
        print ii[1]"\t"i"\t"CDS[i]
    } }' \
  > ${PREFIX}.genes.all.gff3.CDS_lengths.tsv

cat ${PREFIX}.genes.all.gff3.CDS_lengths.tsv \
  | sort -k1,1 -k3,3nr \
  | awk '!seen[$1]++{print $2}' \
  > ${PREFIX}.genes.all.gff3.longest.txt

cat ${PREFIX}.genes.all.gff3 \
  | awk -F'\t' '{split($9,a,"="); print $0"\t"a[2]}' \
  | ~/scripts/grepf_column.py -c 10 -f ${PREFIX}.genes.all.gff3.longest.txt \
  | cut -f 1-9 \
  > ${PREFIX}.genes.gff3

gffread -S \
  -g ${PREFIX}.assembly.fasta \
  -x ${PREFIX}.genes.cds.fna.tmp \
  -y ${PREFIX}.genes.pep.faa.tmp \
  ${PREFIX}.genes.gff3
```

**CDS**

```bash
cat ${PREFIX}.genes.cds.fna.tmp \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.genes.cds.fna
rm ${PREFIX}.genes.cds.fna.tmp
```

**PEP**

```bash
cat ${PREFIX}.genes.pep.faa.tmp \
  | reformat.sh fastawrap=0 iupacToN=f ignorejunk=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.genes.pep.faa
rm ${PREFIX}.genes.pep.faa.tmp
```

**Reads**

Download data from BioProject PRJDB4187. Data is short read Illumina data from a MiSeq platform. 

```bash
SRR="DRR045165"
fasterq-dump --threads 6 --split-3 --skip-technical --progress --outdir . ${SRR}
reformat.sh deleteinput=t verifypaired=t trimreaddescription=t addslash=t spaceslash=f \
  in=${SRR}_1.fastq in2=${SRR}_2.fastq out=${SRR}_R1.fastq.gz out2=${SRR}_R2.fastq.gz
statswrapper.sh ${SRR}_R1.fastq.gz ${SRR}_R2.fastq.gz 1>read_stats.txt
```

## *Porites astreoides* BBBMv1

From Bailey’s Bay Reef Flats (32°2227N, 64°4437W) in Bermuda

Version: 1.0

Paper: In Press

Website: https://osf.io/7s5z4/

```bash
wget https://sra-download.ncbi.nlm.nih.gov/traces/sra25/SRZ/019144/SRR19144705/filtered_p_ctg.fasta
# Gff + pep + cds from https://osf.io/7s5z4/
```

```bash
PREFIX="Porites_astreoides_BBBMv1"
```

**Genome**

```bash
cat filtered_p_ctg.fasta \
  | sed -e "s/>/>${PREFIX}___/" \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.assembly.fasta
```

**GFF**

`gffread` removes/excludes one gene because it is only 2aa long. Add it back in manually to the `${PREFIX}.genes.gff3` file.

> Porites_astreoides_BBBMv1___000176F   maker  transcript   474987 474995 .    +    .    ID=Porites_astreoides_BBBMv1___18076-RA
>
> Porites_astreoides_BBBMv1___000176F   maker  CDS   474987 474995 .    +    0    Parent=Porites_astreoides_BBBMv1___18076-RA;

```bash
cat Pastreoides_all_v1.gff \
  | sed -e "s/^/${PREFIX}___/" \
  | sed -e "s/ID=Pastreoides/ID=${PREFIX}___/" \
  | sed -e "s/Parent=Pastreoides/Parent=${PREFIX}___/" \
  | awk -F'\t' '$3=="CDS"' \
  | gffread \
  | awk '$1!~"^#"' \
  | awk 'BEGIN{OFS=FS="\t"} { if($3=="gene"){$3="transcript"; print} else if ($3=="CDS") {print; $3="exon"; print} else {print} }' \
  | sed -e 's/;geneID=.*//' \
  > ${PREFIX}.genes.gff3

gffread -S \
  -g ${PREFIX}.assembly.fasta \
  -x ${PREFIX}.genes.cds.fna.tmp \
  -y ${PREFIX}.genes.pep.faa.tmp \
  ${PREFIX}.genes.gff3
```

**CDS**

```bash
cat ${PREFIX}.genes.cds.fna.tmp \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.genes.cds.fna
rm ${PREFIX}.genes.cds.fna.tmp
```

**PEP**

```bash
cat ${PREFIX}.genes.pep.faa.tmp \
  | reformat.sh fastawrap=0 iupacToN=f ignorejunk=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.genes.pep.faa
rm ${PREFIX}.genes.pep.faa.tmp
```

The `Pastreoides_transcripts_v1.fasta` file is the axons not CDS (so contains UTRs), which is why we used `gffread` to regenerate these files. The below commands confirm that the original and new protein files are identical and all is consistent between old and new genes.

```bash
grep -v '>' Porites_astreoides_BBBMv1.genes.pep.faa | sort | md5sum
512b5c3eed317e9e0de2d15ca19fbeb9  -
cat Pastreoides_proteins_v1.fasta | reformat.sh fastawrap=0 iupacToN=f ignorejunk=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa | grep -v '>' | sort | md5sum
512b5c3eed317e9e0de2d15ca19fbeb9  -
```



## *Porites compressa* KBHIv1

From Kāne‘ohe Bay, Oʻahu, Hawaii

Version: 1.0

Genome files from JunMo (2021/12/06)

```bash
PREFIX="Porites_compressa_KBHIv1"
```

**Genome**

```bash
cat Porites_compressa_CANU_HM2x2_ref_assembly.fasta \
  | sed -e "s/>/>${PREFIX}___/" \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.assembly.fasta
```

**CDS**

```bash
cat Porites_compressa_CDS.fasta \
  | sed -e "s/Porites_compressa../${PREFIX}___/" \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.genes.cds.fna
```

**PEP**

```bash
cat Porites_compressa_proteins.fasta \
  | sed -e "s/Porites_compressa../${PREFIX}___/" \
  | reformat.sh fastawrap=0 iupacToN=f ignorejunk=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.genes.pep.faa
```

**GFF**

```bash
cat Porites_compressa_genemodels.gff \
  | sed -e "s/^/${PREFIX}___/" \
  | sed -e "s/ID=/ID=${PREFIX}___/" \
  | sed -e "s/Parent=/Parent=${PREFIX}___/" \
  | awk -F'\t' '$3=="gene" || $3=="CDS"' \
  | gffread \
  | awk '$1!~"^#"' \
  | awk 'BEGIN{OFS=FS="\t"} { if($3=="gene"){$3="transcript"; print} else if ($3=="CDS") {print; $3="exon"; print} else {print} }' \
  | sed -e 's/;geneID=.*//' \
  > ${PREFIX}.genes.gff3
```

**Reads**

Download data from BioProject PRJNA663761. Data is short read Illumina data from a HiSeq 2000 platform. 

```bash
SRR="SRR12695158"
fasterq-dump --threads 6 --split-3 --skip-technical --progress --outdir . ${SRR}
reformat.sh deleteinput=t verifypaired=t trimreaddescription=t addslash=t spaceslash=f \
  in=${SRR}_1.fastq in2=${SRR}_2.fastq out=${SRR}_R1.fastq.gz out2=${SRR}_R2.fastq.gz
statswrapper.sh ${SRR}_R1.fastq.gz ${SRR}_R2.fastq.gz 1>read_stats.txt
```

## *Porites lutea* OIAUv1.1

From Orpheus Island, Great Barrier Reef, Australia

Version: 1.1

Paper: https://doi.org/10.1038/s41564-019-0532-4

Website: http://plut.reefgenomics.org

```bash
wget http://plut.reefgenomics.org/download/plut_final_2.1.fasta.gz
wget http://plut.reefgenomics.org/download/plut2v1.1.transcripts.fasta.gz
wget http://plut.reefgenomics.org/download/plut2v1.1.proteins.fasta.gz
wget http://plut.reefgenomics.org/download/plut2v1.1.genes.gff3.gz
```

```bash
PREFIX="Porites_lutea_OIAUv1.1"
```

**Genome**

```bash
zcat plut_final_2.1.fasta.gz \
  | sed -e "s/>/>${PREFIX}___/" \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.assembly.fasta
```

**CDS**

```bash
zcat plut2v1.1.transcripts.fasta.gz \
  | sed -e "s/>plut2./>${PREFIX}___/" \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.genes.cds.fna
```

**PEP**

```bash
zcat plut2v1.1.proteins.fasta.gz \
  | sed -e "s/>plut2./>${PREFIX}___/" \
  | reformat.sh fastawrap=0 iupacToN=f ignorejunk=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.genes.pep.faa
```

**GFF**

```bash
zcat plut2v1.1.genes.gff3.gz \
  | awk '$1!~"^#"' \
  | sed -e "s/^/${PREFIX}___/" \
  | sed -e "s/ID=plut2./ID=${PREFIX}___/" \
  | sed -e "s/Parent=plut2./Parent=${PREFIX}___/" \
  | sed -e "s/Name=plut2./Name=${PREFIX}___/" \
  | awk -F'\t' '$3=="mRNA" || $3=="CDS" || $3=="exon"' \
  | gffread \
  | awk '$1!~"^#"' \
  | awk 'BEGIN{OFS=FS="\t"} { if($3=="mRNA") {$3="transcript"; print} else {print} }' \
  > ${PREFIX}.genes.gff3
```

**Reads**

Download data from BioProject PRJEB6884. Data is short read Illumina data from a HiSeq 2500 platform; 7 Runs from 6 Experiments. 

```bash
for SRR in ERR571463 ERR571462 ERR571461 ERR571460 ERR571459 ERR571457 ERR571458;
do
  fasterq-dump --threads 6 --split-3 --skip-technical --progress --outdir . ${SRR}
  reformat.sh deleteinput=t verifypaired=t trimreaddescription=t addslash=t spaceslash=f \
    in=${SRR}_1.fastq in2=${SRR}_2.fastq out=${SRR}_R1.fastq.gz out2=${SRR}_R2.fastq.gz
  statswrapper.sh ${SRR}_R1.fastq.gz ${SRR}_R2.fastq.gz 1>read_stats.txt
done
```

## *Porites rus* NAIDv1

From: imported from Indonesia in 2007

Version: 1.0

Paper: https://doi.org/10.1093/gigascience/giy075

Website: Genome: (NCBI: GCA_900290455.1)

```bash
# Genome assembly from NCBI
wget https://ftp.cngb.org/pub/gigadb/pub/10.5524/100001_101000/100462/Prus_cds.fas
wget https://ftp.cngb.org/pub/gigadb/pub/10.5524/100001_101000/100462/Prus_cds_aa.fas
wget https://ftp.cngb.org/pub/gigadb/pub/10.5524/100001_101000/100462/Prus_gene_prediction.gff
```

```bash
PREFIX="Porites_rus_NAIDv1"
```

**Genome**

```bash
gunzip -c ncbi-genomes-2022-04-15/GCA_900290455.1_Prus_genomic.fna.gz \
  | awk '{ if($1~"^>") {gsub(",","",$9); print ">"$9} else {print} }' \
  | sed -e "s/>/>${PREFIX}___/" \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.assembly.fasta
```

**GFF**

Only take "CDS" features and let `gffread` add "transcript" features. Have to do this as "transcript" features arent formated correctly and `gffread` ignores them anyway. 

```bash
cat Prus_gene_prediction.gff \
  | awk -F'\t' '$1!~"^#"' \
  | awk -F'\t' '$3=="CDS"' \
  | gffread \
  | awk '$1!~"^#"' \
  | sed -e "s/^/${PREFIX}___/" -e "s/=/=${PREFIX}___/g" \
  | awk 'BEGIN{OFS=FS="\t"} { if($3!="CDS"){print} else {print; $3="exon"; print} }' \
  | sed -e 's/;geneID=.*//' \
  > ${PREFIX}.genes.gff3
```

**CDS**

Remove the CDS what are called "Prus_finalgXXX.t1", these look to be updated?? Versions of existing gene models however they are not in the `gff` so we will ignore them for now. 

```bash
cat Prus_cds.fas \
  | sed -e "s/>/>${PREFIX}___/" \
  | seqkit fx2tab | grep -v 'Prus_finalg' | seqkit tab2fx \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.genes.cds.fna
```

**PEP**

Remove the PEP what are called "Prus_finalgXXX.t1", these look to be updated?? Versions of existing gene models however they are not in the `gff` so we will ignore them for now. 

```bash
cat Prus_cds_aa.fas \
  | sed -e "s/>/>${PREFIX}___/" -e 's/*$//' \
  | seqkit fx2tab | grep -v 'Prus_finalg' | seqkit tab2fx \
  | reformat.sh fastawrap=0 iupacToN=f ignorejunk=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.genes.pep.faa
```

**Reads**

Download data from BioProject PRJEB25136. Data is short read Illumina data but generated from a Illumina HiSeq 1500 machine. 

```bash
SRR="ERR2322283"
fasterq-dump --threads 6 --split-3 --skip-technical --progress --outdir . ${SRR}
reformat.sh deleteinput=t verifypaired=t trimreaddescription=t addslash=t spaceslash=f \
  in=${SRR}_1.fastq in2=${SRR}_2.fastq out=${SRR}_R1.fastq.gz out2=${SRR}_R2.fastq.gz
statswrapper.sh ${SRR}_R1.fastq.gz ${SRR}_R2.fastq.gz 2>read_stats.txt
```

## *Stylophora pistillata* GAJOv1

From: in front of the Marine Science Station, Gulf of Aqaba, Jordan

**NOTE: There are 39453 genes reported in the paper and in the `gff` file but there are 39511 in the supplied protein and CDS `fasta` files. There are 58 proteins (named "Prus_finalg...") which are in the supplied `fasta` files but not the `gff` file.**

Version: 1.0

Paper: https://doi.org/10.1038/s41598-017-17484-x

Website: http://spis.reefgenomics.org

```bash
wget http://spis.reefgenomics.org/download/Spis.genome.scaffold.final.fa.gz
wget http://spis.reefgenomics.org/download/Spis.genome.annotation.CDS.longest.fa.gz
wget http://spis.reefgenomics.org/download/Spis.genome.annotation.pep.longest.fa.gz
wget http://spis.reefgenomics.org/download/Spis.genome.annotation.gff3.gz
```

```bash
PREFIX="Stylophora_pistillata_GAJOv1"
```

**Genome**

```bash
zcat Spis.genome.scaffold.final.fa.gz \
  | sed -e "s/>Spis./>${PREFIX}___/" \
  | sed -e 's/|.*//' \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.assembly.fasta
```

**GFF**

```bash
zcat Spis.genome.annotation.gff3.gz \
  | awk '$1!~"^#"' \
  | sed -e "s/^Spis./${PREFIX}___/" -e 's/|[^\t]*\t/\t/' \
  | sed -e "s/ID=Spis/ID=${PREFIX}___/" \
  | sed -e "s/Parent=Spis/Parent=${PREFIX}___/" \
  | sed -e "s/Name=Spis/Name=${PREFIX}___/" \
  | awk -F'\t' '$3=="mRNA" || $3=="CDS" || $3=="exon"' \
  | gffread \
  | awk '$1!~"^#"' \
  | awk 'BEGIN{OFS=FS="\t"} { if($3=="mRNA") {$3="transcript"; print} else {print} }' \
  > ${PREFIX}.genes.all.gff3

cat ${PREFIX}.genes.all.gff3 \
  | awk -F'\t' '$3=="CDS" {
      split($9,a,"="); 
      CDS[a[2]]=CDS[a[2]]+(($5-$4)+1)
    } END {
      for (i in CDS){
        split(i,ii,"\\.t")
        print ii[1]"\t"i"\t"CDS[i]
    } }' \
  > ${PREFIX}.genes.all.gff3.CDS_lengths.tsv

cat ${PREFIX}.genes.all.gff3.CDS_lengths.tsv \
  | sort -k1,1 -k3,3nr \
  | awk '!seen[$1]++{print $2}' \
  > ${PREFIX}.genes.all.gff3.longest.txt

cat ${PREFIX}.genes.all.gff3 \
  | awk -F'\t' '{split($9,a,"="); print $0"\t"a[2]}' \
  | ~/scripts/grepf_column.py -c 10 -f ${PREFIX}.genes.all.gff3.longest.txt \
  | cut -f 1-9 \
  > ${PREFIX}.genes.gff3
```

**CDS**

```bash
zcat Spis.genome.annotation.CDS.longest.fa.gz \
  | sed -e "s/>Spis/>${PREFIX}___/" \
  | awk '{print $1}' \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.genes.cds.fna
```

**PEP**

```bash
zcat Spis.genome.annotation.pep.longest.fa.gz \
  | sed -e "s/>Spis/>${PREFIX}___/" \
  | awk '{print $1}' \
  | reformat.sh fastawrap=0 iupacToN=f ignorejunk=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  | sed -e 's/*$//' \
  > ${PREFIX}.genes.pep.faa
```

---

---

---

---

---

# Coral Transcripts

## *Astreopora* sp.

Transcripts published in: https://doi.org/10.7554/eLife.13288

Original paper: 

From: Unpublished

Version: 1.0

Website: http://comparative.reefgenomics.org

```bash
wget http://comparative.reefgenomics.org/faa/Coral/Astreopora_sp_peptides_100.final.clstr.faa
wget http://comparative.reefgenomics.org/fna/Coral/Astreopora_sp_cds_100.final.clstr.fna
wget http://comparative.reefgenomics.org/tx_annots/Astreopora_sp.annot.tsv
```

```bash
PREFIX="Astreopora_sp_REEFv1"
```

**CDS**

```bash
cat Astreopora_sp_cds_100.final.clstr.fna \
  | sed -e "s/>Astreopora_sp_/>${PREFIX}___/" \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.transcripts.cds.fna
```

**PEP**

```bash
cat Astreopora_sp_peptides_100.final.clstr.faa \
  | sed -e "s/>Astreopora_sp_/>${PREFIX}___/" -e 's/*$//' \
  | reformat.sh fastawrap=0 iupacToN=f ignorejunk=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.transcripts.pep.faa
```

## *Favia* sp.

Transcripts published in: https://doi.org/10.7554/eLife.13288

Original paper: 

From: Gulf of Eilat, in the northern Red Sea

Version: 1.0

Website: 

```bash
wget http://comparative.reefgenomics.org/faa/Coral/Favia_sp_peptides_100.final.clstr.faa
wget http://comparative.reefgenomics.org/fna/Coral/Favia_sp_cds_100.final.clstr.fna
wget http://comparative.reefgenomics.org/tx_annots/Favia_sp.annot.tsv
```

```bash
PREFIX="Favia_sp_REEFv1"
```

**CDS**

```bash
cat Favia_sp_cds_100.final.clstr.fna \
  | sed -e "s/>Favia_sp_/>${PREFIX}___/" \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.transcripts.cds.fna
```

**PEP**

```bash
cat Favia_sp_peptides_100.final.clstr.faa \
  | sed -e "s/>Favia_sp_/>${PREFIX}___/" -e 's/*$//' \
  | reformat.sh fastawrap=0 iupacToN=f ignorejunk=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.transcripts.pep.faa
```

## *Fungia scutaria*

Transcripts published in: https://doi.org/10.7554/eLife.13288

Original paper: 

From: Coconut Island, Hawaii

Version: 1.0

Website: 

```bash
wget http://comparative.reefgenomics.org/faa/Coral/Fungia_scutaria_peptides_100.final.clstr.faa
wget http://comparative.reefgenomics.org/fna/Coral/Fungia_scutaria_cds_100.final.clstr.fna
wget http://comparative.reefgenomics.org/tx_annots/Fungia_scutaria.annot.tsv
```

```bash
PREFIX="Fungia_scutaria_REEFv1"
```

**CDS**

```bash
cat Fungia_scutaria_cds_100.final.clstr.fna \
  | sed -e "s/>Fungia_scutaria_/>${PREFIX}___/" \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.transcripts.cds.fna
```

**PEP**

```bash
cat Fungia_scutaria_peptides_100.final.clstr.faa \
  | sed -e "s/>Fungia_scutaria_/>${PREFIX}___/" -e 's/*$//' \
  | reformat.sh fastawrap=0 iupacToN=f ignorejunk=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.transcripts.pep.faa
```

## *Madracis auretenra*

Transcripts published in: https://doi.org/10.7554/eLife.13288

Original paper: 

From: 

Version: 1.0

Website: 

```bash
wget http://comparative.reefgenomics.org/faa/Coral/Madracis_auretenra_peptides_100.final.clstr.faa
wget http://comparative.reefgenomics.org/fna/Coral/Madracis_auretenra_cds_100.final.clstr.fna
wget http://comparative.reefgenomics.org/tx_annots/Madracis_auretenra.annot.tsv
```

```bash
PREFIX="Madracis_auretenra_REEFv1"
```

**CDS**

```bash
cat Madracis_auretenra_cds_100.final.clstr.fna \
  | sed -e "s/>Madracis_auretenra_/>${PREFIX}___/" \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.transcripts.cds.fna
```

**PEP**

```bash
cat Madracis_auretenra_peptides_100.final.clstr.faa \
  | sed -e "s/>Madracis_auretenra_/>${PREFIX}___/" -e 's/*$//' \
  | reformat.sh fastawrap=0 iupacToN=f ignorejunk=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.transcripts.pep.faa
```

## *Montastraea cavernosa*

Transcripts published in: https://doi.org/10.7554/eLife.13288

Original paper: 

From: 

Version: 1.0

Website: 

```bash
wget http://comparative.reefgenomics.org/faa/Coral/Montastraea_cavernosa_peptides_100.final.clstr.faa
wget http://comparative.reefgenomics.org/fna/Coral/Montastraea_cavernosa_cds_100.final.clstr.fna
wget http://comparative.reefgenomics.org/tx_annots/Montastraea_cavernosa.annot.tsv
```

```bash
PREFIX="Montastraea_cavernosa_REEFv1"
```

**CDS**

```bash
cat Montastraea_cavernosa_cds_100.final.clstr.fna \
  | sed -e "s/>Montastraea_cavernosa_/>${PREFIX}___/" \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.transcripts.cds.fna
```

**PEP**

```bash
cat Montastraea_cavernosa_peptides_100.final.clstr.faa \
  | sed -e "s/>Montastraea_cavernosa_/>${PREFIX}___/" -e 's/*$//' \
  | reformat.sh fastawrap=0 iupacToN=f ignorejunk=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.transcripts.pep.faa
```

## *Montastraea faveolata*

Transcripts published in: https://doi.org/10.7554/eLife.13288

Original paper: 

From: 

Version: 1.0

Website: 

```bash
wget http://comparative.reefgenomics.org/faa/Coral/Montastraea_faveolata_peptides_100.final.clstr.faa
wget http://comparative.reefgenomics.org/fna/Coral/Montastraea_faveolata_cds_100.final.clstr.fna
wget http://comparative.reefgenomics.org/tx_annots/Montastraea_faveolata.annot.tsv
```

```bash
PREFIX="Montastraea_faveolata_REEFv1"
```

**CDS**

```bash
cat Montastraea_faveolata_cds_100.final.clstr.fna \
  | sed -e "s/>Montastraea_faveolata_/>${PREFIX}___/" -e 's/-/N/g' \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.transcripts.cds.fna
```

**PEP**

```bash
cat Montastraea_faveolata_peptides_100.final.clstr.faa \
  | sed -e "s/>Montastraea_faveolata_/>${PREFIX}___/" -e 's/*$//' \
  | reformat.sh fastawrap=0 iupacToN=f ignorejunk=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.transcripts.pep.faa
```

## *Platygyra carnosus*

Transcripts published in: https://doi.org/10.7554/eLife.13288

Original paper: 

From: 

Version: 1.0

Website: 

```bash
wget http://comparative.reefgenomics.org/faa/Coral/Platygyra_carnosus_peptides_100.final.clstr.faa
wget http://comparative.reefgenomics.org/fna/Coral/Platygyra_carnosus_cds_100.final.clstr.fna
wget http://comparative.reefgenomics.org/tx_annots/Platygyra_carnosus.annot.tsv
```

```bash
PREFIX="Platygyra_carnosus_REEFv1"
```

**CDS**

```bash
cat Platygyra_carnosus_cds_100.final.clstr.fna \
  | sed -e "s/>Platygyra_carnosus_/>${PREFIX}___/" \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.transcripts.cds.fna
```

**PEP**

```bash
cat Platygyra_carnosus_peptides_100.final.clstr.faa \
  | sed -e "s/>Platygyra_carnosus_/>${PREFIX}___/" -e 's/*$//' \
  | reformat.sh fastawrap=0 iupacToN=f ignorejunk=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.transcripts.pep.faa
```

## *Porites astreoides*

Transcripts published in: https://doi.org/10.7554/eLife.13288

Original paper: 

From: 

Version: 1.0

Website: 

```bash
wget http://comparative.reefgenomics.org/faa/Coral/Porites_astreoides_peptides_100.final.clstr.faa
wget http://comparative.reefgenomics.org/fna/Coral/Porites_astreoides_cds_100.final.clstr.fna
wget http://comparative.reefgenomics.org/tx_annots/Porites_astreoides.annot.tsv
```

```bash
PREFIX="Porites_astreoides_REEFv1"
```

**CDS**

```bash
cat Porites_astreoides_cds_100.final.clstr.fna \
  | sed -e "s/>Porites_astreoides_/>${PREFIX}___/" \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.transcripts.cds.fna
```

**PEP**

```bash
cat Porites_astreoides_peptides_100.final.clstr.faa \
  | sed -e "s/>Porites_astreoides_/>${PREFIX}___/" -e 's/*$//' \
  | reformat.sh fastawrap=0 iupacToN=f ignorejunk=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.transcripts.pep.faa
```

## *Porites lobata*

Transcripts published in: https://doi.org/10.7554/eLife.13288

Original paper: 

From: 

Version: 1.0

Website: 

```bash
wget http://comparative.reefgenomics.org/faa/Coral/Porites_lobata_peptides_100.final.clstr.faa
wget http://comparative.reefgenomics.org/fna/Coral/Porites_lobata_cds_100.final.clstr.fna
wget http://comparative.reefgenomics.org/tx_annots/Porites_lobata.annot.tsv
```

```bash
PREFIX="Porites_lobata_REEFv1"
```

**CDS**

```bash
cat Porites_lobata_cds_100.final.clstr.fna \
  | sed -e "s/>Porites_lobata_/>${PREFIX}___/" \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.transcripts.cds.fna
```

**PEP**

```bash
cat Porites_lobata_peptides_100.final.clstr.faa \
  | sed -e "s/>Porites_lobata_/>${PREFIX}___/" -e 's/*$//' \
  | reformat.sh fastawrap=0 iupacToN=f ignorejunk=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.transcripts.pep.faa
```

## *Pseudodiploria strigosa*

Transcripts published in: https://doi.org/10.7554/eLife.13288

Original paper: 

From: 

Version: 1.0

Website: 

```bash
wget http://comparative.reefgenomics.org/faa/Coral/Pseudodiploria_strigosa_peptides_100.final.clstr.faa
wget http://comparative.reefgenomics.org/fna/Coral/Pseudodiploria_strigosa_cds_100.final.clstr.fna
wget http://comparative.reefgenomics.org/tx_annots/Pseudodiploria_strigosa.annot.tsv
```

```bash
PREFIX="Pseudodiploria_strigosa_REEFv1"
```

**CDS**

```bash
cat Pseudodiploria_strigosa_cds_100.final.clstr.fna \
  | sed -e "s/>Pseudodiploria_strigosa_/>${PREFIX}___/" \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.transcripts.cds.fna
```

**PEP**

```bash
cat Pseudodiploria_strigosa_peptides_100.final.clstr.faa \
  | sed -e "s/>Pseudodiploria_strigosa_/>${PREFIX}___/" -e 's/*$//' \
  | reformat.sh fastawrap=0 iupacToN=f ignorejunk=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.transcripts.pep.faa
```

## *Seriatopora hystrix*

Transcripts published in: https://doi.org/10.7554/eLife.13288

Original paper: 

From: 

Version: 1.0

Website: 

```bash
wget http://comparative.reefgenomics.org/faa/Coral/Seriatopora_hystrix_peptides_100.final.clstr.faa
wget http://comparative.reefgenomics.org/fna/Coral/Seriatopora_hystrix_cds_100.final.clstr.fna
wget http://comparative.reefgenomics.org/tx_annots/Seriatopora_hystrix.annot.tsv
```

```bash
PREFIX="Seriatopora_hystrix_REEFv1"
```

**CDS**

```bash
cat Seriatopora_hystrix_cds_100.final.clstr.fna \
  | sed -e "s/>Seriatopora_hystrix_/>${PREFIX}___/" \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.transcripts.cds.fna
```

**PEP**

```bash
cat Seriatopora_hystrix_peptides_100.final.clstr.faa \
  | sed -e "s/>Seriatopora_hystrix_/>${PREFIX}___/" -e 's/*$//' \
  | reformat.sh fastawrap=0 iupacToN=f ignorejunk=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.transcripts.pep.faa
```

---

---

---

---

---

# Out group Genomes

## *Amplexidiscus fenestrafer* 

> Corallimorpharia

From: Mediterranean Sea

Version: 1.0

Paper: https://doi.org/10.1111/1755-0998.12680

Website: http://corallimorpharia.reefgenomics.org/download/

```bash
wget http://corallimorpharia.reefgenomics.org/download/afen.genome.fa.gz
wget http://corallimorpharia.reefgenomics.org/download/afen.cds.fa.gz
wget http://corallimorpharia.reefgenomics.org/download/afen.prot.fa.gz
wget http://corallimorpharia.reefgenomics.org/download/afen.gene_models.gff3.gz
```

```bash
PREFIX="Amplexidiscus_fenestrafer_MSMEv1"
```

**Genome**

```bash
zcat afen.genome.fa.gz \
  | sed -e "s/>/>${PREFIX}___/" \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.assembly.fasta
```

**CDS**

```bash
zcat afen.cds.fa.gz \
  | sed -e "s/>/>${PREFIX}___/" \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.genes.cds.fna
```

**PEP**

```bash
zcat afen.prot.fa.gz \
  | sed -e "s/>/>${PREFIX}___/" \
  | reformat.sh fastawrap=0 iupacToN=f ignorejunk=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.genes.pep.faa
```

**GFF**

```bash
zcat afen.gene_models.gff3.gz \
  | awk '$1!~"^#"' \
  | sed -e "s/^/${PREFIX}___/" \
  | sed -e "s/ID=/ID=${PREFIX}___/" \
  | sed -e "s/Parent=/Parent=${PREFIX}___/" \
  | sed -e "s/Name=/Name=${PREFIX}___/" \
  | awk -F'\t' '$3=="mRNA" || $3=="CDS" || $3=="exon"' \
  | gffread \
  | awk '$1!~"^#"' \
  | awk 'BEGIN{OFS=FS="\t"} { if($3=="mRNA") {$3="transcript"; print} else {print} }' \
  > ${PREFIX}.genes.gff3
```

## *Discosoma* sp. 

> Corallimorpharia

From: Mediterranean Sea

Version: 1.0

Paper: https://doi.org/10.1111/1755-0998.12680

Website: http://corallimorpharia.reefgenomics.org/download/

```bash
wget http://corallimorpharia.reefgenomics.org/download/dspp.genome.fa.gz
wget http://corallimorpharia.reefgenomics.org/download/dspp.cds.fa.gz
wget http://corallimorpharia.reefgenomics.org/download/dspp.prot.fa.gz
wget http://corallimorpharia.reefgenomics.org/download/dspp.gene_models.gff3.gz
```

```bash
PREFIX="Discosoma_sp_MSMEv1"
```

**Genome**

```bash
zcat dspp.genome.fa.gz \
  | sed -e "s/>/>${PREFIX}___/" \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.assembly.fasta
```

**CDS**

```bash
zcat dspp.cds.fa.gz \
  | sed -e "s/>/>${PREFIX}___/" \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.genes.cds.fna
```

**PEP**

```bash
zcat dspp.prot.fa.gz \
  | sed -e "s/>/>${PREFIX}___/" \
  | reformat.sh fastawrap=0 iupacToN=f ignorejunk=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.genes.pep.faa
```

**GFF**

```bash
zcat dspp.gene_models.gff3.gz \
  | awk '$1!~"^#"' \
  | sed -e "s/^/${PREFIX}___/" \
  | sed -e "s/ID=/ID=${PREFIX}___/" \
  | sed -e "s/Parent=/Parent=${PREFIX}___/" \
  | sed -e "s/Name=/Name=${PREFIX}___/" \
  | awk -F'\t' '$3=="mRNA" || $3=="CDS" || $3=="exon"' \
  | gffread \
  | awk '$1!~"^#"' \
  | awk 'BEGIN{OFS=FS="\t"} { if($3=="mRNA") {$3="transcript"; print} else {print} }' \
  > ${PREFIX}.genes.gff3
```























## *Exaiptasia diaphana* (*Aiptasia* strain CC7)

> Actiniaria

From: Mediterranean Sea

Version: 1.0

Paper: https://doi.org/10.1073/pnas.1513318112

Website: http://aiptasia.reefgenomics.org/download/ AND https://www.ncbi.nlm.nih.gov/genome/?term=txid2652724

```bash
wget http://aiptasia.reefgenomics.org/download/aiptasia_genome.scaffolds.fa.gz
wget http://aiptasia.reefgenomics.org/download/aiptasia_genome.mRNA.fa.gz
wget http://aiptasia.reefgenomics.org/download/aiptasia_genome.proteins.fa.gz
wget http://aiptasia.reefgenomics.org/download/aiptasia_genome.gff3.gz
```

```bash
PREFIX="Exaiptasia_diaphana_CC7v1"
```

**Genome**

```bash
zcat aiptasia_genome.scaffolds.fa.gz \
  | sed -e "s/>/>${PREFIX}___/" \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.assembly.fasta
```

**GFF**

```bash
zcat aiptasia_genome.gff3.gz \
  | awk '$1!~"^#"' \
  | sed -e "s/^/${PREFIX}___/" \
  | sed -e "s/ID=/ID=${PREFIX}___/" \
  | sed -e "s/Parent=/Parent=${PREFIX}___/" \
  | sed -e "s/Name=/Name=${PREFIX}___/" \
  | awk -F'\t' '$3=="mRNA" || $3=="CDS" || $3=="exon"' \
  | gffread \
  | awk '$1!~"^#"' \
  | awk 'BEGIN{OFS=FS="\t"} { if($3=="mRNA") {$3="transcript"; print} else {print} }' \
  > ${PREFIX}.genes.gff3.orig
  
gffread -S \
  -g ${PREFIX}.assembly.fasta \
  -x ${PREFIX}.genes.cds.fna.tmp \
  -y ${PREFIX}.genes.pep.faa.tmp \
  ${PREFIX}.genes.gff3
```

**NOTE: There are a number of gene models with duplicate names in the input `gff` file. Remove the duplicates which have only 'exon' features.**

> Exaiptasia_diaphana_CC7v1\_\_\_AIPGENE627
> Exaiptasia_diaphana_CC7v1\_\_\_AIPGENE628
> Exaiptasia_diaphana_CC7v1\_\_\_AIPGENE629
> Exaiptasia_diaphana_CC7v1\_\_\_AIPGENE10730
> Exaiptasia_diaphana_CC7v1\_\_\_AIPGENE19416
> Exaiptasia_diaphana_CC7v1\_\_\_AIPGENE19417
> Exaiptasia_diaphana_CC7v1\_\_\_AIPGENE13859
> Exaiptasia_diaphana_CC7v1\_\_\_AIPGENE13860 (remove first single CDS gene as it only covers a large region of N's)

**CDS**

The mRNA sequences (in `aiptasia_genome.mRNA.fa.gz`), which are available through the project website are the full transcripts (including UTRs), so it does not exactly mirror the proteins. Extract the CDS-only transcripts using the new `gff` that we have reformatted. 

```bash
cat ${PREFIX}.genes.cds.fna.tmp \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.genes.cds.fna
rm ${PREFIX}.genes.cds.fna.tmp
rm ${PREFIX}.genes.pep.faa.tmp
```

**PEP**

**NOTE: There looks to be something weird with the phase of some of the genes in the `gff` that causes a large number of inflame stop codons. Use the original protein file available through the project website - it should mirror up with the CDS file we just generated well enough for any downstream analysis.**

```bash
zcat aiptasia_genome.proteins.fa.gz \
  | sed -e "s/>/>${PREFIX}___/" -e 's/*$//' \
  | reformat.sh fastawrap=0 iupacToN=f ignorejunk=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.genes.pep.faa
```

## *Nematostella vectensis* 

> Actiniaria

From: Rhode River in Maryland, USA (from paper https://link.springer.com/content/pdf/10.1007/s00427-002-0214-7.pdf)

Version: 1.0

Paper: https://doi.org/10.1101/2020.10.30.359448

Website: https://simrbase.stowers.org/starletseaanemone

```bash
wget --no-check-certificate https://simrbase.stowers.org/files/pub/nematostella/Nvec/genomes/Nvec200/Nvec200.fasta
wget --no-check-certificate https://simrbase.stowers.org/files/pub/nematostella/Nvec/genomes/Nvec200/analysis/MAKERS/NVEC200.20200813.transcripts.fasta
wget --no-check-certificate https://simrbase.stowers.org/files/pub/nematostella/Nvec/genomes/Nvec200/analysis/MAKERS/NVEC200.20200813.proteins.fasta
wget --no-check-certificate https://simrbase.stowers.org/files/pub/nematostella/Nvec/genomes/Nvec200/analysis/MAKERS/NVEC200.20200813.gff
```

```bash
PREFIX="Nematostella_vectensis_RRUSv1"
```

**Genome**

```bash
cat Nvec200.fasta \
  | sed -e "s/>/>${PREFIX}___/" \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.assembly.fasta
```

**GFF**

```bash
cat NVEC200.20200813.gff \
  | awk '$3=="CDS"' \
  | sed -e "s/^/${PREFIX}___/" \
  | sed -e "s/ID=/ID=${PREFIX}___/" \
  | sed -e "s/Parent=/Parent=${PREFIX}___/" \
  | sed -e "s/Name=/Name=${PREFIX}___/" \
  | awk -F'\t' '$3=="mRNA" || $3=="CDS" || $3=="exon"' \
  | gffread \
  | awk '$1!~"^#"' \
  > ${PREFIX}.genes.gff3
  
gffread -S \
  -g ${PREFIX}.assembly.fasta \
  -x ${PREFIX}.genes.cds.fna.tmp \
  -y ${PREFIX}.genes.pep.faa.tmp \
  ${PREFIX}.genes.gff3
```

**CDS**

The transcript sequences (in `NVEC200.20200813.transcripts.fasta`), which are available through the project website are the full transcripts (including UTRs), so it does not exactly mirror the proteins. Extract the CDS-only transcripts using the new `gff` that we have reformatted. 

```bash
cat ${PREFIX}.genes.cds.fna.tmp \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.genes.cds.fna
rm ${PREFIX}.genes.cds.fna.tmp
rm ${PREFIX}.genes.pep.faa.tmp
```

**PEP**

**The proteins that we have extracted using `gffread` are identical to the ones from the project website.**

```bash
cat NVEC200.20200813.proteins.fasta \
  | sed -e "s/>/>${PREFIX}___/" -e 's/*$//' \
  | reformat.sh fastawrap=0 iupacToN=f ignorejunk=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.genes.pep.faa
```

## *Dendronephthya gigantea* 

> Actiniaria

From: approximately 20 m underwater near Seogwipo, Jeju Island, South Korea (33° 13′39″ N, 126° 34′03″ E) on May 22, 2015

Version: 1.0

Paper: https://doi.org/10.1093/gbe/evz043

Website: https://www.ncbi.nlm.nih.gov/genome/?term=txid151771[orgn] ( RSEI01000000 )

```bash
# Download genome+gff3+CDS+PEP from NCBI
tar -xvf genome_assemblies_cds_fasta.tar
tar -xvf genome_assemblies_prot_fasta.tar
tar -xvf genome_assemblies_genome_gff.tar
tar -xvf genome_assemblies_genome_fasta.tar
```

```bash
PREFIX="Dendronephthya_gigantea_SJKRv1"
```

**Genome**

```bash
zcat ncbi-genomes-2022-05-23/GCF_004324835.1_DenGig_1.0_genomic.fna.gz \
  | sed -e "s/>/>${PREFIX}___/" \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.assembly.fasta
```

**Get longest ORF per gene**

```bash
zcat ncbi-genomes-2022-05-23/GCF_004324835.1_DenGig_1.0_genomic.gff.gz \
  | awk -F'\t' '$3=="CDS"' \
  | sed -e "s/[^\t]*GeneID:\([^,;]*\).*Name=\([^;]*\).*/${PREFIX}___\2\t\1/" \
  > ncbi-genomes-2022-05-23.CDS_features.tsv
  
cat ncbi-genomes-2022-05-23.CDS_features.tsv \
  | awk -F'\t' '$3=="CDS" {
      CDS[$9]=CDS[$9]+(($5-$4)+1)
      GENE[$9]=$10
    } END {
      for (i in CDS){
        print GENE[i]"\t"i"\t"CDS[i]
    } }' \
  > ncbi-genomes-2022-05-23.CDS_lengths.tsv

cat ncbi-genomes-2022-05-23.CDS_lengths.tsv \
  | sort -k1,1 -k3,3nr \
  | awk '!seen[$1]++{print $2}' \
  > ncbi-genomes-2022-05-23.longest.txt
```

**CDS**

```bash
zcat ncbi-genomes-2022-05-23/GCF_004324835.1_DenGig_1.0_cds_from_genomic.fna.gz \
  | sed -e "s/.*protein_id=\([^]]*\)\].*/>${PREFIX}___\1/" \
  | ~/scripts/grepf_fasta.py -f ncbi-genomes-2022-05-23.longest.txt \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.genes.cds.fna
```

**PEP**

```bash
zcat ncbi-genomes-2022-05-23/GCF_004324835.1_DenGig_1.0_protein.faa.gz \
  | sed -e "s/>/>${PREFIX}___/" \
  | ~/scripts/grepf_fasta.py -f ncbi-genomes-2022-05-23.longest.txt \
  | reformat.sh fastawrap=0 iupacToN=f ignorejunk=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.genes.pep.faa
```

**GFF**

```bash
zcat ncbi-genomes-2022-05-23/GCF_004324835.1_DenGig_1.0_genomic.gff.gz \
  | awk -F'\t' '$3=="CDS"' \
  | sed -e "s/[^\t]*Name=\([^;]*\).*/Parent=${PREFIX}___\1\t${PREFIX}___\1/" \
  | sed -e "s/^/${PREFIX}___/" \
  | ~/scripts/grepf_column.py -c 10 -f ncbi-genomes-2022-05-23.longest.txt \
  | cut -f1-9 \
  | gffread \
  | awk '$1!~"^#"' \
  | awk 'BEGIN{OFS=FS="\t"} { if($3!="CDS"){print} else {print; $3="exon"; print} }' \
  > ${PREFIX}.genes.gff3
```

## *Hydra vulgaris* 

> Anthoathecata

Hydra vulgaris (HVU) strain 105 (formerly Hydra magnipapillata strain 105)

From: Strain 105 was established from a single polyp collected by Dr. Tsutomu Sugiyama in September, 1973 from a swamp adjacent to the National Institute of Genetics in Mishima, Japan and has been propagated asexually in a number of labs ever since. (See https://www.ncbi.nlm.nih.gov/biosample/SAMN18321928/)

Version: 3.0

Paper: https://doi.org/10.1126/sciadv.abi5884

Website: https://www.ncbi.nlm.nih.gov/genome/?term=txid6087 ( JAGKSS000000000.1 )

```bash
# Download genome+gff3+CDS+PEP from NCBI
tar -xvf genome_assemblies_cds_fasta.tar
tar -xvf genome_assemblies_prot_fasta.tar
tar -xvf genome_assemblies_genome_gff.tar
tar -xvf genome_assemblies_genome_fasta.tar
```

```bash
PREFIX="Hydra_vulgaris_MIJPv3"
```

**Genome**

```bash
zcat ncbi-genomes-2022-05-23/GCF_022113875.1_Hydra_105_v3_genomic.fna.gz \
  | sed -e "s/>/>${PREFIX}___/" \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.assembly.fasta
```

**Get longest ORF per gene**

```bash
zcat ncbi-genomes-2022-05-23/GCF_022113875.1_Hydra_105_v3_genomic.gff.gz \
  | awk -F'\t' '$3=="CDS"' \
  | sed -e "s/[^\t]*GeneID:\([^,;]*\).*Name=\([^;]*\).*/${PREFIX}___\2\t\1/" \
  > ncbi-genomes-2022-05-23.CDS_features.tsv
  
cat ncbi-genomes-2022-05-23.CDS_features.tsv \
  | awk -F'\t' '$3=="CDS" {
      CDS[$9]=CDS[$9]+(($5-$4)+1)
      GENE[$9]=$10
    } END {
      for (i in CDS){
        print GENE[i]"\t"i"\t"CDS[i]
    } }' \
  > ncbi-genomes-2022-05-23.CDS_lengths.tsv

cat ncbi-genomes-2022-05-23.CDS_lengths.tsv \
  | sort -k1,1 -k3,3nr \
  | awk '!seen[$1]++{print $2}' \
  > ncbi-genomes-2022-05-23.longest.txt
```

**CDS**

```bash
zcat ncbi-genomes-2022-05-23/GCF_022113875.1_Hydra_105_v3_cds_from_genomic.fna.gz \
  | sed -e "s/.*protein_id=\([^]]*\)\].*/>${PREFIX}___\1/" \
  | ~/scripts/grepf_fasta.py -f ncbi-genomes-2022-05-23.longest.txt \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.genes.cds.fna
```

**PEP**

```bash
zcat ncbi-genomes-2022-05-23/GCF_022113875.1_Hydra_105_v3_protein.faa.gz \
  | sed -e "s/>/>${PREFIX}___/" \
  | ~/scripts/grepf_fasta.py -f ncbi-genomes-2022-05-23.longest.txt \
  | reformat.sh fastawrap=0 iupacToN=f ignorejunk=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.genes.pep.faa
```

**GFF**

```bash
zcat ncbi-genomes-2022-05-23/GCF_022113875.1_Hydra_105_v3_genomic.gff.gz \
  | awk -F'\t' '$3=="CDS"' \
  | sed -e "s/[^\t]*Name=\([^;]*\).*/Parent=${PREFIX}___\1\t${PREFIX}___\1/" \
  | sed -e "s/^/${PREFIX}___/" \
  | ~/scripts/grepf_column.py -c 10 -f ncbi-genomes-2022-05-23.longest.txt \
  | cut -f1-9 \
  | gffread \
  | awk '$1!~"^#"' \
  | awk 'BEGIN{OFS=FS="\t"} { if($3!="CDS"){print} else {print; $3="exon"; print} }' \
  > ${PREFIX}.genes.gff3
```

# Out groups Transcripts

## *Anthopleura elegantissima* 

> Actiniaria

Transcripts published in: https://doi.org/10.7554/eLife.13288

Original paper: 

From: 

Version: 1.0

Website: http://comparative.reefgenomics.org/datasets.html

```bash
wget http://comparative.reefgenomics.org/faa/Cnidaria/Anthopleura_elegantissima_peptides_100.final.clstr.faa
wget http://comparative.reefgenomics.org/fna/Cnidaria/Anthopleura_elegantissima_cds_100.final.clstr.fna
wget http://comparative.reefgenomics.org/tx_annots/Anthopleura_elegantissima.annot.tsv
```

```bash
PREFIX="Anthopleura_elegantissima_REEFv1"
```

**CDS**

```bash
cat Anthopleura_elegantissima_cds_100.final.clstr.fna \
  | sed -e "s/>Anthopleura_elegantissima_/>${PREFIX}___/" \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.transcripts.cds.fna
```

**PEP**

```bash
cat Anthopleura_elegantissima_peptides_100.final.clstr.faa \
  | sed -e "s/>Anthopleura_elegantissima_/>${PREFIX}___/" -e 's/*$//' \
  | reformat.sh fastawrap=0 iupacToN=f ignorejunk=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.transcripts.pep.faa
```





























# QC

Get the number of uniq gene names across the GFF, PEP, and CDS files. Make sure this number matches the GFF file.

```bash
for PREFIX in *; do echo $PREFIX; cd $PREFIX/genome_assembly; awk -F'\t' '$3=="transcript"' $PREFIX.genes.gff3 | sed -e 's/.*ID=//' | wc -l; sort <(awk -F'\t' '$3=="transcript"' $PREFIX.genes.gff3 | sed -e 's/.*ID=//') <(grep '>' $PREFIX.genes.pep.faa | sed -e 's/>//') <(grep '>' $PREFIX.genes.cds.fna | sed -e 's/>//') | uniq | wc -l; cd ../../; done
```



```bash
head -n1 */genome_assembly/*.genes.pep.faa
head -n1 */genome_assembly/*.genes.cds.fna
head -n1 */genome_assembly/*.genes.gff3
```



















