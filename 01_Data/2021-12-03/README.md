# Download and reformat coral genome data

Download coral genome assembly + predicted genes (`gff` file + `CDS` and/or `PEP` sequences if available).

Rename the files, and sequences in the files, so they both use the same informative prefix. Reformat genome+CDS+PEP `fasta` files using unique prefix; also reformat sequence by removing description, non-ATGCN characters (genome+CDS only), making isingle line `fasta`, and all uppercase characters.

Helpful list of data: https://xylotrupesgideon.github.io/Cnidarian-Sequence-resources/

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
- *Ephydatia muelleri*
    - There are 84 proteins in the available protein set without available CDS or off features. They will be removed and ignored in downstream analysis.
- Amphimedon queenslandica* HIAUv1.1
    - Removed 1 protein without accurate gff features (i.e., two gene features with same name that produced sequences that didn't match protein). 





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

> Stony coral

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

> Stony coral

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

> Stony coral

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

> Stony coral

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

> Stony coral

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

> Stony coral

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

> Stony coral

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

> Stony coral

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

> Stony coral

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

> Stony coral

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

> Stony coral

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

> Stony coral

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

> Stony coral

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

> Stony coral

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

> Stony coral

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

> Stony coral

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

> Stony coral

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

> Stony coral

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

> Stony coral

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

> Stony coral

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

> Stony coral

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

> Stony coral

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

> Stony coral

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

**Functional annotations from JunMo**

```bash
cat Montipora_capitata_Conserved_Domain_results.txt \
  | sed -e "s/Montipora_capitata../${PREFIX}___/" \
  > ${PREFIX}.genes.Conserved_Domain_results.txt

cat Montipora_capitata_EggNog_results.txt \
  | sed -e "s/Montipora_capitata../${PREFIX}___/" \
  > ${PREFIX}.genes.EggNog_results.txt

cat Montipora_capitata_KEGG_results.txt \
  | sed -e "s/Montipora_capitata../${PREFIX}___/" \
  > ${PREFIX}.genes.KEGG_results.txt
```

## *Montipora capitata* WTHIv1.1

> Stony coral

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

> Stony coral

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

> Stony coral

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

## *Orbicella faveolata* NAPAv1

> Stony coral

From: Panama

Version: 1.0

Paper: Unpublished

Website: https://www.ncbi.nlm.nih.gov/assembly/GCF_002042975.1

```bash
# Download genome+gff3+CDS+PEP from NCBI
unzip GCF_002042975.1.zip 
```

```bash
PREFIX="Orbicella_faveolata_NAPAv1"
```

**Genome**

```bash
cat ncbi_dataset/data/GCF_002042975.1/GCF_002042975.1_ofav_dov_v1_genomic.fna \
  | sed -e "s/>/>${PREFIX}___/" \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.assembly.fasta
```

**Get longest ORF per gene**

```bash
cat ncbi_dataset/data/GCF_002042975.1/genomic.gff \
  | awk -F'\t' '$3=="CDS"' \
  | sed -e "s/[^\t]*GeneID:\([^,;]*\).*Name=\([^;]*\).*/${PREFIX}___\2\t\1/" \
  > ncbi_dataset.CDS_features.tsv
  
cat ncbi_dataset.CDS_features.tsv \
  | awk -F'\t' '$3=="CDS" {
      CDS[$9]=CDS[$9]+(($5-$4)+1)
      GENE[$9]=$10
    } END {
      for (i in CDS){
        print GENE[i]"\t"i"\t"CDS[i]
    } }' \
  > ncbi_dataset.CDS_lengths.tsv

cat ncbi_dataset.CDS_lengths.tsv \
  | sort -k1,1 -k3,3nr \
  | awk '!seen[$1]++{print $2}' \
  > ncbi_dataset.longest.txt
```

**CDS**

```bash
cat ncbi_dataset/data/GCF_002042975.1/cds_from_genomic.fna \
  | sed -e "s/.*protein_id=\([^]]*\)\].*/>${PREFIX}___\1/" \
  | ~/scripts/grepf_fasta.py -f ncbi_dataset.longest.txt \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.genes.cds.fna
```

**PEP**

```bash
cat ncbi_dataset/data/GCF_002042975.1/protein.faa \
  | sed -e "s/>/>${PREFIX}___/" \
  | ~/scripts/grepf_fasta.py -f ncbi_dataset.longest.txt \
  | reformat.sh fastawrap=0 iupacToN=f ignorejunk=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.genes.pep.faa
```

**GFF**

```bash
cat ncbi_dataset/data/GCF_002042975.1/genomic.gff \
  | awk -F'\t' '$3=="CDS"' \
  | sed -e "s/[^\t]*Name=\([^;]*\).*/Parent=${PREFIX}___\1\t${PREFIX}___\1/" \
  | sed -e "s/^/${PREFIX}___/" \
  | ~/scripts/grepf_column.py -c 10 -f ncbi_dataset.longest.txt \
  | cut -f1-9 \
  | gffread \
  | awk '$1!~"^#"' \
  | awk 'BEGIN{OFS=FS="\t"} { if($3!="CDS"){print} else {print; $3="exon"; print} }' \
  > ${PREFIX}.genes.gff3
```

## *Pachyseris speciosa* OIAUv0.12

> Stony coral

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

> Stony coral

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

> Stony coral

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

**Functional annotations from JunMo**

```bash
cat Pocillopora_acuta_Conserved_Domain_results.txt \
  | sed -e "s/Pocillopora_acuta../${PREFIX}___/" \
  > ${PREFIX}.genes.Conserved_Domain_results.txt

cat Pocillopora_acuta_EggNog_results.txt \
  | sed -e "s/Pocillopora_acuta../${PREFIX}___/" \
  > ${PREFIX}.genes.EggNog_results.txt

cat Pocillopora_acuta_KEGG_results.txt \
  | sed -e "s/Pocillopora_acuta../${PREFIX}___/" \
  > ${PREFIX}.genes.KEGG_results.txt
```

## *Pocillopora acuta* LBIDv1

> Stony coral

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

It appears that TransDecoder v5.5.0, when allowed to predict multiple genes per transcript, can produce gene models with the wrong strand information when multiple gene models are predicted on opposite strands on a single exon gene. Without introns to hard enforce a strand (i.e., to make TransDecoder discard genes on the opposite strand) the program gets confused and makes all ORFs on one of the strands. Leading to gene models that when translated are in the wrong direction and are full if stop codons. 

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
  > ${PREFIX}.genes.tmp.gff3 


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

> Stony coral

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

## *Pocillopora cf. effusa* NAPFv3

> Stony coral

From: French Polynesia

Version: 3.0

Paper: https://doi.org/10.1186/s13059-023-02960-7

Website: https://www.genoscope.cns.fr/corals/genomes.html

```bash
wget https://www.genoscope.cns.fr/corals/data/Pocillopora_effusa_v3.fa
wget https://www.genoscope.cns.fr/corals/data/Pocillopora_effusa_v3.annot.gff
wget https://www.genoscope.cns.fr/corals/data/Pocillopora_effusa_v3.annot.pep.fa
```

```bash
PREFIX="Pocillopora_cf_effusa_NAPFv3"
```

**Genome**

```bash
cat Pocillopora_effusa_v3.fa \
  | sed -e "s/>Pocillopora_effusa_/>${PREFIX}___/" \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.assembly.fasta
```

**GFF**

```bash
cat Pocillopora_effusa_v3.annot.gff \
  | sed -e "s/^Pocillopora_meandrina_/${PREFIX}___/" \
  | sed -e "s/ID=Peff_/ID=${PREFIX}___/" \
  | sed -e "s/Parent=Peff/Parent=${PREFIX}___/" \
  | sed -e "s/Alias=Pmea_/Alias=${PREFIX}___/" \
  | awk -F'\t' '$3=="mRNA" || $3=="CDS" || $3=="exon"' \
  | gffread \
  | awk '$1!~"^#"' \
  | awk 'BEGIN{OFS=FS="\t"} { if($3=="mRNA") {$3="transcript"; print} else {print} }' \
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
cat Pocillopora_effusa_v3.annot.pep.fa \
  | sed -e "s/>Peff_/>${PREFIX}___/" \
  | reformat.sh fastawrap=0 iupacToN=f ignorejunk=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  | sed -e 's/*$//' \
  > ${PREFIX}.genes.pep.faa
```

## *Pocillopora meandrina* KBHIv1

> Stony coral

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

**Functional annotations from JunMo**

```bash
cat Pocillopora_meandrina_Conserved_Domain_results.txt \
  | sed -e "s/Pocillopora_meandrina../${PREFIX}___/" \
  > ${PREFIX}.genes.Conserved_Domain_results.txt

cat Pocillopora_meandrina_EggNog_results.txt \
  | sed -e "s/Pocillopora_meandrina../${PREFIX}___/" \
  > ${PREFIX}.genes.EggNog_results.txt

cat Pocillopora_meandrina_KEGG_results.txt \
  | sed -e "s/Pocillopora_meandrina../${PREFIX}___/" \
  > ${PREFIX}.genes.KEGG_results.txt
```

## *Pocillopora verrucosa* RSSAv1

> Stony coral

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

## *Porites astreoides* BBBMv1

Stony coral

From Bailey’s Bay Reef Flats (32°2227N, 64°4437W) in Bermuda

Version: 1.0

Paper: https://doi.org/10.46471/gigabyte.65

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

The `Pastreoides_transcripts_v1.fasta` file is the exons not CDS (so contains UTRs), which is why we used `gffread` to regenerate these files. The below commands confirm that the original and new protein files are identical and all is consistent between old and new genes.

```bash
grep -v '>' Porites_astreoides_BBBMv1.genes.pep.faa | sort | md5sum
512b5c3eed317e9e0de2d15ca19fbeb9  -
cat Pastreoides_proteins_v1.fasta | reformat.sh fastawrap=0 iupacToN=f ignorejunk=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa | grep -v '>' | sort | md5sum
512b5c3eed317e9e0de2d15ca19fbeb9  -
```

## *Porites australiensis* SLJPv1

> Stony coral

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

## *Porites compressa* KBHIv1

> Stony coral

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

**Functional annotations from JunMo**

```bash
cat Porites_compressa_Conserved_Domain_results.txt \
  | sed -e "s/Porites_compressa../${PREFIX}___/" \
  > ${PREFIX}.genes.Conserved_Domain_results.txt

cat Porites_compressa_EggNog_results.txt \
  | sed -e "s/Porites_compressa../${PREFIX}___/" \
  > ${PREFIX}.genes.EggNog_results.txt

cat Porites_compressa_KEGG_results.txt \
  | sed -e "s/Porites_compressa../${PREFIX}___/" \
  > ${PREFIX}.genes.KEGG_results.txt
```









## *Porites evermanni* NAPFv1

> Stony coral

From: French Polynesia

Version: 1.0

Paper: https://doi.org/10.1186/s13059-023-02960-7

Website: https://www.genoscope.cns.fr/corals/genomes.html

```bash
wget https://www.genoscope.cns.fr/corals/data/Porites_evermanni_v1.fa
wget https://www.genoscope.cns.fr/corals/data/Porites_evermanni_v1.annot.gff
wget https://www.genoscope.cns.fr/corals/data/Porites_evermanni_v1.annot.pep.fa
```

```bash
PREFIX="Porites_evermanni_NAPFv1"
```

**Genome**

```bash
cat Porites_evermanni_v1.fa \
  | sed -e "s/>Porites_evermani_/>${PREFIX}___/" \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.assembly.fasta
```

**GFF**

```bash
cat Porites_evermanni_v1.annot.gff \
  | sed -e "s/^Porites_evermani_/${PREFIX}___/" \
  | sed -e "s/ID=Peve_/ID=${PREFIX}___/" \
  | sed -e "s/Parent=Peve_/Parent=${PREFIX}___/" \
  | sed -e "s/Name=Peve_/Name=${PREFIX}___/" \
  | awk -F'\t' '$3=="mRNA" || $3=="CDS" || $3=="exon"' \
  | gffread \
  | awk '$1!~"^#"' \
  | awk 'BEGIN{OFS=FS="\t"} { if($3=="mRNA") {$3="transcript"; print} else {print} }' \
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
cat Porites_evermanni_v1.annot.pep.fa \
  | sed -e "s/>Peve_/>${PREFIX}___/" \
  | reformat.sh fastawrap=0 iupacToN=f ignorejunk=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  | sed -e 's/*$//' \
  > ${PREFIX}.genes.pep.faa
```

## *Porites lobata* NAPFv3

> Stony coral

From: French Polynesia

Version: 3.0

Paper: https://doi.org/10.1186/s13059-023-02960-7

Website: https://www.genoscope.cns.fr/corals/genomes.html

```bash
wget https://www.genoscope.cns.fr/corals/data/Porites_lobata_v3.fa
wget https://www.genoscope.cns.fr/corals/data/Porites_lobata_v3.annot.gff
wget https://www.genoscope.cns.fr/corals/data/Porites_lobata_v3.annot.pep.fa
```

```bash
PREFIX="Porites_lobata_NAPFv3"
```

**Genome**

```bash
cat Porites_lobata_v3.fa \
  | sed -e "s/>Porites_lobata_/>${PREFIX}___/" \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.assembly.fasta
```

**GFF**

```bash
cat Porites_lobata_v3.annot.gff \
  | sed -e "s/^Porites_lobata_/${PREFIX}___/" \
  | sed -e "s/ID=Plob_/ID=${PREFIX}___/" \
  | sed -e "s/Parent=Plob_/Parent=${PREFIX}___/" \
  | sed -e "s/Alias=Plob_/Alias=${PREFIX}___/" \
  | awk -F'\t' '$3=="mRNA" || $3=="CDS" || $3=="exon"' \
  | gffread \
  | awk '$1!~"^#"' \
  | awk 'BEGIN{OFS=FS="\t"} { if($3=="mRNA") {$3="transcript"; print} else {print} }' \
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
cat Porites_lobata_v3.annot.pep.fa \
  | sed -e "s/>Plob_/>${PREFIX}___/" \
  | reformat.sh fastawrap=0 iupacToN=f ignorejunk=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  | sed -e 's/*$//' \
  > ${PREFIX}.genes.pep.faa
```

## *Porites lutea* OIAUv1.1

> Stony coral

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
  | sed -e "s/>plut2./>${PREFIX}___/" -e "s/>jamg1.model./>${PREFIX}___/" \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.genes.cds.fna
```

**PEP**

```bash
zcat plut2v1.1.proteins.fasta.gz \
  | sed -e "s/>plut2./>${PREFIX}___/"  -e "s/>jamg1.model./>${PREFIX}___/" \
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
  | sed -e "s/ID=jamg1.model./ID=${PREFIX}___/" \
  | sed -e "s/Parent=jamg1.model./Parent=${PREFIX}___/" \
  | sed -e "s/Name=jamg1.model./Name=${PREFIX}___/" \
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

> Stony coral

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

> Stony coral

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

## *Astreopora* sp. REEFv1

> Stony coral

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

## *Acropora hyacinthus* NANAv1

> Stony coral

Original paper: Unpublished

From: Unpublished

Version: 1.0

Website: https://www.dropbox.com/s/yc9y8mxw02tc0ik/ahyacinthus_july2014.zip?dl=0&file_subpath=%2Fahyacinthus_july2014

```bash
# Downloaded manually
unzip ahyacinthus_july2014.zip
```

```bash
PREFIX="Acropora_hyacinthus_NANAv1"
```

**CDS**

```bash
cat ahyacinthus_july2014/ahya_CDS.fas \
  | sed -e "s/>/>${PREFIX}___/" \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.transcripts.cds.fna
```

**PEP**

```bash
cat ahyacinthus_july2014/ahya_PRO.fas \
  | sed -e "s/>/>${PREFIX}___/" -e 's/*$//' \
  | reformat.sh fastawrap=0 iupacToN=f ignorejunk=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.transcripts.pep.faa
```

## *Acropora muricata* XSCSv1

> Stony coral

Paper: https://doi.org/10.3390/ijms231911135

From: Xisha Islands in the South China Sea (latitude 15°40′–17°10′ N, longitude 111°–113° E)

Version: 1.0

Website: https://figshare.com/articles/dataset/Full-length_transcriptome_maps_of_reef-building_coral_illuminate_the_molecular_basis_of_calcification_symbiosis_and_circadian_genes/19403021

```bash
# Downloaded manually
```

```bash
PREFIX="Acropora_muricata_XSCSv1"
```

**CDS**

```bash
cat A.muricata_cds.fasta \
  | sed -e "s/>/>${PREFIX}___/" -e 's@/@_@g' -e 's@|@_@g' \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.transcripts.cds.fna
```

**PEP**

```bash
cat A.muricata_pep.fasta \
  | sed -e "s/>/>${PREFIX}___/" -e 's@/@_@g' -e 's@|@_@g' -e 's/*$//' \
  | reformat.sh fastawrap=0 iupacToN=f ignorejunk=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.transcripts.pep.faa
```

## *Ctenactis echinata* NAJPv1

> Stony coral

Paper: NA

From: Okinawa, Japan

Version: 1.0

Website: NCBI:  GDZV00000000.1

```bash
wget https://sra-download.ncbi.nlm.nih.gov/traces/wgs03/wgs_aux/GD/ZV/GDZV01/GDZV01.1.fsa_nt.gz
```

```bash
PREFIX="Ctenactis_echinata_NAJPv1"
```

**Transcripts**

Extract transcripts and predict ORFs.

```bash
zcat GDZV01.1.fsa_nt.gz \
  | sed -e "s/>/>${PREFIX}___/" \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.transcripts.fasta
```

Get gene to transcript name mappings.

```bash
zcat GDZV01.1.fsa_nt.gz \
  | grep '>' | sed -e "s/>/${PREFIX}___/" \
  | awk '{print $5"\t"$1}' | sed -e 's/_seq[^\t]*//' \
  > ${PREFIX}.transcripts.gene2trans.txt
```

Predict ORFs in transcripts using `transdecoder`.

```bash
./run_Transdecoder.sh
```

Get best ORF per "gene"

```bash
grep '>' ${PREFIX}.transcripts.fasta.transdecoder.pep \
  | sed -e 's/.* \([^~]*\)~~\([^ ]*\) .*score=\([^ ,]*\).*/\1\t\2\t\3/' \
  | sort -k1,1 -k3,3nr \
  | awk -F'\t' '!seen[$1]++{print $2}' \
  > ${PREFIX}.transcripts.fasta.transdecoder.topPerGene.txt
```

**CDS**

```bash
cat ${PREFIX}.transcripts.fasta.transdecoder.cds \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  | ~/scripts/grepf_fasta.py -f ${PREFIX}.transcripts.fasta.transdecoder.topPerGene.txt \
  > ${PREFIX}.transcripts.cds.fna
```

**PEP**

```bash
cat ${PREFIX}.transcripts.fasta.transdecoder.pep \
  | sed -e 's/*$//' \
  | reformat.sh fastawrap=0 iupacToN=f ignorejunk=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  | ~/scripts/grepf_fasta.py -f ${PREFIX}.transcripts.fasta.transdecoder.topPerGene.txt \
  > ${PREFIX}.transcripts.pep.faa
```

## *Dipsastraea lizardensis* NAJPv1 (formally *Favia lizardensis*)

> Stony coral

Paper: NA

From: Okinawa, Japan

Version: 1.0

Website: NCBI: GDZU00000000.1

```bash
wget https://sra-download.ncbi.nlm.nih.gov/traces/wgs03/wgs_aux/GD/ZU/GDZU01/GDZU01.1.fsa_nt.gz
```

```bash
PREFIX="Dipsastraea_lizardensis_NAJPv1"
```

**Transcripts**

Extract transcripts and predict ORFs.

```bash
zcat GDZU01.1.fsa_nt.gz \
  | sed -e "s/>/>${PREFIX}___/" \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.transcripts.fasta
```

Get gene to transcript name mappings.

```bash
zcat GDZU01.1.fsa_nt.gz \
  | grep '>' | sed -e "s/>/${PREFIX}___/" \
  | awk '{print $5"\t"$1}' | sed -e 's/_seq[^\t]*//' \
  > ${PREFIX}.transcripts.gene2trans.txt
```

Predict ORFs in transcripts using `transdecoder`.

```bash
./run_Transdecoder.sh
```

Get best ORF per "gene"

```bash
grep '>' ${PREFIX}.transcripts.fasta.transdecoder.pep \
  | sed -e 's/.* \([^~]*\)~~\([^ ]*\) .*score=\([^ ,]*\).*/\1\t\2\t\3/' \
  | sort -k1,1 -k3,3nr \
  | awk -F'\t' '!seen[$1]++{print $2}' \
  > ${PREFIX}.transcripts.fasta.transdecoder.topPerGene.txt
```

**CDS**

```bash
cat ${PREFIX}.transcripts.fasta.transdecoder.cds \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  | ~/scripts/grepf_fasta.py -f ${PREFIX}.transcripts.fasta.transdecoder.topPerGene.txt \
  > ${PREFIX}.transcripts.cds.fna
```

**PEP**

```bash
cat ${PREFIX}.transcripts.fasta.transdecoder.pep \
  | sed -e 's/*$//' \
  | reformat.sh fastawrap=0 iupacToN=f ignorejunk=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  | ~/scripts/grepf_fasta.py -f ${PREFIX}.transcripts.fasta.transdecoder.topPerGene.txt \
  > ${PREFIX}.transcripts.pep.faa
```

## *Favia* sp. REEFv1

> Stony coral

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

## *Fungia scutaria* REEFv1

> Stony coral

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

## *Madracis auretenra* REEFv1

> Stony coral

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

## *Montastraea cavernosa* REEFv1

> Stony coral

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

## *Montastraea faveolata* REEFv1

> Stony coral

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

## *Montipora aequituberculata* NANAv1

> Stony coral

Original paper: Unpublished

From: Unpublished

Version: 1.0

Website: https://www.dropbox.com/s/mnkqwz7oi3l226c/maequituberculata_transcriptome_july2014.zip?dl=0

```bash
# Downloaded manually
unzip maequituberculata_transcriptome_july2014.zip
```

```bash
PREFIX="Montipora_aequituberculata_NANAv1"
```

**CDS**

```bash
cat maequituberculata_transcriptome_july2014/maeq_coral_CDS.fas \
  | sed -e "s/>/>${PREFIX}___/" \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.transcripts.cds.fna
```

**PEP**

```bash
cat maequituberculata_transcriptome_july2014/maeq_coral_PRO.fas \
  | sed -e "s/>/>${PREFIX}___/" -e 's/*$//' \
  | reformat.sh fastawrap=0 iupacToN=f ignorejunk=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.transcripts.pep.faa
```

## *Montipora foliosa* XSCSv1

> Stony coral

Paper: https://doi.org/10.3390/ijms231911135

From: Xisha Islands in the South China Sea (latitude 15°40′–17°10′ N, longitude 111°–113° E)

Version: 1.0

Website: https://figshare.com/articles/dataset/Full-length_transcriptome_maps_of_reef-building_coral_illuminate_the_molecular_basis_of_calcification_symbiosis_and_circadian_genes/19403021

```bash
# Downloaded manually
```

```bash
PREFIX="Montipora_foliosa_XSCSv1"
```

**CDS**

```bash
cat M.foliosa_cds.fasta \
  | sed -e "s/>/>${PREFIX}___/" -e 's@/@_@g' -e 's@|@_@g' \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.transcripts.cds.fna
```

**PEP**

```bash
cat M.foliosa_pep.fasta \
  | sed -e "s/>/>${PREFIX}___/" -e 's@/@_@g' -e 's@|@_@g' -e 's/*$//' \
  | reformat.sh fastawrap=0 iupacToN=f ignorejunk=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.transcripts.pep.faa
```

## *Platygyra carnosus* REEFv1

> Stony coral

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

## *Pocillopora damicornis* XSCSv1

> Stony coral

Paper: https://doi.org/10.3390/ijms231911135

From: Xisha Islands in the South China Sea (latitude 15°40′–17°10′ N, longitude 111°–113° E)

Version: 1.0

Website: https://figshare.com/articles/dataset/Full-length_transcriptome_maps_of_reef-building_coral_illuminate_the_molecular_basis_of_calcification_symbiosis_and_circadian_genes/19403021

```bash
# Downloaded manually
```

```bash
PREFIX="Pocillopora_damicornis_XSCSv1"
```

**CDS**

```bash
cat P.damicornis_cds.fasta \
  | sed -e "s/>/>${PREFIX}___/" -e 's@/@_@g' -e 's@|@_@g' \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.transcripts.cds.fna
```

**PEP**

```bash
cat P.damicornis_pep.fasta \
  | sed -e "s/>/>${PREFIX}___/" -e 's@/@_@g' -e 's@|@_@g' -e 's/*$//' \
  | reformat.sh fastawrap=0 iupacToN=f ignorejunk=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.transcripts.pep.faa
```

## *Pocillopora verrucosa* XSCSv1

> Stony coral

Paper: https://doi.org/10.3390/ijms231911135

From: Xisha Islands in the South China Sea (latitude 15°40′–17°10′ N, longitude 111°–113° E)

Version: 1.0

Website: https://figshare.com/articles/dataset/Full-length_transcriptome_maps_of_reef-building_coral_illuminate_the_molecular_basis_of_calcification_symbiosis_and_circadian_genes/19403021

```bash
# Downloaded manually
```

```bash
PREFIX="Pocillopora_verrucosa_XSCSv1"
```

**CDS**

```bash
cat P.verrucosa_cds.fasta \
  | sed -e "s/>/>${PREFIX}___/" -e 's@/@_@g' -e 's@|@_@g' \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.transcripts.cds.fna
```

**PEP**

```bash
cat P.verrucosa_pep.fasta \
  | sed -e "s/>/>${PREFIX}___/" -e 's@/@_@g' -e 's@|@_@g' -e 's/*$//' \
  | reformat.sh fastawrap=0 iupacToN=f ignorejunk=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.transcripts.pep.faa
```

## *Porites astreoides* REEFv1

> Stony coral

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

## *Porites lobata* REEFv1

> Stony coral

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

## *Pseudodiploria strigosa* REEFv1

> Stony coral

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

## *Seriatopora hystrix* REEFv1

> Stony coral

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

# Outgroup Genomes

## *Actinia equina* RHUKv1

> Sea anemone

From: United Kingdom: Rhosneigr

Version: 1.0

Paper: https://doi.org/10.1016/j.margen.2020.100753

Website: http://aequ.reefgenomics.org/download/

```bash
wget http://aequ.reefgenomics.org/download/equina_smartden.arrow4.noredun.fa.gz
wget http://aequ.reefgenomics.org/download/equina_smart.rnam-trna.merged.ggf.curated.remredun.proteins.gff3.gz
wget http://aequ.reefgenomics.org/download/equina_smart.rnam-trna.merged.ggf.curated.remredun.aa.fa.gz
wget http://aequ.reefgenomics.org/download/equina_smart.rnam-trna.merged.ggf.curated.remredun.nucl.fa.gz
```

```bash
PREFIX="Actinia_equina_RHUKv1"
```

**Genome**

```bash
zcat equina_smartden.arrow4.noredun.fa.gz \
  | sed -e "s/>/>${PREFIX}___/" \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.assembly.fasta
```

**Get longest ORF per gene**

Ignore the tRNAscan and RNAmmer features as we only want protein coding genes.

```bash
zcat equina_smart.rnam-trna.merged.ggf.curated.remredun.proteins.gff3.gz \
  | dos2unix \
  | awk -F'\t' '$3=="mRNA"' \
  | sed -e "s/.*ID=\([^,;]*\).*;Parent=\([^,;]*\).*/${PREFIX}___\1\t${PREFIX}___\2/" \
  > dataset.gene2mRNA.tsv

zcat equina_smart.rnam-trna.merged.ggf.curated.remredun.proteins.gff3.gz \
  | dos2unix \
  | awk -F'\t' '$3=="CDS"' \
  | awk -F'\t' '$2!="tRNAscan-SE-1.3.1" && $2!="RNAmmer-1.2"' \
  | sed -e "s/[^\t]*ID=\([^,;]*\).*;Parent=\([^,;]*\).*/${PREFIX}___\2\t${PREFIX}___\1/" \
  > dataset.CDS_features.tsv
  
cat dataset.CDS_features.tsv \
  | awk -F'\t' '$3=="CDS" {
      CDS[$9]=CDS[$9]+(($5-$4)+1)
      GENE[$9]=$10
    } END {
      for (i in CDS){
        print GENE[i]"\t"i"\t"CDS[i]
    } }' \
  > dataset.CDS_lengths.tsv

cat dataset.gene2mRNA.tsv \
  | ~/scripts/add_value_to_table.py -a <(cut -f2,3 dataset.CDS_lengths.tsv) \
  | sort -k2,2 -k3,3nr -k1,1 \
  > dataset.gene2mRNA.CDS_lengths.tsv

cat dataset.gene2mRNA.CDS_lengths.tsv \
  | awk '!seen[$2]++{print $1}' \
  > dataset.longest.txt
```

**CDS**

```bash
zcat equina_smart.rnam-trna.merged.ggf.curated.remredun.nucl.fa.gz \
  | sed -e "s/>/>${PREFIX}___/" \
  | ~/scripts/grepf_fasta.py -f dataset.longest.txt \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.genes.cds.fna
```

**PEP**

```bash
zcat equina_smart.rnam-trna.merged.ggf.curated.remredun.aa.fa.gz \
  | sed -e "s/>/>${PREFIX}___/" \
  |  ~/scripts/grepf_fasta.py -f dataset.longest.txt \
  | reformat.sh fastawrap=0 iupacToN=f ignorejunk=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  | sed -e 's/*$//' \
  > ${PREFIX}.genes.pep.faa
```

**GFF**

```bash
zcat equina_smart.rnam-trna.merged.ggf.curated.remredun.proteins.gff3.gz \
  | awk '$1!~"^#"' \
  | sed -e "s/^/${PREFIX}___/" \
  | sed -e "s/ID=/ID=${PREFIX}___/" \
  | sed -e "s/Parent=/Parent=${PREFIX}___/" \
  | sed -e "s/Name=/Name=${PREFIX}___/" \
  | awk -F'\t' '$3=="CDS" || $3=="exon"' \
  | awk -F'\t' '$2!="tRNAscan-SE-1.3.1" && $2!="RNAmmer-1.2"' \
  | dos2unix \
  | sed -e 's/.*Parent=\([^,;]*\).*/\0\t\1/' \
  | ~/scripts/grepf_column.py -c 10 -f dataset.longest.txt \
  | cut -f1-9 \
  | gffread \
  | awk '$1!~"^#"' \
  | awk 'BEGIN{OFS=FS="\t"} { if($3=="mRNA") {$3="transcript"; print} else {print} }' \
  > ${PREFIX}.genes.gff3
```

## *Alatina alata* BNNLv1

> Jellyfish

From: Bonaire, The Netherlands (April 2014, 22:00–01:00)

Version: 1.0

Paper: https://doi.org/10.1093/gigascience/giz069

Website: http://ryanlab.whitney.ufl.edu/genomes/Aala/

```bash
wget http://ryanlab.whitney.ufl.edu/genomes/Aala/downloads/Aala.11.genome.fsa.gz
git clone https://github.com/josephryan/Ohdera_et_al_2018.git
```

```bash
PREFIX="Alatina_alata_BNNLv1"
```

**Genome**

```bash
zcat Aala.11.genome.fsa.gz \
  | sed -e "s/>Aala./>${PREFIX}___/" \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.assembly.fasta
```

**CDS**

```bash
zcat Ohdera_et_al_2018/Condingseq/Aala.augustus.07.codingseq.gz \
  | sed -e "s/>Aala./>${PREFIX}___/" \
  | seqkit rename \
  | reformat.sh fastawrap=0 iupacToN=f ignorejunk=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.genes.cds.fna
```

**PEP**

```bash
zcat Ohdera_et_al_2018/AA_Files/Aala.augustus.07.faa.gz \
  | sed -e "s/>Aala_/>${PREFIX}___/" \
  | seqkit rename \
  | reformat.sh fastawrap=0 iupacToN=f ignorejunk=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  | sed -e 's/*$//' \
  > ${PREFIX}.genes.pep.faa
```

**GFF**

No CDS or GFF files were provided.

## *Amphimedon queenslandica* HIAUv1.1

> Sponge

From: Heron Island Reef, Australia (23.45 S 151.92 E)

Version: 1.1

Paper: https://doi.org/10.1038/nature09201

Website: NCBI: GCA_000090795.2

NOTE: There are two CDS and GFF features that have the same name; there is only one protein with this same. The protein (NP_001266236.1) does not match the full CDs of either gene (the protein is much longer and has a lot of unexplained sequence), so we will remove it from downstream analysis.

```bash
# From NCBI
unzip GCF_000090795.2.zip
```

```bash
PREFIX="Amphimedon_queenslandica_HIAUv1.1"
```

**Genome**

```bash
cat ncbi_dataset/data/GCF_000090795.2/GCF_000090795.2_v1.1_genomic.fna \
  | sed -e "s/>/>${PREFIX}___/" \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.assembly.fasta
```

**Get longest ORF per gene**

```bash
cat ncbi_dataset/data/GCF_000090795.2/genomic.gff \
  | awk -F'\t' '$3=="CDS"' \
  | sed -e "s/[^\t]*GeneID:\([^,;]*\).*Name=\([^;]*\).*/${PREFIX}___\2\t\1/" \
  > ncbi_dataset.CDS_features.tsv
  
cat ncbi_dataset.CDS_features.tsv \
  | awk -F'\t' '$3=="CDS" {
      CDS[$9]=CDS[$9]+(($5-$4)+1)
      GENE[$9]=$10
    } END {
      for (i in CDS){
        print GENE[i]"\t"i"\t"CDS[i]
    } }' \
  > ncbi_dataset.CDS_lengths.tsv

cat ncbi_dataset.CDS_lengths.tsv \
  | grep -v 'Amphimedon_queenslandica_HIAUv1.1___NP_001266236.1' \
  | sort -k1,1 -k3,3nr \
  | awk '!seen[$1]++{print $2}' \
  > ncbi_dataset.longest.txt
```

**CDS**

```bash
cat ncbi_dataset/data/GCF_000090795.2/cds_from_genomic.fna \
  | sed -e "s/.*protein_id=\([^]]*\)\].*/>${PREFIX}___\1/" \
  | ~/scripts/grepf_fasta.py -f ncbi_dataset.longest.txt \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.genes.cds.fna
```

**PEP**

```bash
cat ncbi_dataset/data/GCF_000090795.2/protein.faa \
  | sed -e "s/>/>${PREFIX}___/" \
  | ~/scripts/grepf_fasta.py -f ncbi_dataset.longest.txt \
  | reformat.sh fastawrap=0 iupacToN=f ignorejunk=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.genes.pep.faa
```

**GFF**

```bash
cat ncbi_dataset/data/GCF_000090795.2/genomic.gff \
  | awk -F'\t' '$3=="CDS"' \
  | sed -e "s/[^\t]*Name=\([^;]*\).*/Parent=${PREFIX}___\1\t${PREFIX}___\1/" \
  | sed -e "s/^/${PREFIX}___/" \
  | ~/scripts/grepf_column.py -c 10 -f ncbi_dataset.longest.txt \
  | cut -f 1-9 \
  | gffread \
  | awk '$1!~"^#"' \
  | awk 'BEGIN{OFS=FS="\t"} { if($3!="CDS"){print} else {print; $3="exon"; print} }' \
  > ${PREFIX}.genes.gff3
```

## *Amplexidiscus fenestrafer* MSMEv1

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

## *Aurelia aurita* loc. Atlantic Ocean ALTOv1

> Jellyfish

From: Atlantic Ocean

Version: 1.0

Paper: https://doi.org/10.1038/s41559-019-0853-y

Website: https://marinegenomics.oist.jp/aurelia_aurita/viewer/download?project_id=69

```bash
wget https://marinegenomics.oist.jp/aurelia_aurita/download/AUR21_r04_250316.fa.gz
wget https://marinegenomics.oist.jp/aurelia_aurita/download/AUR21_r04_wref.gff3.gz
wget https://marinegenomics.oist.jp/aurelia_aurita/download/AUR21_r04_proteins.fa.gz
wget https://marinegenomics.oist.jp/aurelia_aurita/download/AUR21_r04_mRNA.fa.gz
```

```bash
PREFIX="Aurelia_aurita_ALTOv1"
```

**Genome**

```bash
zcat AUR21_r04_250316.fa.gz \
  | sed -e "s/>/>${PREFIX}___/" \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.assembly.fasta
```

**Get longest ORF per gene**

```bash
zcat AUR21_r04_wref.gff3.gz \
  | awk -F'\t' '$1!~"#" && $3=="transcript"' \
  | sed -e "s/.*ID=\([^,;]*\).*Parent=\([^,;]*\).*/${PREFIX}___\1\t${PREFIX}___\2/" \
  > dataset.mRNA2gene.tsv
  
zcat AUR21_r04_wref.gff3.gz \
  | awk -F'\t' '$1!~"#" && $3=="CDS" {
      split($9,a,"Parent=");
      split(a[2],b,";");
      ID=b[1]
      CDS[ID]=CDS[ID]+(($5-$4)+1)
    } END {
      for (i in CDS){
        print i"\t"CDS[i]
    } }' \
  | sed -e "s/^/${PREFIX}___/" \
  > dataset.CDS_length_per_mRNA.tsv

~/scripts/add_value_to_table_SQLite3.py -d "NA" \
    -i dataset.mRNA2gene.tsv \
    -a dataset.CDS_length_per_mRNA.tsv \
  | sort -k2,2 -k3,3nr \
  | awk -F'\t' '!seen[$2]++{print $1}' \
  > ncbi_dataset.longest.txt
```

**PEP**

NOTE: One mRNA doesn't have a protein sequence (possible pseudo gene?), so have to use the PEP file to filer CDS and GFF downstream.

```bash
zcat AUR21_r04_proteins.fa.gz \
  | sed -e "s/>/>${PREFIX}___/" \
  | reformat.sh fastawrap=0 iupacToN=f ignorejunk=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  | ~/scripts/grepf_fasta.py -f ncbi_dataset.longest.txt \
  > ${PREFIX}.genes.pep.faa
```

**CDS**

```bash
zcat AUR21_r04_mRNA.fa.gz \
  | sed -e "s/>/>${PREFIX}___/" \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  | ~/scripts/grepf_fasta.py -f <(grep '>' ${PREFIX}.genes.pep.faa | sed -e 's/>//') \
  > ${PREFIX}.genes.cds.fna
```

**GFF**

```bash
zcat AUR21_r04_wref.gff3.gz \
  | awk '$1!~"^#"' \
  | sed -e "s/^/${PREFIX}___/" \
  | sed -e "s/ID=/ID=${PREFIX}___/" \
  | sed -e "s/Parent=/Parent=${PREFIX}___/" \
  | sed -e "s/Name=/Name=${PREFIX}___/" \
  | awk -F'\t' '$3=="mRNA" || $3=="CDS" || $3=="exon"' \
  | gffread \
  | awk '$1!~"^#"' \
  | awk 'BEGIN{OFS=FS="\t"} { if($3=="mRNA") {$3="transcript"; print} else {print} }' \
  | sed -e 's/.*=\([^=]*\)/\0\t\1/' \
  | ~/scripts/grepf_column.py -c 10 -f <(grep '>' ${PREFIX}.genes.pep.faa | sed -e 's/>//') \
  | cut -f1-9 \
  > ${PREFIX}.genes.gff3
```

## *Aurelia aurita* loc. Roscoff RSPOv1

> Jellyfish

From: Atlantic Ocean

Version: 1.0

Paper: https://doi.org/10.1038/s41559-019-0853-y

Website: https://marinegenomics.oist.jp/aurelia_aurita_pacific/viewer/download?project_id=75

```bash
wget https://marinegenomics.oist.jp/aurelia_aurita_pacific/download/ARSv1_genome_assembly.fa.gz
wget https://marinegenomics.oist.jp/aurelia_aurita_pacific/download/ARSv1_genome_assembly.gff3.gz
wget https://marinegenomics.oist.jp/aurelia_aurita_pacific/download/ARSv1_proteins.fa.gz
wget https://marinegenomics.oist.jp/aurelia_aurita_pacific/download/ARSv1_mRNA.fa.gz
```

```bash
PREFIX="Aurelia_aurita_RSPOv1"
```

**Genome**

```bash
zcat ARSv1_genome_assembly.fa.gz \
  | sed -e "s/>/>${PREFIX}___/" \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.assembly.fasta
```

**Get longest ORF per gene**

```bash
zcat ARSv1_genome_assembly.gff3.gz \
  | awk -F'\t' '$1!~"#" && $3=="transcript"' \
  | sed -e "s/.*ID=\([^,;]*\).*Parent=\([^,;]*\).*/${PREFIX}___\1\t${PREFIX}___\2/" \
  > dataset.mRNA2gene.tsv
  
zcat ARSv1_genome_assembly.gff3.gz \
  | awk -F'\t' '$1!~"#" && $3=="CDS" {
      split($9,a,"Parent=");
      split(a[2],b,";");
      ID=b[1]
      CDS[ID]=CDS[ID]+(($5-$4)+1)
    } END {
      for (i in CDS){
        print i"\t"CDS[i]
    } }' \
  | sed -e "s/^/${PREFIX}___/" \
  > dataset.CDS_length_per_mRNA.tsv

~/scripts/add_value_to_table_SQLite3.py -d "NA" \
    -i dataset.mRNA2gene.tsv \
    -a dataset.CDS_length_per_mRNA.tsv \
  | sort -k2,2 -k3,3nr \
  | awk -F'\t' '!seen[$2]++{print $1}' \
  > ncbi_dataset.longest.txt
```

**PEP**

NOTE: 6 mRNAs doesn't have a protein sequences (possible pseudo genes?), so have to use the PEP file to filer CDS and GFF downstream.

```bash
zcat ARSv1_proteins.fa.gz \
  | sed -e "s/>/>${PREFIX}___/" \
  | reformat.sh fastawrap=0 iupacToN=f ignorejunk=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  | ~/scripts/grepf_fasta.py -f ncbi_dataset.longest.txt \
  > ${PREFIX}.genes.pep.faa
```

**CDS**

```bash
zcat ARSv1_mRNA.fa.gz \
  | sed -e "s/>/>${PREFIX}___/" \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  | ~/scripts/grepf_fasta.py -f <(grep '>' ${PREFIX}.genes.pep.faa | sed -e 's/>//') \
  > ${PREFIX}.genes.cds.fna
```

**GFF**

```bash
zcat ARSv1_genome_assembly.gff3.gz \
  | awk '$1!~"^#"' \
  | sed -e "s/^/${PREFIX}___/" \
  | sed -e "s/ID=/ID=${PREFIX}___/" \
  | sed -e "s/Parent=/Parent=${PREFIX}___/" \
  | sed -e "s/Name=/Name=${PREFIX}___/" \
  | awk -F'\t' '$3=="mRNA" || $3=="CDS" || $3=="exon"' \
  | gffread \
  | awk '$1!~"^#"' \
  | awk 'BEGIN{OFS=FS="\t"} { if($3=="mRNA") {$3="transcript"; print} else {print} }' \
  | sed -e 's/.*=\([^=]*\)/\0\t\1/' \
  | ~/scripts/grepf_column.py -c 10 -f <(grep '>' ${PREFIX}.genes.pep.faa | sed -e 's/>//') \
  | cut -f1-9 \
  > ${PREFIX}.genes.gff3
```

## *Calvadosia cruxmelitensis* CRGBv3.2

> Jellyfish

From: Chimney Rock, off the coast of Penzance, Cornwall, England (January 2013)

Version: 3.2

Paper: https://doi.org/10.1093/gigascience/giz069

Website: http://ryanlab.whitney.ufl.edu/genomes/Ccrux/

```bash
wget http://ryanlab.whitney.ufl.edu/genomes/Ccrux/downloads/Ccrux.v3.greater200.fasta.gz
git clone https://github.com/josephryan/Ohdera_et_al_2018.git
```

```bash
PREFIX="Calvadosia_cruxmelitensis_CRGBv3.2"
```

**Genome**

```bash
zcat Ccrux.v3.greater200.fasta.gz \
  | sed -e "s/>Ccrux./>${PREFIX}___/" \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.assembly.fasta
```

**CDS**

```bash
zcat Ohdera_et_al_2018/Condingseq/Ccrux.augustus.codingseq.gz \
  | sed -e "s/>Ccrux./>${PREFIX}___/" \
  | seqkit rename \
  | reformat.sh fastawrap=0 iupacToN=f ignorejunk=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.genes.cds.fna
```

**PEP**

```bash
zcat Ohdera_et_al_2018/AA_Files/Ccrux.augustus.3.3.faa.gz \
  | sed -e "s/>/>${PREFIX}___/" \
  | seqkit rename \
  | reformat.sh fastawrap=0 iupacToN=f ignorejunk=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  | sed -e 's/*$//' \
  > ${PREFIX}.genes.pep.faa
```

**GFF**

No CDS or GFF files were provided.

## *Clytia hemispaerica* NANAv1

> Hydrozoa

From: Ondo Fishing Port, Hiroshima Prefecture, Japan

Version: 1.0

Paper: Unpublished

Website: http://marimba.obs-vlfr.fr/downloads/

```bash
wget http://marimba.obs-vlfr.fr/sites/default/files/download/clytia_hm2.fasta
wget http://marimba.obs-vlfr.fr/sites/default/files/download/merged_transcript_models.gff3
wget http://marimba.obs-vlfr.fr/sites/default/files/download/full_nr_align.fasta
wget http://marimba.obs-vlfr.fr/sites/default/files/download/transcripts.fa
```

```bash
PREFIX="Clytia_hemispaerica_NANAv1"
```

**Genome**

```bash
cat clytia_hm2.fasta \
  | sed -e "s/>/>${PREFIX}___/" \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.assembly.fasta
```

**PEP**

```bash
cat full_nr_align.fasta \
  | sed -e "s/>/>${PREFIX}___/" -e 's/-protein//' \
  | reformat.sh fastawrap=0 iupacToN=f ignorejunk=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.genes.pep.faa
```

**CDS**

```bash
cat transcripts.fa \
  | sed -e "s/>/>${PREFIX}___/" \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  | ~/scripts/grepf_fasta.py -f <(grep '>' ${PREFIX}.genes.pep.faa | sed -e 's/>//') \
  > ${PREFIX}.genes.cds.fna
```

**GFF**

```bash
cat merged_transcript_models.gff3 \
  | awk '$1!~"^#"' \
  | sed -e "s/^/${PREFIX}___/" \
  | sed -e "s/ID=/ID=${PREFIX}___/" \
  | sed -e "s/Parent=/Parent=${PREFIX}___/" \
  | sed -e "s/Name=/Name=${PREFIX}___/" \
  | awk -F'\t' '$3=="transcript" || $3=="exon"' \
  | gffread \
  | awk '$1!~"^#"' \
  | awk 'BEGIN{OFS=FS="\t"} { if($3=="mRNA") {$3="transcript"; print} else {print} }' \
  | sed -e 's/.*=\([^=]*\)/\0\t\1/' \
  | ~/scripts/grepf_column.py -c 10 -f <(grep '>' ${PREFIX}.genes.pep.faa | sed -e 's/>//') \
  | cut -f1-9 \
  > ${PREFIX}.genes.gff3
```

## *Dendronephthya gigantea* SJKRv1

> Soft coral

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

## *Discosoma* sp. MSMEv1

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

## *Ephydatia muelleri* BCCAv1

> Sponge

From: Sooke Reservoir, Victoria, British Columbia

Version: 1.0

Paper: https://doi.org/10.1038/s41467-020-17397-w

Website: https://spaces.facsci.ualberta.ca/ephybase/

```bash
wget https://bitbucket.org/EphydatiaGenome/ephydatiagenome/downloads/Emu_genome_v1.fa.gz
wget https://bitbucket.org/EphydatiaGenome/ephydatiagenome/downloads/Emu_genome_v1.gff.gz
wget https://bitbucket.org/EphydatiaGenome/ephydatiagenome/downloads/Emu_v1_prots.fasta.gz
wget https://bitbucket.org/EphydatiaGenome/ephydatiagenome/downloads/Emu_genome_v1.codingseq.gz
```

```bash
PREFIX="Ephydatia_muelleri_BCCAv1"
```

**Genome**

```bash
zcat Emu_genome_v1.fa.gz \
  | sed -e "s/>/>${PREFIX}___/" \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.assembly.fasta
```

**CDS**

```bash
zcat Emu_genome_v1.codingseq.gz \
  | sed -e "s/>/>${PREFIX}___/" \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.genes.cds.fna
```

**PEP**

```bash
zcat Emu_v1_prots.fasta.gz \
  | sed -e "s/>/>${PREFIX}___/" \
  | awk '{if($1~"^>"){print $0".t1"}else{print}}' \
  | reformat.sh fastawrap=0 iupacToN=f ignorejunk=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  | ~/scripts/grepf_fasta.py -f <(grep '>' ${PREFIX}.genes.cds.fna | sed -e 's/>//') \
  > ${PREFIX}.genes.pep.faa
```

**GFF**

```bash
zcat Emu_genome_v1.gff.gz \
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

## *Exaiptasia diaphana* (*Aiptasia* strain CC7) CC7v1

> Sea anemone

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







## *Hoilungia hongkongensis* HCHKv1

> Placozoan

From: Ho Chung  River close to a small mangrove at Heung Chung village, Hong Kong (22.352728N  114.251733E), on June 6, 2012

Version: 1.0

Paper: https://doi.org/10.1371/journal.pbio.2005359

Website: https://bitbucket.org/molpalmuc/hoilungia-genome/src/master/

```bash
git clone https://bitbucket.org/molpalmuc/hoilungia-genome.git
```

```bash
PREFIX="Hoilungia_hongkongensis_HCHKv1"
```

**Genome**

```bash
zcat hoilungia-genome/sequences/Hhon_final_contigs_unmasked.fasta.gz \
  | sed -e "s/>/>${PREFIX}___/" \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.assembly.fasta
```

**Get longest ORF per gene**

Ignore CDS from pseudogenes. Only 3 pseudo genes in dataset.

```bash
zcat hoilungia-genome/tracks/Hhon_BRAKER1_CDS.gff3.gz \
  | awk -F'\t' '$3=="CDS"' \
  | sed -e "s/[^\t]*ID=\([^,;]*\).*/${PREFIX}___\1\t${PREFIX}___\1/" \
  | sed -e 's/\.t[0-9]*$//' \
  > ncbi_dataset.CDS_features.tsv
  
cat ncbi_dataset.CDS_features.tsv \
  | awk -F'\t' '$3=="CDS" {
      CDS[$9]=CDS[$9]+(($5-$4)+1)
      GENE[$9]=$10
    } END {
      for (i in CDS){
        print GENE[i]"\t"i"\t"CDS[i]
    } }' \
  > ncbi_dataset.CDS_lengths.tsv

cat ncbi_dataset.CDS_lengths.tsv \
  | sort -k1,1 -k3,3nr \
  | awk '!seen[$1]++{print $2}' \
  > ncbi_dataset.longest.txt
```

**CDS**

```bash
zcat hoilungia-genome/sequences/Hhon_BRAKER1_CDS.fasta.gz \
  | sed -e "s/>/>${PREFIX}___/" \
  | ~/scripts/grepf_fasta.py -f ncbi_dataset.longest.txt \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.genes.cds.fna
```

**PEP**

```bash
zcat hoilungia-genome/sequences/Hhon_BRAKER1_proteins.fasta.gz \
  | sed -e "s/>/>${PREFIX}___/" \
  | ~/scripts/grepf_fasta.py -f ncbi_dataset.longest.txt \
  | reformat.sh fastawrap=0 iupacToN=f ignorejunk=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.genes.pep.faa
```

**GFF**

```bash
zcat hoilungia-genome/tracks/Hhon_BRAKER1_CDS.gff3.gz \
  | awk -F'\t' '$3=="CDS"' \
  | sed -e "s/[^\t]*Name=\([^;]*\).*/Parent=${PREFIX}___\1\t${PREFIX}___\1/" \
  | sed -e "s/^/${PREFIX}___/" \
  | ~/scripts/grepf_column.py -c 10 -f ncbi_dataset.longest.txt \
  | cut -f 1-9 \
  | gffread \
  | awk '$1!~"^#"' \
  | awk 'BEGIN{OFS=FS="\t"} { if($3!="CDS"){print} else {print; $3="exon"; print} }' \
  > ${PREFIX}.genes.gff3
```

## *Hormiphora californensis* MBUSv1

> Hydrozoa

From: Monterey Bay, California

Version: 1.0

Paper: https://doi.org/10.1093%2Fg3journal%2Fjkab302

Website: https://github.com/conchoecia/hormiphora

```bash
wget https://github.com/conchoecia/hormiphora/blob/master/annotation/raw_files/UCSC_Hcal_v1.fa.gz
wget https://github.com/conchoecia/hormiphora/releases/download/Hcv1.av93zen/Hcv1av93_release.tar.gz
```

```bash
PREFIX="Hormiphora_californensis_MBUSv1"
```

NOTE: 

- Need to remove the Mitochondrial scaffold from the assembly.
- Don't use phased assembly/genes. Just use the full un-phased datasets.

**Genome**

```bash
zcat UCSC_Hcal_v1.fa.gz \
  | sed -e "s/>/>${PREFIX}___/" \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  | ~/scripts/grepf_fasta.py -v -f <(echo "Hormiphora_californensis_MBUSv1___M") \
  > ${PREFIX}.assembly.fasta
```

**Get longest ORF per gene**

```bash
zcat Hcv1av93_release/Hcv1av93.gff.gz \
  | awk -F'\t' '$1!~"#" && $3=="transcript" && $1!="M"' \
  | sed -e "s/.*ID=\([^,;]*\).*Parent=\([^,;]*\).*/${PREFIX}___\1\t${PREFIX}___\2/" \
  | sed -e 's/Hcv1.av93.//g' \
  > dataset.mRNA2gene.tsv
  
zcat Hcv1av93_release/Hcv1av93.gff.gz \
  | awk -F'\t' '$1!~"#" && $3=="exon" {
      split($9,a,"Parent=");
      split(a[2],b,";");
      ID=b[1]
      CDS[ID]=CDS[ID]+(($5-$4)+1)
    } END {
      for (i in CDS){
        print i"\t"CDS[i]
    } }' \
  | sed -e "s/^/${PREFIX}___/" -e 's/Hcv1.av93.//' \
  > dataset.CDS_length_per_mRNA.tsv

~/scripts/add_value_to_table_SQLite3.py -d "NA" \
    -i dataset.mRNA2gene.tsv \
    -a dataset.CDS_length_per_mRNA.tsv \
  | sort -k2,2 -k3,3nr \
  | awk -F'\t' '!seen[$2]++{print $1}' \
  > ncbi_dataset.longest.txt
```

**PEP**

NOTE: A number of mRNA don't have a protein sequences (possible pseudo genes?), so have to use the PEP file to filer CDS and GFF downstream.

```bash
zcat Hcv1av93_release/Hcv1av93_model_proteins.pep.gz \
  | sed -e "s/>/>${PREFIX}___/" -e 's/Hcv1.av93.//' \
  | reformat.sh fastawrap=0 iupacToN=f ignorejunk=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  | ~/scripts/grepf_fasta.py -f ncbi_dataset.longest.txt \
  > ${PREFIX}.genes.pep.faa
```

**CDS**

```bash
zcat Hcv1av93_release/Hcv1av93_transcripts.fasta.gz \
  | sed -e "s/>/>${PREFIX}___/" -e 's/Hcv1.av93.//' \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  | ~/scripts/grepf_fasta.py -f <(grep '>' ${PREFIX}.genes.pep.faa | sed -e 's/>//') \
  > ${PREFIX}.genes.cds.fna
```

**GFF**

```bash
zcat Hcv1av93_release/Hcv1av93.gff.gz \
  | awk '$1!~"^#"' \
  | sed -e "s/^/${PREFIX}___/" \
  | sed -e "s/ID=/ID=${PREFIX}___/" \
  | sed -e "s/Parent=/Parent=${PREFIX}___/" \
  | sed -e "s/Name=/Name=${PREFIX}___/" \
  | sed -e 's/Hcv1.av93.//g' \
  | awk -F'\t' '$3=="mRNA" || $3=="CDS" || $3=="exon"' \
  | gffread \
  | awk '$1!~"^#"' \
  | awk 'BEGIN{OFS=FS="\t"} { if($3=="mRNA") {$3="transcript"; print} else {print} }' \
  | sed -e 's/.*=\([^=]*\)/\0\t\1/' \
  | ~/scripts/grepf_column.py -c 10 -f <(grep '>' ${PREFIX}.genes.pep.faa | sed -e 's/>//') \
  | cut -f1-9 \
  > ${PREFIX}.genes.gff3
```

## *Hydra viridissima* 99AUv1

> Hydrozoa

From: Australia (isolate A99)

Version: 1.0

Paper: Unpublished

Website: https://marinegenomics.oist.jp/hydra_viridissima_a99/viewer/download?project_id=82

```bash
wget https://marinegenomics.oist.jp/hydra_viridissima_a99/download/hvir_genome_hm2_250116_renamed.fa.gz
wget https://marinegenomics.oist.jp/hydra_viridissima_a99/download/hvir_r06.all.recounted.gff.gz
wget https://marinegenomics.oist.jp/hydra_viridissima_a99/download/hvir_r06.all.recounted.aa.gz
wget https://marinegenomics.oist.jp/hydra_viridissima_a99/download/hvir_r06.all.recounted.mrna.gz
```

```bash
PREFIX="Hydra_viridissima_99AUv1"
```

**Genome**

```bash
zcat hvir_genome_hm2_250116_renamed.fa.gz \
  | sed -e "s/>/>${PREFIX}___/" \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.assembly.fasta
```

**Get longest ORF per gene**

```bash
zcat hvir_r06.all.recounted.gff.gz \
  | awk -F'\t' '$1!~"#" && $3=="transcript"' \
  | sed -e "s/.*ID=\([^,;]*\).*Parent=\([^,;]*\).*/${PREFIX}___\1\t${PREFIX}___\2/" \
  > dataset.mRNA2gene.tsv
  
zcat hvir_r06.all.recounted.gff.gz \
  | awk -F'\t' '$1!~"#" && $3=="CDS" {
      split($9,a,"Parent=");
      split(a[2],b,";");
      ID=b[1]
      CDS[ID]=CDS[ID]+(($5-$4)+1)
    } END {
      for (i in CDS){
        print i"\t"CDS[i]
    } }' \
  | sed -e "s/^/${PREFIX}___/" \
  > dataset.CDS_length_per_mRNA.tsv

~/scripts/add_value_to_table_SQLite3.py -d "NA" \
    -i dataset.mRNA2gene.tsv \
    -a dataset.CDS_length_per_mRNA.tsv \
  | sort -k2,2 -k3,3nr \
  | awk -F'\t' '!seen[$2]++{print $1}' \
  > ncbi_dataset.longest.txt
```

**PEP**

NOTE: One mRNA doesn't have a protein sequence (possible pseudo gene?), so have to use the PEP file to filer CDS and GFF downstream.

```bash
zcat hvir_r06.all.recounted.aa.gz \
  | sed -e "s/>/>${PREFIX}___/" \
  | reformat.sh fastawrap=0 iupacToN=f ignorejunk=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  | ~/scripts/grepf_fasta.py -f ncbi_dataset.longest.txt \
  > ${PREFIX}.genes.pep.faa
```

**CDS**

```bash
zcat hvir_r06.all.recounted.mrna.gz \
  | sed -e "s/>/>${PREFIX}___/" \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  | ~/scripts/grepf_fasta.py -f <(grep '>' ${PREFIX}.genes.pep.faa | sed -e 's/>//') \
  > ${PREFIX}.genes.cds.fna
```

**GFF**

```bash
zcat hvir_r06.all.recounted.gff.gz \
  | awk '$1!~"^#"' \
  | sed -e "s/^/${PREFIX}___/" \
  | sed -e "s/ID=/ID=${PREFIX}___/" \
  | sed -e "s/Parent=/Parent=${PREFIX}___/" \
  | sed -e "s/Name=/Name=${PREFIX}___/" \
  | awk -F'\t' '$3=="mRNA" || $3=="CDS" || $3=="exon"' \
  | gffread \
  | awk '$1!~"^#"' \
  | awk 'BEGIN{OFS=FS="\t"} { if($3=="mRNA") {$3="transcript"; print} else {print} }' \
  | sed -e 's/.*=\([^=]*\)/\0\t\1/' \
  | ~/scripts/grepf_column.py -c 10 -f <(grep '>' ${PREFIX}.genes.pep.faa | sed -e 's/>//') \
  | cut -f1-9 \
  > ${PREFIX}.genes.gff3
```

## *Hydra vulgaris* MIJPv3

> Hydrozoa

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

## *Hydractinia echinata* NWAOv1

> Hydrozoa

From: Northwestern Atlantic Ocean

Version: 1.0

Paper: Unpublished

Website: https://research.nhgri.nih.gov/hydractinia/

```bash
wget https://research.nhgri.nih.gov/hydractinia/download/assembly/echinata/Hech_primary_v1.0.fa.gz
wget https://research.nhgri.nih.gov/hydractinia/download/genemodels_gff3/echinata/Hech_primary_v1.0.gff3.gz
wget https://research.nhgri.nih.gov/hydractinia/download/protein_models/echinata/Hech_primary_v1.0.aa.gz
wget https://research.nhgri.nih.gov/hydractinia/download/gene_models/echinata/Hech_primary_v1.0.nt.gz
```

```bash
PREFIX="Hydractinia_echinata_NWAOv1"
```

**Genome**

```bash
zcat Hech_primary_v1.0.fa.gz \
  | sed -e "s/>HyE/>${PREFIX}___/" \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.assembly.fasta
```

**GFF**

```bash
zcat Hech_primary_v1.0.gff3.gz \
  | awk '$1!~"^#"' \
  | sed -e "s/^HyE/${PREFIX}___/" \
  | sed -e "s/ID=HyE/ID=${PREFIX}___/" \
  | sed -e "s/Parent=HyE/Parent=${PREFIX}___/" \
  | sed -e "s/Name=HyE/Name=${PREFIX}___/" \
  | awk -F'\t' '$3=="mRNA" || $3=="CDS" || $3=="exon"' \
  | gffread \
  | awk '$1!~"^#"' \
  | awk 'BEGIN{OFS=FS="\t"} { if($3=="mRNA") {$3="transcript"; print} else {print} }' \
  > ${PREFIX}.genes.gff3
```

Extract sequence names from `gff` file (more sequences in CDS and PEP files than GFF).

```bash
awk -F'\t' '$3=="transcript"{print $9}' ${PREFIX}.genes.gff3 \
  | sed -e 's/ID=//' \
  > ${PREFIX}.genes.names.txt
```

**PEP**

```bash
zcat Hech_primary_v1.0.aa.gz \
  | sed -e "s/>HyE/>${PREFIX}___/" \
  | reformat.sh fastawrap=0 iupacToN=f ignorejunk=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  | ~/scripts/grepf_fasta.py -f ${PREFIX}.genes.names.txt \
  > ${PREFIX}.genes.pep.faa
```

**CDS**

```bash
zcat Hech_primary_v1.0.nt.gz \
  | sed -e "s/>HyE/>${PREFIX}___/" \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  | ~/scripts/grepf_fasta.py -f ${PREFIX}.genes.names.txt \
  > ${PREFIX}.genes.cds.fna
```

## *Hydractinia symbiolongicarpus* NEAOv1

> Hydrozoa

From: Northeastern Atlantic Ocean

Version: 1.0

Paper: Unpublished

Website: https://research.nhgri.nih.gov/hydractinia/

```bash
wget https://research.nhgri.nih.gov/hydractinia/download/assembly/symbio/Hsym_primary_v1.0.fa.gz
wget https://research.nhgri.nih.gov/hydractinia/download/genemodels_gff3/symbio/Hsym_primary_v1.0.gff3.gz
wget https://research.nhgri.nih.gov/hydractinia/download/protein_models/symbio/Hsym_primary_v1.0.aa.gz
wget https://research.nhgri.nih.gov/hydractinia/download/gene_models/symbio/Hsym_primary_v1.0.nt.gz
```

```bash
PREFIX="Hydractinia_symbiolongicarpus_NEAOv1"
```

**Genome**

```bash
zcat Hsym_primary_v1.0.fa.gz \
  | sed -e "s/>HyS/>${PREFIX}___/" \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.assembly.fasta
```

**GFF**

```bash
zcat Hsym_primary_v1.0.gff3.gz \
  | awk '$1!~"^#"' \
  | sed -e "s/^HyS/${PREFIX}___/" \
  | sed -e "s/ID=HyS/ID=${PREFIX}___/" \
  | sed -e "s/Parent=HyS/Parent=${PREFIX}___/" \
  | sed -e "s/Name=HyS/Name=${PREFIX}___/" \
  | awk -F'\t' '$3=="mRNA" || $3=="CDS" || $3=="exon"' \
  | gffread \
  | awk '$1!~"^#"' \
  | awk 'BEGIN{OFS=FS="\t"} { if($3=="mRNA") {$3="transcript"; print} else {print} }' \
  > ${PREFIX}.genes.gff3
```

Extract sequence names from `gff` file (more sequences in CDS and PEP files than GFF).

```bash
awk -F'\t' '$3=="transcript"{print $9}' ${PREFIX}.genes.gff3 \
  | sed -e 's/ID=//' \
  > ${PREFIX}.genes.names.txt
```

**PEP**

```bash
zcat Hsym_primary_v1.0.aa.gz \
  | sed -e "s/>HyS/>${PREFIX}___/" \
  | reformat.sh fastawrap=0 iupacToN=f ignorejunk=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  | ~/scripts/grepf_fasta.py -f ${PREFIX}.genes.names.txt \
  > ${PREFIX}.genes.pep.faa
```

**CDS**

```bash
zcat Hsym_primary_v1.0.nt.gz \
  | sed -e "s/>HyS/>${PREFIX}___/" \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  | ~/scripts/grepf_fasta.py -f ${PREFIX}.genes.names.txt \
  > ${PREFIX}.genes.cds.fna
```

## *Mnemiopsis leidyi* WHUSv2.2

> Ctenophore

From: Vineyard Sound near Woods Hole, Massachusetts, USA

Version: 2.2

Paper: https://doi.org/10.1126%2Fscience.1242592

Website: https://research.nhgri.nih.gov/mnemiopsis/

```bash
wget https://research.nhgri.nih.gov/mnemiopsis/download/genome/MlScaffold09.nt.gz
wget https://research.nhgri.nih.gov/mnemiopsis/download/proteome/ML2.2.gff3.gz
wget https://research.nhgri.nih.gov/mnemiopsis/download/proteome/ML2.2.aa.gz
wget https://research.nhgri.nih.gov/mnemiopsis/download/transcriptome/ML2.2.nt.gz
```

```bash
PREFIX="Mnemiopsis_leidyi_WHUSv2.2"
```

**Genome**

```bash
zcat MlScaffold09.nt.gz \
  | sed -e "s/>/>${PREFIX}___/" \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.assembly.fasta
```

**CDS**

```bash
zcat ML2.2.nt.gz \
  | sed -e "s/>/>${PREFIX}___/" \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.genes.cds.fna
```

**PEP**

```bash
zcat ML2.2.aa.gz \
  | sed -e "s/>/>${PREFIX}___/" \
  | reformat.sh fastawrap=0 iupacToN=f ignorejunk=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.genes.pep.faa
```

**GFF**

```bash
zcat ML2.2.gff3.gz \
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

## *Monosiga brevicollis* CCBMv1

> Choanoflagellate

From: Seawater from Church Cave, Bemuda

Version: 1.0

Paper:  https://doi.org/10.1038/nature06617

Website: NCBI: GCF_000002865.3

```bash
# tar from NCBI
unzip GCF_000002865.3.zip
```

```bash
PREFIX="Monosiga_brevicollis_CCBMv1"
```

**Genome**

```bash
cat ncbi_dataset/data/GCF_000002865.3/GCF_000002865.3_V1.0_genomic.fna \
  | sed -e "s/>/>${PREFIX}___/" \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.assembly.fasta
```

**Get longest ORF per gene**

Ignore CDS from pseudogenes. Only 3 pseudo genes in dataset.

```bash
cat ncbi_dataset/data/GCF_000002865.3/genomic.gff \
  | awk -F'\t' '$3=="CDS"' \
  | grep -v 'pseudo=true' \
  | sed -e "s/[^\t]*GeneID:\([^,;]*\).*Name=\([^;]*\).*/${PREFIX}___\2\t\1/" \
  > ncbi_dataset.CDS_features.tsv
  
cat ncbi_dataset.CDS_features.tsv \
  | awk -F'\t' '$3=="CDS" {
      CDS[$9]=CDS[$9]+(($5-$4)+1)
      GENE[$9]=$10
    } END {
      for (i in CDS){
        print GENE[i]"\t"i"\t"CDS[i]
    } }' \
  > ncbi_dataset.CDS_lengths.tsv

cat ncbi_dataset.CDS_lengths.tsv \
  | sort -k1,1 -k3,3nr \
  | awk '!seen[$1]++{print $2}' \
  > ncbi_dataset.longest.txt
```

**CDS**

```bash
cat ncbi_dataset/data/GCF_000002865.3/cds_from_genomic.fna \
  | sed -e "s/.*protein_id=\([^]]*\)\].*/>${PREFIX}___\1/" \
  | ~/scripts/grepf_fasta.py -f ncbi_dataset.longest.txt \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.genes.cds.fna
```

**PEP**

```bash
cat ncbi_dataset/data/GCF_000002865.3/protein.faa \
  | sed -e "s/>/>${PREFIX}___/" \
  | ~/scripts/grepf_fasta.py -f ncbi_dataset.longest.txt \
  | reformat.sh fastawrap=0 iupacToN=f ignorejunk=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.genes.pep.faa
```

**GFF**

```bash
cat ncbi_dataset/data/GCF_000002865.3/genomic.gff \
  | awk -F'\t' '$3=="CDS"' \
  | grep -v 'pseudo=true' \
  | sed -e "s/[^\t]*Name=\([^;]*\).*/Parent=${PREFIX}___\1\t${PREFIX}___\1/" \
  | sed -e "s/^/${PREFIX}___/" \
  | ~/scripts/grepf_column.py -c 10 -f ncbi_dataset.longest.txt \
  | cut -f 1-9 \
  | gffread \
  | awk '$1!~"^#"' \
  | awk 'BEGIN{OFS=FS="\t"} { if($3!="CDS"){print} else {print; $3="exon"; print} }' \
  > ${PREFIX}.genes.gff3
```

## *Morbakka virulenta* ODJPv1

> Jellyfish

From: Ondo Fishing Port, Hiroshima Prefecture, Japan

Version: 1.0

Paper: https://doi.org/10.1038/s41559-019-0853-y

Website: https://marinegenomics.oist.jp/morbakka_virulenta/viewer/download?project_id=70

```bash
wget https://marinegenomics.oist.jp/morbakka_virulenta/download/MOR05_r06_genome.fa.gz
wget https://marinegenomics.oist.jp/morbakka_virulenta/download/MOR05_r06_wref.gff3.gz
wget https://marinegenomics.oist.jp/morbakka_virulenta/download/MOR05_r06_proteins.fa.gz
wget https://marinegenomics.oist.jp/morbakka_virulenta/download/MOR05_r06_mRNA.fa.gz
```

```bash
PREFIX="Morbakka_virulenta_ODJPv1"
```

**Genome**

```bash
zcat MOR05_r06_genome.fa.gz \
  | sed -e "s/>/>${PREFIX}___/" \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.assembly.fasta
```

**Get longest ORF per gene**

```bash
zcat MOR05_r06_wref.gff3.gz \
  | awk -F'\t' '$1!~"#" && $3=="transcript"' \
  | sed -e "s/.*ID=\([^,;]*\).*Parent=\([^,;]*\).*/${PREFIX}___\1\t${PREFIX}___\2/" \
  > dataset.mRNA2gene.tsv
  
zcat MOR05_r06_wref.gff3.gz \
  | awk -F'\t' '$1!~"#" && $3=="CDS" {
      split($9,a,"Parent=");
      split(a[2],b,";");
      ID=b[1]
      CDS[ID]=CDS[ID]+(($5-$4)+1)
    } END {
      for (i in CDS){
        print i"\t"CDS[i]
    } }' \
  | sed -e "s/^/${PREFIX}___/" \
  > dataset.CDS_length_per_mRNA.tsv

~/scripts/add_value_to_table_SQLite3.py -d "NA" \
    -i dataset.mRNA2gene.tsv \
    -a dataset.CDS_length_per_mRNA.tsv \
  | sort -k2,2 -k3,3nr \
  | awk -F'\t' '!seen[$2]++{print $1}' \
  > ncbi_dataset.longest.txt
```

**PEP**

```bash
zcat MOR05_r06_proteins.fa.gz \
  | sed -e "s/>/>${PREFIX}___/" \
  | reformat.sh fastawrap=0 iupacToN=f ignorejunk=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  | ~/scripts/grepf_fasta.py -f ncbi_dataset.longest.txt \
  > ${PREFIX}.genes.pep.faa
```

**CDS**

```bash
zcat MOR05_r06_mRNA.fa.gz \
  | sed -e "s/>/>${PREFIX}___/" \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  | ~/scripts/grepf_fasta.py -f ncbi_dataset.longest.txt \
  > ${PREFIX}.genes.cds.fna
```

**GFF**

```bash
zcat MOR05_r06_wref.gff3.gz \
  | awk '$1!~"^#"' \
  | sed -e "s/^/${PREFIX}___/" \
  | sed -e "s/ID=/ID=${PREFIX}___/" \
  | sed -e "s/Parent=/Parent=${PREFIX}___/" \
  | sed -e "s/Name=/Name=${PREFIX}___/" \
  | awk -F'\t' '$3=="mRNA" || $3=="CDS" || $3=="exon"' \
  | gffread \
  | awk '$1!~"^#"' \
  | awk 'BEGIN{OFS=FS="\t"} { if($3=="mRNA") {$3="transcript"; print} else {print} }' \
  | sed -e 's/.*=\([^=]*\)/\0\t\1/' \
  | ~/scripts/grepf_column.py -c 10 -f ncbi_dataset.longest.txt \
  | cut -f1-9 \
  > ${PREFIX}.genes.gff3
```

## *Nematostella vectensis* RRUSv1

> Sea anemone

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

## *Paramuricea clavata* CTESv1

> Soft coral

From: Cova de la Vaca (42°2’52.97’’N; 3°13’34.76’’E) in The Montgrí, Medes Islands and Baix Ter Natural Park (Catalunya, Spain)

Version: 1.0

Paper: https://doi.org/10.1534/g3.120.401371

Website: https://denovo.cnag.cat/pclavata

```bash
# Downloaded manually from site
```

```bash
PREFIX="Paramuricea_clavata_CTESv1"
```

**Genome**

```bash
zcat pcla8.scaffolds.fa.gz \
  | sed -e "s/>pcla8_/>${PREFIX}___/" \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.assembly.fasta
```

**PEP**

```bash
zcat Pacla8A.longest_peptide.fa.gz \
  | sed -e "s/>Pacla8A/>${PREFIX}___/" \
  | reformat.sh fastawrap=0 iupacToN=f ignorejunk=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  | sed -e 's/*$//' \
  > ${PREFIX}.genes.pep.faa
```

**CDS**

```bash
zcat Pacla8A.cds.fa.gz  \
  | sed -e "s/>Pacla8A/>${PREFIX}___/" \
  | ~/scripts/grepf_fasta.py -f <(grep '>' ${PREFIX}.genes.pep.faa | sed -e 's/>//') \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.genes.cds.fna
```

**GFF**

```bash
zcat Pacla8A.gff3.gz \
  | awk '$1!~"^#"' \
  | sed -e "s/^pcla8_/${PREFIX}___/" \
  | sed -e "s/ID=Pacla8A/ID=${PREFIX}___/" \
  | sed -e "s/Parent=Pacla8A/Parent=${PREFIX}___/" \
  | sed -e "s/Name=Pacla8A/Name=${PREFIX}___/" \
  | awk -F'\t' '$3=="CDS" || $3=="exon"' \
  | dos2unix \
  | sed -e 's/.*Parent=\([^,;]*\).*/\0\t\1/' \
  | sed -e 's/T\([0-9]\)/P\1/g' \
  | ~/scripts/grepf_column.py -c 10 -f <(grep '>' ${PREFIX}.genes.pep.faa | sed -e 's/>//') \
  | cut -f1-9 \
  | gffread \
  | awk '$1!~"^#"' \
  | awk 'BEGIN{OFS=FS="\t"} { if($3=="mRNA") {$3="transcript"; print} else {print} }' \
  > ${PREFIX}.genes.gff3
```

## *Renilla muelleri* FLUSv1

> Sea pansy

From: Gulf Specimen Marine Lab (Panacea, FL, USA), which collects specimens off the panhandle of Florida in the Gulf of Mexico

Version: 1.0

Paper: https://doi.org/10.1093/gigascience/giz026

Website: http://rmue.reefgenomics.org/

```bash
wget http://rmue.reefgenomics.org/download/final_renilla_genome.fa.gz
wget http://rmue.reefgenomics.org/download/renilla_gene_model.gff.gz
wget http://rmue.reefgenomics.org/download/renilla_predicted_proteins.fa.gz
```

```bash
PREFIX="Renilla_muelleri_FLUSv1"
```

**Genome**

```bash
zcat final_renilla_genome.fa.gz \
  | sed -e "s/>/>${PREFIX}___/" \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.assembly.fasta
```

**PEP**

```bash
zcat renilla_predicted_proteins.fa.gz \
  | sed -e "s/>/>${PREFIX}___/" \
  | reformat.sh fastawrap=0 iupacToN=f ignorejunk=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.genes.pep.faa
```

**GFF**

```bash
zcat renilla_gene_model.gff.gz \
  | sed -e "s/^/${PREFIX}___/" \
        -e "s/transcript_id \"/transcript_id \"${PREFIX}___/g" \
        -e "s/gene_id \"/gene_id \"${PREFIX}___/g" \
  | awk 'BEGIN{OFS=FS="\t"} { if($3=="CDS"){print; $3="exon"; print} }' \
  | gffread \
  | awk '$1!~"^#"' \
  | sed -e 's/;geneID=.*//' \
  > ${PREFIX}.genes.gff3

gffread -S \
  -g ${PREFIX}.assembly.fasta \
  -x ${PREFIX}.genes.cds.fna.tmp \
  -y ${PREFIX}.genes.pep.faa.tmp \
  ${PREFIX}.genes.gff3
rm ${PREFIX}.genes.pep.faa.tmp
```

**CDS**

```bash
cat ${PREFIX}.genes.cds.fna.tmp \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.genes.cds.fna
rm ${PREFIX}.genes.cds.fna.tmp
```

## *Renilla reniformis* FLUSv1

> Sea pansy

From: Fort George Inlet, Jacksonville, Florida, USA

Version: 1.0

Paper: https://doi.org/10.1186/s12862-018-1142-0

Website: http://ryanlab.whitney.ufl.edu/genomes/Renilla_reniformis/

```bash
wget http://ryanlab.whitney.ufl.edu/genomes/Renilla_reniformis/downloads/Renilla_reniformis.v1.fa.gz
wget http://ryanlab.whitney.ufl.edu/genomes/Renilla_reniformis/downloads/Reni_reni.v1.gff.gz
wget http://ryanlab.whitney.ufl.edu/genomes/Renilla_reniformis/downloads/Reni_reni.v1.aa.gz
wget http://ryanlab.whitney.ufl.edu/genomes/Renilla_reniformis/downloads/Reni_reni.v1.cds.gz
```

```bash
PREFIX="Renilla_reniformis_FLUSv1"
```

**Genome**

```bash
zcat Renilla_reniformis.v1.fa.gz \
  | sed -e "s/>Renilla_reniformis\.v1\.[^1-9]*/>${PREFIX}___/" \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.assembly.fasta
```

**PEP**

```bash
zcat Reni_reni.v1.aa.gz \
  | sed -e "s/>/>${PREFIX}___/" \
  | reformat.sh fastawrap=0 iupacToN=f ignorejunk=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.genes.pep.faa
```

**CDS**

```bash
zcat Reni_reni.v1.cds.gz \
  | sed -e "s/>/>${PREFIX}___/" \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.genes.cds.fna
```

**GFF**

```bash
zcat Reni_reni.v1.gff.gz \
  | awk '$1!~"^#"' \
  | awk -F'\t' '$3=="transcript" || $3=="CDS"' \
  | awk -F'\t' -vP="${PREFIX}___" 'BEGIN{OFS=FS="\t"}{if($3=="transcript"){$9="ID="P$9; print}else{print}}' \
  | sed -e 's/"//g' \
        -e "s/transcript_id /ID=${PREFIX}___/" \
        -e "s/ gene_id /Parent=${PREFIX}___/" \
  | gffread \
  | awk '$1!~"^#"' \
  | awk 'BEGIN{OFS=FS="\t"} { if($3=="mRNA") {$3="transcript"; print} else {print} }' \
  | sed -e "s/^scaffold[^1-9]*/${PREFIX}___/" -e 's/_cov[^\t]*//' \
  > ${PREFIX}.genes.gff3
```

## *Salpingoeca rosetta* HIUSv1

> Choanoflagellate

From: Hog Island, Virginia, USA

Version: 1.0

Paper: https://doi.org/10.1186/gb-2013-14-2-r15

Website: NCBI:  GCF_000188695.1

```bash
# Download from NCBI
unzip GCF_000188695.1.zip
```

```bash
PREFIX="Salpingoeca_rosetta_HIUSv1"
```

**Genome**

```bash
cat ncbi_dataset/data/GCF_000188695.1/GCF_000188695.1_Proterospongia_sp_ATCC50818_genomic.fna \
  | sed -e "s/>/>${PREFIX}___/" \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.assembly.fasta
```

**Get longest ORF per gene**

Ignore CDS from pseudogenes. Only 3 pseudo genes in dataset.

```bash
cat ncbi_dataset/data/GCF_000188695.1/genomic.gff \
  | awk -F'\t' '$3=="CDS"' \
  | sed -e "s/[^\t]*GeneID:\([^,;]*\).*Name=\([^;]*\).*/${PREFIX}___\2\t\1/" \
  > ncbi_dataset.CDS_features.tsv
  
cat ncbi_dataset.CDS_features.tsv \
  | awk -F'\t' '$3=="CDS" {
      CDS[$9]=CDS[$9]+(($5-$4)+1)
      GENE[$9]=$10
    } END {
      for (i in CDS){
        print GENE[i]"\t"i"\t"CDS[i]
    } }' \
  > ncbi_dataset.CDS_lengths.tsv

cat ncbi_dataset.CDS_lengths.tsv \
  | sort -k1,1 -k3,3nr \
  | awk '!seen[$1]++{print $2}' \
  > ncbi_dataset.longest.txt
```

**CDS**

```bash
cat ncbi_dataset/data/GCF_000188695.1/cds_from_genomic.fna \
  | sed -e "s/.*protein_id=\([^]]*\)\].*/>${PREFIX}___\1/" \
  | ~/scripts/grepf_fasta.py -f ncbi_dataset.longest.txt \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.genes.cds.fna
```

**PEP**

```bash
cat ncbi_dataset/data/GCF_000188695.1/protein.faa \
  | sed -e "s/>/>${PREFIX}___/" \
  | ~/scripts/grepf_fasta.py -f ncbi_dataset.longest.txt \
  | reformat.sh fastawrap=0 iupacToN=f ignorejunk=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.genes.pep.faa
```

**GFF**

```bash
cat ncbi_dataset/data/GCF_000188695.1/genomic.gff \
  | awk -F'\t' '$3=="CDS"' \
  | sed -e "s/[^\t]*Name=\([^;]*\).*/Parent=${PREFIX}___\1\t${PREFIX}___\1/" \
  | sed -e "s/^/${PREFIX}___/" \
  | ~/scripts/grepf_column.py -c 10 -f ncbi_dataset.longest.txt \
  | cut -f 1-9 \
  | gffread \
  | awk '$1!~"^#"' \
  | awk 'BEGIN{OFS=FS="\t"} { if($3!="CDS"){print} else {print; $3="exon"; print} }' \
  > ${PREFIX}.genes.gff3
```

## *Scolanthus callimorphus* CRFRv1

> Sea anemone

From: lle Callot, Carantec, France

Version: 1.0

Paper: https://doi.org/10.1101/2020.10.30.359448

Website: https://simrbase.stowers.org/wormanemone

```bash
wget --no-check-certificate https://simrbase.stowers.org/files/pub/nematostella/Scal/genomes/Scal100/Scal100.fasta
wget --no-check-certificate https://simrbase.stowers.org/files/pub/seaanemone/Scal/genomes/Scal100/aligned/NY_Scal100_v1/NY_Scal100_v1.20200813.gff
wget --no-check-certificate https://simrbase.stowers.org/files/pub/nematostella/Scal/genomes/Scal100/aligned/NY_Scal100_v1/NY_Scal100_v1.20200813.proteins.fasta
wget --no-check-certificate https://simrbase.stowers.org/files/pub/nematostella/Scal/genomes/Scal100/aligned/NY_Scal100_v1/NY_Scal100_v1.20200813.transcript.fasta
```

```bash
PREFIX="Scolanthus_callimorphus_CRFRv1"
```

**Genome**

```bash
cat Scal100.fasta \
  | sed -e "s/>/>${PREFIX}___/" \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.assembly.fasta
```

**CDS**

```bash
cat NY_Scal100_v1.20200813.transcript.fasta \
  | sed -e "s/>NY_Scal100_v1./>${PREFIX}___/" \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.genes.cds.fna
```

**PEP**

```bash
cat NY_Scal100_v1.20200813.proteins.fasta \
  | sed -e "s/>NY_Scal100_v1./>${PREFIX}___/" \
  | reformat.sh fastawrap=0 iupacToN=f ignorejunk=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.genes.pep.faa
```

**GFF**

```bash
cat NY_Scal100_v1.20200813.gff \
  | awk '$1!~"^#"' \
  | sed -e "s/^/${PREFIX}___/" \
  | sed -e "s/ID=NY_Scal100_v1./ID=${PREFIX}___/" \
  | sed -e "s/Parent=NY_Scal100_v1./Parent=${PREFIX}___/" \
  | sed -e "s/Name=NY_Scal100_v1./Name=${PREFIX}___/" \
  | awk -F'\t' '$3=="CDS" || $3=="exon"' \
  | dos2unix \
  | sed -e 's/.*Parent=\([^,;]*\).*/\0\t\1/' \
  | gffread \
  | awk '$1!~"^#"' \
  | awk 'BEGIN{OFS=FS="\t"} { if($3=="mRNA") {$3="transcript"; print} else {print} }' \
  > ${PREFIX}.genes.gff3
```

## *Stylissa carteri* RSSAv1

> Sponge

From: Fsar Reef (22.228408 N, 39.028187E) on the Red Sea coast of Saudi Arabia

Version: 1.0

Paper: https://doi.org/10.1186/s12864-016-2501-0

Website: http://sc.reefgenomics.org/

```bash
wget http://sc.reefgenomics.org/download/SC_genome_800.fa.gz
wget http://sc.reefgenomics.org/download/SC_final.gff.gz
wget http://sc.reefgenomics.org/download/SC_proteins_ltAED075.fa.gz
wget http://sc.reefgenomics.org/download/SC_transcripts_ltAED075.fa.gz
```

```bash
PREFIX="Stylissa_carteri_RSSAv1"
```

**Genome**

```bash
zcat SC_genome_800.fa.gz \
  | sed -e "s/>SC_/>${PREFIX}___/" \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.assembly.fasta
```

**Get longest ORF per gene**

```bash
zcat SC_final.gff.gz \
  | awk '$1!~"^#"' \
  | awk -F'\t' '$2=="maker" && $3=="mRNA"' \
  | sed -e "s/.*ID=\([^,;]*\).*Parent=\([^,;]*\).*/${PREFIX}___\2\t${PREFIX}___\1/" \
  | sort -k1,1 -k2,2 \
  > dataset.gene2mRNA.tsv

zcat SC_final.gff.gz \
  | awk '$1!~"^#"' \
  | awk -F'\t' '$2=="maker" && $3=="CDS"' \
  | sed -e "s/.*Parent=\([^,;]*\).*/\0\t${PREFIX}___\1/" \
  | awk -F'\t' '{
      CDS[$10]=CDS[$10]+(($5-$4)+1)
    } END {
      for (i in CDS){
        print i"\t"CDS[i]
    } }' \
  > dataset.CDS_lengths.tsv

~/scripts/add_value_to_table.py -c 2 -d "NA" \
  -i dataset.gene2mRNA.tsv \
  -a dataset.CDS_lengths.tsv \
  -o dataset.gene2mRNA2CDS_lengths.tsv

cat dataset.gene2mRNA2CDS_lengths.tsv \
  | ~/scripts/grepf_column.py -c 2 \
      -f <(zcat SC_proteins_ltAED075.fa.gz | grep '>' | sed -e "s/>/${PREFIX}___/" -e 's/ .*//') \
  | sort -k1,1 -k3,3nr \
  | awk '!seen[$1]++{print $2}' \
  > dataset.longest.txt
```

**GFF**

```bash
zcat SC_final.gff.gz \
  | awk '$1!~"^#"' \
  | awk -F'\t' '$2=="maker" && ($3=="exon" || $3=="CDS")' \
  | sed -e "s/^SC_/${PREFIX}___/" \
  | sed -e "s/ID=/ID=${PREFIX}___/" \
  | sed -e "s/Parent=/Parent=${PREFIX}___/" \
  | sed -e "s/Name=/Name=${PREFIX}___/" \
  | dos2unix \
  | sed -e 's/.*Parent=\([^,;]*\).*/\0\t\1/' \
  | ~/scripts/grepf_column.py -c 10 -f dataset.longest.txt \
  | cut -f1-9 \
  | gffread \
  | awk '$1!~"^#"' \
  | awk 'BEGIN{OFS=FS="\t"} { if($3=="mRNA") {$3="transcript"; print} else {print} }' \
  > ${PREFIX}.genes.gff3
```

4 genes don't have features returned by `gffread` because they are misinformed (multiple overlapping exons) and very short (the shortest is 1 aa). Extract the IDs of these genes so we can use them for downstream analysis.

```bash
grep 'transcript' ${PREFIX}.genes.gff3 \
  | cut -f9 | sed -e 's/ID=//' | sort \
  > dataset.genes_in_gff.txt
```

**CDS**

```bash
zcat SC_transcripts_ltAED075.fa.gz \
  | sed -e "s/>/>${PREFIX}___/" \
  | ~/scripts/grepf_fasta.py -f <(cat dataset.longest.txt dataset.genes_in_gff.txt \
                                    | sort | uniq -d) \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.genes.cds.fna
```

**PEP**

```bash
zcat SC_proteins_ltAED075.fa.gz \
  | sed -e "s/>/>${PREFIX}___/" \
  | ~/scripts/grepf_fasta.py -f <(cat dataset.longest.txt dataset.genes_in_gff.txt \
                                    | sort | uniq -d) \
  | reformat.sh fastawrap=0 iupacToN=f ignorejunk=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.genes.pep.faa
```

## *Thelohanellus kitauei* TJCNv1

> Myxozoa

From: Wuqing, Tianjin, China in July 2007

Version: 1.0

Paper: https://doi.org/10.1093/gbe/evu247

Website: JWZT00000000.1 

```bash
# Download genome+gff3+PEP+CDS from NCBI
tar -xvf genome_assemblies_cds_fasta.tar
tar -xvf genome_assemblies_genome_fasta.tar
tar -xvf genome_assemblies_genome_gff.tar
tar -xvf genome_assemblies_prot_fasta.tar
```

```bash
PREFIX="Thelohanellus_kitauei_TJCNv1"
```

**Genome**

```bash
zcat ncbi-genomes-2023-06-06/GCA_000827895.1_ASM82789v1_genomic.fna.gz \
  | sed -e "s/>/>${PREFIX}___/" \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.assembly.fasta
```

**Get longest ORF per gene**

```bash
zcat ncbi-genomes-2023-06-06/GCA_000827895.1_ASM82789v1_genomic.gff.gz \
  | awk -F'\t' '$3=="CDS"' \
  | sed -e "s/[^\t]*ID=cds-\([^,;]*\).*Name=\([^;]*\).*/${PREFIX}___\2\t\1/" \
  > ncbi-genomes-2023-06-06.CDS_features.tsv
  
cat ncbi-genomes-2023-06-06.CDS_features.tsv \
  | awk -F'\t' '$3=="CDS" {
      CDS[$9]=CDS[$9]+(($5-$4)+1)
      GENE[$9]=$10
    } END {
      for (i in CDS){
        print GENE[i]"\t"i"\t"CDS[i]
    } }' \
  > ncbi-genomes-2023-06-06.CDS_lengths.tsv

cat ncbi-genomes-2023-06-06.CDS_lengths.tsv \
  | sort -k1,1 -k3,3nr \
  | awk '!seen[$1]++{print $2}' \
  > ncbi-genomes-2023-06-06.longest.txt
```

**CDS**

```bash
zcat ncbi-genomes-2023-06-06/GCA_000827895.1_ASM82789v1_cds_from_genomic.fna.gz \
  | sed -e "s/.*protein_id=\([^]]*\)\].*/>${PREFIX}___\1/" \
  | ~/scripts/grepf_fasta.py -f ncbi-genomes-2023-06-06.longest.txt \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.genes.cds.fna
```

**PEP**

```bash
zcat ncbi-genomes-2023-06-06/GCA_000827895.1_ASM82789v1_protein.faa.gz \
  | sed -e "s/>/>${PREFIX}___/" \
  | ~/scripts/grepf_fasta.py -f ncbi-genomes-2023-06-06.longest.txt \
  | reformat.sh fastawrap=0 iupacToN=f ignorejunk=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.genes.pep.faa
```

**GFF**

```bash
zcat ncbi-genomes-2023-06-06/GCA_000827895.1_ASM82789v1_genomic.gff.gz \
  | awk -F'\t' '$3=="CDS"' \
  | sed -e "s/[^\t]*Name=\([^;]*\).*/Parent=${PREFIX}___\1\t${PREFIX}___\1/" \
  | sed -e "s/^/${PREFIX}___/" \
  | ~/scripts/grepf_column.py -c 10 -f ncbi-genomes-2023-06-06.longest.txt \
  | cut -f1-9 \
  | gffread \
  | awk '$1!~"^#"' \
  | awk 'BEGIN{OFS=FS="\t"} { if($3!="CDS"){print} else {print; $3="exon"; print} }' \
  > ${PREFIX}.genes.gff3
```

## *Trichoplax adhaerens* RSSAv1

> Placozoan

From: Red Sea

Version: 1.0

Paper: https://doi.org/10.1038/nature07191

Website: https://mycocosm.jgi.doe.gov/Triad1/Triad1.home.html

```bash
# Download manually from JGI
# Use the "best" model per locus that JGI provides
```

```bash
PREFIX="Trichoplax_adhaerens_RSSAv1"
```

**Genome**

```bash
zcat Triad1_genomic_scaffolds.fasta.gz \
  | sed -e "s/>/>${PREFIX}___/" \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.assembly.fasta
```

**CDS**

```bash
zcat Triad1_best_transcripts.fasta.gz \
  | sed -e "s/>.*|/>${PREFIX}___/" \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.genes.cds.fna
```

**PEP**

```bash
zcat Triad1_best_proteins.fasta.gz \
  | sed -e "s/>.*|/>${PREFIX}___/" \
  | reformat.sh fastawrap=0 iupacToN=f ignorejunk=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  | sed -e 's/*$//' \
  > ${PREFIX}.genes.pep.faa
```

**GFF**

```bash
zcat Triad1_best_genes.gff.gz \
  | awk '$1!~"^#"' \
  | awk -F'\t' '$3=="exon" || $3=="CDS"' \
  | sed -e "s/^/${PREFIX}___/" \
  | sed -e 's/"//g' -e "s/name \([^,;]*\)/ID=${PREFIX}___\1;Parent=${PREFIX}___\1/" \
  | dos2unix \
  | gffread \
  | awk '$1!~"^#"' \
  | awk 'BEGIN{OFS=FS="\t"} { if($3=="mRNA") {$3="transcript"; print} else {print} }' \
  > ${PREFIX}.genes.gff3
```

## *Trichoplax* sp. H2 CRPAv1

> Placozoan

From: Caribbean, Bocas del Toro, Panama in 2003

Version: 1.0

Paper: https://doi.org/10.1038%2Fs41598-018-29400-y

Website: NCBI:  GCA_003344405.1

```bash
# Download from NCBI
unzip GCA_003344405.1.zip
```

```bash
PREFIX="Trichoplax_sp_H2_CRPAv1"
```

**Genome**

```bash
cat ncbi_dataset/data/GCA_003344405.1/GCA_003344405.1_TrispH2_1.0_genomic.fna \
  | sed -e "s/>/>${PREFIX}___/" \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.assembly.fasta
```

**Get longest ORF per gene**

Ignore CDS from pseudogenes. Only 3 pseudo genes in dataset.

```bash
cat ncbi_dataset/data/GCA_003344405.1/genomic.gff \
  | awk -F'\t' '$3=="CDS"' \
  | sed -e "s/[^\t]*Parent=\([^,;]*\).*Name=\([^;]*\).*/${PREFIX}___\2\t\1/" \
  > ncbi_dataset.CDS_features.tsv
  
cat ncbi_dataset.CDS_features.tsv \
  | awk -F'\t' '$3=="CDS" {
      CDS[$9]=CDS[$9]+(($5-$4)+1)
      GENE[$9]=$10
    } END {
      for (i in CDS){
        print GENE[i]"\t"i"\t"CDS[i]
    } }' \
  > ncbi_dataset.CDS_lengths.tsv

cat ncbi_dataset.CDS_lengths.tsv \
  | sort -k1,1 -k3,3nr \
  | awk '!seen[$1]++{print $2}' \
  > ncbi_dataset.longest.txt
```

**CDS**

```bash
cat ncbi_dataset/data/GCA_003344405.1/cds_from_genomic.fna \
  | sed -e "s/.*protein_id=\([^]]*\)\].*/>${PREFIX}___\1/" \
  | ~/scripts/grepf_fasta.py -f ncbi_dataset.longest.txt \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.genes.cds.fna
```

**PEP**

```bash
cat ncbi_dataset/data/GCA_003344405.1/protein.faa \
  | sed -e "s/>/>${PREFIX}___/" \
  | ~/scripts/grepf_fasta.py -f ncbi_dataset.longest.txt \
  | reformat.sh fastawrap=0 iupacToN=f ignorejunk=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.genes.pep.faa
```

**GFF**

```bash
cat ncbi_dataset/data/GCA_003344405.1/genomic.gff \
  | awk -F'\t' '$3=="CDS"' \
  | sed -e "s/[^\t]*Name=\([^;]*\).*/Parent=${PREFIX}___\1\t${PREFIX}___\1/" \
  | sed -e "s/^/${PREFIX}___/" \
  | ~/scripts/grepf_column.py -c 10 -f ncbi_dataset.longest.txt \
  | cut -f 1-9 \
  | gffread \
  | awk '$1!~"^#"' \
  | awk 'BEGIN{OFS=FS="\t"} { if($3!="CDS"){print} else {print; $3="exon"; print} }' \
  > ${PREFIX}.genes.gff3
```

## *Xenia sp.* CTEAv1

> Soft coral

From: obtained from a local coral aquarium shop called CTE Aquatics

Version: 1.0

Paper: https://doi.org/10.1038/s41586-020-2385-7

Website: https://cmo.carnegiescience.edu/data

```bash
wget https://cmo.carnegiescience.edu/endosymbiosis/genome/xenSp1.scaffolds.fa
wget https://cmo.carnegiescience.edu/endosymbiosis/genome/xenSp1.gff3
wget https://cmo.carnegiescience.edu/endosymbiosis/genome/xenSp1.proteins.fa
wget https://cmo.carnegiescience.edu/endosymbiosis/genome/xenSp1.transcripts.fa
```

```bash
PREFIX="Xenia_sp_CTEAv1"
```

**Genome**

```bash
cat xenSp1.scaffolds.fa \
  | sed -e "s/>/>${PREFIX}___/" \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.assembly.fasta
```

**Get longest ORF per gene**

```bash
cat xenSp1.gff3 \
  | awk '$1!~"^#"' \
  | awk -F'\t' '$3=="mRNA"' \
  | sed -e "s/.*ID=\([^,;]*\).*Parent=\([^,;]*\).*/${PREFIX}___\2\t${PREFIX}___\1/" \
  | sort -k1,1 -k2,2 \
  > dataset.gene2mRNA.tsv

cat xenSp1.gff3 \
  | awk '$1!~"^#"' \
  | awk -F'\t' '$3=="CDS"' \
  | sed -e "s/.*Parent=\([^,;]*\).*/\0\t${PREFIX}___\1/" \
  | awk -F'\t' '{
      CDS[$10]=CDS[$10]+(($5-$4)+1)
    } END {
      for (i in CDS){
        print i"\t"CDS[i]
    } }' \
  > dataset.CDS_lengths.tsv

~/scripts/add_value_to_table.py -c 2 -d "NA" \
  -i dataset.gene2mRNA.tsv \
  -a dataset.CDS_lengths.tsv \
  -o dataset.gene2mRNA2CDS_lengths.tsv

cat dataset.gene2mRNA2CDS_lengths.tsv \
  | sort -k1,1 -k3,3nr \
  | awk '!seen[$1]++{print $2}' \
  > dataset.longest.txt
```

**GFF**

```bash
cat xenSp1.gff3 \
  | awk '$1!~"^#"' \
  | awk -F'\t' '$3=="exon" || $3=="CDS"' \
  | sed -e "s/^/${PREFIX}___/" \
  | sed -e "s/ID=/ID=${PREFIX}___/" \
  | sed -e "s/Parent=/Parent=${PREFIX}___/" \
  | dos2unix \
  | sed -e 's/.*Parent=\([^,;]*\).*/\0\t\1/' \
  | ~/scripts/grepf_column.py -c 10 -f dataset.longest.txt \
  | cut -f1-9 \
  | gffread \
  | awk '$1!~"^#"' \
  | awk 'BEGIN{OFS=FS="\t"} { if($3=="mRNA") {$3="transcript"; print} else {print} }' \
  > ${PREFIX}.genes.gff3
```

**CDS**

```bash
cat xenSp1.transcripts.fa \
  | sed -e "s/>/>${PREFIX}___/" \
  | ~/scripts/grepf_fasta.py -f dataset.longest.txt \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.genes.cds.fna
```

**PEP**

```bash
cat xenSp1.proteins.fa \
  | sed -e "s/>/>${PREFIX}___/" \
  | ~/scripts/grepf_fasta.py -f dataset.longest.txt \
  | reformat.sh fastawrap=0 iupacToN=f ignorejunk=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.genes.pep.faa
```

## *Xestospongia testudinaria* RSSAv1

> Sponge

From: Fsar Reef (22.228408 N, 39.028187E) on the Red Sea coast of Saudi Arabia

Version: 1.0

Paper: https://doi.org/10.1186/s12864-016-2501-0

Website: http://xt.reefgenomics.org/

```bash
wget http://xt.reefgenomics.org/download/XT_genome_800.fa.gz
wget http://xt.reefgenomics.org/download/XT_final.gff.gz
wget http://xt.reefgenomics.org/download/XT_proteins_ltAED075.fa.gz
wget http://xt.reefgenomics.org/download/XT_transcripts_ltAED075.fa.gz
```

```bash
PREFIX="Xestospongia_testudinaria_RSSAv1"
```

**Genome**

```bash
zcat XT_genome_800.fa.gz \
  | sed -e "s/>XT_/>${PREFIX}___/" \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.assembly.fasta
```

**Get longest ORF per gene**

```bash
zcat XT_final.gff.gz \
  | awk '$1!~"^#"' \
  | awk -F'\t' '$2=="maker" && $3=="mRNA"' \
  | sed -e "s/.*ID=\([^,;]*\).*Parent=\([^,;]*\).*/${PREFIX}___\2\t${PREFIX}___\1/" \
  | sort -k1,1 -k2,2 \
  > dataset.gene2mRNA.tsv

zcat XT_final.gff.gz \
  | awk '$1!~"^#"' \
  | awk -F'\t' '$2=="maker" && $3=="CDS"' \
  | sed -e "s/.*Parent=\([^,;]*\).*/\0\t${PREFIX}___\1/" \
  | awk -F'\t' '{
      CDS[$10]=CDS[$10]+(($5-$4)+1)
    } END {
      for (i in CDS){
        print i"\t"CDS[i]
    } }' \
  > dataset.CDS_lengths.tsv

~/scripts/add_value_to_table.py -c 2 -d "NA" \
  -i dataset.gene2mRNA.tsv \
  -a dataset.CDS_lengths.tsv \
  -o dataset.gene2mRNA2CDS_lengths.tsv

cat dataset.gene2mRNA2CDS_lengths.tsv \
  | ~/scripts/grepf_column.py -c 2 \
      -f <(zcat XT_proteins_ltAED075.fa.gz | grep '>' | sed -e "s/>/${PREFIX}___/" -e 's/ .*//') \
  | sort -k1,1 -k3,3nr \
  | awk '!seen[$1]++{print $2}' \
  > dataset.longest.txt
```

**GFF**

```bash
zcat XT_final.gff.gz \
  | awk '$1!~"^#"' \
  | awk -F'\t' '$2=="maker" && ($3=="exon" || $3=="CDS")' \
  | sed -e "s/^XT_/${PREFIX}___/" \
  | sed -e "s/ID=/ID=${PREFIX}___/" \
  | sed -e "s/Parent=/Parent=${PREFIX}___/" \
  | dos2unix \
  | sed -e 's/.*Parent=\([^,;]*\).*/\0\t\1/' \
  | ~/scripts/grepf_column.py -c 10 -f dataset.longest.txt \
  | cut -f1-9 \
  | gffread \
  | awk '$1!~"^#"' \
  | awk 'BEGIN{OFS=FS="\t"} { if($3=="mRNA") {$3="transcript"; print} else {print} }' \
  > ${PREFIX}.genes.gff3
```

1 gene don't have features returned by `gffread` because it is misinformed (multiple overlapping exons) and very short (2 aa). Extract IDs of valid genes so we can use them for downstream analysis.

```bash
grep 'transcript' ${PREFIX}.genes.gff3 \
  | cut -f9 | sed -e 's/ID=//' | sort \
  > dataset.genes_in_gff.txt
```

**CDS**

```bash
zcat XT_transcripts_ltAED075.fa.gz \
  | sed -e "s/>/>${PREFIX}___/" \
  | ~/scripts/grepf_fasta.py -f dataset.genes_in_gff.txt \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.genes.cds.fna
```

**PEP**

```bash
zcat XT_proteins_ltAED075.fa.gz \
  | sed -e "s/>/>${PREFIX}___/" \
  | ~/scripts/grepf_fasta.py -f <(cat dataset.longest.txt dataset.genes_in_gff.txt \
                                    | sort | uniq -d) \
  | reformat.sh fastawrap=0 iupacToN=f ignorejunk=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.genes.pep.faa
```

# Out groups Transcripts

## *Alatina alata* BNNLv1-TRANS

> Jellyfish

From: Bonaire, The Netherlands (April 2014, 22:00–01:00)

Version: 1.0

Paper: https://doi.org/10.1093/gigascience/giz069

Website: http://ryanlab.whitney.ufl.edu/genomes/Aala/

```bash
git clone https://github.com/josephryan/Ohdera_et_al_2018.git
```

```bash
PREFIX="Alatina_alata_BNNLv1-TRANS"
```

**CDS**

CDS not available.

**PEP**

```bash
zcat Ohdera_et_al_2018/AA_Files/Alatina_alata.transdecoder.faa.gz \
  | sed -e "s/>Alatina_/>${PREFIX}___/" \
  | reformat.sh fastawrap=0 iupacToN=f ignorejunk=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  | sed -e 's/*$//' \
  > ${PREFIX}.transcripts.pep.faa
```

## *Anthopleura elegantissima* REEFv1

> Sea anemone

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







## *Calvadosia cruxmelitensis* NANAv1

> Stauromedusae

Paper: NA

From: NA

Version: 1.0

Website: NCBI:  HAHC00000000.1

```bash
wget https://sra-download.ncbi.nlm.nih.gov/traces/wgs03/wgs_aux/HA/HC/HAHC01/HAHC01.1.fsa_nt.gz
```

```bash
PREFIX="Calvadosia_cruxmelitensis_NANAv1"
```

**Transcripts**

Extract transcripts and predict ORFs.

```bash
zcat HAHC01.1.fsa_nt.gz \
  | sed -e "s/>/>${PREFIX}___/" \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.transcripts.fasta
```

Get gene to transcript name mappings.

```bash
zcat HAHC01.1.fsa_nt.gz \
  | grep '>' | sed -e "s/>/${PREFIX}___/" \
  | awk '{print $6"\t"$1}' | sed -e 's/_i[^\t]*//' \
  > ${PREFIX}.transcripts.gene2trans.txt
```

Predict ORFs in transcripts using `transdecoder`.

```bash
./run_Transdecoder.sh
```

Get best ORF per "gene"

```bash
grep '>' ${PREFIX}.transcripts.fasta.transdecoder.pep \
  | sed -e 's/.* \([^~]*\)~~\([^ ]*\) .*score=\([^ ,]*\).*/\1\t\2\t\3/' \
  | sort -k1,1 -k3,3nr \
  | awk -F'\t' '!seen[$1]++{print $2}' \
  > ${PREFIX}.transcripts.fasta.transdecoder.topPerGene.txt
```

**CDS**

```bash
cat ${PREFIX}.transcripts.fasta.transdecoder.cds \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  | ~/scripts/grepf_fasta.py -f ${PREFIX}.transcripts.fasta.transdecoder.topPerGene.txt \
  > ${PREFIX}.transcripts.cds.fna
```

**PEP**

```bash
cat ${PREFIX}.transcripts.fasta.transdecoder.pep \
  | sed -e 's/*$//' \
  | reformat.sh fastawrap=0 iupacToN=f ignorejunk=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  | ~/scripts/grepf_fasta.py -f ${PREFIX}.transcripts.fasta.transdecoder.topPerGene.txt \
  > ${PREFIX}.transcripts.pep.faa
```

## *Lucernaria quadricornis* AIUSv1

> Stauromedusae

Paper: NA

From: Appledore  Island, Maine, USA

Version: 1.0

Website: NCBI: HAHD00000000.1

```bash
wget https://sra-download.ncbi.nlm.nih.gov/traces/wgs03/wgs_aux/HA/HD/HAHD01/HAHD01.1.fsa_nt.gz
```

```bash
PREFIX="Lucernaria_quadricornis_AIUSv1"
```

**Transcripts**

Extract transcripts and predict ORFs.

```bash
zcat HAHD01.1.fsa_nt.gz \
  | sed -e "s/>/>${PREFIX}___/" \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.transcripts.fasta
```

Get gene to transcript name mappings.

```bash
zcat HAHD01.1.fsa_nt.gz \
  | grep '>' | sed -e "s/>/${PREFIX}___/" \
  | awk '{print $6"\t"$1}' | sed -e 's/_i[^\t]*//' \
  > ${PREFIX}.transcripts.gene2trans.txt
```

Predict ORFs in transcripts using `transdecoder`.

```bash
./run_Transdecoder.sh
```

Get best ORF per "gene"

```bash
grep '>' ${PREFIX}.transcripts.fasta.transdecoder.pep \
  | sed -e 's/.* \([^~]*\)~~\([^ ]*\) .*score=\([^ ,]*\).*/\1\t\2\t\3/' \
  | sort -k1,1 -k3,3nr \
  | awk -F'\t' '!seen[$1]++{print $2}' \
  > ${PREFIX}.transcripts.fasta.transdecoder.topPerGene.txt
```

**CDS**

```bash
cat ${PREFIX}.transcripts.fasta.transdecoder.cds \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  | ~/scripts/grepf_fasta.py -f ${PREFIX}.transcripts.fasta.transdecoder.topPerGene.txt \
  > ${PREFIX}.transcripts.cds.fna
```

**PEP**

```bash
cat ${PREFIX}.transcripts.fasta.transdecoder.pep \
  | sed -e 's/*$//' \
  | reformat.sh fastawrap=0 iupacToN=f ignorejunk=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  | ~/scripts/grepf_fasta.py -f ${PREFIX}.transcripts.fasta.transdecoder.topPerGene.txt \
  > ${PREFIX}.transcripts.pep.faa
```

## *Pachycerianthus borealis* NANAv1 (AKA *Cerianthus borealis*)

> Stauromedusae

Paper: NA

From: NA

Version: 1.0

Website: NCBI: HAGY00000000.1

```bash
wget https://sra-download.ncbi.nlm.nih.gov/traces/wgs03/wgs_aux/HA/GY/HAGY01/HAGY01.1.fsa_nt.gz
```

```bash
PREFIX="Pachycerianthus_borealis_NANAv1"
```

**Transcripts**

Extract transcripts and predict ORFs.

```bash
zcat HAGY01.1.fsa_nt.gz \
  | sed -e "s/>/>${PREFIX}___/" \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.transcripts.fasta
```

Get gene to transcript name mappings.

```bash
zcat HAGY01.1.fsa_nt.gz \
  | grep '>' | sed -e "s/>/${PREFIX}___/" \
  | awk '{print $6"\t"$1}' | sed -e 's/_i[^\t]*//' \
  > ${PREFIX}.transcripts.gene2trans.txt
```

Predict ORFs in transcripts using `transdecoder`.

```bash
./run_Transdecoder.sh
```

Get best ORF per "gene"

```bash
grep '>' ${PREFIX}.transcripts.fasta.transdecoder.pep \
  | sed -e 's/.* \([^~]*\)~~\([^ ]*\) .*score=\([^ ,]*\).*/\1\t\2\t\3/' \
  | sort -k1,1 -k3,3nr \
  | awk -F'\t' '!seen[$1]++{print $2}' \
  > ${PREFIX}.transcripts.fasta.transdecoder.topPerGene.txt
```

**CDS**

```bash
cat ${PREFIX}.transcripts.fasta.transdecoder.cds \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  | ~/scripts/grepf_fasta.py -f ${PREFIX}.transcripts.fasta.transdecoder.topPerGene.txt \
  > ${PREFIX}.transcripts.cds.fna
```

**PEP**

```bash
cat ${PREFIX}.transcripts.fasta.transdecoder.pep \
  | sed -e 's/*$//' \
  | reformat.sh fastawrap=0 iupacToN=f ignorejunk=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  | ~/scripts/grepf_fasta.py -f ${PREFIX}.transcripts.fasta.transdecoder.topPerGene.txt \
  > ${PREFIX}.transcripts.pep.faa
```

## *Haliclystus sanjuanensis* SJUSv1

> Stauromedusae

Paper: NA

From: San Juan Island, Washington, USA

Version: 1.0

Website: NCBI: HAHB00000000.1

```bash
wget https://sra-download.ncbi.nlm.nih.gov/traces/wgs03/wgs_aux/HA/HB/HAHB01/HAHB01.1.fsa_nt.gz
```

```bash
PREFIX="Haliclystus_sanjuanensis_SJUSv1"
```

**Transcripts**

Extract transcripts and predict ORFs.

```bash
zcat HAHB01.1.fsa_nt.gz \
  | sed -e "s/>/>${PREFIX}___/" \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.transcripts.fasta
```

Get gene to transcript name mappings.

```bash
zcat HAHB01.1.fsa_nt.gz \
  | grep '>' | sed -e "s/>/${PREFIX}___/" \
  | awk '{print $6"\t"$1}' \
  | sed -e 's/,//' \
  > ${PREFIX}.transcripts.gene2trans.txt
```

Predict ORFs in transcripts using `transdecoder`.

```bash
./run_Transdecoder.sh
```

Get best ORF per "gene"

```bash
grep '>' ${PREFIX}.transcripts.fasta.transdecoder.pep \
  | sed -e 's/.* \([^~]*\)~~\([^ ]*\) .*score=\([^ ,]*\).*/\1\t\2\t\3/' \
  | sort -k1,1 -k3,3nr \
  | awk -F'\t' '!seen[$1]++{print $2}' \
  > ${PREFIX}.transcripts.fasta.transdecoder.topPerGene.txt
```

**CDS**

```bash
cat ${PREFIX}.transcripts.fasta.transdecoder.cds \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  | ~/scripts/grepf_fasta.py -f ${PREFIX}.transcripts.fasta.transdecoder.topPerGene.txt \
  > ${PREFIX}.transcripts.cds.fna
```

**PEP**

```bash
cat ${PREFIX}.transcripts.fasta.transdecoder.pep \
  | sed -e 's/*$//' \
  | reformat.sh fastawrap=0 iupacToN=f ignorejunk=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  | ~/scripts/grepf_fasta.py -f ${PREFIX}.transcripts.fasta.transdecoder.topPerGene.txt \
  > ${PREFIX}.transcripts.pep.faa
```

## *Haliclystus auricula* NANAv1

> Stauromedusae

Paper: NA

From: NA

Version: 1.0

Website: NCBI: HAHA00000000.1

```bash
wget https://sra-download.ncbi.nlm.nih.gov/traces/wgs03/wgs_aux/HA/HA/HAHA01/HAHA01.1.fsa_nt.gz
```

```bash
PREFIX="Haliclystus_auricula_NANAv1"
```

**Transcripts**

Extract transcripts and predict ORFs.

```bash
zcat HAHA01.1.fsa_nt.gz \
  | sed -e "s/>/>${PREFIX}___/" \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.transcripts.fasta
```

Get gene to transcript name mappings.

```bash
zcat HAHA01.1.fsa_nt.gz \
  | grep '>' | sed -e "s/>/${PREFIX}___/" \
  | awk '{print $6"\t"$1}' | sed -e 's/_i[^\t]*//' \
  > ${PREFIX}.transcripts.gene2trans.txt
```

Predict ORFs in transcripts using `transdecoder`.

```bash
./run_Transdecoder.sh
```

Get best ORF per "gene"

```bash
grep '>' ${PREFIX}.transcripts.fasta.transdecoder.pep \
  | sed -e 's/.* \([^~]*\)~~\([^ ]*\) .*score=\([^ ,]*\).*/\1\t\2\t\3/' \
  | sort -k1,1 -k3,3nr \
  | awk -F'\t' '!seen[$1]++{print $2}' \
  > ${PREFIX}.transcripts.fasta.transdecoder.topPerGene.txt
```

**CDS**

```bash
cat ${PREFIX}.transcripts.fasta.transdecoder.cds \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  | ~/scripts/grepf_fasta.py -f ${PREFIX}.transcripts.fasta.transdecoder.topPerGene.txt \
  > ${PREFIX}.transcripts.cds.fna
```

**PEP**

```bash
cat ${PREFIX}.transcripts.fasta.transdecoder.pep \
  | sed -e 's/*$//' \
  | reformat.sh fastawrap=0 iupacToN=f ignorejunk=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  | ~/scripts/grepf_fasta.py -f ${PREFIX}.transcripts.fasta.transdecoder.topPerGene.txt \
  > ${PREFIX}.transcripts.pep.faa
```

## *Craterolophus convolvulus* RYUSv1

> Stauromedusae

Paper: NA

From: Rye Harbor State Park, New Hampshire, USA

Version: 1.0

Website: NCBI: HAGZ00000000.1

```bash
wget https://sra-download.ncbi.nlm.nih.gov/traces/wgs03/wgs_aux/HA/GZ/HAGZ01/HAGZ01.1.fsa_nt.gz
```

```bash
PREFIX="Craterolophus_convolvulus_RYUSv1"
```

**Transcripts**

Extract transcripts and predict ORFs.

```bash
zcat HAGZ01.1.fsa_nt.gz \
  | sed -e "s/>/>${PREFIX}___/" \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.transcripts.fasta
```

Get gene to transcript name mappings.

```bash
zcat HAGZ01.1.fsa_nt.gz \
  | grep '>' | sed -e "s/>/${PREFIX}___/" \
  | awk '{print $6"\t"$1}' | sed -e 's/_i[^\t]*//' \
  > ${PREFIX}.transcripts.gene2trans.txt
```

Predict ORFs in transcripts using `transdecoder`.

```bash
./run_Transdecoder.sh
```

Get best ORF per "gene"

```bash
grep '>' ${PREFIX}.transcripts.fasta.transdecoder.pep \
  | sed -e 's/.* \([^~]*\)~~\([^ ]*\) .*score=\([^ ,]*\).*/\1\t\2\t\3/' \
  | sort -k1,1 -k3,3nr \
  | awk -F'\t' '!seen[$1]++{print $2}' \
  > ${PREFIX}.transcripts.fasta.transdecoder.topPerGene.txt
```

**CDS**

```bash
cat ${PREFIX}.transcripts.fasta.transdecoder.cds \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  | ~/scripts/grepf_fasta.py -f ${PREFIX}.transcripts.fasta.transdecoder.topPerGene.txt \
  > ${PREFIX}.transcripts.cds.fna
```

**PEP**

```bash
cat ${PREFIX}.transcripts.fasta.transdecoder.pep \
  | sed -e 's/*$//' \
  | reformat.sh fastawrap=0 iupacToN=f ignorejunk=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  | ~/scripts/grepf_fasta.py -f ${PREFIX}.transcripts.fasta.transdecoder.topPerGene.txt \
  > ${PREFIX}.transcripts.pep.faa
```

## *Hydra oligactis* GNCHv1

> Hydrozoa

Paper: NA

From: Genev, Switzerland

Version: 1.0

Website: NCBI: GHUC00000000.1

```bash
wget https://sra-download.ncbi.nlm.nih.gov/traces/wgs01/wgs_aux/GH/UC/GHUC01/GHUC01.1.fsa_nt.gz
```

```bash
PREFIX="Hydra_oligactis_GNCHv1"
```

**Transcripts**

Extract transcripts and predict ORFs.

```bash
zcat GHUC01.1.fsa_nt.gz \
  | sed -e "s/>/>${PREFIX}___/" \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.transcripts.fasta
```

Get gene to transcript name mappings.

```bash
zcat GHUC01.1.fsa_nt.gz \
  | grep '>' | sed -e "s/>/${PREFIX}___/" \
  | awk '{print $5"\t"$1}' | sed -e 's/_i[^\t]*//' \
  > ${PREFIX}.transcripts.gene2trans.txt
```

Predict ORFs in transcripts using `transdecoder`.

```bash
./run_Transdecoder.sh
```

Get best ORF per "gene"

```bash
grep '>' ${PREFIX}.transcripts.fasta.transdecoder.pep \
  | sed -e 's/.* \([^~]*\)~~\([^ ]*\) .*score=\([^ ,]*\).*/\1\t\2\t\3/' \
  | sort -k1,1 -k3,3nr \
  | awk -F'\t' '!seen[$1]++{print $2}' \
  > ${PREFIX}.transcripts.fasta.transdecoder.topPerGene.txt
```

**CDS**

```bash
cat ${PREFIX}.transcripts.fasta.transdecoder.cds \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  | ~/scripts/grepf_fasta.py -f ${PREFIX}.transcripts.fasta.transdecoder.topPerGene.txt \
  > ${PREFIX}.transcripts.cds.fna
```

**PEP**

```bash
cat ${PREFIX}.transcripts.fasta.transdecoder.pep \
  | sed -e 's/*$//' \
  | reformat.sh fastawrap=0 iupacToN=f ignorejunk=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  | ~/scripts/grepf_fasta.py -f ${PREFIX}.transcripts.fasta.transdecoder.topPerGene.txt \
  > ${PREFIX}.transcripts.pep.faa
```

## *Palythoa variabilis* PCBRv1 (AKA *Protopalythoa variabilis*)

> Zoantharia

Paper: NA

From: Pernambuco coast, Brazil

Version: 1.0

Website: NCBI: GCVI00000000.1

```bash
wget https://sra-download.ncbi.nlm.nih.gov/traces/wgs03/wgs_aux/GC/VI/GCVI01/GCVI01.1.fsa_nt.gz
```

```bash
PREFIX="Palythoa_variabilis_PCBRv1"
```

**Transcripts**

Extract transcripts and predict ORFs.

```bash
zcat GCVI01.1.fsa_nt.gz \
  | sed -e "s/>/>${PREFIX}___/" \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.transcripts.fasta
```

Get gene to transcript name mappings.

```bash
zcat GCVI01.1.fsa_nt.gz \
  | grep '>' | sed -e "s/>/${PREFIX}___/" \
  | awk '{print $5"\t"$1}' \
  > ${PREFIX}.transcripts.gene2trans.txt
```

Predict ORFs in transcripts using `transdecoder`.

```bash
./run_Transdecoder.sh
```

Get best ORF per "gene"

```bash
grep '>' ${PREFIX}.transcripts.fasta.transdecoder.pep \
  | sed -e 's/.* \([^~]*\)~~\([^ ]*\) .*score=\([^ ,]*\).*/\1\t\2\t\3/' \
  | sort -k1,1 -k3,3nr \
  | awk -F'\t' '!seen[$1]++{print $2}' \
  > ${PREFIX}.transcripts.fasta.transdecoder.topPerGene.txt
```

**CDS**

```bash
cat ${PREFIX}.transcripts.fasta.transdecoder.cds \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  | ~/scripts/grepf_fasta.py -f ${PREFIX}.transcripts.fasta.transdecoder.topPerGene.txt \
  > ${PREFIX}.transcripts.cds.fna
```

**PEP**

```bash
cat ${PREFIX}.transcripts.fasta.transdecoder.pep \
  | sed -e 's/*$//' \
  | reformat.sh fastawrap=0 iupacToN=f ignorejunk=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  | ~/scripts/grepf_fasta.py -f ${PREFIX}.transcripts.fasta.transdecoder.topPerGene.txt \
  > ${PREFIX}.transcripts.pep.faa
```

## *Tripedalia cystophora* (Conant 1897) 1897v1

> Cubozoa

Paper: NA

From: NA

Version: 1.0

Website: NCBI: GGWE00000000.1

```bash
wget https://sra-download.ncbi.nlm.nih.gov/traces/wgs03/wgs_aux/GG/WE/GGWE01/GGWE01.1.fsa_nt.gz
```

```bash
PREFIX="Tripedalia_cystophora_Conant1897_1897v1"
```

**Transcripts**

Extract transcripts and predict ORFs.

```bash
zcat GGWE01.1.fsa_nt.gz \
  | sed -e "s/>/>${PREFIX}___/" \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.transcripts.fasta
```

Get gene to transcript name mappings.

```bash
zcat GGWE01.1.fsa_nt.gz \
  | grep '>' | sed -e "s/>/${PREFIX}___/" \
  | awk '{print $8"\t"$1}' \
  | sed -e 's/,//' \
  > ${PREFIX}.transcripts.gene2trans.txt
```

Predict ORFs in transcripts using `transdecoder`.

```bash
./run_Transdecoder.sh
```

Get best ORF per "gene"

```bash
grep '>' ${PREFIX}.transcripts.fasta.transdecoder.pep \
  | sed -e 's/.* \([^~]*\)~~\([^ ]*\) .*score=\([^ ,]*\).*/\1\t\2\t\3/' \
  | sort -k1,1 -k3,3nr \
  | awk -F'\t' '!seen[$1]++{print $2}' \
  > ${PREFIX}.transcripts.fasta.transdecoder.topPerGene.txt
```

**CDS**

```bash
cat ${PREFIX}.transcripts.fasta.transdecoder.cds \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  | ~/scripts/grepf_fasta.py -f ${PREFIX}.transcripts.fasta.transdecoder.topPerGene.txt \
  > ${PREFIX}.transcripts.cds.fna
```

**PEP**

```bash
cat ${PREFIX}.transcripts.fasta.transdecoder.pep \
  | sed -e 's/*$//' \
  | reformat.sh fastawrap=0 iupacToN=f ignorejunk=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  | ~/scripts/grepf_fasta.py -f ${PREFIX}.transcripts.fasta.transdecoder.topPerGene.txt \
  > ${PREFIX}.transcripts.pep.faa
```

## *Corallium rubrum* MRFRv1

> Alcyonacea

Paper: https://doi.org/10.1111/1755-0998.12383

From: Marseilles,  France

Version: 1.0

Website: https://datadryad.org/stash/dataset/doi:10.5061/dryad.31f77

```bash
unzip doi_10.5061_dryad.31f77__v1.zip
```

```bash
PREFIX="Corallium_rubrum_MRFRv1"
```

**Transcripts**

Extract transcripts and predict ORFs.

```bash
cat CoralliumRubrum_Contigs.tfa \
  | sed -e "s/>/>${PREFIX}___/" \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.transcripts.fasta
```

Get gene to transcript name mappings.

```bash
cat CoralliumRubrum_Contigs.tfa \
  | grep '>' | sed -e "s/>/${PREFIX}___/" \
  | awk '{print $1"\t"$1}' \
  > ${PREFIX}.transcripts.gene2trans.txt
```

Predict ORFs in transcripts using `transdecoder`.

```bash
./run_Transdecoder.sh
```

Get best ORF per "gene"

```bash
grep '>' ${PREFIX}.transcripts.fasta.transdecoder.pep \
  | sed -e 's/.* \([^~]*\)~~\([^ ]*\) .*score=\([^ ,]*\).*/\1\t\2\t\3/' \
  | sort -k1,1 -k3,3nr \
  | awk -F'\t' '!seen[$1]++{print $2}' \
  > ${PREFIX}.transcripts.fasta.transdecoder.topPerGene.txt
```

**CDS**

```bash
cat ${PREFIX}.transcripts.fasta.transdecoder.cds \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  | ~/scripts/grepf_fasta.py -f ${PREFIX}.transcripts.fasta.transdecoder.topPerGene.txt \
  > ${PREFIX}.transcripts.cds.fna
```

**PEP**

```bash
cat ${PREFIX}.transcripts.fasta.transdecoder.pep \
  | sed -e 's/*$//' \
  | reformat.sh fastawrap=0 iupacToN=f ignorejunk=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  | ~/scripts/grepf_fasta.py -f ${PREFIX}.transcripts.fasta.transdecoder.topPerGene.txt \
  > ${PREFIX}.transcripts.pep.faa
```

## *Gorgonia ventalina* LRPRv1

> Alcyonacea

Paper: https://doi.org/10.3389/fphys.2013.00180

From: Laurel Patch Reef, La Parguera, Puerto Rico

Version: 1.0

Website: https://figshare.com/articles/dataset/Gorgonia_ventalina_transcriptome/94326/2

```bash
# Download from figshare
```

```bash
PREFIX="Gorgonia_ventalina_LRPRv1"
```

**Transcripts**

Extract transcripts and predict ORFs.

```bash
cat CombinedGV.fa \
  | sed -e "s/>/>${PREFIX}___/" \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.transcripts.fasta
```

Get gene to transcript name mappings.

```bash
cat CombinedGV.fa \
  | grep '>' | sed -e "s/>/${PREFIX}___/" \
  | awk '{print $1"\t"$1}' \
  > ${PREFIX}.transcripts.gene2trans.txt
```

Predict ORFs in transcripts using `transdecoder`.

```bash
./run_Transdecoder.sh
```

Get best ORF per "gene"

```bash
grep '>' ${PREFIX}.transcripts.fasta.transdecoder.pep \
  | sed -e 's/.* \([^~]*\)~~\([^ ]*\) .*score=\([^ ,]*\).*/\1\t\2\t\3/' \
  | sort -k1,1 -k3,3nr \
  | awk -F'\t' '!seen[$1]++{print $2}' \
  > ${PREFIX}.transcripts.fasta.transdecoder.topPerGene.txt
```

**CDS**

```bash
cat ${PREFIX}.transcripts.fasta.transdecoder.cds \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  | ~/scripts/grepf_fasta.py -f ${PREFIX}.transcripts.fasta.transdecoder.topPerGene.txt \
  > ${PREFIX}.transcripts.cds.fna
```

**PEP**

```bash
cat ${PREFIX}.transcripts.fasta.transdecoder.pep \
  | sed -e 's/*$//' \
  | reformat.sh fastawrap=0 iupacToN=f ignorejunk=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  | ~/scripts/grepf_fasta.py -f ${PREFIX}.transcripts.fasta.transdecoder.topPerGene.txt \
  > ${PREFIX}.transcripts.pep.faa
```

## *Kudoa iwatai* RSILv1

> Myxozoa

Paper: NA

From: Eilat, Red sea, Israel

Version: 1.0

Website: NCBI: GBGI00000000.1

```bash
wget https://sra-download.ncbi.nlm.nih.gov/traces/wgs03/wgs_aux/GB/GI/GBGI01/GBGI01.1.fsa_nt.gz
```

```bash
PREFIX="Kudoa_iwatai_RSILv1"
```

**Transcripts**

Extract transcripts and predict ORFs.

```bash
zcat GBGI01.1.fsa_nt.gz \
  | sed -e "s/>/>${PREFIX}___/" \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.transcripts.fasta
```

Get gene to transcript name mappings.

```bash
zcat GBGI01.1.fsa_nt.gz \
  | grep '>' | sed -e "s/>/${PREFIX}___/" \
  | awk '{print $5"\t"$1}' | sed -e 's/_seq[^\t]*//' \
  > ${PREFIX}.transcripts.gene2trans.txt
```

Predict ORFs in transcripts using `transdecoder`.

```bash
./run_Transdecoder.sh
```

Get best ORF per "gene"

```bash
grep '>' ${PREFIX}.transcripts.fasta.transdecoder.pep \
  | sed -e 's/.* \([^~]*\)~~\([^ ]*\) .*score=\([^ ,]*\).*/\1\t\2\t\3/' \
  | sort -k1,1 -k3,3nr \
  | awk -F'\t' '!seen[$1]++{print $2}' \
  > ${PREFIX}.transcripts.fasta.transdecoder.topPerGene.txt
```

**CDS**

```bash
cat ${PREFIX}.transcripts.fasta.transdecoder.cds \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  | ~/scripts/grepf_fasta.py -f ${PREFIX}.transcripts.fasta.transdecoder.topPerGene.txt \
  > ${PREFIX}.transcripts.cds.fna
```

**PEP**

```bash
cat ${PREFIX}.transcripts.fasta.transdecoder.pep \
  | sed -e 's/*$//' \
  | reformat.sh fastawrap=0 iupacToN=f ignorejunk=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  | ~/scripts/grepf_fasta.py -f ${PREFIX}.transcripts.fasta.transdecoder.topPerGene.txt \
  > ${PREFIX}.transcripts.pep.faa
```

## *Myxobolus cerebralis* NANAv1

> Myxozoa

Paper: NA

From: NA

Version: 1.0

Website: NCBI: GBKL00000000.1

```bash
wget https://sra-download.ncbi.nlm.nih.gov/traces/wgs03/wgs_aux/GB/KL/GBKL01/GBKL01.1.fsa_nt.gz
```

```bash
PREFIX="Myxobolus_cerebralis_NANAv1"
```

**Transcripts**

Extract transcripts and predict ORFs.

```bash
zcat GBKL01.1.fsa_nt.gz \
  | sed -e "s/>/>${PREFIX}___/" \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.transcripts.fasta
```

Get gene to transcript name mappings.

```bash
zcat GBKL01.1.fsa_nt.gz \
  | grep '>' | sed -e "s/>/${PREFIX}___/" \
  | awk '{print $5"\t"$1}' | sed -e 's/_seq[^\t]*//' \
  > ${PREFIX}.transcripts.gene2trans.txt
```

Predict ORFs in transcripts using `transdecoder`.

```bash
./run_Transdecoder.sh
```

Get best ORF per "gene"

```bash
grep '>' ${PREFIX}.transcripts.fasta.transdecoder.pep \
  | sed -e 's/.* \([^~]*\)~~\([^ ]*\) .*score=\([^ ,]*\).*/\1\t\2\t\3/' \
  | sort -k1,1 -k3,3nr \
  | awk -F'\t' '!seen[$1]++{print $2}' \
  > ${PREFIX}.transcripts.fasta.transdecoder.topPerGene.txt
```

**CDS**

```bash
cat ${PREFIX}.transcripts.fasta.transdecoder.cds \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  | ~/scripts/grepf_fasta.py -f ${PREFIX}.transcripts.fasta.transdecoder.topPerGene.txt \
  > ${PREFIX}.transcripts.cds.fna
```

**PEP**

```bash
cat ${PREFIX}.transcripts.fasta.transdecoder.pep \
  | sed -e 's/*$//' \
  | reformat.sh fastawrap=0 iupacToN=f ignorejunk=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  | ~/scripts/grepf_fasta.py -f ${PREFIX}.transcripts.fasta.transdecoder.topPerGene.txt \
  > ${PREFIX}.transcripts.pep.faa
```

## *Myxobolus pendula* ONCAv1

> Myxozoa

Paper: https://doi.org/10.1186%2Fs12864-015-2039-6

From: NA

Version: 1.0

Website: Additional  file 7 from Manuscript

```bash
# 12864_2015_2039_MOESM7_ESM.fa from Manuscript
```

```bash
PREFIX="Myxobolus_pendula_ONCAv1"
```

**Transcripts**

Extract transcripts and predict ORFs.

```bash
cat 12864_2015_2039_MOESM7_ESM.fa \
  | sed -e "s/>/>${PREFIX}___/" \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.transcripts.fasta
```

Get gene to transcript name mappings.

```bash
cat 12864_2015_2039_MOESM7_ESM.fa \
  | grep '>' | sed -e "s/>/${PREFIX}___/" \
  | awk '{print $1"\t"$1}' \
  | sed -e "s/${PREFIX}___//" -e 's/_seq[^\t]*//' \
  > ${PREFIX}.transcripts.gene2trans.txt
```

Predict ORFs in transcripts using `transdecoder`.

```bash
./run_Transdecoder.sh
```

Get best ORF per "gene"

```bash
grep '>' ${PREFIX}.transcripts.fasta.transdecoder.pep \
  | sed -e 's/.* \([^~]*\)~~\([^ ]*\) .*score=\([^ ,]*\).*/\1\t\2\t\3/' \
  | sort -k1,1 -k3,3nr \
  | awk -F'\t' '!seen[$1]++{print $2}' \
  > ${PREFIX}.transcripts.fasta.transdecoder.topPerGene.txt
```

**CDS**

```bash
cat ${PREFIX}.transcripts.fasta.transdecoder.cds \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  | ~/scripts/grepf_fasta.py -f ${PREFIX}.transcripts.fasta.transdecoder.topPerGene.txt \
  > ${PREFIX}.transcripts.cds.fna
```

**PEP**

```bash
cat ${PREFIX}.transcripts.fasta.transdecoder.pep \
  | sed -e 's/*$//' \
  | reformat.sh fastawrap=0 iupacToN=f ignorejunk=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  | ~/scripts/grepf_fasta.py -f ${PREFIX}.transcripts.fasta.transdecoder.topPerGene.txt \
  > ${PREFIX}.transcripts.pep.faa
```

## *Polypodium hydriforme* OKUSv1

> Polypodiozoa

Paper: NA

From: NA

Version: 1.0

Website: NCBI: GBKL00000000.1

```bash
# NCBI Nucleotide
```

```bash
PREFIX="Polypodium_hydriforme_OKUSv1"
```

**Transcripts**

Extract transcripts and predict ORFs.

```bash
cat sequence.fasta \
  | sed -e "s/>/>${PREFIX}___/" \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.transcripts.fasta
```

Get gene to transcript name mappings.

```bash
cat sequence.fasta \
  | grep '>' | sed -e "s/>/${PREFIX}___/" \
  | awk '{print $5"\t"$1}' | sed -e 's/_i[^\t]*//' \
  > ${PREFIX}.transcripts.gene2trans.txt
```

Predict ORFs in transcripts using `transdecoder`.

```bash
./run_Transdecoder.sh
```

Get best ORF per "gene"

```bash
grep '>' ${PREFIX}.transcripts.fasta.transdecoder.pep \
  | sed -e 's/.* \([^~]*\)~~\([^ ]*\) .*score=\([^ ,]*\).*/\1\t\2\t\3/' \
  | sort -k1,1 -k3,3nr \
  | awk -F'\t' '!seen[$1]++{print $2}' \
  > ${PREFIX}.transcripts.fasta.transdecoder.topPerGene.txt
```

**CDS**

```bash
cat ${PREFIX}.transcripts.fasta.transdecoder.cds \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  | ~/scripts/grepf_fasta.py -f ${PREFIX}.transcripts.fasta.transdecoder.topPerGene.txt \
  > ${PREFIX}.transcripts.cds.fna
```

**PEP**

```bash
cat ${PREFIX}.transcripts.fasta.transdecoder.pep \
  | sed -e 's/*$//' \
  | reformat.sh fastawrap=0 iupacToN=f ignorejunk=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  | ~/scripts/grepf_fasta.py -f ${PREFIX}.transcripts.fasta.transdecoder.topPerGene.txt \
  > ${PREFIX}.transcripts.pep.faa
```

## *Calvadosia cruxmelitensis* CRGBv3.2-TRANS

> Jellyfish

From: Chimney Rock, off the coast of Penzance, Cornwall, England (January 2013)

Version: 3.2

Paper: https://doi.org/10.1093/gigascience/giz069

Website: http://ryanlab.whitney.ufl.edu/genomes/Ccrux/

```bash
git clone https://github.com/josephryan/Ohdera_et_al_2018.git
```

```bash
PREFIX="Calvadosia_cruxmelitensis_CRGBv3.2-TRANS"
```

**CDS**

No CDS available.

**PEP**

```bash
zcat Ohdera_et_al_2018/AA_Files/Ccrux.Trinity.faa.gz \
  | sed -e "s/>Ccrux./>${PREFIX}___/" \
  | reformat.sh fastawrap=0 iupacToN=f ignorejunk=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  | sed -e 's/*$//' \
  > ${PREFIX}.transcripts.pep.faa
```

## *Cassiopea xamachana* T1-Av2

> Jellyfish

From: Unknown (from line maintained in culture [T1-A])

Version: 2.0

Paper: https://doi.org/10.1093/gigascience/giz069

Website: http://ryanlab.whitney.ufl.edu/genomes/Cxam/

```bash
git clone https://github.com/josephryan/Ohdera_et_al_2018.git
```

```bash
PREFIX="Cassiopea_xamachana_T1-Av2"
```

**Genome**

Not available.

**CDS**

```bash
zcat Ohdera_et_al_2018/Condingseq/Cxam.augustus.codingseq.gz \
  | sed -e "s/>/>${PREFIX}___/" \
  | reformat.sh fastawrap=0 iupacToN=f ignorejunk=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.transcripts.cds.fna
```

**PEP**

```bash
zcat Ohdera_et_al_2018/AA_Files/Cxam.augustus.1.3.faa.gz \
  | sed -e "s/>Cxam_/>${PREFIX}___/" \
  | reformat.sh fastawrap=0 iupacToN=f ignorejunk=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  | sed -e 's/*$//' \
  > ${PREFIX}.transcripts.pep.faa
```

**GFF**

No CDS or GFF files were provided.

## *Cassiopea xamachana* T1-Av2-TRANS

> Jellyfish

From: Unknown (from line maintained in culture [T1-A])

Version: 2.0

Paper: https://doi.org/10.1093/gigascience/giz069

Website: http://ryanlab.whitney.ufl.edu/genomes/Cxam/

```bash
git clone https://github.com/josephryan/Ohdera_et_al_2018.git
```

```bash
PREFIX="Cassiopea_xamachana_T1-Av2-TRANS"
```

**CDS**

Not available.

**PEP**

```bash
zcat Ohdera_et_al_2018/AA_Files/Cxam.Scyph.transdecoder.faa.gz \
  | sed -e "s/>CxamS_/>${PREFIX}___/" -e 's/::/_/g' \
  | reformat.sh fastawrap=0 iupacToN=f ignorejunk=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  | sed -e 's/*$//' \
  > ${PREFIX}.transcripts.pep.faa
```

## *Corynactis australis* JBAUv1

> Corallimorpharia

Paper: NA

From: Jervis Bay, NSW, Australia

Version: 1.0

Website: NCBI:  GELM00000000.1

```bash
wget https://sra-download.ncbi.nlm.nih.gov/traces/wgs03/wgs_aux/GE/LM/GELM01/GELM01.1.fsa_nt.gz
```

```bash
PREFIX="Corynactis_australis_JBAUv1"
```

**Transcripts**

Extract transcripts and predict ORFs.

```bash
zcat GELM01.1.fsa_nt.gz \
  | sed -e "s/>/>${PREFIX}___/" -e 's/gb|//' -e 's/|//' \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.transcripts.fasta
```

Get gene to transcript name mappings.

```bash
zcat GELM01.1.fsa_nt.gz \
  | grep '>' | sed -e "s/>/${PREFIX}___/" -e 's/gb|//' -e 's/|//' \
  | awk '{print $5"\t"$1}' | sed -e 's/_seq[^\t]*//' \
  > ${PREFIX}.transcripts.gene2trans.txt
```

Predict ORFs in transcripts using `transdecoder`.

```bash
./run_Transdecoder.sh
```

Get best ORF per "gene"

```bash
grep '>' ${PREFIX}.transcripts.fasta.transdecoder.pep \
  | sed -e 's/.* \([^~]*\)~~\([^ ]*\) .*score=\([^ ,]*\).*/\1\t\2\t\3/' \
  | sort -k1,1 -k3,3nr \
  | awk -F'\t' '!seen[$1]++{print $2}' \
  > ${PREFIX}.transcripts.fasta.transdecoder.topPerGene.txt
```

**CDS**

```bash
cat ${PREFIX}.transcripts.fasta.transdecoder.cds \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  | ~/scripts/grepf_fasta.py -f ${PREFIX}.transcripts.fasta.transdecoder.topPerGene.txt \
  > ${PREFIX}.transcripts.cds.fna
```

**PEP**

```bash
cat ${PREFIX}.transcripts.fasta.transdecoder.pep \
  | sed -e 's/*$//' \
  | reformat.sh fastawrap=0 iupacToN=f ignorejunk=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  | ~/scripts/grepf_fasta.py -f ${PREFIX}.transcripts.fasta.transdecoder.topPerGene.txt \
  > ${PREFIX}.transcripts.pep.faa
```

## *Plexaura homomalla* PRUSv1

> Soft coral

Paper: https://doi.org/10.1007/s00338-021-02175-x

From: Puerto Rico: Caribbean Sea (17.9367 N 67.046 W)

Version: 1.0

Website: https://github.com/jessiepelosi/plexaurids

As far as I can tell the assembly was generated from this run: https://www.ncbi.nlm.nih.gov/sra/SRX9516166[accn]

```bash
git clone https://github.com/jessiepelosi/plexaurids.git
```

```bash
PREFIX="Plexaura_homomalla_PRUSv1"
```

**Transcripts**

Extract transcripts and predict ORFs.

```bash
zcat plexaurids/data/HF4_Trinity.annotated.fasta.gz \
  | sed -e "s/>/>${PREFIX}___/" \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.transcripts.fasta
```

Get gene to transcript name mappings.

```bash
grep '>' ${PREFIX}.transcripts.fasta | sed -e 's/>//' | awk '{print $1"\t"$1}' | sed -e 's/_i[^_]*$//' > ${PREFIX}.transcripts.gene2trans.txt
```

Predict ORFs in transcripts using `transdecoder`.

```bash
./run_Transdecoder.sh
```

**CDS**

```bash
cat ${PREFIX}.transcripts.fasta.transdecoder.cds \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.transcripts.cds.fna
```

**PEP**

```bash
cat ${PREFIX}.transcripts.fasta.transdecoder.pep \
  | sed -e 's/*$//' \
  | reformat.sh fastawrap=0 iupacToN=f ignorejunk=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.transcripts.pep.faa
```

## *Plexaura kukenthali* PRUSv1

> Soft coral

Paper: https://doi.org/10.1007/s00338-021-02175-x

From: Puerto Rico: Caribbean Sea (17.9367 N 67.046 W)

Version: 1.0

Website: https://github.com/jessiepelosi/plexaurids

As far as I can tell the assembly was generated from this run: https://www.ncbi.nlm.nih.gov/sra/SRX9516168[accn]

```bash
git clone https://github.com/jessiepelosi/plexaurids.git
```

```bash
PREFIX="Plexaura_kukenthali_PRUSv1"
```

**Transcripts**

Extract transcripts and predict ORFs.

```bash
zcat plexaurids/data/KG1_Trinity.annotated.fasta.gz \
  | sed -e "s/>/>${PREFIX}___/" \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.transcripts.fasta
```

Get gene to transcript name mappings.

```bash
grep '>' ${PREFIX}.transcripts.fasta | sed -e 's/>//' | awk '{print $1"\t"$1}' | sed -e 's/_i[^_]*$//' > ${PREFIX}.transcripts.gene2trans.txt
```

Predict ORFs in transcripts using `transdecoder`.

```bash
./run_Transdecoder.sh
```

**CDS**

```bash
cat ${PREFIX}.transcripts.fasta.transdecoder.cds \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.transcripts.cds.fna
```

**PEP**

```bash
cat ${PREFIX}.transcripts.fasta.transdecoder.pep \
  | sed -e 's/*$//' \
  | reformat.sh fastawrap=0 iupacToN=f ignorejunk=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.transcripts.pep.faa
```

## *Dynamena pumila* WSRUv1

> Hydrozoa

Paper: NA

From: Russia:  White Sea Biological Station, Republic Karelia

Version: 1.0

Website: NCBI: GHMC00000000.1

```bash
wget https://sra-download.ncbi.nlm.nih.gov/traces/wgs01/wgs_aux/GH/MC/GHMC01/GHMC01.1.fsa_nt.gz
```

```bash
PREFIX="Dynamena_pumila_WSRUv1"
```

**Transcripts**

Extract transcripts and predict ORFs.

```bash
zcat GHMC01.1.fsa_nt.gz \
  | sed -e "s/>/>${PREFIX}___/" \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.transcripts.fasta
```

Get gene to transcript name mappings.

```bash
zcat GHMC01.1.fsa_nt.gz \
  | grep '>' | sed -e "s/>/${PREFIX}___/" \
  | awk '{print $1"\t"$1}' \
  > ${PREFIX}.transcripts.gene2trans.txt
```

Predict ORFs in transcripts using `transdecoder`.

```bash
./run_Transdecoder.sh
```

Get best ORF per "gene"

```bash
grep '>' ${PREFIX}.transcripts.fasta.transdecoder.pep \
  | sed -e 's/.* \([^~]*\)~~\([^ ]*\) .*score=\([^ ,]*\).*/\1\t\2\t\3/' \
  | sort -k1,1 -k3,3nr \
  | awk -F'\t' '!seen[$1]++{print $2}' \
  > ${PREFIX}.transcripts.fasta.transdecoder.topPerGene.txt
```

**CDS**

```bash
cat ${PREFIX}.transcripts.fasta.transdecoder.cds \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  | ~/scripts/grepf_fasta.py -f ${PREFIX}.transcripts.fasta.transdecoder.topPerGene.txt \
  > ${PREFIX}.transcripts.cds.fna
```

**PEP**

```bash
cat ${PREFIX}.transcripts.fasta.transdecoder.pep \
  | sed -e 's/*$//' \
  | reformat.sh fastawrap=0 iupacToN=f ignorejunk=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  | ~/scripts/grepf_fasta.py -f ${PREFIX}.transcripts.fasta.transdecoder.topPerGene.txt \
  > ${PREFIX}.transcripts.pep.faa
```

## *Millepora alcicornis* PRUSv1

> Hydrozoa

Paper: NA

From: Puerto Rico: Lajas, Enrique Reef

Version: 1.0

Website: NCBI:  GFAS00000000.1

```bash
wget https://sra-download.ncbi.nlm.nih.gov/traces/wgs03/wgs_aux/GF/AS/GFAS01/GFAS01.1.fsa_nt.gz
wget https://sra-download.ncbi.nlm.nih.gov/traces/wgs03/wgs_aux/GF/AS/GFAS01/GFAS01.2.fsa_nt.gz
```

```bash
PREFIX="Millepora_alcicornis_PRUSv1"
```

**Transcripts**

Extract transcripts and predict ORFs.

```bash
zcat GFAS01.1.fsa_nt.gz GFAS01.2.fsa_nt.gz \
  | sed -e "s/>/>${PREFIX}___/" -e 's/gb|//' -e 's/|//' \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.transcripts.fasta
```

Get gene to transcript name mappings.

```bash
zcat GFAS01.1.fsa_nt.gz GFAS01.2.fsa_nt.gz \
  | grep '>' | sed -e "s/>/${PREFIX}___/" -e 's/gb|//' -e 's/|//' \
  | awk '{print $5"\t"$1}' | sed -e 's/_i[^\t]*//' \
  > ${PREFIX}.transcripts.gene2trans.txt
```

Predict ORFs in transcripts using `transdecoder`.

```bash
./run_Transdecoder.sh
```

Get best ORF per "gene"

```bash
grep '>' ${PREFIX}.transcripts.fasta.transdecoder.pep \
  | sed -e 's/.* \([^~]*\)~~\([^ ]*\) .*score=\([^ ,]*\).*/\1\t\2\t\3/' \
  | sort -k1,1 -k3,3nr \
  | awk -F'\t' '!seen[$1]++{print $2}' \
  > ${PREFIX}.transcripts.fasta.transdecoder.topPerGene.txt
```

**CDS**

```bash
cat ${PREFIX}.transcripts.fasta.transdecoder.cds \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  | ~/scripts/grepf_fasta.py -f ${PREFIX}.transcripts.fasta.transdecoder.topPerGene.txt \
  > ${PREFIX}.transcripts.cds.fna
```

**PEP**

```bash
cat ${PREFIX}.transcripts.fasta.transdecoder.pep \
  | sed -e 's/*$//' \
  | reformat.sh fastawrap=0 iupacToN=f ignorejunk=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  | ~/scripts/grepf_fasta.py -f ${PREFIX}.transcripts.fasta.transdecoder.topPerGene.txt \
  > ${PREFIX}.transcripts.pep.faa
```

## *Millepora squarrosa* PRUSv1

> Hydrozoa

Paper: NA

From: NA

Version: 1.0

Website: NCBI:  GFGU00000000.1

```bash
wget https://sra-download.ncbi.nlm.nih.gov/traces/wgs03/wgs_aux/GF/GU/GFGU01/GFGU01.1.fsa_nt.gz
```

```bash
PREFIX="Millepora_squarrosa_PRUSv1"
```

**Transcripts**

Extract transcripts and predict ORFs.

```bash
zcat GFGU01.1.fsa_nt.gz \
  | sed -e "s/>/>${PREFIX}___/" \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.transcripts.fasta
```

Get gene to transcript name mappings.

```bash
zcat GFGU01.1.fsa_nt.gz \
  | grep '>' | sed -e "s/>/${PREFIX}___/" \
  | awk '{print $5"\t"$1}' | sed -e 's/_i[^\t]*//' \
  > ${PREFIX}.transcripts.gene2trans.txt
```

Predict ORFs in transcripts using `transdecoder`.

```bash
./run_Transdecoder.sh
```

Get best ORF per "gene"

```bash
grep '>' ${PREFIX}.transcripts.fasta.transdecoder.pep \
  | sed -e 's/.* \([^~]*\)~~\([^ ]*\) .*score=\([^ ,]*\).*/\1\t\2\t\3/' \
  | sort -k1,1 -k3,3nr \
  | awk -F'\t' '!seen[$1]++{print $2}' \
  > ${PREFIX}.transcripts.fasta.transdecoder.topPerGene.txt
```

**CDS**

```bash
cat ${PREFIX}.transcripts.fasta.transdecoder.cds \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  | ~/scripts/grepf_fasta.py -f ${PREFIX}.transcripts.fasta.transdecoder.topPerGene.txt \
  > ${PREFIX}.transcripts.cds.fna
```

**PEP**

```bash
cat ${PREFIX}.transcripts.fasta.transdecoder.pep \
  | sed -e 's/*$//' \
  | reformat.sh fastawrap=0 iupacToN=f ignorejunk=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  | ~/scripts/grepf_fasta.py -f ${PREFIX}.transcripts.fasta.transdecoder.topPerGene.txt \
  > ${PREFIX}.transcripts.pep.faa
```

## *Millepora complanata* PRUSv1

> Hydrozoa

Paper: NA

From: Puerto Rico: Enrique Key, La Parguera, Lajas

Version: 1.0

Website: NCBI:  GFGT00000000.1

```bash
wget https://sra-download.ncbi.nlm.nih.gov/traces/wgs03/wgs_aux/GF/GT/GFGT01/GFGT01.1.fsa_nt.gz
```

```bash
PREFIX="Millepora_complanata_PRUSv1"
```

**Transcripts**

Extract transcripts and predict ORFs.

```bash
zcat GFGT01.1.fsa_nt.gz \
  | sed -e "s/>/>${PREFIX}___/" \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.transcripts.fasta
```

Get gene to transcript name mappings.

```bash
zcat GFGT01.1.fsa_nt.gz \
  | grep '>' | sed -e "s/>/${PREFIX}___/" \
  | awk '{print $5"\t"$1}' | sed -e 's/_i[^\t]*//' \
  > ${PREFIX}.transcripts.gene2trans.txt
```

Predict ORFs in transcripts using `transdecoder`.

```bash
./run_Transdecoder.sh
```

Get best ORF per "gene"

```bash
grep '>' ${PREFIX}.transcripts.fasta.transdecoder.pep \
  | sed -e 's/.* \([^~]*\)~~\([^ ]*\) .*score=\([^ ,]*\).*/\1\t\2\t\3/' \
  | sort -k1,1 -k3,3nr \
  | awk -F'\t' '!seen[$1]++{print $2}' \
  > ${PREFIX}.transcripts.fasta.transdecoder.topPerGene.txt
```

**CDS**

```bash
cat ${PREFIX}.transcripts.fasta.transdecoder.cds \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  | ~/scripts/grepf_fasta.py -f ${PREFIX}.transcripts.fasta.transdecoder.topPerGene.txt \
  > ${PREFIX}.transcripts.cds.fna
```

**PEP**

```bash
cat ${PREFIX}.transcripts.fasta.transdecoder.pep \
  | sed -e 's/*$//' \
  | reformat.sh fastawrap=0 iupacToN=f ignorejunk=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  | ~/scripts/grepf_fasta.py -f ${PREFIX}.transcripts.fasta.transdecoder.topPerGene.txt \
  > ${PREFIX}.transcripts.pep.faa
```

## *Millepora* sp. RR-2016 PRUSv1

> Hydrozoa

Paper: NA

From: Puerto Rico: Mata Seca Key, La Parguera, Lajas

Version: 1.0

Website: NCBI: GFGV00000000.1

```bash
wget https://sra-download.ncbi.nlm.nih.gov/traces/wgs03/wgs_aux/GF/GV/GFGV01/GFGV01.1.fsa_nt.gz
```

```bash
PREFIX="Millepora_sp_RR-2016_PRUSv1"
```

**Transcripts**

Extract transcripts and predict ORFs.

```bash
zcat GFGV01.1.fsa_nt.gz \
  | sed -e "s/>/>${PREFIX}___/" \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.transcripts.fasta
```

Get gene to transcript name mappings.

```bash
zcat GFGV01.1.fsa_nt.gz \
  | grep '>' | sed -e "s/>/${PREFIX}___/" \
  | awk '{print $6"\t"$1}' | sed -e 's/_i[^\t]*//' \
  > ${PREFIX}.transcripts.gene2trans.txt
```

Predict ORFs in transcripts using `transdecoder`.

```bash
./run_Transdecoder.sh
```

Get best ORF per "gene"

```bash
grep '>' ${PREFIX}.transcripts.fasta.transdecoder.pep \
  | sed -e 's/.* \([^~]*\)~~\([^ ]*\) .*score=\([^ ,]*\).*/\1\t\2\t\3/' \
  | sort -k1,1 -k3,3nr \
  | awk -F'\t' '!seen[$1]++{print $2}' \
  > ${PREFIX}.transcripts.fasta.transdecoder.topPerGene.txt
```

**CDS**

```bash
cat ${PREFIX}.transcripts.fasta.transdecoder.cds \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  | ~/scripts/grepf_fasta.py -f ${PREFIX}.transcripts.fasta.transdecoder.topPerGene.txt \
  > ${PREFIX}.transcripts.cds.fna
```

**PEP**

```bash
cat ${PREFIX}.transcripts.fasta.transdecoder.pep \
  | sed -e 's/*$//' \
  | reformat.sh fastawrap=0 iupacToN=f ignorejunk=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  | ~/scripts/grepf_fasta.py -f ${PREFIX}.transcripts.fasta.transdecoder.topPerGene.txt \
  > ${PREFIX}.transcripts.pep.faa
```

## *Podocoryna carnea* NANAv1

> Hydrozoa

Paper: NA

From: NA

Version: 1.0

Website: NCBI: GBEH00000000.1

```bash
wget https://sra-download.ncbi.nlm.nih.gov/traces/wgs03/wgs_aux/GB/EH/GBEH01/GBEH01.1.fsa_nt.gz
```

```bash
PREFIX="Podocoryna_carnea_NANAv1"
```

**Transcripts**

Extract transcripts and predict ORFs.

```bash
zcat GBEH01.1.fsa_nt.gz \
  | sed -e "s/>/>${PREFIX}___/" \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.transcripts.fasta
```

Get gene to transcript name mappings.

```bash
zcat GBEH01.1.fsa_nt.gz \
  | grep '>' | sed -e "s/>/${PREFIX}___/" \
  | awk '{print $5"\t"$1}' | sed -e 's/_i[^\t]*//' \
  > ${PREFIX}.transcripts.gene2trans.txt
```

Predict ORFs in transcripts using `transdecoder`.

```bash
./run_Transdecoder.sh
```

Get best ORF per "gene"

```bash
grep '>' ${PREFIX}.transcripts.fasta.transdecoder.pep \
  | sed -e 's/.* \([^~]*\)~~\([^ ]*\) .*score=\([^ ,]*\).*/\1\t\2\t\3/' \
  | awk -F'\t' '{split($1,a,"_"); print a[1]"_"a[2]"\t"$2"\t"$3}' \
  | sort -k1,1 -k3,3nr \
  | awk -F'\t' '!seen[$1]++{print $2}' \
  > ${PREFIX}.transcripts.fasta.transdecoder.topPerGene.txt
```

**CDS**

```bash
cat ${PREFIX}.transcripts.fasta.transdecoder.cds \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  | ~/scripts/grepf_fasta.py -f ${PREFIX}.transcripts.fasta.transdecoder.topPerGene.txt \
  > ${PREFIX}.transcripts.cds.fna
```

**PEP**

```bash
cat ${PREFIX}.transcripts.fasta.transdecoder.pep \
  | sed -e 's/*$//' \
  | reformat.sh fastawrap=0 iupacToN=f ignorejunk=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  | ~/scripts/grepf_fasta.py -f ${PREFIX}.transcripts.fasta.transdecoder.topPerGene.txt \
  > ${PREFIX}.transcripts.pep.faa
```

## *Porpita porpita* YOJPv1

> Hydrozoa

Paper: NA

From: Okinawa,Yomitan,Nagahama

Version: 1.0

Website: NCBI: GHBA00000000.1

```bash
wget https://sra-download.ncbi.nlm.nih.gov/traces/wgs03/wgs_aux/GH/BA/GHBA01/GHBA01.1.fsa_nt.gz
```

```bash
PREFIX="Porpita_porpita_YOJPv1"
```

**Transcripts**

Extract transcripts and predict ORFs.

```bash
zcat GHBA01.1.fsa_nt.gz \
  | sed -e "s/>/>${PREFIX}___/" \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.transcripts.fasta
```

Get gene to transcript name mappings.

```bash
zcat GHBA01.1.fsa_nt.gz \
  | grep '>' | sed -e "s/>/${PREFIX}___/" \
  | awk '{print $7"\t"$1}' | sed -e 's/,//' \
  > ${PREFIX}.transcripts.gene2trans.txt
```

Predict ORFs in transcripts using `transdecoder`.

```bash
./run_Transdecoder.sh
```

Get best ORF per "gene"

```bash
grep '>' ${PREFIX}.transcripts.fasta.transdecoder.pep \
  | sed -e 's/.* \([^~]*\)~~\([^ ]*\) .*score=\([^ ,]*\).*/\1\t\2\t\3/' \
  | sort -k1,1 -k3,3nr \
  | awk -F'\t' '!seen[$1]++{print $2}' \
  > ${PREFIX}.transcripts.fasta.transdecoder.topPerGene.txt
```

**CDS**

```bash
cat ${PREFIX}.transcripts.fasta.transdecoder.cds \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  | ~/scripts/grepf_fasta.py -f ${PREFIX}.transcripts.fasta.transdecoder.topPerGene.txt \
  > ${PREFIX}.transcripts.cds.fna
```

**PEP**

```bash
cat ${PREFIX}.transcripts.fasta.transdecoder.pep \
  | sed -e 's/*$//' \
  | reformat.sh fastawrap=0 iupacToN=f ignorejunk=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  | ~/scripts/grepf_fasta.py -f ${PREFIX}.transcripts.fasta.transdecoder.topPerGene.txt \
  > ${PREFIX}.transcripts.pep.faa
```

## *Rhodactis indosinensis* NATWv1

> Corallimorpharia

Paper: NA

From: Taiwan

Version: 1.0

Website: NCBI: GELO00000000.1

```bash
wget https://sra-download.ncbi.nlm.nih.gov/traces/wgs03/wgs_aux/GE/LO/GELO01/GELO01.1.fsa_nt.gz
```

```bash
PREFIX="Rhodactis_indosinensis_NATWv1"
```

**Transcripts**

Extract transcripts and predict ORFs.

```bash
zcat GELO01.1.fsa_nt.gz \
  | sed -e "s/>/>${PREFIX}___/" -e 's/gb|//' -e 's/|//' \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.transcripts.fasta
```

Get gene to transcript name mappings.

```bash
zcat GELO01.1.fsa_nt.gz \
  | grep '>' | sed -e "s/>/${PREFIX}___/" -e 's/gb|//' -e 's/|//' \
  | awk '{print $5"\t"$1}' | sed -e 's/_seq[^\t]*//' \
  > ${PREFIX}.transcripts.gene2trans.txt
```

Predict ORFs in transcripts using `transdecoder`.

```bash
./run_Transdecoder.sh
```

Get best ORF per "gene"

```bash
grep '>' ${PREFIX}.transcripts.fasta.transdecoder.pep \
  | sed -e 's/.* \([^~]*\)~~\([^ ]*\) .*score=\([^ ,]*\).*/\1\t\2\t\3/' \
  | sort -k1,1 -k3,3nr \
  | awk -F'\t' '!seen[$1]++{print $2}' \
  > ${PREFIX}.transcripts.fasta.transdecoder.topPerGene.txt
```

**CDS**

```bash
cat ${PREFIX}.transcripts.fasta.transdecoder.cds \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  | ~/scripts/grepf_fasta.py -f ${PREFIX}.transcripts.fasta.transdecoder.topPerGene.txt \
  > ${PREFIX}.transcripts.cds.fna
```

**PEP**

```bash
cat ${PREFIX}.transcripts.fasta.transdecoder.pep \
  | sed -e 's/*$//' \
  | reformat.sh fastawrap=0 iupacToN=f ignorejunk=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  | ~/scripts/grepf_fasta.py -f ${PREFIX}.transcripts.fasta.transdecoder.topPerGene.txt \
  > ${PREFIX}.transcripts.pep.faa
```

## *Ricordea yuma* NAAUv1

> Corallimorpharia

Paper: NA

From: Great Barrier Reef, Australia

Version: 1.0

Website: NCBI: GELN00000000.1

```bash
wget https://sra-download.ncbi.nlm.nih.gov/traces/wgs03/wgs_aux/GE/LN/GELN01/GELN01.1.fsa_nt.gz
```

```bash
PREFIX="Ricordea_yuma_NAAUv1"
```

**Transcripts**

Extract transcripts and predict ORFs.

```bash
zcat GELN01.1.fsa_nt.gz \
  | sed -e "s/>/>${PREFIX}___/" -e 's/gb|//' -e 's/|//' \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.transcripts.fasta
```

Get gene to transcript name mappings.

```bash
zcat GELN01.1.fsa_nt.gz \
  | grep '>' | sed -e "s/>/${PREFIX}___/" -e 's/gb|//' -e 's/|//' \
  | awk '{print $5"\t"$1}' | sed -e 's/_seq[^\t]*//' \
  > ${PREFIX}.transcripts.gene2trans.txt
```

Predict ORFs in transcripts using `transdecoder`.

```bash
./run_Transdecoder.sh
```

Get best ORF per "gene"

```bash
grep '>' ${PREFIX}.transcripts.fasta.transdecoder.pep \
  | sed -e 's/.* \([^~]*\)~~\([^ ]*\) .*score=\([^ ,]*\).*/\1\t\2\t\3/' \
  | sort -k1,1 -k3,3nr \
  | awk -F'\t' '!seen[$1]++{print $2}' \
  > ${PREFIX}.transcripts.fasta.transdecoder.topPerGene.txt
```

**CDS**

```bash
cat ${PREFIX}.transcripts.fasta.transdecoder.cds \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  | ~/scripts/grepf_fasta.py -f ${PREFIX}.transcripts.fasta.transdecoder.topPerGene.txt \
  > ${PREFIX}.transcripts.cds.fna
```

**PEP**

```bash
cat ${PREFIX}.transcripts.fasta.transdecoder.pep \
  | sed -e 's/*$//' \
  | reformat.sh fastawrap=0 iupacToN=f ignorejunk=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  | ~/scripts/grepf_fasta.py -f ${PREFIX}.transcripts.fasta.transdecoder.topPerGene.txt \
  > ${PREFIX}.transcripts.pep.faa
```

## *Turritopsis* sp. SK-2016 Xv1

> Hydrozoa

Paper: https://doi.org/10.2108/zs150186

From: Tanabe Bay, Japan (33.68′E, 135.36′N)

Version: 1.0

Website: NCBI:  IAAF00000000.1

```bash
wget https://sra-download.ncbi.nlm.nih.gov/traces/wgs03/wgs_aux/IA/AF/IAAF01/IAAF01.1.fsa_nt.gz
```

```bash
PREFIX="Turritopsis_sp_SK-2016_TBJPv1"
```

**Transcripts**

Extract transcripts and predict ORFs.

```bash
zcat IAAF01.1.fsa_nt.gz \
  | sed -e "s/>/>${PREFIX}___/" \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.transcripts.fasta
```

Get gene to transcript name mappings.

```bash
zcat IAAF01.1.fsa_nt.gz \
  | grep '>' | sed -e "s/>/${PREFIX}___/" \
  | awk '{print $8"\t"$1}' | sed -e 's/,//' \
  > ${PREFIX}.transcripts.gene2trans.txt
```

Predict ORFs in transcripts using `transdecoder`.

```bash
./run_Transdecoder.sh
```

Get best ORF per "gene"

```bash
grep '>' ${PREFIX}.transcripts.fasta.transdecoder.pep \
  | sed -e 's/.* \([^~]*\)~~\([^ ]*\) .*score=\([^ ,]*\).*/\1\t\2\t\3/' \
  | awk -F'\t' '{split($1,a,"_"); print a[1]"_"a[2]"\t"$2"\t"$3}' \
  | sort -k1,1 -k3,3nr \
  | awk -F'\t' '!seen[$1]++{print $2}' \
  > ${PREFIX}.transcripts.fasta.transdecoder.topPerGene.txt
```

**CDS**

```bash
cat ${PREFIX}.transcripts.fasta.transdecoder.cds \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  | ~/scripts/grepf_fasta.py -f ${PREFIX}.transcripts.fasta.transdecoder.topPerGene.txt \
  > ${PREFIX}.transcripts.cds.fna
```

**PEP**

```bash
cat ${PREFIX}.transcripts.fasta.transdecoder.pep \
  | sed -e 's/*$//' \
  | reformat.sh fastawrap=0 iupacToN=f ignorejunk=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  | ~/scripts/grepf_fasta.py -f ${PREFIX}.transcripts.fasta.transdecoder.topPerGene.txt \
  > ${PREFIX}.transcripts.pep.faa
```

## *Velella velella* Xv1

> Hydrozoa

Paper: NA

From: NA

Version: 1.0

Website: NCBI: GHAZ00000000.1

```bash
wget https://sra-download.ncbi.nlm.nih.gov/traces/wgs03/wgs_aux/GH/AZ/GHAZ01/GHAZ01.1.fsa_nt.gz
```

```bash
PREFIX="Velella_velella_OKJPv1"
```

**Transcripts**

Extract transcripts and predict ORFs.

```bash
zcat GHAZ01.1.fsa_nt.gz \
  | sed -e "s/>/>${PREFIX}___/" \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.transcripts.fasta
```

Get gene to transcript name mappings.

```bash
zcat GHAZ01.1.fsa_nt.gz \
  | grep '>' | sed -e "s/>/${PREFIX}___/" \
  | awk '{print $7"\t"$1}' | sed -e 's/,//' \
  > ${PREFIX}.transcripts.gene2trans.txt
```

Predict ORFs in transcripts using `transdecoder`.

```bash
./run_Transdecoder.sh
```

Get best ORF per "gene"

```bash
grep '>' ${PREFIX}.transcripts.fasta.transdecoder.pep \
  | sed -e 's/.* \([^~]*\)~~\([^ ]*\) .*score=\([^ ,]*\).*/\1\t\2\t\3/' \
  | sort -k1,1 -k3,3nr \
  | awk -F'\t' '!seen[$1]++{print $2}' \
  > ${PREFIX}.transcripts.fasta.transdecoder.topPerGene.txt
```

**CDS**

```bash
cat ${PREFIX}.transcripts.fasta.transdecoder.cds \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  | ~/scripts/grepf_fasta.py -f ${PREFIX}.transcripts.fasta.transdecoder.topPerGene.txt \
  > ${PREFIX}.transcripts.cds.fna
```

**PEP**

```bash
cat ${PREFIX}.transcripts.fasta.transdecoder.pep \
  | sed -e 's/*$//' \
  | reformat.sh fastawrap=0 iupacToN=f ignorejunk=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  | ~/scripts/grepf_fasta.py -f ${PREFIX}.transcripts.fasta.transdecoder.topPerGene.txt \
  > ${PREFIX}.transcripts.pep.faa
```

## *Physalia physalis* OKJPv1

> Hydrozoa

Paper: NA

From: Japan: Okinawa,Yomitan,Nagahama

Version: 1.0

Website: NCBI: GHBB00000000.1

```bash
wget https://sra-download.ncbi.nlm.nih.gov/traces/wgs03/wgs_aux/GH/BB/GHBB01/GHBB01.1.fsa_nt.gz
```

```bash
PREFIX="Physalia_physalis_OKJPv1"
```

**Transcripts**

Extract transcripts and predict ORFs.

```bash
zcat GHBB01.1.fsa_nt.gz \
  | sed -e "s/>/>${PREFIX}___/" \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.transcripts.fasta
```

Get gene to transcript name mappings.

```bash
zcat GHBB01.1.fsa_nt.gz \
  | grep '>' | sed -e "s/>/${PREFIX}___/" \
  | awk '{print $7"\t"$1}' | sed -e 's/,//' \
  > ${PREFIX}.transcripts.gene2trans.txt
```

Predict ORFs in transcripts using `transdecoder`.

```bash
./run_Transdecoder.sh
```

Get best ORF per "gene"

```bash
grep '>' ${PREFIX}.transcripts.fasta.transdecoder.pep \
  | sed -e 's/.* \([^~]*\)~~\([^ ]*\) .*score=\([^ ,]*\).*/\1\t\2\t\3/' \
  | sort -k1,1 -k3,3nr \
  | awk -F'\t' '!seen[$1]++{print $2}' \
  > ${PREFIX}.transcripts.fasta.transdecoder.topPerGene.txt
```

**CDS**

```bash
cat ${PREFIX}.transcripts.fasta.transdecoder.cds \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  | ~/scripts/grepf_fasta.py -f ${PREFIX}.transcripts.fasta.transdecoder.topPerGene.txt \
  > ${PREFIX}.transcripts.cds.fna
```

**PEP**

```bash
cat ${PREFIX}.transcripts.fasta.transdecoder.pep \
  | sed -e 's/*$//' \
  | reformat.sh fastawrap=0 iupacToN=f ignorejunk=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  | ~/scripts/grepf_fasta.py -f ${PREFIX}.transcripts.fasta.transdecoder.topPerGene.txt \
  > ${PREFIX}.transcripts.pep.faa
```

# QC

### `run_name_count.sh`

Get the number of uniq gene names across the GFF, PEP, and CDS files. Make sure this number matches the GFF file.



### `run_seq_names.sh`

Check that the `PREFIX` used for each genome dataset is consistent.



---

---

---

---

---

### `run_extract_proteins_01.sh` and `run_extract_proteins_02.sh`

```bash
export PATH="$PATH:/home/timothy/programs/seqkit_v0.15.0"
```

---

**Acropora_millepora_JSIDv2.1**

```bash
PREFIX="Acropora_millepora_JSIDv2.1"
~/scripts/add_value_to_table.py \
    -i <(seqkit fx2tab ${PREFIX}.genes.pep.faa.OLD \
           | cut -f1,2 | sed -e 's/*/X/g') \
    -a <(seqkit fx2tab ${PREFIX}.genes.pep.faa.NEW \
           | cut -f1,2 | sed -e 's/*/X/g') \
  | awk -F'\t' '$2!=$3' | wc -l
# 543 

diff <(seqkit fx2tab ${PREFIX}.genes.pep.faa.OLD \
           | cut -f1,2 | sed -e 's/*/X/g') \
     <(seqkit fx2tab ${PREFIX}.genes.pep.faa.NEW \
           | cut -f1,2 | sed -e 's/*/X/g')
```

Seems to be a lot (543 ) of proteins that differ by a single additional amino acid (internal insertion in the NEW proteins). Does not change the frame so the majority of the protein match between the two sets.

---

**Montastrea_cavernosa_FGUSv2**

```bash
PREFIX="Montastrea_cavernosa_FGUSv2"
~/scripts/add_value_to_table.py \
    -i <(seqkit fx2tab ${PREFIX}.genes.pep.faa.OLD \
           | cut -f1,2 | sed -e 's/*/X/g') \
    -a <(seqkit fx2tab ${PREFIX}.genes.pep.faa.NEW \
           | cut -f1,2 | sed -e 's/*/X/g') \
  | awk -F'\t' '$2!=$3' | wc -l
# 2

diff <(seqkit fx2tab ${PREFIX}.genes.pep.faa.OLD \
           | cut -f1,2 | sed -e 's/*/X/g') \
     <(seqkit fx2tab ${PREFIX}.genes.pep.faa.NEW \
           | cut -f1,2 | sed -e 's/*/X/g')
```

Differs only by two proteins. Looks like one of the exons (towards the ends of the genes) have weird coords that cause there ends of the proteins to have different frames and thus different sequences. 

---

**Montipora_capitata_KBHIv3**

```bash
PREFIX="Montipora_capitata_KBHIv3"
~/scripts/add_value_to_table.py \
    -i <(seqkit fx2tab ${PREFIX}.genes.pep.faa.OLD \
           | cut -f1,2 | sed -e 's/*/X/g') \
    -a <(seqkit fx2tab ${PREFIX}.genes.pep.faa.NEW \
           | cut -f1,2 | sed -e 's/*/X/g') \
  | awk -F'\t' '$2!=$3' | wc -l
# 826

~/scripts/add_value_to_table.py \
    -i <(seqkit fx2tab ${PREFIX}.genes.pep.faa.OLD \
           | cut -f1,2 | sed -e 's/*/X/g' -e 's/\t[ML]/\t/') \
    -a <(seqkit fx2tab ${PREFIX}.genes.pep.faa.NEW \
           | cut -f1,2 | sed -e 's/*/X/g' -e 's/\t[ML]/\t/') \
  | awk -F'\t' '$2!=$3' | wc -l
# 0
# Differences are from the use of alternative start codons ("L") in the original OLD set vs. the NEW set. 
```

---

**Montipora_capitata_WTHIv1.1**

```bash
PREFIX="Montipora_capitata_WTHIv1.1"
GFF=$(ls -1 ../*_genomes/${PREFIX}/genome_assembly/${PREFIX}.genes.gff3)
GEN=$(ls -1 ../*_genomes/${PREFIX}/genome_assembly/${PREFIX}.assembly.fasta)
gffread -g ${GEN} <(grep -v 'xfSc0014868' ${GFF}) -S \
        -x >(reformat.sh fastawrap=0 iupacToN=f ignorejunk=t \
                         touppercase=t trimreaddescription=t \
                         in=stdin.fa out=stdout.fa \
               | seqkit sort --line-width 0 \
               > ${PREFIX}.genes.cds.fna.NEW) \
        -y >(reformat.sh fastawrap=0 iupacToN=f ignorejunk=t \
                         touppercase=t trimreaddescription=t \
                         in=stdin.fa out=stdout.fa \
               | seqkit sort --line-width 0 \
               > ${PREFIX}.genes.pep.faa.NEW)


~/scripts/add_value_to_table.py \
    -i <(seqkit fx2tab ${PREFIX}.genes.pep.faa.OLD \
           | cut -f1,2 | sed -e 's/*/X/g') \
    -a <(seqkit fx2tab ${PREFIX}.genes.pep.faa.NEW \
           | cut -f1,2 | sed -e 's/*/X/g') \
  | awk -F'\t' '$2!=$3' | wc -l
# 49

diff <(seqkit fx2tab ${PREFIX}.genes.pep.faa.OLD \
           | cut -f1,2 | sed -e 's/*/X/g') \
     <(seqkit fx2tab ${PREFIX}.genes.pep.faa.NEW \
           | cut -f1,2 | sed -e 's/*/X/g')
```

After removing the one scaffold from the GFF file that caused `gffread` to terminate early, we only have 49 genes that differ between the protein sets. This appears (visually) to be because some ambiguous residues (`X`) are not ambiguous in the OLD set. 

There is also one extra protein in the OLD set that is not in the NEW set - likely a result of one of the scaffolds being missing from the assembly that is listed in the GFF file.

---

**Montipora_sp1_aff_capitata_ULFMv1**

```bash
PREFIX="Montipora_sp1_aff_capitata_ULFMv1"
~/scripts/add_value_to_table.py \
    -i <(seqkit fx2tab ${PREFIX}.genes.pep.faa.OLD \
           | cut -f1,2 | sed -e 's/*/X/g') \
    -a <(seqkit fx2tab ${PREFIX}.genes.pep.faa.NEW \
           | cut -f1,2 | sed -e 's/*/X/g') \
  | awk -F'\t' '$2!=$3' | wc -l
# 5093

~/scripts/add_value_to_table.py \
    -i <(seqkit fx2tab ${PREFIX}.genes.pep.faa.OLD \
           | cut -f1,2 | sed -e 's/*/X/g' -e 's/\tX/\t/' -e 's/X$//') \
    -a <(seqkit fx2tab ${PREFIX}.genes.pep.faa.NEW \
           | cut -f1,2 | sed -e 's/*/X/g' -e 's/\tX/\t/' -e 's/X$//') \
  | awk -F'\t' '$2!=$3' | wc -l
# 220 - After removing leading and trailing 'X' characters

diff <(seqkit fx2tab ${PREFIX}.genes.pep.faa.OLD \
           | cut -f1,2 | sed -e 's/*/X/g' -e 's/\tX/\t/' -e 's/X$//') \
     <(seqkit fx2tab ${PREFIX}.genes.pep.faa.NEW \
           | cut -f1,2 | sed -e 's/*/X/g' -e 's/\tX/\t/' -e 's/X$//')
```

Looks to be a problem with how terminal residues (in partially predicted proteins) are treated. Looks like an extra residue is being predicted by `gffread` that is not included in the OLD set - often these are `X` residues, however there are 220 cases where it is an actual amino acid. 

---

**Orbicella_faveolata_FLUSv1**

```bash
PREFIX="Orbicella_faveolata_FLUSv1"
~/scripts/add_value_to_table.py \
    -i <(seqkit fx2tab ${PREFIX}.genes.pep.faa.OLD \
           | cut -f1,2 | sed -e 's/*/X/g') \
    -a <(seqkit fx2tab ${PREFIX}.genes.pep.faa.NEW \
           | cut -f1,2 | sed -e 's/*/X/g') \
  | awk -F'\t' '$2!=$3' | wc -l
# 93

diff <(seqkit fx2tab ${PREFIX}.genes.pep.faa.OLD \
           | cut -f1,2 | sed -e 's/*/X/g') \
     <(seqkit fx2tab ${PREFIX}.genes.pep.faa.NEW \
           | cut -f1,2 | sed -e 's/*/X/g')
```

Looks to 93 cases (visually estimated) where an ambiguous residue `X` has been predicted in the OLD set and a non-ambiguous residue has been predicted in the NEW set. Only one protein that was <5% shorter (but still in same frame) as protein in other set. 

---

**Platygyra_sinensis_SDTHv1**

```bash
PREFIX="Platygyra_sinensis_SDTHv1"
~/scripts/add_value_to_table.py \
    -i <(seqkit fx2tab ${PREFIX}.genes.pep.faa.OLD \
           | cut -f1,2 | sed -e 's/*/X/g') \
    -a <(seqkit fx2tab ${PREFIX}.genes.pep.faa.NEW \
           | cut -f1,2 | sed -e 's/*/X/g') \
  | awk -F'\t' '$2!=$3' | wc -l
# 136

~/scripts/add_value_to_table.py \
    -i <(seqkit fx2tab ${PREFIX}.genes.pep.faa.OLD \
           | cut -f1,2 | sed -e 's/*/X/g' -e 's/\t[ML]/\t/') \
    -a <(seqkit fx2tab ${PREFIX}.genes.pep.faa.NEW \
           | cut -f1,2 | sed -e 's/*/X/g' -e 's/\t[ML]/\t/') \
  | awk -F'\t' '$2!=$3' | wc -l
# 0
# Differences are from the use of alternative start codons ("L") in the original OLD set vs. the NEW set.
```

---

**Pocillopora_acuta_KBHIv2**

```bash
PREFIX="Pocillopora_acuta_KBHIv2"
~/scripts/add_value_to_table.py \
    -i <(seqkit fx2tab ${PREFIX}.genes.pep.faa.OLD \
           | cut -f1,2 | sed -e 's/*/X/g') \
    -a <(seqkit fx2tab ${PREFIX}.genes.pep.faa.NEW \
           | cut -f1,2 | sed -e 's/*/X/g') \
  | awk -F'\t' '$2!=$3' | wc -l
# 500

~/scripts/add_value_to_table.py \
    -i <(seqkit fx2tab ${PREFIX}.genes.pep.faa.OLD \
           | cut -f1,2 | sed -e 's/*/X/g' -e 's/\t[ML]/\t/') \
    -a <(seqkit fx2tab ${PREFIX}.genes.pep.faa.NEW \
           | cut -f1,2 | sed -e 's/*/X/g' -e 's/\t[ML]/\t/') \
  | awk -F'\t' '$2!=$3' | wc -l
# 0
# Differences are from the use of alternative start codons ("L") in the original OLD set vs. the NEW set.
```

---

**Pocillopora_acuta_LBIDv1**

```bash
PREFIX="Pocillopora_acuta_LBIDv1"
~/scripts/add_value_to_table.py \
    -i <(seqkit fx2tab ${PREFIX}.genes.pep.faa.OLD \
           | cut -f1,2 | sed -e 's/*/X/g') \
    -a <(seqkit fx2tab ${PREFIX}.genes.pep.faa.NEW \
           | cut -f1,2 | sed -e 's/*/X/g') \
  | awk -F'\t' '$2!=$3' | wc -l
# 93

diff <(seqkit fx2tab ${PREFIX}.genes.pep.faa.OLD \
           | cut -f1,2 | sed -e 's/*/X/g') \
     <(seqkit fx2tab ${PREFIX}.genes.pep.faa.NEW \
           | cut -f1,2 | sed -e 's/*/X/g')
```

Looks to be differences in how internal or terminal ambiguous characters (`X`) are handled. OLD set has predicted non-ambiguous residue where as NEW has ambiguous character. 

---

**Pocillopora_meandrina_KBHIv1**

```bash
PREFIX="Pocillopora_meandrina_KBHIv1"
~/scripts/add_value_to_table.py \
    -i <(seqkit fx2tab ${PREFIX}.genes.pep.faa.OLD \
           | cut -f1,2 | sed -e 's/*/X/g') \
    -a <(seqkit fx2tab ${PREFIX}.genes.pep.faa.NEW \
           | cut -f1,2 | sed -e 's/*/X/g') \
  | awk -F'\t' '$2!=$3' | wc -l
# 334

~/scripts/add_value_to_table.py \
    -i <(seqkit fx2tab ${PREFIX}.genes.pep.faa.OLD \
           | cut -f1,2 | sed -e 's/*/X/g' -e 's/\t[ML]/\t/') \
    -a <(seqkit fx2tab ${PREFIX}.genes.pep.faa.NEW \
           | cut -f1,2 | sed -e 's/*/X/g' -e 's/\t[ML]/\t/') \
  | awk -F'\t' '$2!=$3' | wc -l
# 0
# Differences are from the use of alternative start codons ("L") in the original OLD set vs. the NEW set.
```

---

**Pocillopora_verrucosa_RSSAv1**

```bash
PREFIX="Pocillopora_verrucosa_RSSAv1"
~/scripts/add_value_to_table.py \
    -i <(seqkit fx2tab ${PREFIX}.genes.pep.faa.OLD \
           | cut -f1,2 | sed -e 's/*/X/g') \
    -a <(seqkit fx2tab ${PREFIX}.genes.pep.faa.NEW \
           | cut -f1,2 | sed -e 's/*/X/g') \
  | awk -F'\t' '$2!=$3' | wc -l
# 2353

diff <(seqkit fx2tab ${PREFIX}.genes.pep.faa.OLD \
           | cut -f1,2 | sed -e 's/*/X/g') \
     <(seqkit fx2tab ${PREFIX}.genes.pep.faa.NEW \
           | cut -f1,2 | sed -e 's/*/X/g')
```

Looks to be differences in how internal or terminal ambiguous characters (`X`) are handled. Some of the proteins in the OLD set have extra internal ambiguous residues compared with the NEW set. There are also a limited number of proteins have appear to be completely different between the two sets. Likely a problem with the frame of the GFF record. 

**Large number of proteins that are different between sets**



---

**Porites_compressa_KBHIv1**

```bash
PREFIX="Porites_compressa_KBHIv1"
~/scripts/add_value_to_table.py \
    -i <(seqkit fx2tab ${PREFIX}.genes.pep.faa.OLD \
           | cut -f1,2 | sed -e 's/*/X/g') \
    -a <(seqkit fx2tab ${PREFIX}.genes.pep.faa.NEW \
           | cut -f1,2 | sed -e 's/*/X/g') \
  | awk -F'\t' '$2!=$3' | wc -l
# 611

~/scripts/add_value_to_table.py \
    -i <(seqkit fx2tab ${PREFIX}.genes.pep.faa.OLD \
           | cut -f1,2 | sed -e 's/*/X/g' -e 's/\t[ML]/\t/') \
    -a <(seqkit fx2tab ${PREFIX}.genes.pep.faa.NEW \
           | cut -f1,2 | sed -e 's/*/X/g' -e 's/\t[ML]/\t/') \
  | awk -F'\t' '$2!=$3' | wc -l
# 0
# Differences are from the use of alternative start codons ("L") in the original OLD set vs. the NEW set.
```

---

**Porites_lutea_OIAUv1.1**

```bash
PREFIX="Porites_lutea_OIAUv1.1"
~/scripts/add_value_to_table.py \
    -i <(seqkit fx2tab ${PREFIX}.genes.pep.faa.OLD \
           | cut -f1,2 | sed -e 's/*/X/g') \
    -a <(seqkit fx2tab ${PREFIX}.genes.pep.faa.NEW \
           | cut -f1,2 | sed -e 's/*/X/g') \
  | awk -F'\t' '$2!=$3' | wc -l
# 64

diff <(seqkit fx2tab ${PREFIX}.genes.pep.faa.OLD \
           | cut -f1,2 | sed -e 's/*/X/g') \
     <(seqkit fx2tab ${PREFIX}.genes.pep.faa.NEW \
           | cut -f1,2 | sed -e 's/*/X/g')
```

Looks to be differences in how internal or terminal ambiguous characters (`X`) are handled. OLD set has predicted non-ambiguous residue where as NEW has ambiguous character. 

---

**Stylophora_pistillata_GAJOv1**

```bash
PREFIX="Stylophora_pistillata_GAJOv1"
~/scripts/add_value_to_table.py \
    -i <(seqkit fx2tab ${PREFIX}.genes.pep.faa.OLD \
           | cut -f1,2 | sed -e 's/*/X/g') \
    -a <(seqkit fx2tab ${PREFIX}.genes.pep.faa.NEW \
           | cut -f1,2 | sed -e 's/*/X/g') \
  | awk -F'\t' '$2!=$3' | wc -l
# 115

diff <(seqkit fx2tab ${PREFIX}.genes.pep.faa.OLD \
           | cut -f1,2 | sed -e 's/*/X/g') \
     <(seqkit fx2tab ${PREFIX}.genes.pep.faa.NEW \
           | cut -f1,2 | sed -e 's/*/X/g')
```

Looks to be differences in how internal or terminal ambiguous characters (`X`) are handled. OLD set has predicted non-ambiguous residue where as NEW has ambiguous character. Only seven (short) proteins that are massively different between the two sets (i.e., in a different frame). 

---

**Alatina_alata_BNNLv1**

Missing `gff` file.

---

**Amphimedon_queenslandica_HIAUv1.1**

```bash
PREFIX="Amphimedon_queenslandica_HIAUv1.1"
~/scripts/add_value_to_table.py \
    -i <(seqkit fx2tab ${PREFIX}.genes.pep.faa.OLD \
           | cut -f1,2 | sed -e 's/*$//' -e 's/*/X/g') \
    -a <(seqkit fx2tab ${PREFIX}.genes.pep.faa.NEW \
           | cut -f1,2 | sed -e 's/*$//' -e 's/*/X/g') \
  | awk -F'\t' '$2!=$3' | wc -l
# 121

diff <(seqkit fx2tab ${PREFIX}.genes.pep.faa.OLD \
           | cut -f1,2 | sed -e 's/*/X/g') \
     <(seqkit fx2tab ${PREFIX}.genes.pep.faa.NEW \
           | cut -f1,2 | sed -e 's/*/X/g')
```

Looks to be differences in how internal ambiguous characters (`X`) are handled or very minor (one or two per protein) mismatch differences. **Didnt extensively manually check proteins**

---

**Amplexidiscus_fenestrafer_MSMEv1**

```bash
PREFIX="Amplexidiscus_fenestrafer_MSMEv1"
~/scripts/add_value_to_table.py \
    -i <(seqkit fx2tab ${PREFIX}.genes.pep.faa.OLD \
           | cut -f1,2 | sed -e 's/*/X/g') \
    -a <(seqkit fx2tab ${PREFIX}.genes.pep.faa.NEW \
           | cut -f1,2 | sed -e 's/*/X/g') \
  | awk -F'\t' '$2!=$3' | wc -l
# 166

diff <(seqkit fx2tab ${PREFIX}.genes.pep.faa.OLD \
           | cut -f1,2 | sed -e 's/*/X/g') \
     <(seqkit fx2tab ${PREFIX}.genes.pep.faa.NEW \
           | cut -f1,2 | sed -e 's/*/X/g')
```

Looks to be differences in how internal or terminal ambiguous characters (`X`) are handled. Some of the proteins in the OLD set have extra internal ambiguous residues compared with the NEW set. There is also a significant number of (mostly small) proteins that appear to be completely different between the two sets. Likely a problem with the frame of the GFF record. 

**Large number of proteins that are different between sets**

---

**Aurelia_aurita_ALTOv1**





**Aurelia_aurita_RSPOv1**







**Calvadosia_cruxmelitensis_CRGBv3.2**

Missing `off` file.





---

**Clytia_hemispaerica_NANAv1**







---

**Dendronephthya_gigantea_SJKRv1**

```bash
PREFIX="Dendronephthya_gigantea_SJKRv1"
~/scripts/add_value_to_table.py \
    -i <(seqkit fx2tab ${PREFIX}.genes.pep.faa.OLD \
           | cut -f1,2 | sed -e 's/*/X/g') \
    -a <(seqkit fx2tab ${PREFIX}.genes.pep.faa.NEW \
           | cut -f1,2 | sed -e 's/*/X/g') \
  | awk -F'\t' '$2!=$3' | wc -l
# 191

diff <(seqkit fx2tab ${PREFIX}.genes.pep.faa.OLD \
           | cut -f1,2 | sed -e 's/*/X/g') \
     <(seqkit fx2tab ${PREFIX}.genes.pep.faa.NEW \
           | cut -f1,2 | sed -e 's/*/X/g')
```

Looks to be differences in how internal or terminal ambiguous characters (`X`) are handled. Some of the proteins in the OLD set have extra internal ambiguous residues compared with the NEW set. **Didnt extensively manually check proteins**

---

**Discosoma_sp_MSMEv1**

```bash
PREFIX="Discosoma_sp_MSMEv1"
~/scripts/add_value_to_table.py \
    -i <(seqkit fx2tab ${PREFIX}.genes.pep.faa.OLD \
           | cut -f1,2 | sed -e 's/*/X/g') \
    -a <(seqkit fx2tab ${PREFIX}.genes.pep.faa.NEW \
           | cut -f1,2 | sed -e 's/*/X/g') \
  | awk -F'\t' '$2!=$3' | wc -l
# 250

diff <(seqkit fx2tab ${PREFIX}.genes.pep.faa.OLD \
           | cut -f1,2 | sed -e 's/*/X/g') \
     <(seqkit fx2tab ${PREFIX}.genes.pep.faa.NEW \
           | cut -f1,2 | sed -e 's/*/X/g')
```

Looks to be differences in how internal or terminal ambiguous characters (`X`) are handled. Some of the proteins in the OLD set have extra internal ambiguous residues compared with the NEW set. There is also a significant number of (mostly small) proteins that appear to be completely different between the two sets. Likely a problem with the frame of the GFF record. There are also a limited number of proteins that have deletions of single amino acids. 

**Large number of proteins that are different between sets**



---

**Ephydatia_muelleri_BCCAv1**

```bash
PREFIX="Ephydatia_muelleri_BCCAv1"
~/scripts/add_value_to_table.py \
    -i <(seqkit fx2tab ${PREFIX}.genes.pep.faa.OLD \
           | cut -f1,2 | sed -e 's/*/X/g') \
    -a <(seqkit fx2tab ${PREFIX}.genes.pep.faa.NEW \
           | cut -f1,2 | sed -e 's/*/X/g') \
  | awk -F'\t' '$2!=$3' | wc -l
# 3

diff <(seqkit fx2tab ${PREFIX}.genes.pep.faa.OLD \
           | cut -f1,2 | sed -e 's/*/X/g') \
     <(seqkit fx2tab ${PREFIX}.genes.pep.faa.NEW \
           | cut -f1,2 | sed -e 's/*/X/g')
```

Differences are because of `*` characters being predicted next to long shared homo-*-regions.

No major changes to proteins or frame.

---

**Exaiptasia_diaphana_CC7v1**

```bash
PREFIX="Exaiptasia_diaphana_CC7v1"
~/scripts/add_value_to_table.py \
    -i <(seqkit fx2tab ${PREFIX}.genes.pep.faa.OLD \
           | cut -f1,2 | sed -e 's/*/X/g') \
    -a <(seqkit fx2tab ${PREFIX}.genes.pep.faa.NEW \
           | cut -f1,2 | sed -e 's/*/X/g') \
  | awk -F'\t' '$2!=$3' | wc -l
# 15760

diff <(seqkit fx2tab ${PREFIX}.genes.pep.faa.OLD \
           | cut -f1,2 | sed -e 's/*/X/g') \
     <(seqkit fx2tab ${PREFIX}.genes.pep.faa.NEW \
           | cut -f1,2 | sed -e 's/*/X/g')
```

**Lots of indwells between the two sets.** 

Generally this does not significantly change the sequence of the encoded protein.



---

**Hydra_vulgaris_MIJPv3**

```bash
PREFIX="Hydra_vulgaris_MIJPv3"
~/scripts/add_value_to_table.py \
    -i <(seqkit fx2tab ${PREFIX}.genes.pep.faa.OLD \
           | cut -f1,2 | sed -e 's/*/X/g') \
    -a <(seqkit fx2tab ${PREFIX}.genes.pep.faa.NEW \
           | cut -f1,2 | sed -e 's/*/X/g') \
  | awk -F'\t' '$2!=$3' | wc -l
# 455

diff <(seqkit fx2tab ${PREFIX}.genes.pep.faa.OLD \
           | cut -f1,2 | sed -e 's/*/X/g') \
     <(seqkit fx2tab ${PREFIX}.genes.pep.faa.NEW \
           | cut -f1,2 | sed -e 's/*/X/g')
```

Looks to be differences in how internal or terminal stop characters (`*`) are handled. Some of the proteins in the OLD set have extra internal ambiguous residues compared with the NEW set. There are also a limited number of proteins that have deletions of single amino acids. 

**Large number of proteins that are different between sets**



---

**Mnemiopsis_leidyi_WHUSv2.2**

```bash
PREFIX="Mnemiopsis_leidyi_WHUSv2.2"
~/scripts/add_value_to_table.py \
    -i <(seqkit fx2tab ${PREFIX}.genes.pep.faa.OLD \
           | cut -f1,2 | sed -e 's/*/X/g') \
    -a <(seqkit fx2tab ${PREFIX}.genes.pep.faa.NEW \
           | cut -f1,2 | sed -e 's/*/X/g') \
  | awk -F'\t' '$2!=$3' | wc -l
# 502

diff <(seqkit fx2tab ${PREFIX}.genes.pep.faa.OLD \
           | cut -f1,2 | sed -e 's/*/X/g') \
     <(seqkit fx2tab ${PREFIX}.genes.pep.faa.NEW \
           | cut -f1,2 | sed -e 's/*/X/g')
```

Proteins are completely different. Looks to be major issues with frame of proteins for this subset.

**Didnt extensively manually check proteins**



---

**Monosiga_brevicollis_CCBMv1**

```bash
PREFIX="Monosiga_brevicollis_CCBMv1"
~/scripts/add_value_to_table.py \
    -i <(seqkit fx2tab ${PREFIX}.genes.pep.faa.OLD \
           | cut -f1,2 | sed -e 's/*/X/g') \
    -a <(seqkit fx2tab ${PREFIX}.genes.pep.faa.NEW \
           | cut -f1,2 | sed -e 's/*/X/g') \
  | awk -F'\t' '$2!=$3' | wc -l
# 32

diff <(seqkit fx2tab ${PREFIX}.genes.pep.faa.OLD \
           | cut -f1,2 | sed -e 's/*/X/g') \
     <(seqkit fx2tab ${PREFIX}.genes.pep.faa.NEW \
           | cut -f1,2 | sed -e 's/*/X/g')
```

Difference is because of internal stop codons that in the OLD set are 'W' characters. Also appears to be come terminal characters that are in the OLD but not NEW proteins.

**Didnt extensively manually check proteins**



---

**Paramuricea_clavata_CTESv1**

```bash
PREFIX="Paramuricea_clavata_CTESv1"
~/scripts/add_value_to_table.py \
    -i <(seqkit fx2tab ${PREFIX}.genes.pep.faa.OLD \
           | cut -f1,2 | sed -e 's/*/X/g') \
    -a <(seqkit fx2tab ${PREFIX}.genes.pep.faa.NEW \
           | cut -f1,2 | sed -e 's/*/X/g') \
  | awk -F'\t' '$2!=$3' | wc -l
# 932

diff <(seqkit fx2tab ${PREFIX}.genes.pep.faa.OLD \
           | cut -f1,2 | sed -e 's/*/X/g') \
     <(seqkit fx2tab ${PREFIX}.genes.pep.faa.NEW \
           | cut -f1,2 | sed -e 's/*/X/g')
```

Proteins generally appear to be the same (i.e., no major frame shifts) but with some differences such as missing terminal  and start amino acids. Sometimes the difference in length can be large but still results in highly similar proteins. The OLD always appear to be longer than the NEW.

**Didnt extensively manually check proteins**



---

**Renilla_muelleri_FLUSv1**

```bash
PREFIX="Renilla_muelleri_FLUSv1"
~/scripts/add_value_to_table.py \
    -i <(seqkit fx2tab ${PREFIX}.genes.pep.faa.OLD \
           | cut -f1,2 | sed -e 's/*/X/g') \
    -a <(seqkit fx2tab ${PREFIX}.genes.pep.faa.NEW \
           | cut -f1,2 | sed -e 's/*/X/g') \
  | awk -F'\t' '$2!=$3' | wc -l
# 2

diff <(seqkit fx2tab ${PREFIX}.genes.pep.faa.OLD \
           | cut -f1,2 | sed -e 's/*/X/g') \
     <(seqkit fx2tab ${PREFIX}.genes.pep.faa.NEW \
           | cut -f1,2 | sed -e 's/*/X/g')
```

Difference in the two proteins is because of mismatches around homo-X regions. Is 'X' in OLD but characterized ('A', 'R', etc.) in NEW.



---

**Scolanthus_callimorphus_CRFRv1**

```bash
PREFIX="Scolanthus_callimorphus_CRFRv1"
~/scripts/add_value_to_table.py \
    -i <(seqkit fx2tab ${PREFIX}.genes.pep.faa.OLD \
           | cut -f1,2 | sed -e 's/*/X/g') \
    -a <(seqkit fx2tab ${PREFIX}.genes.pep.faa.NEW \
           | cut -f1,2 | sed -e 's/*/X/g') \
  | awk -F'\t' '$2!=$3' | wc -l
# 2029

diff <(seqkit fx2tab ${PREFIX}.genes.pep.faa.OLD \
           | cut -f1,2 | sed -e 's/*/X/g') \
     <(seqkit fx2tab ${PREFIX}.genes.pep.faa.NEW \
           | cut -f1,2 | sed -e 's/*/X/g')
```

Proteins generally appear to be the same (i.e., no major frame shifts) but with some differences such as missing terminal  and start amino acids. Sometimes the difference in length can be large but still results in highly similar proteins. The OLD always appear to be longer than the NEW.

**Didnt extensively manually check proteins**



---

**Stylissa_carteri_RSSAv1**

```bash
PREFIX="Stylissa_carteri_RSSAv1"
~/scripts/add_value_to_table.py \
    -i <(seqkit fx2tab ${PREFIX}.genes.pep.faa.OLD \
           | cut -f1,2 | sed -e 's/*/X/g') \
    -a <(seqkit fx2tab ${PREFIX}.genes.pep.faa.NEW \
           | cut -f1,2 | sed -e 's/*/X/g') \
  | awk -F'\t' '$2!=$3' | wc -l
# 44723

diff <(seqkit fx2tab ${PREFIX}.genes.pep.faa.OLD \
           | cut -f1,2 | sed -e 's/*/X/g') \
     <(seqkit fx2tab ${PREFIX}.genes.pep.faa.NEW \
           | cut -f1,2 | sed -e 's/*/X/g')
```

**`gffread` appears to fail with an unhelpful error in this dataset!** Not possible to really assess differences.



---

**Trichoplax_adhaerens_RSSAv1**

```bash
PREFIX="Trichoplax_adhaerens_RSSAv1"
~/scripts/add_value_to_table.py \
    -i <(seqkit fx2tab ${PREFIX}.genes.pep.faa.OLD \
           | cut -f1,2 | sed -e 's/*/X/g') \
    -a <(seqkit fx2tab ${PREFIX}.genes.pep.faa.NEW \
           | cut -f1,2 | sed -e 's/*/X/g') \
  | awk -F'\t' '$2!=$3' | wc -l
# 2829

diff <(seqkit fx2tab ${PREFIX}.genes.pep.faa.OLD \
           | cut -f1,2 | sed -e 's/*/X/g') \
     <(seqkit fx2tab ${PREFIX}.genes.pep.faa.NEW \
           | cut -f1,2 | sed -e 's/*/X/g')
```

Proteins are completely different. Looks to be major issues with frame of proteins for this subset.

**Didnt extensively manually check proteins**



---

**Xenia_sp_CTEAv1**

```bash
PREFIX="Xenia_sp_CTEAv1"
~/scripts/add_value_to_table.py \
    -i <(seqkit fx2tab ${PREFIX}.genes.pep.faa.OLD \
           | cut -f1,2 | sed -e 's/*/X/g') \
    -a <(seqkit fx2tab ${PREFIX}.genes.pep.faa.NEW \
           | cut -f1,2 | sed -e 's/*/X/g') \
  | awk -F'\t' '$2!=$3' | wc -l
# 35

diff <(seqkit fx2tab ${PREFIX}.genes.pep.faa.OLD \
           | cut -f1,2 | sed -e 's/*/X/g') \
     <(seqkit fx2tab ${PREFIX}.genes.pep.faa.NEW \
           | cut -f1,2 | sed -e 's/*/X/g')
```

Proteins generally appear to be the same (i.e., no major frame shifts) but with some differences such as missing terminal  and start amino acids. The OLD always appear to be longer than the NEW. Some of the old also have `-` characters as their initial residue. 



---

**Xestospongia_testudinaria_RSSAv1**

```bash
PREFIX="Xestospongia_testudinaria_RSSAv1"
~/scripts/add_value_to_table.py \
    -i <(seqkit fx2tab ${PREFIX}.genes.pep.faa.OLD \
           | cut -f1,2 | sed -e 's/*/X/g') \
    -a <(seqkit fx2tab ${PREFIX}.genes.pep.faa.NEW \
           | cut -f1,2 | sed -e 's/*/X/g') \
  | awk -F'\t' '$2!=$3' | wc -l
# 52575

diff <(seqkit fx2tab ${PREFIX}.genes.pep.faa.OLD \
           | cut -f1,2 | sed -e 's/*/X/g') \
     <(seqkit fx2tab ${PREFIX}.genes.pep.faa.NEW \
           | cut -f1,2 | sed -e 's/*/X/g')
```

**`gffread` appears to fail with an unhelpful error in this dataset!** Not possible to really assess differences.



---





















```bash
./run_name_count-transcripts.sh
./run_seq_names-transcripts.sh
```

































## *----Millepora squarrosa* PRUSv1

> Hydrozoa

Paper: NA

From: NA

Version: 1.0

Website: NCBI:  GFGU00000000.1

```bash
wget 
```

```bash
PREFIX="Millepora_squarrosa_PRUSv1"
```

**Transcripts**

Extract transcripts and predict ORFs.

```bash
zcat GFGU01.1.fsa_nt.gz \
  | sed -e "s/>/>${PREFIX}___/" \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  > ${PREFIX}.transcripts.fasta
```

Get gene to transcript name mappings.

```bash
zcat GFGU01.1.fsa_nt.gz \
  | grep '>' | sed -e "s/>/${PREFIX}___/" \
  | awk '{print $5"\t"$1}' | sed -e 's/_i[^\t]*//' \
  > ${PREFIX}.transcripts.gene2trans.txt
```

Predict ORFs in transcripts using `transdecoder`.

```bash
./run_Transdecoder.sh
```

Get best ORF per "gene"

```bash
grep '>' ${PREFIX}.transcripts.fasta.transdecoder.pep \
  | sed -e 's/.* \([^~]*\)~~\([^ ]*\) .*score=\([^ ,]*\).*/\1\t\2\t\3/' \
  | sort -k1,1 -k3,3nr \
  | awk -F'\t' '!seen[$1]++{print $2}' \
  > ${PREFIX}.transcripts.fasta.transdecoder.topPerGene.txt
```

**CDS**

```bash
cat ${PREFIX}.transcripts.fasta.transdecoder.cds \
  | reformat.sh fastawrap=0 iupacToN=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  | ~/scripts/grepf_fasta.py -f ${PREFIX}.transcripts.fasta.transdecoder.topPerGene.txt \
  > ${PREFIX}.transcripts.cds.fna
```

**PEP**

```bash
cat ${PREFIX}.transcripts.fasta.transdecoder.pep \
  | sed -e 's/*$//' \
  | reformat.sh fastawrap=0 iupacToN=f ignorejunk=t touppercase=t trimreaddescription=t in=stdin.fa out=stdout.fa \
  | ~/scripts/grepf_fasta.py -f ${PREFIX}.transcripts.fasta.transdecoder.topPerGene.txt \
  > ${PREFIX}.transcripts.pep.faa
```

## 













