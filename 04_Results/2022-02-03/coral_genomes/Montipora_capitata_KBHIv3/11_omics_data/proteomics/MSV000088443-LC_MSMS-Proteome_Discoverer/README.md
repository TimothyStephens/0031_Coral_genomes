# Reformat proteomic data

Reformat proteomic data to use a standard sample naming scheme and so that missing (zero) values are replaced with actual zeros (`0`). Also remove any columns that we won't really need downstream (i.e., only protein name and normalized abundance).

## Setup analysis directory

Setup bash environment.

```bash
conda activate py27
```



## Protein results

### Get list of column names

List the column names so that we can check their order.

```bash
cut -f4,21-34 Montipora_capitata_KBHIv3.protein_abundance.tsv | head -n 1 | sed -e 's/\t/\n/g'
```

> Accession
>
> "Abundances (Normalized): F1: 128N, Sample, T1A"
>
> "Abundances (Normalized): F1: 129N, Sample, T1A"
>
> "Abundances (Normalized): F1: 129C, Sample, T1H"
>
> "Abundances (Normalized): F1: 130N, Sample, T1H"
>
> "Abundances (Normalized): F1: 130C, Sample, T3A"
>
> "Abundances (Normalized): F1: 131N, Sample, T3A"
>
> "Abundances (Normalized): F1: 131C, Sample, T3H"
>
> "Abundances (Normalized): F1: 132N, Sample, T3H"
>
> "Abundances (Normalized): F1: 132C, Sample, T5A"
>
> "Abundances (Normalized): F1: 133N, Sample, T5A"
>
> "Abundances (Normalized): F1: 133C, Sample, T5H"
>
> "Abundances (Normalized): F1: 134N, Sample, T5H"
>
> "Abundances (Normalized): F1: 126, Sample, Wild"
>
> "Abundances (Normalized): F1: 127N, Sample, Wild"

### Rename columns

| OLD                                                 | NEW                 |
| --------------------------------------------------- | ------------------- |
| Abundances (Normalized): F1: 128N,  Sample, T1A, 1  | MC-289_TP1-Amb_1603 |
| Abundances (Normalized): F1: 129N,  Sample, T1A, 2  | MC-289_TP1-Amb_1609 |
| Abundances (Normalized): F1: 129C,  Sample, T1H, 1  | MC-289_TP1-HiT_2878 |
| Abundances (Normalized): F1: 130N,  Sample, T1H, 2  | MC-289_TP1-HiT_2023 |
| Abundances (Normalized): F1: 130C,  Sample, T3A, 1  | MC-289_TP3-Amb_1595 |
| Abundances (Normalized): F1: 131N,  Sample, T3A, 2  | MC-289_TP3-Amb_2741 |
| Abundances (Normalized): F1: 131C,  Sample, T3H, 1  | MC-289_TP3-HiT_2183 |
| Abundances (Normalized): F1: 132N,  Sample, T3H, 2  | MC-289_TP3-HiT_2058 |
| Abundances (Normalized): F1: 132C,  Sample, T5A, 1  | MC-289_TP5-Amb_1721 |
| Abundances (Normalized): F1: 133N,  Sample, T5A, 2  | MC-289_TP5-Amb_2874 |
| Abundances (Normalized): F1: 133C,  Sample, T5H, 1  | MC-289_TP5-HiT_1341 |
| Abundances (Normalized): F1: 134N,  Sample, T5H, 2  | MC-289_TP5-HiT_2998 |
| Abundances (Normalized): F1: 126, Sample,  Wild, 1  | MC-289_Field_289-1  |
| Abundances (Normalized): F1: 127N,  Sample, Wild, 2 | MC-289_Field_289-2  |

Add new column headers to file.

```bash
echo -e 'Name\tMC-289_TP1-Amb_1603\tMC-289_TP1-Amb_1609\tMC-289_TP1-HiT_2878\tMC-289_TP1-HiT_2023\tMC-289_TP3-Amb_1595\tMC-289_TP3-Amb_2741\tMC-289_TP3-HiT_2183\tMC-289_TP3-HiT_2058\tMC-289_TP5-Amb_1721\tMC-289_TP5-Amb_2874\tMC-289_TP5-HiT_1341\tMC-289_TP5-HiT_2998\tMC-289_Field_289-1\tMC-289_Field_289-2' > Montipora_capitata_KBHIv3.protein_abundance.normalized.tsv
```

Add normalized abundance values to file.

```bash
cut -f4,21-34 Montipora_capitata_KBHIv3.protein_abundance.tsv \
  | awk 'BEGIN{OFS=FS="\t"} { for(i=2; i<=NF; i++){ if($i==""){$i=0} }; print }' \
  | awk 'NR>1' \
  >> Montipora_capitata_KBHIv3.protein_abundance.normalized.tsv
```



## Peptide results

### Get list of column names

List the column names so that we can check their order.

```bash
cut -f3,13,20-33 Montipora_capitata_KBHIv3.peptide_abundance.tsv | head -n 1 | sed -e 's/\t/\n/g'
```

> Sequence
>
> Positions in Master Proteins
>
> "Abundances (Normalized): F1: 128N, Sample, T1A"
>
> "Abundances (Normalized): F1: 129N, Sample, T1A"
>
> "Abundances (Normalized): F1: 129C, Sample, T1H"
>
> "Abundances (Normalized): F1: 130N, Sample, T1H"
>
> "Abundances (Normalized): F1: 130C, Sample, T3A"
>
> "Abundances (Normalized): F1: 131N, Sample, T3A"
>
> "Abundances (Normalized): F1: 131C, Sample, T3H"
>
> "Abundances (Normalized): F1: 132N, Sample, T3H"
>
> "Abundances (Normalized): F1: 132C, Sample, T5A"
>
> "Abundances (Normalized): F1: 133N, Sample, T5A"
>
> "Abundances (Normalized): F1: 133C, Sample, T5H"
>
> "Abundances (Normalized): F1: 134N, Sample, T5H"
>
> "Abundances (Normalized): F1: 126, Sample, Wild"
>
> "Abundances (Normalized): F1: 127N, Sample, Wild"

### Rename columns

| OLD                                                 | NEW                 |
| --------------------------------------------------- | ------------------- |
| Abundances (Normalized): F1: 128N,  Sample, T1A, 1  | MC-289_TP1-Amb_1603 |
| Abundances (Normalized): F1: 129N,  Sample, T1A, 2  | MC-289_TP1-Amb_1609 |
| Abundances (Normalized): F1: 129C,  Sample, T1H, 1  | MC-289_TP1-HiT_2878 |
| Abundances (Normalized): F1: 130N,  Sample, T1H, 2  | MC-289_TP1-HiT_2023 |
| Abundances (Normalized): F1: 130C,  Sample, T3A, 1  | MC-289_TP3-Amb_1595 |
| Abundances (Normalized): F1: 131N,  Sample, T3A, 2  | MC-289_TP3-Amb_2741 |
| Abundances (Normalized): F1: 131C,  Sample, T3H, 1  | MC-289_TP3-HiT_2183 |
| Abundances (Normalized): F1: 132N,  Sample, T3H, 2  | MC-289_TP3-HiT_2058 |
| Abundances (Normalized): F1: 132C,  Sample, T5A, 1  | MC-289_TP5-Amb_1721 |
| Abundances (Normalized): F1: 133N,  Sample, T5A, 2  | MC-289_TP5-Amb_2874 |
| Abundances (Normalized): F1: 133C,  Sample, T5H, 1  | MC-289_TP5-HiT_1341 |
| Abundances (Normalized): F1: 134N,  Sample, T5H, 2  | MC-289_TP5-HiT_2998 |
| Abundances (Normalized): F1: 126, Sample,  Wild, 1  | MC-289_Field_289-1  |
| Abundances (Normalized): F1: 127N,  Sample, Wild, 2 | MC-289_Field_289-2  |

Add new column headers to file.

```bash
echo -e 'Peptide_sequence\tPositions_in_master_proteins\tMC-289_TP1-Amb_1603\tMC-289_TP1-Amb_1609\tMC-289_TP1-HiT_2878\tMC-289_TP1-HiT_2023\tMC-289_TP3-Amb_1595\tMC-289_TP3-Amb_2741\tMC-289_TP3-HiT_2183\tMC-289_TP3-HiT_2058\tMC-289_TP5-Amb_1721\tMC-289_TP5-Amb_2874\tMC-289_TP5-HiT_1341\tMC-289_TP5-HiT_2998\tMC-289_Field_289-1\tMC-289_Field_289-2' > Montipora_capitata_KBHIv3.peptide_abundance.normalized.tsv
```

Add normalized abundance values to file.

```bash
cut -f3,13,20-33 Montipora_capitata_KBHIv3.peptide_abundance.tsv \
  | awk 'BEGIN{OFS=FS="\t"} { gsub("; ", ";", $2); gsub(" \\[", ":", $2); gsub("\\]", "", $2); for(i=3; i<=NF; i++){ if($i==""){$i=0} }; print }' \
  | awk 'NR>1' \
  >> Montipora_capitata_KBHIv3.peptide_abundance.normalized.tsv
```



NOTE: There can be multiple instances (rows) for the same peptide sequence if that particular fragment was detected with different modifications. ~2000 instances of this happening.

> Total peptides: 23,793
>
> Unique peptide sequences + abundances: 23,410
>
> Unique peptide sequences: 22,430
