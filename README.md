# L1 expresssion for NGS RNA-seq

## 1. mapping with hisat2

this program is used to map the reads to the reference genome and get the uniq reads.

```bash
sh 01mapping_R2.sh
```

## 2. calaulate the reads count with featureCounts

```bash
sbatch 02fc.sh
```

## 3. calculate the total reads count with samtools

```bash
cd hista2
sbatch total_rc.sh

cd ../hisat2_ERCC
sbatch total_rc.sh
```

## 4. calculate the repeatability of the replicates

```bash
cd fc_gene

Rscript repeatability.R # sbatch
```

## 5. calculate the L1 expression levels

```bash
cd fc_L1

Rscript L1_exp_level.R # sbatch
```
