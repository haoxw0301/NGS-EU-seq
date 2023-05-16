#!/bin/bash
#SBATCH --job-name=fc_L1
#SBATCH -p shenxhlab
#SBATCH --ntasks=5
#SBATCH -N 1

## pwd: /NFSdata01/haoxw/cellcylce_EU/
source /WORK/Samples/bio.sh

mkdir fc_L1
mkdir fc_gene

featureCounts -t exon -g gene_id \
    -a ~/haoxw/reference/mm10/L1xn/L1xn_FC_all.gtf \
    -o fc_L1/fc_L1.txt \
    -s 1 hisat2/*uniq.bam

featureCounts -t exon -g gene_id \
    -a ~/haoxw/reference/mm10/coding_gene/mm10.ncbiRefSeq.gtf \
    -o fc_gene/fc_gene.txt \
    -s 1 hisat2/*uniq.bam