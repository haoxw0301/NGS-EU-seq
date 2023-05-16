#!/bin/bash
#SBATCH --job-name=total_rc
#SBATCH -p shenxhlab
#SBATCH --ntasks=5
#SBATCH -N 1
source /WORK/Samples/bio.sh

for i in `ls *unique.bam`
do 
    c=$(samtools view -c ${i})
    echo ${i} ${c} >> total_rc.txt
done 