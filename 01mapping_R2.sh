#!/bin/bash
for i in `ls ~/xlq/euSeq20201207fq/trim_galore/*eu*2.fq.gz`
do
    j=$(basename $i "_R2_val_2.fq.gz")
    u="NH:i:1"
    echo "#!/bin/bash" > ${j}.sh
    echo "#SBATCH --job-name=${j}" >> ${j}.sh
    echo "#SBATCH -p shenxhlab" >> ${j}.sh
    echo "#SBATCH --ntasks=5" >> ${j}.sh
    echo "#SBATCH -N 1" >> ${j}.sh
    echo "source /WORK/Samples/bio.sh " >> ${j}.sh
    echo "
    ## mm10
    hisat2 -x ~/haoxw/reference/mm10/hisat-index/mm10 \
        -U ${i} -S hisat2/${j}.sam
    samtools view -bS -F 4 hisat2/${j}.sam > hisat2/${j}.bam
    samtools sort hisat2/${j}.bam -o hisat2/${j}.bam
    samtools index hisat2/${j}.bam

    grep ${u} hisat2/${j}.sam > hisat2/${j}_Uniq.sam
    samtools view -H -S hisat2/${j}.sam > hisat2/${j}_header.sam
    cat hisat2/${j}_header.sam \
        hisat2/${j}_Uniq.sam > hisat2/${j}_uniq.sam
    samtools view -bS -F 4 hisat2/${j}_uniq.sam > hisat2/${j}_uniq.bam
    samtools sort hisat2/${j}_uniq.bam -o hisat2/${j}_uniq.bam
    samtools index hisat2/${j}_uniq.bam
    
    ## ERCC
    hisat2 -x ~/haoxw/reference/ERCC/hisat-index/ERCC \
        -U ${i} -S hisat2_ERCC/${j}.sam
    samtools view -bS -F 4 hisat2_ERCC/${j}.sam > hisat2_ERCC/${j}.bam
    samtools sort hisat2_ERCC/${j}.bam -o hisat2_ERCC/${j}.bam
    samtools index hisat2_ERCC/${j}.bam

    grep ${u} hisat2_ERCC/${j}.sam > hisat2_ERCC/${j}_Uniq.sam
    samtools view -H -S hisat2_ERCC/${j}.sam > hisat2_ERCC/${j}_header.sam
    cat hisat2_ERCC/${j}_header.sam \
        hisat2_ERCC/${j}_Uniq.sam > hisat2_ERCC/${j}_uniq.sam
    samtools view -bS -F 4 hisat2_ERCC/${j}_uniq.sam > hisat2_ERCC/${j}_uniq.bam
    samtools sort hisat2_ERCC/${j}_uniq.bam -o hisat2_ERCC/${j}_uniq.bam
    samtools index hisat2_ERCC/${j}_uniq.bam    
    " >> ${j}.sh

    sbatch ${j}.sh
done