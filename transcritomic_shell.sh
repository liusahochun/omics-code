#!/bin/bash
set -euo pipefail

# 并行任务数（建议设为 CPU 核心数的一半，比如16核机设为8）
THREADS=8

###################### 文件准备 ######################
ls -1 *1.fastq.gz > 1.txt
ls -1 *2.fastq.gz > 2.txt
paste 1.txt 2.txt > 3.txt
awk '{print $1"\t"$2"\t"substr($1,1,9)}' 3.txt > 4.txt

###################### Step1: fastp #################
echo ">>> Step1: fastp QC" $(date)

cat 4.txt | parallel -j $THREADS --colsep '\t' '
  fq1={1}; fq2={2}; sample={3};
  echo "start fastp for $sample" $(date);
  fastp -w '$THREADS' \
    -i $fq1 -I $fq2 \
    -o ${fq1%.fastq.gz}.clean.fastq.gz \
    -O ${fq2%.fastq.gz}.clean.fastq.gz \
    > ${sample}_fastp.log 2>&1;
  echo "end fastp for $sample" $(date);
'

###################### Step2: HISAT2 ################
echo ">>> Step2: HISAT2 alignment" $(date)

hisat2-build -p $THREADS SY14参考基因组.fa SY14-index

cat 4.txt | parallel -j $THREADS --colsep '\t' '
  fq1={1}; fq2={2}; sample={3};
  echo "start alignment for $sample" $(date);
  hisat2 -p '$THREADS' -q -x SY14-index \
    -1 ${fq1%.fastq.gz}.clean.fastq.gz \
    -2 ${fq2%.fastq.gz}.clean.fastq.gz \
    -S $sample.sam \
    > ${sample}_hisat2.log 2>&1;
  echo "end alignment for $sample" $(date);
'

###################### Step3: BAM ###################
echo ">>> Step3: SAM -> BAM & sort" $(date)

cat 4.txt | parallel -j $THREADS --colsep '\t' '
  sample={3};
  echo "start bam for $sample" $(date);
  samtools view -@ '$THREADS' -bS $sample.sam > $sample.bam;
  samtools sort -@ '$THREADS' -o $sample.sorted.bam $sample.bam;
  samtools index $sample.sorted.bam;
  rm $sample.sam $sample.bam;
  echo "end bam for $sample" $(date);
'

###################### Step4: featureCounts #########
echo ">>> Step4: featureCounts" $(date)

featureCounts -T $THREADS -t exon -g gene_id \
  -a genomic-budding.gtf -p \
  -o gene_counts.txt *.sorted.bam

echo ">>> ALL DONE" $(date)

