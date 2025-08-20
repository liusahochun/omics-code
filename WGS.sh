#!/bin/bash
# 全基因组一键分析脚本
# 使用方法：./auto_wgs_analysis.sh

################### 配置区域 ###################
ref="pombe.fa"                    # 参考基因组
gtf="SP666.gtf"                   # 基因注释文件
threads=16                        # 并行线程数
output_dir="wgs_analysis_results" # 输出目录
##############################################

echo "🚀 开始全基因组一键分析..."
start_time=$(date +%s)

# 创建输出目录
mkdir -p $output_dir
cd $output_dir

# 检查输入文件
echo "📋 检查输入文件..."
[ ! -f "../$ref" ] && echo "❌ 参考基因组不存在" && exit 1
[ ! -f "../$gtf" ] && echo "❌ GTF文件不存在" && exit 1

# 链接输入文件
ln -sf ../*.fq.gz . 2>/dev/null
ln -sf ../$ref .
ln -sf ../$gtf .

################### 分析流程 ###################

# 1. 质控清洗
echo "🔍 步骤1: FASTQ数据质控清洗"
for fq1 in *_1.fq.gz; do
    sample=$(basename "$fq1" _1.fq.gz)
    echo "  处理样本: $sample"
    
    fastp --thread $threads \
        -i "$fq1" -o "${sample}_1.clean.fq.gz" \
        -I "${sample}_2.fq.gz" -O "${sample}_2.clean.fq.gz" \
        --html "${sample}.fastp.html" --json "${sample}.fastp.json"
done

# 2. 建立索引（如需要）
echo "📚 步骤2: 建立参考基因组索引"
[ ! -f "${ref}.bwt" ] && bwa index $ref
[ ! -f "${ref}.fai" ] && samtools faidx $ref
[ ! -f "${ref%.fa}.dict" ] && gatk CreateSequenceDictionary -R $ref

# 3. 序列比对
echo "🧬 步骤3: 序列比对与排序"
for fq1 in *_1.clean.fq.gz; do
    sample=$(basename "$fq1" _1.clean.fq.gz)
    echo "  比对样本: $sample"
    
    bwa mem -t $threads -R "@RG\tID:$sample\tSM:$sample\tPL:Illumina" \
        $ref "${sample}_1.clean.fq.gz" "${sample}_2.clean.fq.gz" \
        | samtools sort -@ $threads -o "${sample}.sorted.bam" -
    
    samtools index "${sample}.sorted.bam"
done

# 4. 标记重复序列
echo "🔖 步骤4: 标记PCR重复序列"
for bam in *.sorted.bam; do
    sample=$(basename "$bam" .sorted.bam)
    gatk MarkDuplicates -I "$bam" -O "${sample}.markdup.bam" -M "${sample}.metrics.txt"
    samtools index "${sample}.markdup.bam"
done

# 5. 变异检测
echo "🔎 步骤5: 变异检测"
for bam in *.markdup.bam; do
    sample=$(basename "$bam" .markdup.bam)
    echo "  检测样本: $sample"
    
    # 生成gVCF
    gatk HaplotypeCaller -R $ref -I "$bam" --emit-ref-confidence GVCF -O "${sample}.g.vcf"
    
    # 基因分型
    gatk GenotypeGVCFs -R $ref -V "${sample}.g.vcf" -O "${sample}.vcf"
    
    # 压缩索引
    bgzip -f "${sample}.vcf"
    tabix -p vcf "${sample}.vcf.gz"
done

# 6. 变异过滤
echo "🧹 步骤6: 变异过滤"
for vcf in *.vcf.gz; do
    sample=$(basename "$vcf" .vcf.gz)
    
    # SNP过滤
    gatk SelectVariants -V "$vcf" -select-type SNP -O "${sample}.snp.vcf.gz"
    gatk VariantFiltration -V "${sample}.snp.vcf.gz" \
        --filter-expression "QD < 2.0 || MQ < 40.0 || FS > 60.0 || SOR > 3.0" \
        --filter-name "SNP_FILTER" -O "${sample}.snp.filtered.vcf.gz"
    
    # INDEL过滤
    gatk SelectVariants -V "$vcf" -select-type INDEL -O "${sample}.indel.vcf.gz"
    gatk VariantFiltration -V "${sample}.indel.vcf.gz" \
        --filter-expression "QD < 2.0 || FS > 200.0 || SOR > 10.0" \
        --filter-name "INDEL_FILTER" -O "${sample}.indel.filtered.vcf.gz"
    
    # 合并结果
    gatk MergeVcfs -I "${sample}.snp.filtered.vcf.gz" -I "${sample}.indel.filtered.vcf.gz" \
        -O "${sample}.final.vcf.gz"
done

# 7. 变异注释
echo "🏷️  步骤7: 变异注释"
for vcf in *.final.vcf.gz; do
    sample=$(basename "$vcf" .final.vcf.gz)
    echo "  注释样本: $sample"
    
    # 使用ANNOVAR或VEP进行注释
    # 这里使用简化版的注释流程
    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%FILTER\n' "$vcf" \
        > "${sample}.annotated.tsv"
done

# 生成分析报告
echo "📊 生成分析报告..."
{
    echo "# 全基因组分析报告"
    echo "生成时间: $(date)"
    echo "分析样本数: $(ls *.fq.gz | wc -l)"
    echo "总变异数: $(bcftools view *.final.vcf.gz | grep -v '^#' | wc -l)"
} > "analysis_report.md"

end_time=$(date +%s)
echo "✅ 分析完成！耗时: $((end_time - start_time)) 秒"
echo "📁 结果保存在: $output_dir"