#!/bin/bash
set -euo pipefail

STRAIN_NAME="GCA_000328405.1_ASM32840v1_genomic"
BAKTA_DIR="/mnt/d/WSL/disk/projects/VP1/bakta_annotations"
GENOME_DIR="/mnt/d/WSL/disk/projects/VP1/clean_genome"
VFDB="/mnt/g/wsl/database/vfdb/VFDB_setA_nt.fas"
TEST_OUTPUT="test_virulence_by_coord.csv"

rm -f tmp_*.tsv tmp_*.fna $TEST_OUTPUT

# 1. 从Bakta注释表提取：基因ID+坐标+链方向（关键！用坐标提序列）
echo "🔍 筛选毒力基因并提取坐标..."
grep -E "virulence|toxin|tdh|trh|hemolysin" "$BAKTA_DIR/$STRAIN_NAME/$STRAIN_NAME.tsv" | \
awk -F "\t" '$1 != "" && $3 != "" && $4 != "" {
    print $1 "\t" $3 "\t" $4 "\t" $6 "\t" $2 "\t" $5  # 基因ID→start→end→链→基因名→功能
}' > tmp_cand_coord.tsv
cand_count=$(wc -l tmp_cand_coord.tsv | awk '{print $1}')
echo "筛选到带坐标的候选基因数：$cand_count"

# 2. 用坐标从基因组提取序列（bedtools更精准，自动处理正负链）
echo -e "\n📑 用坐标提取基因序列..."
# 先安装bedtools（bakta环境没有的话自动安装）
if ! command -v bedtools &> /dev/null; then
    echo "⚠️  安装bedtools（坐标提取序列工具）..."
    conda install -c bioconda bedtools -y
fi

# 转换Bakta坐标为bed格式（chrom→start→end→geneID→score→strand）
awk -F "\t" '{print "chr" "\t" $2-1 "\t" $3 "\t" $1 "\t" "0" "\t" $4}' tmp_cand_coord.tsv > tmp_cand.bed

# 提取序列（bedtools getfasta，100%匹配坐标）
bedtools getfasta -fi "$GENOME_DIR/$STRAIN_NAME.fna" \
                  -bed tmp_cand.bed \
                  -fo tmp_seq.fna \
                  -name  # 序列名用基因ID，方便后续匹配

seq_count=$(grep -c "^>" tmp_seq.fna)
echo "成功提取序列数：$seq_count"

# 3. BLAST比对VFDB验证
echo -e "\n🧬 BLAST比对VFDB..."
blastn -query tmp_seq.fna \
       -db "$VFDB" \
       -outfmt "6 qseqid pident evalue" \
       -evalue 1e-10 \
       -perc_identity 70 \
       -num_threads 4 > tmp_blast.tsv
blast_count=$(wc -l tmp_blast.tsv | awk '{print $1}')
echo "BLAST有效命中数：$blast_count"

# 4. 整合结果
echo -e "\n📊 生成最终结果..."
echo "菌株名,基因ID,基因名,毒力功能,BLAST相似度,E值" > $TEST_OUTPUT
join -t $'\t' <(awk -F "\t" '{print $1 "\t" $5 "\t" $6}' tmp_cand_coord.tsv | sort) \
              <(sort tmp_blast.tsv) | \
awk -v strain="$STRAIN_NAME" -F "\t" '{print strain "," $1 "," $2 "," $3 "," $5 "," $6}' >> $TEST_OUTPUT

echo -e "\n🎉 测试完成！"
echo "最终筛选到真毒力基因数：$(($(wc -l $TEST_OUTPUT | awk '{print $1}') - 1))"
cat $TEST_OUTPUT
