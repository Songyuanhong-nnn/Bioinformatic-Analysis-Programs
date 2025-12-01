#!/bin/bash
set -euo pipefail

STRAIN_NAME="GCA_000328405.1_ASM32840v1_genomic"
VFDB="/mnt/g/wsl/database/vfdb/VFDB_setA_nt.fas"
BAKTA_DIR="/mnt/d/WSL/disk/projects/VP1/bakta_annotations"
GENOME_DIR="/mnt/d/WSL/disk/projects/VP1/clean_genome"
TEST_OUTPUT="test_virulence_fixed.csv"

rm -f tmp_*.tsv tmp_*.fna $TEST_OUTPUT

# 1. 确认文件存在
echo "✅ 确认菌株文件..."
ls -l "$BAKTA_DIR/$STRAIN_NAME/$STRAIN_NAME.tsv"
ls -l "$GENOME_DIR/$STRAIN_NAME.fna"

# 2. 筛选毒力候选基因（只保留有基因ID的）
echo -e "\n🔍 筛选毒力基因..."
grep -E "virulence|toxin|tdh|trh|hemolysin" "$BAKTA_DIR/$STRAIN_NAME/$STRAIN_NAME.tsv" | \
awk -F "\t" '$1 != "" {print $1 "\t" $2 "\t" $5}' > tmp_cand.tsv
cand_count=$(wc -l tmp_cand.tsv | awk '{print $1}')
echo "筛选到候选基因数：$cand_count"

# 3. 优化序列提取（确保只提取有实际序列的基因）
echo -e "\n📑 提取候选基因序列..."
awk 'NR==FNR {
    if ($1 != "") gene_ids[$1] = 1;
    next
} /^>/ {
    split(substr($0,2), id_arr, " ");
    current_id = id_arr[1];
    flag = (current_id in gene_ids) ? 1 : 0;
    print;
    next
} flag' tmp_cand.tsv "$GENOME_DIR/$STRAIN_NAME.fna" > tmp_seq.fna

# 清理空序列（关键！解决BLAST报错）
sed -i '/^>/!{/^$/d}' tmp_seq.fna  # 删除纯空行
seq_count=$(grep -c "^>" tmp_seq.fna)
echo "有效序列数：$seq_count"

# 4. BLAST比对（若有有效序列才运行）
echo -e "\n🧬 BLAST比对VFDB..."
if [ $seq_count -gt 0 ]; then
    blastn -query tmp_seq.fna -db "$VFDB" -outfmt "6 qseqid pident evalue" -evalue 1e-10 -num_threads 4 > tmp_blast.tsv
    blast_count=$(wc -l tmp_blast.tsv | awk '{print $1}')
    echo "BLAST有效命中数：$blast_count"
else
    echo "⚠️  无有效序列，跳过BLAST"
    blast_count=0
fi

# 5. 生成结果（有命中才写入）
echo -e "\n📊 生成测试结果..."
echo "菌株名,基因ID,基因名,毒力功能,BLAST相似度,E值" > $TEST_OUTPUT
if [ $blast_count -gt 0 ]; then
    join -t $'\t' tmp_cand.tsv <(sort tmp_blast.tsv) | awk -v s="$STRAIN_NAME" -F "\t" '{print s "," $1 "," $2 "," $3 "," $5 "," $6}' >> $TEST_OUTPUT
    echo -e "\n🎉 测试成功！结果如下："
    cat $TEST_OUTPUT
else
    echo -e "\n⚠️  该菌株无有效毒力基因命中，换另一株测试即可"
fi
