#!/bin/bash
set -euo pipefail

# 👉 替换成你刚才用ls命令确认存在的菌株名（必须完全一致！）
STRAIN_NAME="GCA_000328405.1_ASM32840v1_genomic"
VFDB="/mnt/g/wsl/database/vfdb/VFDB_setA_nt.fas"
BAKTA_DIR="/mnt/d/WSL/disk/projects/VP1/bakta_annotations"
GENOME_DIR="/mnt/d/WSL/disk/projects/VP1/clean_genome"
TEST_OUTPUT="test_virulence.csv"

# 清理旧文件
rm -f tmp_*.tsv tmp_*.fna $TEST_OUTPUT

# 1. 确认菌株文件存在（这步过了才往下走）
echo "✅ 确认菌株文件..."
ls -l "$BAKTA_DIR/$STRAIN_NAME/$STRAIN_NAME.tsv"
ls -l "$GENOME_DIR/$STRAIN_NAME.fna"

# 2. 筛选毒力候选基因
echo -e "\n🔍 筛选毒力基因..."
grep -E "virulence|toxin|tdh|trh|hemolysin" "$BAKTA_DIR/$STRAIN_NAME/$STRAIN_NAME.tsv" | awk -F "\t" '{print $1 "\t" $2 "\t" $5}' > tmp_cand.tsv
cand_count=$(wc -l tmp_cand.tsv | awk '{print $1}')
echo "筛选到候选基因数：$cand_count"

# 3. 提取基因序列
echo -e "\n📑 提取候选基因序列..."
awk 'NR==FNR{a[$1]=1;next}/^>/{split(substr($0,2),i," ");f=i[1] in a?1:0;print;next}f' tmp_cand.tsv "$GENOME_DIR/$STRAIN_NAME.fna" > tmp_seq.fna
seq_count=$(grep -c "^>" tmp_seq.fna)
echo "提取到序列数：$seq_count"

# 4. BLAST比对验证
echo -e "\n🧬 BLAST比对VFDB..."
blastn -query tmp_seq.fna -db "$VFDB" -outfmt "6 qseqid pident evalue" -evalue 1e-10 -num_threads 4 > tmp_blast.tsv
blast_count=$(wc -l tmp_blast.tsv | awk '{print $1}')
echo "BLAST有效命中数：$blast_count"

# 5. 生成结果文件
echo -e "\n📊 生成测试结果..."
echo "菌株名,基因ID,基因名,毒力功能,BLAST相似度,E值" > $TEST_OUTPUT
join -t $'\t' tmp_cand.tsv <(sort tmp_blast.tsv) | awk -v s="$STRAIN_NAME" -F "\t" '{print s "," $1 "," $2 "," $3 "," $5 "," $6}' >> $TEST_OUTPUT

echo -e "\n🎉 测试成功！结果如下："
cat $TEST_OUTPUT
