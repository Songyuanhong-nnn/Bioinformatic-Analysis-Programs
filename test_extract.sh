#!/bin/bash
set -euo pipefail

# 选1株已注释的菌株（从bakta_annotations目录选1个，比如GCA_000328405.1_ASM32840v1_genomic）
STRAIN_NAME="GCA_000430405.1_ASM430405v1_genomic"  # 替换成你bakta_annotations目录下的任意菌株名
VFDB="/mnt/g/wsl/database/vfdb/VFDB_setA_nt.fas"
BAKTA_DIR="/mnt/d/WSL/disk/projects/VP1/bakta_annotations"
GENOME_DIR="/mnt/d/WSL/disk/projects/VP1/clean_genome"
TEST_OUTPUT="test_virulence.csv"

# 清理之前的临时文件
rm -f tmp_*.tsv tmp_*.fna $TEST_OUTPUT

# 1. 确认该菌株的文件存在
echo "确认菌株文件..."
ls -l "$BAKTA_DIR/$STRAIN_NAME/$STRAIN_NAME.tsv"
ls -l "$GENOME_DIR/$STRAIN_NAME.fna"

# 2. 筛选毒力基因
echo -e "\n筛选毒力基因..."
grep -E "virulence|toxin|tdh|trh" "$BAKTA_DIR/$STRAIN_NAME/$STRAIN_NAME.tsv" | awk -F "\t" '{print $1 "\t" $2 "\t" $5}' > tmp_cand.tsv
echo "筛选到候选基因数：$(wc -l tmp_cand.tsv | awk '{print $1}')"

# 3. 提取序列
echo "提取基因序列..."
awk 'NR==FNR{a[$1]=1;next}/^>/{split(substr($0,2),i," ");f=i[1] in a?1:0;print;next}f' tmp_cand.tsv "$GENOME_DIR/$STRAIN_NAME.fna" > tmp_seq.fna
echo "提取到序列数：$(grep -c "^>" tmp_seq.fna)"

# 4. BLAST验证
echo "BLAST比对..."
blastn -query tmp_seq.fna -db "$VFDB" -outfmt "6 qseqid pident evalue" -evalue 1e-10 -num_threads 4 > tmp_blast.tsv
echo "BLAST命中数：$(wc -l tmp_blast.tsv | awk '{print $1}')"

# 5. 写入结果
echo "菌株名,基因ID,基因名,毒力功能,BLAST相似度,E值" > $TEST_OUTPUT
join -t $'\t' tmp_cand.tsv <(sort tmp_blast.tsv) | awk -v s="$STRAIN_NAME" -F "\t" '{print s "," $1 "," $2 "," $3 "," $5 "," $6}' >> $TEST_OUTPUT

echo -e "\n✅ 测试完成！结果文件：$TEST_OUTPUT"
cat $TEST_OUTPUT  # 直接显示结果
