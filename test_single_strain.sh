#!/bin/bash
set -euo pipefail

PROJECT_DIR="/mnt/d/WSL/disk/projects/VP1"
STRAIN="GCA_000328405.1_ASM32840v1_genomic"
BAKTA_TSV="$PROJECT_DIR/bakta_annotations/$STRAIN/$STRAIN.tsv"
GENOME_FNA="$PROJECT_DIR/clean_genome/$STRAIN.fna"
VFDB_PROT="$PROJECT_DIR/vfdb_online/VFDB_setA_pro.fas"
OUTPUT="$PROJECT_DIR/TEST_SINGLE_VIRULENCE.csv"

echo "===== 测试单个菌株：$STRAIN ====="

# 1. 提取FASTA ID
DEFAULT_FASTA_ID=$(grep "^>" "$GENOME_FNA" | sed 's/^>//' | awk '{print $1}' | head -1)
echo "FASTA ID：$DEFAULT_FASTA_ID"

# 2. 快速提取前20个有效基因（避免基因过多卡住）
grep -v "^#" "$BAKTA_TSV" | awk -F "\t" -v def_id="$DEFAULT_FASTA_ID" '
NF >=7 && $3 ~ /^[0-9]+$/ && $4 ~ /^[0-9]+$/ && $3 < $4 {
    print def_id "\t" $((3-1)) "\t" $4 "\t" $1 "\t" $2 "\t" $5
}' | head -20 > "$PROJECT_DIR/test.tmp"

echo "提取到20个有效基因，开始提取序列..."

# 3. 提取序列
bedtools getfasta -fi "$GENOME_FNA" -bed <(awk '{print $1 "\t" $2 "\t" $3 "\t" $4}' "$PROJECT_DIR/test.tmp") -fo "$PROJECT_DIR/test.seq.fna"

# 4. 快速比对VFDB
echo "开始比对VFDB..."
tblastn -query "$PROJECT_DIR/test.seq.fna" \
        -db "$VFDB_PROT" \
        -outfmt "6 qseqid sseqid pident evalue" \
        -evalue 1e-2 \
        -num_threads 8 \
        -max_target_seqs 1 \
        -quiet > "$PROJECT_DIR/test.blast"

# 5. 解析结果
echo "菌株名,基因ID,基因名,功能描述,坐标,VFDB基因,VFDB描述,相似度(%),E值" > "$OUTPUT"
while read -r BLAST_LINE; do
    QSEQID=$(echo "$BLAST_LINE" | awk '{print $1}')
    VFDB_GENE=$(echo "$BLAST_LINE" | awk '{print $2}')
    PID=$(echo "$BLAST_LINE" | awk '{print $3}')
    EVAL=$(echo "$BLAST_LINE" | awk '{print $4}')
    
    # 匹配基因信息
    INFO=$(grep -w "$QSEQID" "$PROJECT_DIR/test.tmp")
    SEQ_ID=$(echo "$INFO" | awk '{print $1}')
    START=$(echo "$INFO" | awk '{print $2+1}')
    END=$(echo "$INFO" | awk '{print $3}')
    GENE_NAME=$(echo "$INFO" | awk '{print $5}')
    FUNC=$(echo "$INFO" | awk '{print $6}' | sed 's/,/;/g')
    COORD="$SEQ_ID:$START-$END"
    
    # VFDB描述
    VFDB_DESC=$(grep "^$VFDB_GENE|" <(grep "^>" "$VFDB_PROT" | sed 's/^>//; s/\t/|/g') | sed 's/|/\t/g' | awk '{print $2}' || "无")
    
    echo "$STRAIN,$QSEQID,$GENE_NAME,$FUNC,$COORD,$VFDB_GENE,$VFDB_DESC,$PID,$EVAL" >> "$OUTPUT"
done < "$PROJECT_DIR/test.blast"

# 清理临时文件
rm -f "$PROJECT_DIR/test.tmp" "$PROJECT_DIR/test.seq.fna" "$PROJECT_DIR/test.blast"

echo -e "\n✅ 单个菌株测试完成！"
echo "结果文件：$OUTPUT"
echo "找到的毒力基因数：$(grep -v "^菌株名" "$OUTPUT" | wc -l)"
