#!/bin/bash
set -euo pipefail

PROJECT_DIR="/mnt/d/WSL/disk/projects/VP1"
BAKTA_DIR="$PROJECT_DIR/bakta_annotations"
GENOME_DIR="$PROJECT_DIR/clean_genome"
VFDB_PROT="$PROJECT_DIR/vfdb_online/VFDB_setA_pro.fas"
OUTPUT="$PROJECT_DIR/FORCE_VIRULENCE_OFFICIAL.csv"

echo "===== 强制提取：VFDB直接比对所有基因 ====="

# 1. 确保VFDB索引存在
[ ! -f "${VFDB_PROT}.phr" ] && makeblastdb -in "$VFDB_PROT" -dbtype prot -out "$VFDB_PROT" 2>/dev/null

# 2. 提取VFDB官方基因信息
grep "^>" "$VFDB_PROT" | sed 's/^>//; s/\t/|/g' > "$PROJECT_DIR/vfdb_official.tmp"

# 3. 表头
echo "菌株名,基因ID,基因名,功能描述,基因坐标,序列ID,VFDB官方基因名,VFDB官方描述,VFDB相似度(%),VFDB_E值,官方来源" > "$OUTPUT"

# 4. 遍历所有菌株，提取所有基因并比对VFDB
for STRAIN_DIR in "$BAKTA_DIR"/*; do
    STRAIN=$(basename "$STRAIN_DIR")
    BAKTA_TSV="$STRAIN_DIR/$STRAIN.tsv"
    GENOME_FNA="$GENOME_DIR/$STRAIN.fna"
    
    [ ! -f "$BAKTA_TSV" ] && continue
    [ ! -f "$GENOME_FNA" ] && continue
    
    echo "🔧 处理菌株：$STRAIN"
    
    # 提取该菌株所有基因（不筛选关键词，全部提取）
    grep -v "^#" "$BAKTA_TSV" | awk -F "\t" 'NF >=7 {print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7}' > "$PROJECT_DIR/tmp_$STRAIN.all.tmp"
    
    # 临时文件存储所有基因序列
    > "$PROJECT_DIR/tmp_$STRAIN.seq.fna"
    while read -r GENE_ID GENE_NAME START END FUNC STRAND SEQ_ID; do
        # 提取核苷酸序列
        SEQUENCE=$(bedtools getfasta -fi "$GENOME_FNA" -bed <(echo -e "$SEQ_ID\t$((START-1))\t$END\t$GENE_ID") -fo - | tail -n +2)
        [ -z "$SEQUENCE" ] && continue
        # 写入临时序列文件
        echo -e ">$GENE_ID|$STRAIN|$GENE_NAME" >> "$PROJECT_DIR/tmp_$STRAIN.seq.fna"
        echo "$SEQUENCE" >> "$PROJECT_DIR/tmp_$STRAIN.seq.fna"
    done < "$PROJECT_DIR/tmp_$STRAIN.all.tmp"
    
    # 用tblastn比对VFDB官方库（所有基因批量比对）
    tblastn -query "$PROJECT_DIR/tmp_$STRAIN.seq.fna" \
            -db "$VFDB_PROT" \
            -outfmt "6 qseqid sseqid pident evalue" \
            -evalue 1e-3 \
            -perc_identity 50 \
            -num_threads 8 -quiet > "$PROJECT_DIR/tmp_$STRAIN.blast.tmp"
    
    # 解析比对结果，筛选毒力基因
    GENE_CNT=0
    while read -r BLAST_LINE; do
        QSEQID=$(echo "$BLAST_LINE" | awk '{print $1}')
        VFDB_GENE=$(echo "$BLAST_LINE" | awk '{print $2}')
        PID=$(echo "$BLAST_LINE" | awk '{print $3}')
        EVAL=$(echo "$BLAST_LINE" | awk '{print $4}')
        
        # 解析基因信息（QSEQID格式：GENE_ID|STRAIN|GENE_NAME）
        GENE_ID=$(echo "$QSEQID" | cut -d '|' -f1)
        GENE_NAME=$(echo "$QSEQID" | cut -d '|' -f3)
        
        # 从Bakta注释中获取功能描述和坐标
        GENE_INFO=$(grep -w "$GENE_ID" "$PROJECT_DIR/tmp_$STRAIN.all.tmp" | head -1)
        [ -z "$GENE_INFO" ] && continue
        
        START=$(echo "$GENE_INFO" | awk '{print $3}')
        END=$(echo "$GENE_INFO" | awk '{print $4}')
        FUNC=$(echo "$GENE_INFO" | awk '{print $5}' | sed 's/,/;/g')
        STRAND=$(echo "$GENE_INFO" | awk '{print $6}')
        SEQ_ID=$(echo "$GENE_INFO" | awk '{print $7}')
        COORD="$SEQ_ID:$START-$END"
        
        # 从VFDB官方库获取描述
        VFDB_DESC=$(grep "^$VFDB_GENE|" "$PROJECT_DIR/vfdb_official.tmp" | sed 's/|/\t/g' | awk '{print $2}' || echo "无")
        
        # 写入结果
        echo "$STRAIN,$GENE_ID,$GENE_NAME,$FUNC,$COORD,$SEQ_ID,$VFDB_GENE,$VFDB_DESC,$PID,$EVAL,VFDB官方库（华大基因）" >> "$OUTPUT"
        GENE_CNT=$((GENE_CNT+1))
    done < "$PROJECT_DIR/tmp_$STRAIN.blast.tmp"
    
    echo "   ✅ 找到 $GENE_CNT 个VFDB匹配的毒力相关基因"
    
    # 清理临时文件
    rm -f "$PROJECT_DIR/tmp_$STRAIN.all.tmp" "$PROJECT_DIR/tmp_$STRAIN.seq.fna" "$PROJECT_DIR/tmp_$STRAIN.blast.tmp"
done

# 5. 清理临时文件
rm -f "$PROJECT_DIR/vfdb_official.tmp"

echo -e "\n✅ 强制提取完成！"
echo "📁 官方毒力基因表：$OUTPUT"

# 6. 统计成果
TOTAL=$(grep -v "^菌株名" "$OUTPUT" | wc -l)
HIGH_CONFIDENCE=$(grep -E ",[7-9][0-9]\.|,100\." "$OUTPUT" | wc -l)

echo -e "\n===== 成果统计 ====="
echo "1. 总毒力基因数（VFDB匹配）：$TOTAL"
echo "2. 高可信度基因数（相似度≥70%）：$HIGH_CONFIDENCE"
echo "3. 覆盖菌株数：$(grep -v "^菌株名" "$OUTPUT" | awk -F "," '{print $1}' | sort -u | wc -l) 株"
