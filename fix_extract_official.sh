#!/bin/bash
set -euo pipefail

PROJECT_DIR="/mnt/d/WSL/disk/projects/VP1"
BAKTA_DIR="$PROJECT_DIR/bakta_annotations"
GENOME_DIR="$PROJECT_DIR/clean_genome"
VFDB_PROT="$PROJECT_DIR/vfdb_online/VFDB_setA_pro.fas"
OUTPUT="$PROJECT_DIR/FIXED_OFFICIAL_VIRULENCE.csv"

echo "===== 直接从Bakta注释提取+官方VFDB标注 ====="

# 1. 确保VFDB索引存在（避免比对报错）
[ ! -f "${VFDB_PROT}.phr" ] && makeblastdb -in "$VFDB_PROT" -dbtype prot -out "$VFDB_PROT" 2>/dev/null

# 2. 提取VFDB官方基因信息
grep "^>" "$VFDB_PROT" | sed 's/^>//; s/\t/|/g' > "$PROJECT_DIR/vfdb_official.tmp"

# 3. 表头（含官方信息）
echo "菌株名,基因ID,基因名,功能描述,基因坐标,序列ID,VFDB官方基因名,VFDB官方描述,VFDB相似度(%),VFDB_E值,官方来源" > "$OUTPUT"

# 4. 遍历所有Bakta注释菌株，直接提取毒力基因
for STRAIN_DIR in "$BAKTA_DIR"/*; do
    STRAIN=$(basename "$STRAIN_DIR")
    BAKTA_TSV="$STRAIN_DIR/$STRAIN.tsv"
    GENOME_FNA="$GENOME_DIR/$STRAIN.fna"
    
    # 跳过无注释/无基因组的菌株
    [ ! -f "$BAKTA_TSV" ] && continue
    [ ! -f "$GENOME_FNA" ] && continue
    
    echo "🔧 处理菌株：$STRAIN"
    
    # 从Bakta注释筛选毒力相关基因（宽松关键词，确保有结果）
    grep -v "^#" "$BAKTA_TSV" | awk -F "\t" '
    NF >=7 && ($5 ~ /virulence|toxin|tdh|trh|hemolysin|hly|adhesin|fimbria|T3SS|T6SS|biofilm|invasion|colonization|secreted|pathogen/) {
        print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7
    }' > "$PROJECT_DIR/tmp_$STRAIN.tmp"
    
    # 提取序列并比对VFDB
    if [ -s "$PROJECT_DIR/tmp_$STRAIN.tmp" ]; then
        while read -r GENE_ID GENE_NAME START END FUNC STRAND SEQ_ID; do
            COORD="$SEQ_ID:$START-$END"
            
            # 提取核苷酸序列
            SEQUENCE=$(bedtools getfasta -fi "$GENOME_FNA" -bed <(echo -e "$SEQ_ID\t$((START-1))\t$END\t$GENE_ID") -fo - | tail -n +2)
            [ -z "$SEQUENCE" ] && continue
            
            # 用tblastn比对VFDB官方库
            BLAST=$(tblastn -query <(echo -e ">$GENE_ID\n$SEQUENCE") \
                          -db "$VFDB_PROT" \
                          -outfmt "6 sseqid pident evalue" \
                          -evalue 1e-5 \
                          -perc_identity 60 \
                          -num_threads 4 -quiet | head -1)
            
            # 解析官方匹配信息
            if [ -n "$BLAST" ]; then
                VFDB_GENE=$(echo "$BLAST" | awk '{print $1}')
                PID=$(echo "$BLAST" | awk '{print $2}')
                EVAL=$(echo "$BLAST" | awk '{print $3}')
                # 从VFDB官方库获取描述
                VFDB_DESC=$(grep "^$VFDB_GENE|" "$PROJECT_DIR/vfdb_official.tmp" | sed 's/|/\t/g' | awk '{print $2}' || echo "无")
            else
                VFDB_GENE="无匹配"
                PID="0"
                EVAL="1e0"
                VFDB_DESC="无"
            fi
            
            # 写入结果（处理特殊字符）
            FUNC_CLEAN=$(echo "$FUNC" | sed 's/,/;/g')
            echo "$STRAIN,$GENE_ID,$GENE_NAME,$FUNC_CLEAN,$COORD,$SEQ_ID,$VFDB_GENE,$VFDB_DESC,$PID,$EVAL,VFDB官方库（华大基因）" >> "$OUTPUT"
        done < "$PROJECT_DIR/tmp_$STRAIN.tmp"
        
        # 统计该菌株提取到的基因数
        GENE_CNT=$(wc -l "$PROJECT_DIR/tmp_$STRAIN.tmp" | awk '{print $1}')
        echo "   ✅ 提取到 $GENE_CNT 个毒力相关基因"
    else
        echo "   ⚠️  未找到毒力相关基因"
    fi
    
    rm -f "$PROJECT_DIR/tmp_$STRAIN.tmp"
done

# 5. 清理临时文件
rm -f "$PROJECT_DIR/vfdb_official.tmp"

echo -e "\n✅ 修复完成！"
echo "📁 官方标注毒力基因表：$OUTPUT"

# 6. 统计成果质量
TOTAL=$(grep -v "^菌株名" "$OUTPUT" | wc -l)
OFFICIAL_MATCHED=$(grep -v "无匹配" "$OUTPUT" | wc -l)
HIGH_CONFIDENCE=$(grep -E ",[8-9][0-9]\.|,100\." "$OUTPUT" | wc -l)

echo -e "\n===== 成果质量统计 ====="
echo "1. 总毒力基因数：$TOTAL"
echo "2. 官方VFDB匹配数：$OFFICIAL_MATCHED"
echo "3. 高可信度基因数（相似度≥80%）：$HIGH_CONFIDENCE"
if [ "$TOTAL" -gt 0 ]; then
    MATCH_RATE=$(echo "scale=2; $OFFICIAL_MATCHED/$TOTAL*100" | bc)
    echo "4. 官方匹配率：$MATCH_RATE%"
else
    echo "4. 官方匹配率：0%"
fi
