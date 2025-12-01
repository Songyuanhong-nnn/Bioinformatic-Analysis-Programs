#!/bin/bash
set -euo pipefail

PROJECT_DIR="/mnt/d/WSL/disk/projects/VP1"
BAKTA_DIR="$PROJECT_DIR/bakta_annotations"  # Bakta输出的蛋白序列（*.faa）
VFDB_PROT="$PROJECT_DIR/vfdb_online/VFDB_setA_pro.fas"
OUTPUT="$PROJECT_DIR/ULTIMATE_VIRULENCE_OFFICIAL.csv"

echo "===== 极简模式：VFDB直接比对菌株蛋白组 ====="

# 1. 确保VFDB索引（只做1次）
[ ! -f "${VFDB_PROT}.phr" ] && makeblastdb -in "$VFDB_PROT" -dbtype prot -out "$VFDB_PROT" 2>/dev/null

# 2. 提取VFDB官方基因信息（名称+描述）
grep "^>" "$VFDB_PROT" | sed 's/^>//; s/\t/|/g' > "$PROJECT_DIR/vfdb_tmp.txt"

# 3. 表头（含官方标注）
echo "菌株名,菌株蛋白ID,毒力基因名,VFDB官方描述,相似度(%),E值,官方来源" > "$OUTPUT"

# 4. 批量比对所有菌株（只比对蛋白组，跳过所有中间步骤）
for STRAIN_DIR in "$BAKTA_DIR"/*; do
    STRAIN=$(basename "$STRAIN_DIR")
    STRAIN_FAA="$STRAIN_DIR/$STRAIN.faa"  # Bakta自带的蛋白质序列文件（100%存在）
    
    [ ! -f "$STRAIN_FAA" ] && { echo "⚠️  $STRAIN：无蛋白文件，跳过"; continue; }
    
    echo -e "\n🔧 处理菌株：$STRAIN"
    
    # 核心：blastp直接比对（蛋白→蛋白，最快最准）
    blastp -query "$STRAIN_FAA" \
           -db "$VFDB_PROT" \
           -outfmt "6 qseqid sseqid pident evalue" \
           -evalue 1e-5 \
           -perc_identity 60 \
           -num_threads 8 \
           -max_target_seqs 1 \
           -quiet > "$PROJECT_DIR/blast_$STRAIN.tmp"
    
    # 解析结果，只保留有效匹配
    GENE_CNT=0
    while read -r BLAST_LINE; do
        STRAIN_PROT_ID=$(echo "$BLAST_LINE" | awk '{print $1}')  # 菌株蛋白ID
        VFDB_GENE=$(echo "$BLAST_LINE" | awk '{print $2}')       # VFDB官方毒力基因名
        PID=$(echo "$BLAST_LINE" | awk '{print $3}')             # 相似度
        EVAL=$(echo "$BLAST_LINE" | awk '{print $4}')            # E值
        
        # 获取VFDB官方描述
        VFDB_DESC=$(grep "^$VFDB_GENE|" "$PROJECT_DIR/vfdb_tmp.txt" | cut -d '|' -f2 || "无")
        
        # 写入结果
        echo "$STRAIN,$STRAIN_PROT_ID,$VFDB_GENE,$VFDB_DESC,$PID,$EVAL,VFDB官方库（华大基因）" >> "$OUTPUT"
        GENE_CNT=$((GENE_CNT+1))
    done < "$PROJECT_DIR/blast_$STRAIN.tmp"
    
    echo "   ✅ 找到 $GENE_CNT 个毒力相关基因"
    rm -f "$PROJECT_DIR/blast_$STRAIN.tmp"
done

# 5. 清理临时文件
rm -f "$PROJECT_DIR/vfdb_tmp.txt"

echo -e "\n🎉 所有菌株处理完成！"
echo "📁 最终结果：$OUTPUT"

# 6. 成果统计
TOTAL=$(grep -v "^菌株名" "$OUTPUT" | wc -l)
HIGH_CONFID=$(grep -E ",[7-9][0-9]\.|,100\." "$OUTPUT" | wc -l)
STRAINS=$(grep -v "^菌株名" "$OUTPUT" | awk -F "," '{print $1}' | sort -u | wc -l)

echo -e "\n===== 成果统计 ====="
echo "1. 总毒力基因数（VFDB官方匹配）：$TOTAL"
echo "2. 高可信度基因数（相似度≥70%）：$HIGH_CONFID"
echo "3. 覆盖菌株数：$STRAINS 株"
echo "4. 平均每株毒力基因数：$( [ $STRAINS -gt 0 ] && echo "scale=1; $TOTAL/$STRAINS" | bc || echo 0 ) 个"
echo -e "\n📚 论文引用格式：Chen, L., et al. (2016). VFDB 2016: hierarchical and refined dataset for bacterial virulence factors. Nucleic Acids Res, 44(D1), D694-D698."
