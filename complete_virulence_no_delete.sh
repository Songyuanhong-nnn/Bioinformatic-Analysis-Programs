#!/bin/bash
set -euo pipefail

# -------------------------- 核心配置（不修改）--------------------------
PROJECT_DIR="/mnt/d/WSL/disk/projects/VP1"
BAKTA_DIR="$PROJECT_DIR/bakta_annotations"
VFDB_A="$PROJECT_DIR/vfdb_online/VFDB_setA_pro.fas"
VFDB_B="$PROJECT_DIR/vfdb_online/VFDB_setB_pro.fas"
VFDB_ANNO="$PROJECT_DIR/vfdb_online/VFs.xls"
OUTPUT="$PROJECT_DIR/COMPLETE_ALL_FILES_VIRULENCE.csv"
TEMP_MERGE_DIR="$PROJECT_DIR/temp_merged_proteins"
# ----------------------------------------------------------------------

# 1. 创建临时合并目录
mkdir -p "$TEMP_MERGE_DIR"
echo "📁 创建临时合并目录：$TEMP_MERGE_DIR"

# 2. 初始化结果文件（含来源文件追溯字段）
echo "菌株名,蛋白ID,毒力因子名,VFDB来源,相似度(%),E值,功能描述,VFDB分类,参考文献,来源蛋白文件" > "$OUTPUT"

# 3. 收集所有 .faa 文件（包括 hypotheticals.faa）
ALL_FAA_FILES=$(find "$BAKTA_DIR" -name "*.faa" | sort)
echo "📊 共检测到 $(echo "$ALL_FAA_FILES" | wc -l) 个蛋白文件"

# 4. 按菌株分组合并（修正变量拼接语法：$STRAIN_merged → ${STRAIN}_merged）
echo -e "\n🔄 按菌株分组并合并蛋白文件..."
while IFS= read -r FAA_FILE; do
    STRAIN=$(basename "$(dirname "$FAA_FILE")")
    # 修正：变量拼接必须用 ${STRAIN} 明确边界，避免识别为 STAIN_merged 变量
    MERGED_FAA="$TEMP_MERGE_DIR/${STRAIN}_merged.faa"
    
    # 合并当前菌株的所有 .faa 文件
    cat "$FAA_FILE" >> "$MERGED_FAA"
    # 记录蛋白文件来源
    echo "$STRAIN,$FAA_FILE" >> "$TEMP_MERGE_DIR/strains_file_map.csv"
done <<< "$ALL_FAA_FILES"

# 5. 获取去重后的菌株列表
MERGED_STRAINS=$(ls "$TEMP_MERGE_DIR"/*_merged.faa | sort)
TOTAL_STRAINS=$(echo "$MERGED_STRAINS" | wc -l)
echo -e "\n✅ 合并完成！共得到 $TOTAL_STRAINS 个菌株的完整蛋白文件"
echo "------------------------------------------------------------------------"

# 6. 逐个菌株比对
while IFS= read -r MERGED_FAA; do
    STRAIN=$(basename "$MERGED_FAA" "_merged.faa")
    echo -e "\n🔧 开始处理菌株：$STRAIN"

    # 7. SetA库比对（实验验证因子）
    blastp -query "$MERGED_FAA" \
           -db "$VFDB_A" \
           -outfmt "6 qseqid sseqid pident evalue" \
           -evalue 1e-5 \
           -num_threads 6 \
           -max_target_seqs 1 2>/dev/null > "$TEMP_MERGE_DIR/${STRAIN}_setA.blast"

    # 8. SetB库比对（预测因子）
    blastp -query "$MERGED_FAA" \
           -db "$VFDB_B" \
           -outfmt "6 qseqid sseqid pident evalue" \
           -evalue 1e-5 \
           -num_threads 6 \
           -max_target_seqs 1 2>/dev/null > "$TEMP_MERGE_DIR/${STRAIN}_setB.blast"

    # 9. 合并去重
    cat "$TEMP_MERGE_DIR/${STRAIN}_setA.blast" > "$TEMP_MERGE_DIR/${STRAIN}_all.blast"
    grep -v -f <(cut -f1 "$TEMP_MERGE_DIR/${STRAIN}_setA.blast") "$TEMP_MERGE_DIR/${STRAIN}_setB.blast" >> "$TEMP_MERGE_DIR/${STRAIN}_all.blast"

    # 10. 解析结果写入
    VF_COUNT=0
    while IFS= read -r BLAST_LINE; do
        PROT_ID=$(echo "$BLAST_LINE" | cut -f1)
        VF_NAME=$(echo "$BLAST_LINE" | cut -f2)
        PID=$(printf "%.2f" $(echo "$BLAST_LINE" | cut -f3))
        EVAL=$(echo "$BLAST_LINE" | cut -f4)

        # 判断来源
        [ $(grep -w "^$VF_NAME" "$VFDB_A" 2>/dev/null | wc -l) -gt 0 ] && SOURCE="SetA（实验验证）" || SOURCE="SetB（预测补充）"

        # 关联注释
        ANNO=$(grep -w "^$VF_NAME" "$VFDB_ANNO" 2>/dev/null | head -1)
        FUNC=${ANNO:+"$(echo "$ANNO" | cut -f4 | sed 's/,/;/g')"} || FUNC="无官方注释"
        CAT=${ANNO:+"$(echo "$ANNO" | cut -f5 | sed 's/,/;/g')"} || CAT="无"
        REF=${ANNO:+"$(echo "$ANNO" | cut -f7 | sed 's/,/;/g')"} || REF="无"

        # 追溯来源文件
        SOURCE_FILE=$(grep "^$STRAIN," "$TEMP_MERGE_DIR/strains_file_map.csv" | cut -d ',' -f2 | head -1)

        # 写入结果
        echo "$STRAIN,$PROT_ID,$VF_NAME,$SOURCE,$PID,$EVAL,$FUNC,$CAT,$REF,$SOURCE_FILE" >> "$OUTPUT"
        ((VF_COUNT++))
    done < "$TEMP_MERGE_DIR/${STRAIN}_all.blast"

    echo "   ✅ 完成！找到 $VF_COUNT 个毒力因子"
    echo "------------------------------------------------------------------------"

done <<< "$MERGED_STRAINS"

# 11. 清理临时文件
rm -rf "$TEMP_MERGE_DIR"
echo -e "\n🗑️  清理临时文件完成"

# 12. 最终统计
PROCESSED_STRAINS=$(grep -v "^菌株名" "$OUTPUT" | cut -d ',' -f1 | sort -u | wc -l)
TOTAL_VF=$(grep -v "^菌株名" "$OUTPUT" | wc -l)
UNIQUE_VF=$(grep -v "^菌株名" "$OUTPUT" | cut -d ',' -f3 | sort -u | wc -l)

echo -e "\n🎉 所有菌株完整处理完毕！"
echo "📈 最终统计："
echo "   - 总菌株数：$TOTAL_STRAINS"
echo "   - 成功处理菌株数：$PROCESSED_STRAINS"
echo "   - 共鉴定毒力因子数：$TOTAL_VF"
echo "   - 独特毒力因子数：$UNIQUE_VF"
echo "   - 完整结果文件：$OUTPUT"
