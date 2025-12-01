#!/bin/bash
# 关闭严格模式（避免变量空值报错，只针对兜底版）
set -eo pipefail

# -------------------------- 核心配置 --------------------------
PROJECT_DIR="/mnt/d/WSL/disk/projects/VP1"
BAKTA_DIR="$PROJECT_DIR/bakta_annotations"
VFDB_A="$PROJECT_DIR/vfdb_online/VFDB_setA_pro.fas"
VFDB_B="$PROJECT_DIR/vfdb_online/VFDB_setB_pro.fas"
VFDB_ANNO="$PROJECT_DIR/vfdb_online/VFs.xls"
OUTPUT="$PROJECT_DIR/ULTIMATE_VIRULENCE_FINAL.csv"
TEMP_DIR="$PROJECT_DIR/temp_ultimate"
# ----------------------------------------------------------------

# 1. 强制清理所有残留（兜底保障）
rm -rf "$TEMP_DIR" && mkdir -p "$TEMP_DIR"
echo -e "\n📋 初始化完成！终极兜底版开始处理（必自动推进）..."
echo "========================================================================"

# 2. 初始化结果文件
echo "菌株名,蛋白ID,毒力因子名,VFDB来源,相似度(%),E值,功能描述,VFDB分类,参考文献,来源蛋白文件" > "$OUTPUT"

# 3. 用数组存储所有菌株目录（最稳妥的方式）
STRAIN_DIRS=()
while IFS= read -r dir; do
    STRAIN_DIRS+=("$dir")
done < <(find "$BAKTA_DIR" -maxdepth 1 -type d ! -name "bakta_annotations")

TOTAL_STRAINS=${#STRAIN_DIRS[@]}
TOTAL_FAA=$(find "$BAKTA_DIR" -name "*.faa" | wc -l)
echo "📊 全局统计：共 $TOTAL_STRAINS 株菌 | $TOTAL_FAA 个蛋白文件"
echo "========================================================================"

# 4. 用数组索引逐个处理（彻底避开循环兼容问题）
for ((i=0; i<TOTAL_STRAINS; i++)); do
    STRAIN_DIR=${STRAIN_DIRS[$i]}
    STRAIN=$(basename "$STRAIN_DIR")
    PROCESSED_STRAINS=$((i+1))
    STRAIN_VF_COUNT=0
    GLOBAL_VF_COUNT=0

    # 前端实时显示（这次必出来）
    echo -e "\n🔴 当前进度：$PROCESSED_STRAINS/$TOTAL_STRAINS 株菌"
    echo "🟢 当前处理：$STRAIN"
    echo "------------------------------------------------------------------------"

    # 合并蛋白文件（兜底：即使无文件也不报错）
    FAA_FILES=$(find "$STRAIN_DIR" -name "*.faa")
    STRAIN_FAA_COUNT=$(echo "$FAA_FILES" | wc -l)
    echo "📥 正在合并 $STRAIN_FAA_COUNT 个蛋白文件..."
    MERGED_FAA="$TEMP_DIR/${STRAIN}_merged.faa"
    cat $FAA_FILES > "$MERGED_FAA" 2>/dev/null
    echo "✅ 合并完成！蛋白文件大小：$(du -sh "$MERGED_FAA" 2>/dev/null | cut -f1)"

    # 比对VFDB库（显示详细状态）
    echo -e "\n⚡ 开始比对VFDB SetA（实验验证）库..."
    blastp -query "$MERGED_FAA" -db "$VFDB_A" -outfmt "6 qseqid sseqid pident evalue" -evalue 1e-5 -num_threads 4 -max_target_seqs 1 2>/dev/null > "$TEMP_DIR/${STRAIN}_setA.blast"
    echo "⚡ 开始比对VFDB SetB（预测补充）库..."
    blastp -query "$MERGED_FAA" -db "$VFDB_B" -outfmt "6 qseqid sseqid pident evalue" -evalue 1e-5 -num_threads 4 -max_target_seqs 1 2>/dev/null > "$TEMP_DIR/${STRAIN}_setB.blast"

    # 合并去重
    cat "$TEMP_DIR/${STRAIN}_setA.blast" 2>/dev/null > "$TEMP_DIR/${STRAIN}_all.blast"
    grep -v -f <(cut -f1 "$TEMP_DIR/${STRAIN}_setA.blast" 2>/dev/null) "$TEMP_DIR/${STRAIN}_setB.blast" 2>/dev/null >> "$TEMP_DIR/${STRAIN}_all.blast"
    TOTAL_HITS=$(wc -l < "$TEMP_DIR/${STRAIN}_all.blast" 2>/dev/null)
    echo "📊 比对完成！共找到 $TOTAL_HITS 个潜在毒力因子"

    # 解析结果（实时刷新）
    echo -e "\n📝 正在解析结果（实时计数）："
    while IFS= read -r BLAST_LINE; do
        PROT_ID=$(echo "$BLAST_LINE" | cut -f1)
        VF_NAME=$(echo "$BLAST_LINE" | cut -f2)
        PID=$(printf "%.2f" $(echo "$BLAST_LINE" | cut -f3))
        EVAL=$(echo "$BLAST_LINE" | cut -f4)
        
        if grep -w "^$VF_NAME" "$VFDB_A" 2>/dev/null; then
            SOURCE="SetA（实验验证）"
        else
            SOURCE="SetB（预测补充）"
        fi

        ANNO=$(grep -w "^$VF_NAME" "$VFDB_ANNO" 2>/dev/null | head -1)
        FUNC=${ANNO:+"$(echo "$ANNO" | cut -f4 | sed 's/,/;/g')"} || FUNC="无"
        CAT=${ANNO:+"$(echo "$ANNO" | cut -f5 | sed 's/,/;/g')"} || CAT="无"
        REF=${ANNO:+"$(echo "$ANNO" | cut -f7 | sed 's/,/;/g')"} || REF="无"
        SOURCE_FILE=$(find "$STRAIN_DIR" -name "*.faa" | head -1)

        echo "$STRAIN,$PROT_ID,$VF_NAME,$SOURCE,$PID,$EVAL,$FUNC,$CAT,$REF,$SOURCE_FILE" >> "$OUTPUT"
        ((STRAIN_VF_COUNT++))
        ((GLOBAL_VF_COUNT++))
        echo -ne "   当前菌株：$STRAIN_VF_COUNT 个 | 全局累计：$GLOBAL_VF_COUNT 个\r"
        sleep 0.1
    done < "$TEMP_DIR/${STRAIN}_all.blast" 2>/dev/null

    # 菌株处理总结
    echo -e "\n\n✅ 【$STRAIN】处理完成！最终识别 $STRAIN_VF_COUNT 个毒力因子"
    echo "------------------------------------------------------------------------"
    echo "⏳ 剩余菌株：$((TOTAL_STRAINS - PROCESSED_STRAINS)) 株 | 预计剩余：$(( (TOTAL_STRAINS - PROCESSED_STRAINS) * 2 )) 分钟"
done

# 全局总结
rm -rf "$TEMP_DIR"
echo -e "\n========================================================================"
echo "🎉 所有菌株处理完毕！"
echo "📈 最终统计：总菌株 $TOTAL_STRAINS 株 | 总毒力因子 $GLOBAL_VF_COUNT 个"
echo "📄 结果文件：$OUTPUT"
echo "========================================================================"
