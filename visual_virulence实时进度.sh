#!/bin/bash
set -euo pipefail

# -------------------------- 核心配置 --------------------------
PROJECT_DIR="/mnt/d/WSL/disk/projects/VP1"
BAKTA_DIR="$PROJECT_DIR/bakta_annotations"
VFDB_A="$PROJECT_DIR/vfdb_online/VFDB_setA_pro.fas"
VFDB_B="$PROJECT_DIR/vfdb_online/VFDB_setB_pro.fas"
VFDB_ANNO="$PROJECT_DIR/vfdb_online/VFs.xls"
OUTPUT="$PROJECT_DIR/COMPLETE_VISUAL_VIRULENCE.csv"
TEMP_DIR="$PROJECT_DIR/temp_visual"
# ----------------------------------------------------------------

# 1. 初始化环境（清理旧数据，避免冲突）
rm -rf "$TEMP_DIR" && mkdir -p "$TEMP_DIR"
echo -e "\n📋 初始化完成！开始实时可视化毒力因子提取..."
echo "========================================================================"

# 2. 初始化结果文件（含来源追溯）
echo "菌株名,蛋白ID,毒力因子名,VFDB来源,相似度(%),E值,功能描述,VFDB分类,参考文献,来源蛋白文件" > "$OUTPUT"

# 3. 收集所有菌株和蛋白文件（前端显示统计）
ALL_FAA=$(find "$BAKTA_DIR" -name "*.faa" | sort)
TOTAL_FAA=$(echo "$ALL_FAA" | wc -l)
STRAIN_LIST=$(echo "$ALL_FAA" | xargs -I {} dirname {} | xargs -I {} basename {} | sort -u)
TOTAL_STRAINS=$(echo "$STRAIN_LIST" | wc -l)
echo "📊 全局统计：共 $TOTAL_STRAINS 株菌 | $TOTAL_FAA 个蛋白文件"
echo "========================================================================"

# 4. 初始化进度变量
PROCESSED_STRAINS=0
GLOBAL_VF_COUNT=0

# 5. 逐个菌株可视化处理（核心前端逻辑）
while IFS= read -r STRAIN; do
    ((PROCESSED_STRAINS++))
    STRAIN_VF_COUNT=0
    echo -e "\n🔴 当前进度：$PROCESSED_STRAINS/$TOTAL_STRAINS 株菌"
    echo "🟢 当前处理：$STRAIN"
    echo "------------------------------------------------------------------------"

    # 5.1 合并当前菌株的所有蛋白文件（前端显示合并进度）
    FAA_FILES=$(find "$BAKTA_DIR/$STRAIN" -name "*.faa")
    STRAIN_FAA_COUNT=$(echo "$FAA_FILES" | wc -l)
    echo "📥 正在合并 $STRAIN_FAA_COUNT 个蛋白文件..."
    MERGED_FAA="$TEMP_DIR/${STRAIN}_merged.faa"
    cat $FAA_FILES > "$MERGED_FAA"
    echo "✅ 合并完成！蛋白文件大小：$(du -sh "$MERGED_FAA" | cut -f1)"

    # 5.2 前端显示比对开始（明确告知用户正在做什么）
    echo -e "\n⚡ 开始与VFDB库比对（SetA实验验证 + SetB预测）..."
    echo "   线程数：6 | E-value阈值：1e-5 | 最大匹配数：1"

    # 5.3 执行比对（前端实时显示比对状态）
    blastp -query "$MERGED_FAA" -db "$VFDB_A" -outfmt "6 qseqid sseqid pident evalue" -evalue 1e-5 -num_threads 6 -max_target_seqs 1 2>/dev/null > "$TEMP_DIR/${STRAIN}_setA.blast"
    blastp -query "$MERGED_FAA" -db "$VFDB_B" -outfmt "6 qseqid sseqid pident evalue" -evalue 1e-5 -num_threads 6 -max_target_seqs 1 2>/dev/null > "$TEMP_DIR/${STRAIN}_setB.blast"

    # 5.4 合并去重（前端显示结果统计）
    cat "$TEMP_DIR/${STRAIN}_setA.blast" > "$TEMP_DIR/${STRAIN}_all.blast"
    grep -v -f <(cut -f1 "$TEMP_DIR/${STRAIN}_setA.blast") "$TEMP_DIR/${STRAIN}_setB.blast" >> "$TEMP_DIR/${STRAIN}_all.blast"
    TOTAL_HITS=$(wc -l < "$TEMP_DIR/${STRAIN}_all.blast")
    echo "📊 比对完成！共找到 $TOTAL_HITS 个潜在毒力因子"

    # 5.5 解析结果（前端实时刷新当前菌株的毒力因子数）
    echo -e "\n📝 正在解析并写入结果（实时更新）："
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
        SOURCE_FILE=$(find "$BAKTA_DIR/$STRAIN" -name "*.faa" | head -1)

        # 写入结果
        echo "$STRAIN,$PROT_ID,$VF_NAME,$SOURCE,$PID,$EVAL,$FUNC,$CAT,$REF,$SOURCE_FILE" >> "$OUTPUT"

        # 前端实时刷新计数
        ((STRAIN_VF_COUNT++))
        ((GLOBAL_VF_COUNT++))
        echo -ne "   当前菌株已识别：$STRAIN_VF_COUNT 个 | 全局累计：$GLOBAL_VF_COUNT 个\r"
        sleep 0.1  # 控制刷新速度，避免刷屏
    done < "$TEMP_DIR/${STRAIN}_all.blast"

    # 5.6 当前菌株处理完成（前端显示总结）
    echo -e "\n\n✅ 【$STRAIN】处理完成！"
    echo "   最终识别毒力因子：$STRAIN_VF_COUNT 个"
    echo "   其中实验验证（SetA）：$(grep "^$STRAIN," "$OUTPUT" | grep "SetA" | wc -l) 个"
    echo "   其中预测补充（SetB）：$(grep "^$STRAIN," "$OUTPUT" | grep "SetB" | wc -l) 个"
    echo "------------------------------------------------------------------------"
    echo "⏳ 剩余菌株：$((TOTAL_STRAINS - PROCESSED_STRAINS)) 株 | 预计剩余时间：$(( (TOTAL_STRAINS - PROCESSED_STRAINS) * 2 )) 分钟"

done <<< "$STRAIN_LIST"

# 6. 全局处理完成（前端显示最终统计）
rm -rf "$TEMP_DIR"
echo -e "\n========================================================================"
echo "🎉 所有菌株处理完毕！"
echo "📈 最终全局统计："
echo "   - 总处理菌株数：$PROCESSED_STRAINS/$TOTAL_STRAINS"
echo "   - 总识别毒力因子数：$GLOBAL_VF_COUNT 个"
echo "   - 实验验证因子（SetA）：$(grep "SetA" "$OUTPUT" | wc -l) 个"
echo "   - 预测补充因子（SetB）：$(grep "SetB" "$OUTPUT" | wc -l) 个"
echo "   - 完整结果文件：$OUTPUT"
echo "========================================================================"
