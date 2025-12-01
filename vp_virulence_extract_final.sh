#!/bin/bash
set -euo pipefail

# 路径100%匹配你的项目结构
PROJECT_DIR="/mnt/d/WSL/disk/projects/VP1"
BAKTA_FAA_DIR="$PROJECT_DIR/bakta_annotations"  # 菌株蛋白序列
VFDB_A="$PROJECT_DIR/vfdb_online/VFDB_setA_pro.fas"  # 已复制的SetA库
VFDB_B="$PROJECT_DIR/vfdb_online/VFDB_setB_pro.fas"  # 已复制的SetB库
VFDB_ANNO="$PROJECT_DIR/vfdb_online/VFs.xls"  # 已复制的官方注释表
OUTPUT="$PROJECT_DIR/VP_VIRULENCE_FACTORS_FINAL.csv"  # 最终结果

# 初始化结果表头（含官方溯源）
echo "菌株名,蛋白ID,毒力因子名,VFDB来源,相似度(%),E值,功能描述,VFDB分类,参考文献" > "$OUTPUT"

# 遍历所有菌株，批量比对
for STRAIN_DIR in "$BAKTA_FAA_DIR"/*; do
    STRAIN=$(basename "$STRAIN_DIR")
    STRAIN_FAA="$STRAIN_DIR/$STRAIN.faa"
    
    # 容错：跳过无蛋白文件的菌株
    [ ! -f "$STRAIN_FAA" ] && { echo "⚠️  跳过$STRAIN：无蛋白文件"; continue; }
    echo "🔧 处理菌株：$STRAIN"

    # 1. 与SetA库比对（权威实验验证因子）
    blastp -query "$STRAIN_FAA" \
           -db "$VFDB_A" \
           -outfmt "6 qseqid sseqid pident evalue" \
           -evalue 1e-5 \
           -num_threads 4 \
           -max_target_seqs 1 2>/dev/null > /tmp/setA.blast

    # 2. 与SetB库比对（补充预测因子）
    blastp -query "$STRAIN_FAA" \
           -db "$VFDB_B" \
           -outfmt "6 qseqid sseqid pident evalue" \
           -evalue 1e-5 \
           -num_threads 4 \
           -max_target_seqs 1 2>/dev/null > /tmp/setB.blast

    # 3. 合并去重（SetA优先，避免重复）
    cat /tmp/setA.blast > /tmp/all.blast
    grep -v -f <(cut -f1 /tmp/setA.blast) /tmp/setB.blast >> /tmp/all.blast

    # 4. 解析结果+关联官方注释
    while read -r BLAST_LINE; do
        PROT_ID=$(echo "$BLAST_LINE" | cut -f1)
        VF_NAME=$(echo "$BLAST_LINE" | cut -f2)
        PID=$(printf "%.2f" $(echo "$BLAST_LINE" | cut -f3))
        EVAL=$(echo "$BLAST_LINE" | cut -f4)

        # 判断毒力因子来源（SetA/SetB）
        if grep -w "^$VF_NAME" "$VFDB_A" > /dev/null 2>&1; then
            VF_SOURCE="SetA（实验验证）"
        else
            VF_SOURCE="SetB（预测补充）"
        fi

        # 从官方注释表提取功能、分类、参考文献
        ANNO=$(grep -w "^$VF_NAME" "$VFDB_ANNO" | head -1)
        if [ -n "$ANNO" ]; then
            FUNCTION=$(echo "$ANNO" | cut -f4 | sed 's/,/;/g')  # 避免CSV格式错乱
            CATEGORY=$(echo "$ANNO" | cut -f5)
            REF=$(echo "$ANNO" | cut -f7 | sed 's/,/;/g')
        else
            FUNCTION="无"
            CATEGORY="无"
            REF="无"
        fi

        # 写入最终结果
        echo "$STRAIN,$PROT_ID,$VF_NAME,$VF_SOURCE,$PID,$EVAL,$FUNCTION,$CATEGORY,$REF" >> "$OUTPUT"
    done < /tmp/all.blast

    # 输出当前菌株处理结果
    VF_COUNT=$(grep "^$STRAIN," "$OUTPUT" | wc -l)
    echo "   ✅ 完成：$STRAIN 找到 $VF_COUNT 个毒力因子"
done

# 清理临时文件
rm -f /tmp/setA.blast /tmp/setB.blast /tmp/all.blast

echo -e "\n🎉 所有菌株处理完成！最终结果文件：$OUTPUT"
echo "✅ 结果特点：基于VFDB官方库，含实验验证+预测因子，带完整功能/参考文献"
