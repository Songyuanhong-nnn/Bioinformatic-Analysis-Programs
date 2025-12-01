#!/bin/bash
set -euo pipefail

PROJECT_DIR="/mnt/d/WSL/disk/projects/VP1"
SEROTYPE_O="$PROJECT_DIR/kaptive_o_serotype_results.tsv"
SEROTYPE_K="$PROJECT_DIR/kaptive_k_serotype_results.tsv"
OUTPUT="$PROJECT_DIR/strain_serotype_hazard_map.csv"

# 清理旧文件
> "$OUTPUT"

# 表头
echo "菌株名,O血清型,K血清型,危害等级" >> "$OUTPUT"

# 提取所有有效菌株（有基因组的），用while循环避免文件名空格问题
find "$PROJECT_DIR/clean_genome" -name "*.fna" | while read -r GENOME_FILE; do
    # 稳健提取菌株名（去掉路径和.fna后缀）
    STRAIN=$(basename -- "$GENOME_FILE" .fna)
    [ -z "$STRAIN" ] && continue

    # 匹配O血清型（精确匹配菌株名，避免部分匹配）
    O_TYPE=$(grep -w "^$STRAIN" "$SEROTYPE_O" | awk '{print $2}' || echo "未分型")
    # 匹配K血清型（精确匹配菌株名）
    K_TYPE=$(grep -w "^$STRAIN" "$SEROTYPE_K" | awk '{print $2}' || echo "未分型")

    # 危害等级判定（可按你的实际标准修改！）
    # 示例规则：O3/O4/O1 + K6 = 高危害；O5/O6/O7 = 中危害；其他=低危害
    if [[ ( "$O_TYPE" == "O3" || "$O_TYPE" == "O4" || "$O_TYPE" == "O1" ) && "$K_TYPE" == "K6" ]]; then
        HAZARD="高危害"
    elif [[ "$O_TYPE" == "O5" || "$O_TYPE" == "O6" || "$O_TYPE" == "O7" ]]; then
        HAZARD="中危害"
    else
        HAZARD="低危害"
    fi

    # 输出到文件（用引号包裹，避免逗号干扰）
    echo "\"$STRAIN\",\"$O_TYPE\",\"$K_TYPE\",\"$HAZARD\"" >> "$OUTPUT"
done

# 统计结果
TOTAL_STRAINS=$(( $(wc -l "$OUTPUT" | awk '{print $1}') - 1 ))
echo "✅ 血清型对照表生成完成：$OUTPUT"
echo "📊 共匹配 $TOTAL_STRAINS 株菌株"

# 预览各O血清型分布（去重统计）
echo -e "\n各O血清型分布（菌株数从多到少）："
tail -n +2 "$OUTPUT" | awk -F '","' '{print $2}' | sort | uniq -c | sort -nr

# 预览各危害等级分布
echo -e "\n各危害等级分布："
tail -n +2 "$OUTPUT" | awk -F '","' '{print $4}' | sort | uniq -c | sort -nr
