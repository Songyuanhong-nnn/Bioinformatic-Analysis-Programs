#!/bin/bash
set -euo pipefail

PROJECT_DIR="/mnt/d/WSL/disk/projects/VP1"
SEROTYPE_O="$PROJECT_DIR/kaptive_o_serotype_results.tsv"
SEROTYPE_K="$PROJECT_DIR/kaptive_k_serotype_results.tsv"
OUTPUT="$PROJECT_DIR/strain_serotype_hazard_map.csv"

# 表头
echo "菌株名,O血清型,K血清型,危害等级" > "$OUTPUT"

# 提取所有有效菌株（有基因组的）
ALL_STRAINS=$(find "$PROJECT_DIR/clean_genome" -name "*.fna" | xargs basename | sed 's/\.fna$//')

for STRAIN in $ALL_STRAINS; do
    # 匹配O血清型
    O_TYPE=$(grep -w "$STRAIN" "$SEROTYPE_O" | awk '{print $2}' || echo "未分型")
    # 匹配K血清型
    K_TYPE=$(grep -w "$STRAIN" "$SEROTYPE_K" | awk '{print $2}' || echo "未分型")
    # 危害等级（可按你的标准修改，比如O3/O4+K6为高危害，其他为中/低）
    if [[ $O_TYPE =~ "O3|O4|O1" && $K_TYPE == "K6" ]]; then
        HAZARD="高危害"
    elif [[ $O_TYPE =~ "O5|O6|O7" ]]; then
        HAZARD="中危害"
    else
        HAZARD="低危害"
    fi
    echo "$STRAIN,$O_TYPE,$K_TYPE,$HAZARD" >> "$OUTPUT"
done

echo "✅ 血清型对照表生成完成：$OUTPUT"
echo "📊 共匹配 $(wc -l "$OUTPUT" | awk '{print $1-1}') 株菌株"
# 预览各血清型分布
echo -e "\n各O血清型分布："
tail -n +2 "$OUTPUT" | awk -F "," '{print $2}' | sort | uniq -c | sort -nr
