#!/bin/bash
set -euo pipefail

PROJECT_DIR="/mnt/d/WSL/disk/projects/VP1"
ROARY_CSV="$PROJECT_DIR/roary_output_final_1762658265/gene_presence_absence.csv"
SERO_O="$PROJECT_DIR/kaptive_o_serotype_results.tsv"
VFDB_PROT="$PROJECT_DIR/vfdb_online/VFDB_setA_pro.fas"
OUTPUT_CSV="$PROJECT_DIR/FAST_serotype_specific_virulence.csv"

echo "===== 简化版：快速筛选特异性毒力基因 ====="

# 1. 快速提取关键数据
grep "^>" "$VFDB_PROT" | sed 's/^>//; s/ .*//' > "$PROJECT_DIR/vfdb_genes.tmp"
tail -n +2 "$SERO_O" | awk '{print $1 "\t" $2}' > "$PROJECT_DIR/strain_sero.tmp"

# 2. 生成血清型-菌株映射字典（简化匹配逻辑）
declare -A SERO_STRAINS
while read -r STRAIN SERO; do
    SERO_STRAINS[$SERO]+="$STRAIN|"
done < "$PROJECT_DIR/strain_sero.tmp"

# 3. 快速筛选：毒力相关+仅单一血清型存在
echo "血清型,特异性毒力基因ID,基因名,功能描述,VFDB匹配状态" > "$OUTPUT_CSV"
tail -n +2 "$ROARY_CSV" | awk -F "," -v vfdb="$PROJECT_DIR/vfdb_genes.tmp" '
BEGIN {
    while ((getline < vfdb) > 0) vfdb_set[$1]=1;
    # 读取血清型-菌株映射
    while ((getline < "'"$PROJECT_DIR/strain_sero.tmp"'") > 0) {
        strain_to_sero[$1]=$2;
        sero_count[$2]++
    }
}
($14 ~ /virulence|toxin|tdh|trh|hemolysin|adhesin|T3SS|T6SS/ || $2 in vfdb_set) {
    # 统计该基因存在的血清型
    delete sero_set;
    for (i=15;i<=NF;i++) {
        if ($i==1) {
            strain=gensub(/"/,"","g",$(i+14));  # 去除菌株名引号
            sero=strain_to_sero[strain];
            sero_set[sero]=1;
        }
    }
    # 仅单一血清型存在
    if (length(sero_set)==1) {
        sero=first_key(sero_set);
        vfdb_match=($2 in vfdb_set)?"是":"否（功能匹配）";
        print sero "," $1 "," $2 "," $14 "," vfdb_match;
    }
}
function first_key(arr,    k) {
    for (k in arr) return k;
}' >> "$OUTPUT_CSV"

# 4. 统计结果
TOTAL=$(grep -v "血清型" "$OUTPUT_CSV" | wc -l)
echo -e "\n✅ 快速筛选完成！"
echo "   - 共找到 $TOTAL 个血清型特异性毒力基因"
echo "   - 结果文件：$OUTPUT_CSV"

# 5. 清理临时文件
rm -f "$PROJECT_DIR/vfdb_genes.tmp" "$PROJECT_DIR/strain_sero.tmp"
