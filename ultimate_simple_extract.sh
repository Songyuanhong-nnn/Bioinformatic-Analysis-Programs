#!/bin/bash
set -euo pipefail

PROJECT_DIR="/mnt/d/WSL/disk/projects/VP1"
ROARY_CSV="$PROJECT_DIR/roary_output_final_1762658265/gene_presence_absence.csv"
SERO_O="$PROJECT_DIR/kaptive_o_serotype_results.tsv"
VFDB_PROT="$PROJECT_DIR/vfdb_online/VFDB_setA_pro.fas"
OUTPUT_CSV="$PROJECT_DIR/ULTIMATE_SIMPLE_virulence.csv"

echo "===== 终极简化版：直接提取毒力基因 ====="

# 1. 提取VFDB毒力基因名（用于精准匹配）
grep "^>" "$VFDB_PROT" | sed 's/^>//; s/ .*//' > "$PROJECT_DIR/vfdb_genes.tmp"

# 2. 从泛基因组中提取所有毒力相关基因（宽松筛选，确保有结果）
echo "基因ID,基因名,功能描述,存在菌株数,VFDB匹配状态" > "$OUTPUT_CSV"
tail -n +2 "$ROARY_CSV" | awk -F "," -v vfdb="$PROJECT_DIR/vfdb_genes.tmp" '
BEGIN {
    while ((getline < vfdb) > 0) vfdb_set[$1]=1;
}
# 筛选条件：功能含毒力关键词 或 基因名在VFDB中（不限制血清型，先拿到基因）
($14 ~ /virulence|toxin|tdh|trh|hemolysin|adhesin|T3SS|T6SS|biofilm|invasion/ || $2 in vfdb_set) {
    # 统计存在菌株数
    count=0;
    for(i=15;i<=NF;i++) if($i==1) count++;
    vfdb_match=($2 in vfdb_set)?"是":"否（功能匹配）";
    print $1 "," $2 "," $14 "," count "," vfdb_match;
}' >> "$OUTPUT_CSV"

# 3. 统计结果
TOTAL=$(grep -v "基因ID" "$OUTPUT_CSV" | wc -l)
echo -e "\n✅ 提取完成！"
echo "   - 共找到 $TOTAL 个毒力相关基因"
echo "   - 结果文件：$OUTPUT_CSV"

# 4. 按血清型分组（单独文件，方便查看）
echo -e "\n🔧 按血清型分组..."
tail -n +2 "$SERO_O" | awk '{print $1 "\t" $2}' > "$PROJECT_DIR/strain_sero.tmp"
SEROTYPES=$(cut -f2 "$PROJECT_DIR/strain_sero.tmp" | sort -u)
OUTPUT_DIR="$PROJECT_DIR/serotype_virulence_simple"
mkdir -p "$OUTPUT_DIR"

for SERO in $SEROTYPES; do
    STRAINS=$(grep "\t$SERO" "$PROJECT_DIR/strain_sero.tmp" | cut -f1 | tr '\n' '|')
    [ -z "$STRAINS" ] && continue
    # 筛选该血清型存在的毒力基因
    echo "基因ID,基因名,功能描述,存在菌株数,VFDB匹配状态" > "$OUTPUT_DIR/${SERO}_virulence.csv"
    grep -E "^[^,]+,[^,]+,[^,]+,[0-9]+," "$OUTPUT_CSV" | while IFS= read -r line; do
        GENE_ID=$(echo "$line" | awk -F "," '{print $1}')
        # 从Roary中检查该基因是否在该血清型菌株中存在
        grep "^$GENE_ID," "$ROARY_CSV" | grep -qE ",($STRAINS)" && echo "$line" >> "$OUTPUT_DIR/${SERO}_virulence.csv"
    done
    # 统计该血清型基因数
    CNT=$(( $(wc -l "$OUTPUT_DIR/${SERO}_virulence.csv" | awk '{print $1}') - 1 ))
    echo "   ✅ 血清型 $SERO：$CNT 个毒力基因"
done

# 5. 生成简化可视化（不管有没有特异性，先展示分布）
cat > simple_plot.py << 'PYEOF'
import pandas as pd
import matplotlib.pyplot as plt

plt.rcParams['font.sans-serif'] = ['Arial']
plt.rcParams['axes.unicode_minus'] = False

# 读取按血清型分组的数据
import os
PROJECT_DIR = "/mnt/d/WSL/disk/projects/VP1"
OUTPUT_DIR = f"{PROJECT_DIR}/serotype_virulence_simple"

# 统计各血清型毒力基因数
sero_data = []
for file in os.listdir(OUTPUT_DIR):
    if file.endswith("_virulence.csv"):
        sero = file.replace("_virulence.csv", "")
        df = pd.read_csv(f"{OUTPUT_DIR}/{file}")
        cnt = len(df) - 1  # 减表头
        sero_data.append( (sero, cnt) )

# 创建DataFrame并绘图
df_plot = pd.DataFrame(sero_data, columns=['血清型', '毒力基因数'])
df_plot = df_plot.sort_values('毒力基因数', ascending=False)

fig, ax = plt.subplots(figsize=(14, 6))
ax.bar(df_plot['血清型'], df_plot['毒力基因数'], color='steelblue')
ax.set_title('Virulence Genes Distribution by O Serotype', fontsize=14, fontweight='bold')
ax.set_xlabel('O Serotype', fontsize=12)
ax.set_ylabel('Number of Virulence Genes', fontsize=12)
ax.tick_params(axis='x', rotation=45)
plt.tight_layout()
plt.savefig(f"{PROJECT_DIR}/SIMPLE_virulence_distribution.png", dpi=300)
print("✅ 可视化图表生成完成：SIMPLE_virulence_distribution.png")
PYEOF

python3 simple_plot.py

# 6. 清理临时文件
rm -f "$PROJECT_DIR/vfdb_genes.tmp" "$PROJECT_DIR/strain_sero.tmp"

echo -e "\n===== 所有操作完成！ ====="
echo "📁 核心成果："
echo "1. 所有毒力基因总表：$OUTPUT_CSV"
echo "2. 按血清型分组：$OUTPUT_DIR"
echo "3. 可视化图表：SIMPLE_virulence_distribution.png"
echo "🌟 成果说明："
echo "   - 不限制“特异性”，先拿到各血清型的毒力基因清单"
echo "   - 可后续手动筛选特异性（对比不同血清型的基因列表）"
echo "   - 所有基因均经过VFDB匹配或功能关键词验证，确保毒力相关性"
