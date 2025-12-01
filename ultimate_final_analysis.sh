#!/bin/bash
set -euo pipefail

PROJECT_DIR="/mnt/d/WSL/disk/projects/VP1"
echo "===== 开始最终分析：筛选血清型特异性毒力基因 ====="

# 1. 定义已有的核心文件（直接调用，不重新生成）
ROARY_CSV="$PROJECT_DIR/roary_output_final_1762658265/gene_presence_absence.csv"
SERO_O="$PROJECT_DIR/kaptive_o_serotype_results.tsv"
VFDB_PROT="$PROJECT_DIR/vfdb_online/VFDB_setA_pro.fas"
OUTPUT_CSV="$PROJECT_DIR/ULTIMATE_serotype_specific_virulence.csv"
PLOT_PY="$PROJECT_DIR/ultimate_plot.py"

# 2. 提取VFDB毒力基因名+菌株-血清型对应关系
echo "🔧 准备基础数据..."
grep "^>" "$VFDB_PROT" | sed 's/^>//; s/ .*//' > "$PROJECT_DIR/vfdb_genes.tmp"
tail -n +2 "$SERO_O" | awk '{print $1 "\t" $2}' > "$PROJECT_DIR/strain_sero.tmp"
SEROTYPES=$(cut -f2 "$PROJECT_DIR/strain_sero.tmp" | sort -u)
head -1 "$ROARY_CSV" | awk -F "," '{for(i=15;i<=NF;i++) print $i}' > "$PROJECT_DIR/roary_strains.tmp"

# 3. 筛选血清型特异性毒力基因（核心逻辑）
echo "🔍 筛选特异性毒力基因..."
echo "血清型,特异性毒力基因ID,基因名,功能描述,VFDB匹配状态,该血清型存在菌株数,总存在菌株数" > "$OUTPUT_CSV"

for SERO in $SEROTYPES; do
    # 该血清型菌株列表
    SERO_STRAINS=$(grep "\t$SERO" "$PROJECT_DIR/strain_sero.tmp" | cut -f1 | tr '\n' '|')
    SERO_COUNT=$(echo "$SERO_STRAINS" | tr '|' '\n' | wc -l)
    [ "$SERO_COUNT" -eq 0 ] && continue

    # 其他血清型菌株列表
    OTHER_STRAINS=$(grep -v "\t$SERO" "$PROJECT_DIR/strain_sero.tmp" | cut -f1 | tr '\n' '|')
    [ -z "$OTHER_STRAINS" ] && OTHER_STRAINS="无"

    # 筛选满足条件的基因：毒力相关+仅该血清型存在
    tail -n +2 "$ROARY_CSV" | awk -F "," -v s="$SERO_STRAINS" -v o="$OTHER_STRAINS" -v vfdb="$PROJECT_DIR/vfdb_genes.tmp" '
    BEGIN {
        while ((getline < vfdb) > 0) vfdb_set[$1]=1;
        split(s, sero_strains, "|");
        split(o, other_strains, "|");
    }
    ($14 ~ /virulence|toxin|tdh|trh|hemolysin|adhesin|T3SS|T6SS|biofilm/ || $2 in vfdb_set) {
        # 统计该血清型存在数
        sero_cnt=0;
        for (i in sero_strains) {
            if (sero_strains[i] != "") {
                for (j=15;j<=NF;j++) {
                    if ($j == 1 && $0 ~ sero_strains[i]) {sero_cnt++; break;}
                }
            }
        }
        # 统计其他血清型存在数
        other_cnt=0;
        if (o != "无") {
            for (i in other_strains) {
                if (other_strains[i] != "") {
                    for (j=15;j<=NF;j++) {
                        if ($j == 1 && $0 ~ other_strains[i]) {other_cnt++; break;}
                    }
                }
            }
        }
        # 满足特异性条件
        if (sero_cnt > 0 && other_cnt == 0) {
            vfdb_match = ($2 in vfdb_set) ? "是" : "否（功能匹配）";
            total_cnt=sero_cnt+other_cnt;
            print $1 "," $2 "," $14 "," vfdb_match "," sero_cnt "," total_cnt;
        }
    }' >> "$OUTPUT_CSV"

    # 统计该血清型结果
    SPEC_CNT=$(( $(grep "^$SERO," "$OUTPUT_CSV" | wc -l) ))
    echo "   ✅ 血清型 $SERO：找到 $SPEC_CNT 个特异性毒力基因"
done

# 4. 生成验证报告
echo -e "\n===== 结果验证 ====="
echo "📊 全局统计："
TOTAL_SPEC=$(wc -l "$OUTPUT_CSV" | awk '{print $1-1}')
SERO_WITH_SPEC=$(grep -v "血清型" "$OUTPUT_CSV" | awk -F "," '{print $1}' | sort -u | wc -l)
VFDB_MATCH_CNT=$(grep "是" "$OUTPUT_CSV" | wc -l)
VFDB_RATE=$(echo "scale=2; $VFDB_MATCH_CNT/$TOTAL_SPEC*100" | bc 2>/dev/null || echo 0)

echo "   - 总特异性毒力基因数：$TOTAL_SPEC"
echo "   - 含特异性基因的血清型数：$SERO_WITH_SPEC"
echo "   - VFDB直接匹配数：$VFDB_MATCH_CNT"
echo "   - VFDB匹配率：$VFDB_RATE%"

# 5. 生成可视化脚本
cat > "$PLOT_PY" << 'PYEOF'
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

plt.rcParams['font.sans-serif'] = ['Arial']
plt.rcParams['axes.unicode_minus'] = False

# 读取数据
df = pd.read_csv("/mnt/d/WSL/disk/projects/VP1/ULTIMATE_serotype_specific_virulence.csv")

# 统计各血清型特异性基因数
sero_count = df['血清型'].value_counts()

# 创建图表（2个子图：柱状图+饼图）
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))

# 柱状图：各血清型特异性基因数
sero_count.plot(kind='bar', ax=ax1, color='steelblue', edgecolor='black')
ax1.set_title('Number of Serotype-Specific Virulence Genes', fontsize=14, fontweight='bold')
ax1.set_xlabel('O Serotype', fontsize=12)
ax1.set_ylabel('Number of Specific Virulence Genes', fontsize=12)
ax1.tick_params(axis='x', rotation=45, labelsize=10)
ax1.tick_params(axis='y', labelsize=10)

# 饼图：各血清型特异性基因占比（仅显示有基因的血清型）
sero_count.plot(kind='pie', ax=ax2, autopct='%1.1f%%', startangle=90, 
                colors=sns.color_palette('Set2'), textprops={'fontsize': 10})
ax2.set_title('Proportion of Serotype-Specific Virulence Genes', fontsize=14, fontweight='bold')
ax2.set_ylabel('')

# 调整布局
plt.tight_layout()
plt.savefig("/mnt/d/WSL/disk/projects/VP1/ULTIMATE_serotype_specific_virulence_plot.png", 
            dpi=300, bbox_inches='tight')
plt.close()

# 生成统计摘要
summary = f"""
===== 血清型特异性毒力基因分析摘要 =====
1. 分析基础：57株副溶血性弧菌，{len(sero_count)}种O血清型
2. 核心成果：共筛选到 {len(df)} 个血清型特异性毒力基因
3. 分布特征：
   - 特异性基因最多的血清型：{sero_count.index[0]}（{sero_count.iloc[0]} 个）
   - 所有基因均满足“仅目标血清型存在”，特异性100%
4. 可靠性：
   - VFDB毒力库直接匹配率：{df['VFDB匹配状态'].value_counts()['是']/len(df)*100:.1f}%
   - 其余基因为功能关键词匹配（毒力相关）
5. 应用价值：可直接作为血清型分型诊断标志物、毒力差异分析靶点
"""

with open("/mnt/d/WSL/disk/projects/VP1/ULTIMATE_analysis_summary.txt", "w") as f:
    f.write(summary)

print("✅ 可视化报告生成完成：ULTIMATE_serotype_specific_virulence_plot.png")
print("✅ 分析摘要生成完成：ULTIMATE_analysis_summary.txt")
PYEOF

# 6. 运行可视化脚本（安装依赖+执行）
echo -e "\n🔧 生成可视化报告..."
conda activate bakta
conda install -c conda-forge pandas matplotlib seaborn -y > /dev/null 2>&1
python3 "$PLOT_PY"

# 7. 清理临时文件
rm -f "$PROJECT_DIR/vfdb_genes.tmp" "$PROJECT_DIR/strain_sero.tmp" "$PROJECT_DIR/roary_strains.tmp" "$PLOT_PY"

echo -e "\n===== 最终分析完成！ ====="
echo "📁 核心成果文件清单："
echo "1. 特异性毒力基因清单：$OUTPUT_CSV"
echo "2. 可视化图表：ULTIMATE_serotype_specific_virulence_plot.png"
echo "3. 分析摘要：ULTIMATE_analysis_summary.txt"
echo -e "\n🌟 所有成果可直接用于："
echo "   - 血清型分型诊断标志物开发"
echo "   - 论文数据支撑与图表展示"
echo "   - 后续实验验证靶点选择"
