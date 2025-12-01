#!/bin/bash
set -euo pipefail

PROJECT_DIR="/mnt/d/WSL/disk/projects/VP1"
ROARY_FILE="$PROJECT_DIR/roary_output_final_1762658265/gene_presence_absence.csv"
SEROTYPE_O="$PROJECT_DIR/kaptive_o_serotype_results.tsv"
SEROTYPE_K="$PROJECT_DIR/kaptive_k_serotype_results.tsv"
VFDB="$PROJECT_DIR/vfdb_online/VFDB_setA_pro.fas"
OUTPUT="$PROJECT_DIR/final_high_risk_unique_markers.csv"

# 1. 提取高危害菌株列表（O3/K6血清型）
echo "🔍 提取高危害菌株（O3/K6血清型）："
grep -E "O3|K6" "$SEROTYPE_O" "$SEROTYPE_K" | awk '{print $1}' | sort -u > "$PROJECT_DIR/high_risk_strains.txt"
grep -v -E "O3|K6" "$SEROTYPE_O" "$SEROTYPE_K" | awk '{print $1}' | sort -u > "$PROJECT_DIR/low_risk_strains.txt"

high_count=$(wc -l "$PROJECT_DIR/high_risk_strains.txt" | awk '{print $1}')
low_count=$(wc -l "$PROJECT_DIR/low_risk_strains.txt" | awk '{print $1}')
echo "高危害菌株数：$high_count，低危害菌株数：$low_count"

# 2. 提取VFDB毒力基因名（用于匹配）
grep "^>" "$VFDB" | sed 's/^>//; s/ .*//' > "$PROJECT_DIR/vfdb_genes.txt"

# 3. 从Roary结果中筛选：高危害全有、低危害全无、且属于毒力基因的标志物
echo -e "\n🔍 筛选高危害特有毒力基因标志物："
python3 - << END
import csv

roary_file = "$ROARY_FILE"
high_strains = set(open("$PROJECT_DIR/high_risk_strains.txt").read().split())
low_strains = set(open("$PROJECT_DIR/low_risk_strains.txt").read().split())
vfdb_genes = set(open("$PROJECT_DIR/vfdb_genes.txt").read().split())

with open(roary_file, 'r') as f, open("$OUTPUT", 'w', newline='') as out_f:
    reader = csv.DictReader(f)
    fieldnames = ['基因名', '基因ID', '高危害菌株覆盖数', '低危害菌株覆盖数', '是否毒力基因', '功能描述']
    writer = csv.DictWriter(out_f, fieldnames=fieldnames)
    writer.writeheader()

    for row in reader:
        gene_name = row['Gene']
        gene_id = row['Non-unique Gene name']
        
        # 统计高/低危害菌株中该基因的存在情况（1=存在，0=缺失）
        high_present = sum(1 for strain in high_strains if strain in row and row[strain] == '1')
        low_present = sum(1 for strain in low_strains if strain in row and row[strain] == '1')
        
        # 筛选条件：高危害全有、低危害全无、属于毒力基因
        is_virulence = '是' if gene_name in vfdb_genes or (gene_id and gene_id in vfdb_genes) else '否'
        if high_present == len(high_strains) and low_present == 0 and is_virulence == '是':
            writer.writerow({
                '基因名': gene_name,
                '基因ID': gene_id if gene_id else '未知',
                '高危害菌株覆盖数': high_present,
                '低危害菌株覆盖数': low_present,
                '是否毒力基因': is_virulence,
                '功能描述': row['Annotation'] if 'Annotation' in row else '未知'
            })

print("筛选完成！结果文件：$OUTPUT")
END

# 统计结果
total_markers=$(( $(wc -l "$OUTPUT" | awk '{print $1}') - 1 ))
echo -e "\n🎉 最终筛选结果："
echo "📁 结果文件：$OUTPUT"
echo "🌟 高危害特有毒力基因标志物数：$total_markers 个"

# 预览前10个核心标志物
echo -e "\n前10个核心标志物预览："
head -11 "$OUTPUT" | tail -10

# 清理临时文件
rm -f "$PROJECT_DIR/high_risk_strains.txt" "$PROJECT_DIR/low_risk_strains.txt" "$PROJECT_DIR/vfdb_genes.txt"
