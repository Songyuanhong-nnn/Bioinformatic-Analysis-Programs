#!/bin/bash
# 终极版批量毒力因子分析脚本
# 核心：避免文件误执行 + 适配所有环境 + 详细日志
set -euo pipefail  # 遇到错误立即退出，避免连锁报错

# ===================== 配置参数（无需修改）=====================
VFDB_DB="/mnt/d/WSL/disk/databases/VFDB/VFDB_setA"  # 数据库绝对路径（最稳妥）
CLEAN_GENOME_DIR="clean_genome"  # 基因组文件夹
THREADS=8
EVALUE=1e-10
MIN_ID=40  # 相似度≥40%
MIN_COVER=70  # 覆盖率≥70%
# ==================================================================

# 1. 检查必要文件是否存在
if [ ! -d "$CLEAN_GENOME_DIR" ]; then
  echo "❌ 错误：找不到 $CLEAN_GENOME_DIR 文件夹，请确认基因组文件放在该目录下"
  exit 1
fi

if [ ! -f "$VFDB_DB.dmnd" ]; then
  echo "❌ 错误：找不到 VFDB 数据库，请检查路径是否正确"
  exit 1
fi

# 2. 创建输出文件夹（强制清空旧文件）
rm -rf prokka_annotations vf_blast_results final_results
mkdir -p prokka_annotations vf_blast_results final_results

# 3. 安全遍历基因组文件（只处理 FNA/Fasta，绝对不执行文件）
genome_files=($CLEAN_GENOME_DIR/*.fna $CLEAN_GENOME_DIR/*.fasta)
if [ ${#genome_files[@]} -eq 0 ]; then
  echo "❌ 错误：$CLEAN_GENOME_DIR 下没有找到 .fna 或 .fasta 文件"
  exit 1
fi

# 4. 批量处理每个基因组
for genome_file in "${genome_files[@]}"; do
  # 跳过不存在的文件（避免通配符匹配失败）
  if [ ! -f "$genome_file" ]; then
    continue
  fi

  # 提取纯净菌株名（移除路径和后缀，替换特殊字符）
  strain_name=$(basename "$genome_file" | sed -e 's/\.[^.]*$//' -e 's/[^a-zA-Z0-9_]/_/g')
  echo -e "\n=================================================="
  echo "📌 开始处理菌株：$strain_name"
  echo "📂 基因组文件：$genome_file"
  echo "=================================================="

  # 5. Prokka 注释（强制覆盖，避免文件夹冲突）
  echo "🔧 正在进行 Prokka 基因注释..."
  prokka --outdir "prokka_annotations/$strain_name" \
         --prefix "$strain_name" \
         --kingdom Bacteria \
         --genus Vibrio \
         --species parahaemolyticus \
         --cpus $THREADS \
         --force \
         "$genome_file" 2>&1 | tee "prokka_annotations/$strain_name/annotation.log"

  # 检查 Prokka 是否生成了蛋白序列文件（关键！）
  protein_file="prokka_annotations/$strain_name/$strain_name.faa"
  if [ ! -f "$protein_file" ] || [ ! -s "$protein_file" ]; then
    echo "⚠️  警告：Prokka 未生成有效蛋白序列文件，跳过该菌株"
    continue
  fi

  # 6. Diamond 毒力因子比对（绝对路径 + 旧版本兼容参数）
  echo "🔍 正在进行毒力因子比对..."
  diamond blastp \
          --db "$VFDB_DB" \
          --query "$protein_file" \
          --out "vf_blast_results/${strain_name}_vf_results.tsv" \
          --outfmt 6 qseqid sseqid pident length evalue qcovhsp stitle \
          --evalue $EVALUE \
          --id $MIN_ID \
          --query-cover $MIN_COVER \
          --threads $THREADS 2>&1 | tee "vf_blast_results/${strain_name}_blast.log"

  # 7. 整理结果（处理有无匹配的情况）
  result_file="final_results/${strain_name}_virulence_factors.tsv"
  if [ -s "vf_blast_results/${strain_name}_vf_results.tsv" ]; then
    # 有匹配结果：整理字段
    awk -v strain="$strain_name" '{
      print strain "\t" $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7
    }' "vf_blast_results/${strain_name}_vf_results.tsv" > "$result_file"
    vf_count=$(wc -l < "$result_file")
    echo "✅ 处理完成！检测到 $vf_count 个毒力因子"
  else
    # 无匹配结果：生成占位信息
    echo -e "$strain_name\t无\t无\t无\t无\t无\t无\t未检测到已知毒力因子" > "$result_file"
    echo "✅ 处理完成！未检测到已知毒力因子"
  fi
done

# 8. 生成最终汇总表（带清晰表头，可直接用 Excel 打开）
echo -e "\n=================================================="
echo "📊 正在生成所有菌株汇总表..."
echo "=================================================="
summary_file="final_results/all_strains_virulence_summary.tsv"
echo -e "菌株名\t菌株蛋白ID\tVFDB毒力因子ID\t相似度(%)\t序列长度\t置信度(evalue)\t覆盖率(%)\t毒力因子功能描述" > "$summary_file"

# 合并所有菌株结果到汇总表
for result_file in final_results/*_virulence_factors.tsv; do
  if [ -f "$result_file" ]; then
    cat "$result_file" >> "$summary_file"
  fi
done

echo -e "\n🎉 所有菌株分析完成！"
echo "📁 结果存放目录：final_results/"
echo "📋 汇总表：$summary_file（推荐用 Excel 打开）"
echo "🔍 单个菌株详情：final_results/[菌株名]_virulence_factors.tsv"
echo "📝 注释日志：prokka_annotations/[菌株名]/annotation.log"
