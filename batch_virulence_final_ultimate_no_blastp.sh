#!/bin/bash
# 终极版批量毒力因子分析（跳过 blastp，绕开版本检测）
set -euo pipefail

# ===================== 配置参数（无需修改）=====================
VFDB_DB="/mnt/d/WSL/disk/databases/VFDB/VFDB_setA"
CLEAN_GENOME_DIR="clean_genome"
THREADS=8
EVALUE=1e-10
MIN_ID=40
MIN_COVER=70
# ==================================================================

# 1. 检查必要文件
if [ ! -d "$CLEAN_GENOME_DIR" ]; then
  echo "❌ 错误：找不到 $CLEAN_GENOME_DIR 文件夹"
  exit 1
fi
if [ ! -f "$VFDB_DB.dmnd" ]; then
  echo "❌ 错误：找不到 VFDB 数据库"
  exit 1
fi

# 2. 创建输出文件夹
rm -rf prokka_annotations vf_blast_results final_results
mkdir -p prokka_annotations vf_blast_results final_results

# 3. 遍历基因组文件
genome_files=($CLEAN_GENOME_DIR/*.fna $CLEAN_GENOME_DIR/*.fasta)
if [ ${#genome_files[@]} -eq 0 ]; then
  echo "❌ 错误：没有找到 .fna 或 .fasta 文件"
  exit 1
fi

# 4. 批量处理
for genome_file in "${genome_files[@]}"; do
  if [ ! -f "$genome_file" ]; then
    continue
  fi
  strain_name=$(basename "$genome_file" | sed -e 's/\.[^.]*$//' -e 's/[^a-zA-Z0-9_]/_/g')
  echo -e "\n=================================================="
  echo "📌 开始处理菌株：$strain_name"
  echo "📂 基因组文件：$genome_file"
  echo "=================================================="

  # 5. Prokka 注释（关键修改：--noanno 跳过 blastp 功能注释，只做基因预测）
  echo "🔧 正在进行 Prokka 基因预测（跳过 blastp）..."
  prokka --outdir "prokka_annotations/$strain_name" \
         --prefix "$strain_name" \
         --kingdom Bacteria \
         --genus Vibrio \
         --species parahaemolyticus \
         --cpus $THREADS \
         --force \
         --noanno \  # 核心参数：跳过 blastp 比对，绕开版本检测
         "$genome_file" 2>&1 | tee "prokka_annotations/$strain_name/annotation.log"

  # 6. 检查蛋白序列文件
  protein_file="prokka_annotations/$strain_name/$strain_name.faa"
  if [ ! -f "$protein_file" ] || [ ! -s "$protein_file" ]; then
    echo "⚠️  警告：未生成有效蛋白序列文件，跳过该菌株"
    continue
  fi

  # 7. Diamond 毒力因子比对（不受影响）
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

  # 8. 整理结果
  result_file="final_results/${strain_name}_virulence_factors.tsv"
  if [ -s "vf_blast_results/${strain_name}_vf_results.tsv" ]; then
    awk -v strain="$strain_name" '{
      print strain "\t" $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7
    }' "vf_blast_results/${strain_name}_vf_results.tsv" > "$result_file"
    vf_count=$(wc -l < "$result_file")
    echo "✅ 处理完成！检测到 $vf_count 个毒力因子"
  else
    echo -e "$strain_name\t无\t无\t无\t无\t无\t无\t未检测到已知毒力因子" > "$result_file"
    echo "✅ 处理完成！未检测到已知毒力因子"
  fi
done

# 9. 生成汇总表
summary_file="final_results/all_strains_virulence_summary.tsv"
echo -e "菌株名\t菌株蛋白ID\tVFDB毒力因子ID\t相似度(%)\t序列长度\t置信度(evalue)\t覆盖率(%)\t毒力因子功能描述" > "$summary_file"
for result_file in final_results/*_virulence_factors.tsv; do
  if [ -f "$result_file" ]; then
    cat "$result_file" >> "$summary_file"
  fi
done

echo -e "\n🎉 所有菌株分析完成！"
echo "📁 结果存放目录：final_results/"
echo "📋 汇总表：$summary_file"
