#!/bin/bash
# 最终版：Prokka+VFDB双在线数据库 + 修复文件路径解析问题
set -euo pipefail

# ===================== 强制指定绝对路径（关键修复）=====================
# 避免WSL路径解析异常，直接写死绝对路径（根据你的实际路径修改，当前已适配你的环境）
PROKKA_ONLINE_DB="/mnt/d/WSL/disk/projects/VP1/prokka_online_db"
VFDB_ONLINE_DIR="/mnt/d/WSL/disk/projects/VP1/vfdb_online"
CLEAN_GENOME_DIR="/mnt/d/WSL/disk/projects/VP1/clean_genome"
VFDB_DB="$VFDB_ONLINE_DIR/VFDB_setA"

# 其他配置参数
THREADS=8
EVALUE=1e-10
MIN_ID=40
MIN_COVER=70
# ==================================================================

# --------------------------
# 第一步：在线更新Prokka注释数据库
# --------------------------
echo "🔄 正在更新Prokka在线数据库（首次运行较慢）..."
mkdir -p "$PROKKA_ONLINE_DB"
export PROKKA_DB="$PROKKA_ONLINE_DB"  # 强制Prokka使用该路径数据库
prokka --setupdb  # 从官方源拉取最新注释数据库

# --------------------------
# 第二步：在线下载/更新VFDB毒力因子数据库（Diamond用）
# --------------------------
echo -e "\n🔄 正在下载/更新VFDB毒力因子数据库（在线获取最新版）..."
mkdir -p "$VFDB_ONLINE_DIR" && cd "$VFDB_ONLINE_DIR"

# 检查VFDB是否需要更新（通过文件大小校验）
VFDB_URL="https://www.mgc.ac.cn/VFs/Down/VFDB_setA_pro.fas.gz"
LOCAL_FILE="VFDB_setA_pro.fas.gz"
REMOTE_SIZE=$(curl -sI "$VFDB_URL" | grep -i Content-Length | awk '{print $2}' | tr -d '\r' || echo 0)
LOCAL_SIZE=$(if [ -f "$LOCAL_FILE" ]; then wc -c < "$LOCAL_FILE"; else echo 0; fi)

if [ "$LOCAL_SIZE" -ne "$REMOTE_SIZE" ] || [ ! -f "$VFDB_DB.dmnd" ]; then
  echo "📥 正在下载最新VFDB数据库..."
  rm -f VFDB_setA_pro.fas*  # 删除旧文件
  curl -O "$VFDB_URL" || wget "$VFDB_URL"  # 双下载方式容错
  gunzip -f VFDB_setA_pro.fas.gz  # 解压
  diamond makedb --in VFDB_setA_pro.fas --db VFDB_setA --quiet  # 建立Diamond索引
fi
cd /mnt/d/WSL/disk/projects/VP1  # 强制返回VP1目录

# --------------------------
# 第三步：严格检查必要文件和目录
# --------------------------
echo -e "\n✅ 开始校验必要文件..."
if [ ! -d "$CLEAN_GENOME_DIR" ]; then
  echo "❌ 错误：基因组目录不存在！路径：$CLEAN_GENOME_DIR"
  exit 1
fi

# 筛选有效.fna文件（修复空文件、非FASTA文件问题）
valid_genomes=()
for file in "$CLEAN_GENOME_DIR"/*.fna; do
  # 条件：1.是文件 2.非空 3.是FASTA格式（首行含>）
  if [ -f "$file" ] && [ -s "$file" ] && grep -q "^>" "$file" 2>/dev/null; then
    valid_genomes+=("$file")
  fi
done

# 检查是否有有效基因组文件
if [ ${#valid_genomes[@]} -eq 0 ]; then
  echo "❌ 错误：$CLEAN_GENOME_DIR 中未找到有效.fna文件！"
  echo "请检查：1.文件后缀是否为.fna 2.文件是否非空 3.是否为标准FASTA格式（首行以>开头）"
  exit 1
else
  echo "✅ 找到 ${#valid_genomes[@]} 个有效基因组文件，开始分析..."
fi

# --------------------------
# 第四步：创建输出目录（清空旧结果）
# --------------------------
rm -rf prokka_annotations vf_blast_results final_results
mkdir -p prokka_annotations vf_blast_results final_results

# --------------------------
# 第五步：批量处理每个菌株（Prokka注释 + Diamond比对）
# --------------------------
for genome_file in "${valid_genomes[@]}"; do
  # 提取菌株名（移除路径和后缀）
  strain_name=$(basename "$genome_file" .fna | sed 's/[^a-zA-Z0-9_]/_/g')
  echo -e "\n=================================================="
  echo "📌 正在处理菌株：$strain_name"
  echo "📂 基因组文件：$genome_file"
  echo "=================================================="

  # 5.1 Prokka在线注释（使用最新数据库）
  echo "🔧 运行Prokka注释..."
  prokka --outdir "prokka_annotations/$strain_name" \
         --prefix "$strain_name" \
         --kingdom Bacteria \
         --genus Vibrio \
         --species parahaemolyticus \
         --cpus $THREADS \
         --force \
         --usegenus \
         "$genome_file" 2>&1 | tee "prokka_annotations/$strain_name/annotation.log"

  # 5.2 检查注释生成的蛋白文件
  protein_file="prokka_annotations/$strain_name/$strain_name.faa"
  if [ ! -f "$protein_file" ] || [ ! -s "$protein_file" ]; then
    echo "⚠️  警告：未生成有效蛋白文件，跳过该菌株"
    continue
  fi

  # 5.3 Diamond毒力因子比对（使用在线更新的VFDB）
  echo "🔍 运行毒力因子比对..."
  diamond blastp \
          --db "$VFDB_DB" \
          --query "$protein_file" \
          --out "vf_blast_results/${strain_name}_vf_results.tsv" \
          --outfmt 6 qseqid sseqid pident length evalue qcovhsp stitle \
          --evalue $EVALUE \
          --id $MIN_ID \
          --query-cover $MIN_COVER \
          --threads $THREADS 2>&1 | tee "vf_blast_results/${strain_name}_blast.log"

  # 5.4 整理结果
  result_file="final_results/${strain_name}_virulence_factors.tsv"
  if [ -s "vf_blast_results/${strain_name}_vf_results.tsv" ]; then
    # 添加菌株名列，整理格式
    awk -v strain="$strain_name" '{
      printf("%s\t%s\t%s\t%.2f\t%d\t%.1e\t%.2f\t%s\n", 
             strain, $1, $2, $3, $4, $5, $6, $7)
    }' "vf_blast_results/${strain_name}_vf_results.tsv" > "$result_file"
    vf_count=$(wc -l < "$result_file")
    echo "✅ 完成！检测到 $vf_count 个毒力因子"
  else
    echo -e "$strain_name\t无\t无\t无\t无\t无\t无\t未检测到已知毒力因子" > "$result_file"
    echo "✅ 完成！未检测到已知毒力因子"
  fi
done

# --------------------------
# 第六步：生成汇总表（Excel可直接打开）
# --------------------------
summary_file="final_results/all_strains_virulence_summary.tsv"
echo -e "菌株名\t菌株蛋白ID\tVFDB毒力因子ID\t相似度(%)\t序列长度\t置信度(evalue)\t覆盖率(%)\t毒力因子功能描述" > "$summary_file"
cat final_results/*_virulence_factors.tsv >> "$summary_file"

# --------------------------
# 分析完成提示
# --------------------------
echo -e "\n🎉 所有分析流程完成！"
echo "📋 结果文件位置："
echo "  - Prokka注释结果：prokka_annotations/"
echo "  - 比对原始结果：vf_blast_results/"
echo "  - 毒力因子汇总表：$summary_file"
echo "💡 汇总表可直接用Excel打开，包含所有菌株的毒力因子信息"
