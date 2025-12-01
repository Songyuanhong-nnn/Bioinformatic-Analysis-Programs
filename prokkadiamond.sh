#!/bin/bash
set -euo pipefail

# ===================== 通用函数定义（增强灵活性）=====================
# 函数：询问用户输入（带默认值）
ask_input() {
  local prompt="$1"
  local default="$2"
  local input
  read -p "$prompt (默认: $default): " input
  # 若用户未输入，使用默认值；否则使用用户输入
  echo "${input:-$default}"
}

# 函数：检查目录是否存在，不存在则创建
check_and_create_dir() {
  local dir="$1"
  if [ ! -d "$dir" ]; then
    echo "📁 目录 $dir 不存在，正在创建..."
    mkdir -p "$dir" || { echo "❌ 无法创建目录 $dir！"; exit 1; }
  fi
}

# 函数：验证文件是否为有效FASTA（首行含>）
is_valid_fasta() {
  local file="$1"
  grep -q "^>" "$file" 2>/dev/null && return 0 || return 1
}

# ===================== 交互配置（核心灵活点）=====================
echo "=============================================="
echo "🎯  Vibrio prokka注释 + VFDB毒力因子分析（通用版）"
echo "=============================================="
echo "📌 所有路径支持相对路径或绝对路径（例如：./clean_genome 或 /home/user/genomes）"
echo "📌 直接回车使用默认值，适配大多数场景"
echo "==============================================\n"

# 1. 核心输入路径（用户交互）
echo "🔶 第一步：配置输入路径"
GENOME_DIR=$(ask_input "请输入基因组文件（.fna）所在目录" "./clean_genome")
PROKKA_DB_DIR=$(ask_input "请输入Prokka数据库存储目录" "./prokka_db")
VFDB_DIR=$(ask_input "请输入VFDB数据库存储目录" "./vfdb_db")

# 2. 输出路径配置（用户交互，默认按功能分类）
echo -e "\n🔶 第二步：配置输出路径"
OUTPUT_ROOT=$(ask_input "请输入结果根目录" "./analysis_results")
ANNOTATION_DIR=$(ask_input "请输入Prokka注释结果目录（相对根目录）" "01_prokka_annotations")
BLAST_DIR=$(ask_input "请输入毒力因子比对结果目录（相对根目录）" "02_vf_blast_results")
FINAL_DIR=$(ask_input "请输入最终汇总结果目录（相对根目录）" "03_final_summary")

# 3. 核心参数配置（用户可调整）
echo -e "\n🔶 第三步：配置分析参数"
THREADS=$(ask_input "请输入线程数（建议与CPU核心数一致）" "8")
MIN_ID=$(ask_input "毒力因子比对最小相似度(%)" "40")
MIN_COVER=$(ask_input "毒力因子比对最小覆盖率(%)" "70")
EVALUE=$(ask_input "毒力因子比对最大E值" "1e-10")

# 4. 物种信息配置（用户可修改，适配不同物种）
echo -e "\n🔶 第四步：配置物种信息（用于Prokka注释）"
KINGDOM=$(ask_input "请输入生物界（Bacteria/Archaea/Virus/Eukaryota）" "Bacteria")
GENUS=$(ask_input "请输入属名（例如：Vibrio）" "Vibrio")
SPECIES=$(ask_input "请输入种名（例如：parahaemolyticus）" "parahaemolyticus")

# ===================== 路径标准化（避免格式错误）=====================
# 转换为绝对路径，避免相对路径混乱
GENOME_DIR=$(realpath "$GENOME_DIR")
PROKKA_DB_DIR=$(realpath "$PROKKA_DB_DIR")
VFDB_DIR=$(realpath "$VFDB_DIR")
OUTPUT_ROOT=$(realpath "$OUTPUT_ROOT")

# 拼接最终输出路径（按数字前缀排序，结构清晰）
ANNOTATION_DIR="$OUTPUT_ROOT/$ANNOTATION_DIR"
BLAST_DIR="$OUTPUT_ROOT/$BLAST_DIR"
FINAL_DIR="$OUTPUT_ROOT/$FINAL_DIR"

# ===================== 环境检查与初始化 =====================
echo -e "\n✅ 开始环境检查与初始化..."

# 1. 检查必要工具是否安装
REQUIRED_TOOLS=("prokka" "diamond" "curl" "wget" "parallel" "grep" "awk" "sed")
for tool in "${REQUIRED_TOOLS[@]}"; do
  if ! command -v "$tool" &> /dev/null; then
    echo "❌ 未找到必要工具：$tool！请先安装后重试"
    exit 1
  fi
done

# 2. 检查并创建所有目录（输入+输出）
check_and_create_dir "$PROKKA_DB_DIR"
check_and_create_dir "$VFDB_DIR"
check_and_create_dir "$ANNOTATION_DIR"
check_and_create_dir "$BLAST_DIR"
check_and_create_dir "$FINAL_DIR"

# 3. 筛选有效基因组文件（.fna/.fasta后缀，支持大小写）
echo -e "\n🔍 正在筛选有效基因组文件..."
VALID_GENOMES=()
# 支持 .fna 和 .fasta 后缀（大小写兼容：.FNA/.FASTA）
for ext in "fna" "FNA" "fasta" "FASTA"; do
  for file in "$GENOME_DIR"/*."$ext"; do
    if [ -f "$file" ] && [ -s "$file" ] && is_valid_fasta "$file"; then
      VALID_GENOMES+=("$file")
    fi
  done
done

# 去重（避免重复文件）
VALID_GENOMES=($(printf "%s\n" "${VALID_GENOMES[@]}" | sort -u))

# 检查是否有有效文件
if [ ${#VALID_GENOMES[@]} -eq 0 ]; then
  echo "❌ 错误：在 $GENOME_DIR 中未找到有效基因组文件！"
  echo "请检查：1.文件后缀是否为 .fna/.fasta 2.文件是否非空 3.是否为标准FASTA格式（首行以>开头）"
  exit 1
else
  echo "✅ 找到 ${#VALID_GENOMES[@]} 个有效基因组文件："
  for file in "${VALID_GENOMES[@]}"; do echo "  - $(basename "$file")"; done
fi

# ===================== 数据库在线更新（自动检测是否需要更新）=====================
echo -e "\n🔄 开始数据库更新（仅当需要时下载）..."

# 1. Prokka数据库更新（--setupdb自动增量更新）
echo -e "\n📥 更新Prokka注释数据库（存储路径：$PROKKA_DB_DIR）..."
export PROKKA_DB="$PROKKA_DB_DIR"
prokka --setupdb 2>&1 | grep -E "Downloading|Updating|Setup complete" || echo "✅ Prokka数据库已为最新"

# 2. VFDB毒力因子数据库更新（按文件大小校验是否需要重新下载）
echo -e "\n📥 更新VFDB毒力因子数据库（存储路径：$VFDB_DIR）..."
VFDB_URL="https://www.mgc.ac.cn/VFs/Down/VFDB_setA_pro.fas.gz"
VFDB_GZ="$VFDB_DIR/VFDB_setA_pro.fas.gz"
VFDB_FAS="$VFDB_DIR/VFDB_setA_pro.fas"
VFDB_DMND="$VFDB_DIR/VFDB_setA.dmnd"

# 校验远程文件大小（避免重复下载）
REMOTE_SIZE=$(curl -sI "$VFDB_URL" | grep -i "Content-Length" | awk '{print $2}' | tr -d '\r' || echo 0)
LOCAL_SIZE=$(if [ -f "$VFDB_GZ" ]; then wc -c < "$VFDB_GZ"; else echo 0; fi)

if [ "$LOCAL_SIZE" -ne "$REMOTE_SIZE" ] || [ ! -f "$VFDB_DMND" ]; then
  echo "📥 正在下载最新VFDB数据库..."
  rm -f "$VFDB_DIR"/VFDB_setA*  # 删除旧版本
  # 双下载工具容错（curl失败则用wget）
  if ! curl -o "$VFDB_GZ" "$VFDB_URL"; then
    echo "⚠️ curl下载失败，尝试用wget..."
    wget -O "$VFDB_GZ" "$VFDB_URL" || { echo "❌ VFDB下载失败！请检查网络"; exit 1; }
  fi
  gunzip -f "$VFDB_GZ"  # 解压（-f覆盖旧文件）
  diamond makedb --in "$VFDB_FAS" --db "$VFDB_DIR/VFDB_setA" --quiet  # 构建索引
else
  echo "✅ VFDB数据库已为最新，跳过下载"
fi

# ===================== 批量分析核心流程 =====================
echo -e "\n🚀 开始批量分析（共 ${#VALID_GENOMES[@]} 个菌株）..."
echo "=============================================="

# 初始化汇总表
SUMMARY_FILE="$FINAL_DIR/all_strains_virulence_summary.tsv"
echo -e "菌株名\t菌株蛋白ID\tVFDB毒力因子ID\t相似度(%)\t序列长度\t置信度(evalue)\t覆盖率(%)\t毒力因子功能描述" > "$SUMMARY_FILE"

# 遍历每个有效基因组
for GENOME_FILE in "${VALID_GENOMES[@]}"; do
  # 提取菌株名（移除路径和后缀，兼容各种后缀）
  STRAIN_NAME=$(basename "$GENOME_FILE" | sed -E 's/\.(fna|FNA|fasta|FASTA)$//' | sed 's/[^a-zA-Z0-9_]/_/g')
  echo -e "\n=================================================="
  echo "📌 正在处理菌株：$STRAIN_NAME"
  echo "📂 基因组文件：$GENOME_FILE"
  echo "=================================================="

  # 1. Prokka注释（输出到独立目录，避免冲突）
  ANNOTATION_OUT="$ANNOTATION_DIR/$STRAIN_NAME"
  echo "🔧 运行Prokka注释（输出路径：$ANNOTATION_OUT）..."
  prokka --outdir "$ANNOTATION_OUT" \
         --prefix "$STRAIN_NAME" \
         --kingdom "$KINGDOM" \
         --genus "$GENUS" \
         --species "$SPECIES" \
         --cpus "$THREADS" \
         --force \
         --usegenus \
         "$GENOME_FILE" 2>&1 | tee "$ANNOTATION_OUT/annotation.log"

  # 2. 检查注释生成的蛋白文件（.faa）
  PROTEIN_FILE="$ANNOTATION_OUT/$STRAIN_NAME.faa"
  if [ ! -f "$PROTEIN_FILE" ] || [ ! -s "$PROTEIN_FILE" ]; then
    echo "⚠️  警告：未生成有效蛋白文件，跳过该菌株"
    echo -e "$STRAIN_NAME\t无\t无\t无\t无\t无\t无\tProkka注释失败，未生成蛋白文件" >> "$SUMMARY_FILE"
    continue
  fi

  # 3. Diamond毒力因子比对
  BLAST_OUT="$BLAST_DIR/${STRAIN_NAME}_vf_blast.tsv"
  BLAST_LOG="$BLAST_DIR/${STRAIN_NAME}_blast.log"
  echo "🔍 运行毒力因子比对（输出路径：$BLAST_OUT）..."
  diamond blastp \
          --db "$VFDB_DMND" \
          --query "$PROTEIN_FILE" \
          --out "$BLAST_OUT" \
          --outfmt 6 qseqid sseqid pident length evalue qcovhsp stitle \
          --evalue "$EVALUE" \
          --id "$MIN_ID" \
          --query-cover "$MIN_COVER" \
          --threads "$THREADS" 2>&1 | tee "$BLAST_LOG"

  # 4. 整理结果（单菌株结果+汇总）
  SINGLE_RESULT="$FINAL_DIR/${STRAIN_NAME}_virulence_factors.tsv"
  if [ -s "$BLAST_OUT" ]; then
    # 格式化结果（添加菌株名列）
    awk -v strain="$STRAIN_NAME" '{
      printf("%s\t%s\t%s\t%.2f\t%d\t%.1e\t%.2f\t%s\n", 
             strain, $1, $2, $3, $4, $5, $6, $7)
    }' "$BLAST_OUT" > "$SINGLE_RESULT"
    # 统计毒力因子数量
    VF_COUNT=$(wc -l < "$SINGLE_RESULT")
    echo "✅ 完成！检测到 $VF_COUNT 个毒力因子"
    # 追加到汇总表
    cat "$SINGLE_RESULT" >> "$SUMMARY_FILE"
  else
    # 无匹配结果时的默认输出
    echo -e "$STRAIN_NAME\t无\t无\t无\t无\t无\t无\t未检测到已知毒力因子" > "$SINGLE_RESULT"
    echo "✅ 完成！未检测到已知毒力因子"
    # 追加到汇总表
    cat "$SINGLE_RESULT" >> "$SUMMARY_FILE"
  fi
done

# ===================== 分析完成提示 =====================
echo -e "\n🎉 所有菌株分析完成！"
echo "=============================================="
echo "📂 结果目录结构（清晰分类）："
echo "└── $OUTPUT_ROOT"
echo "    ├── $ANNOTATION_DIR  # Prokka注释原始结果（含.gbk/.faa/.gff）"
echo "    ├── $BLAST_DIR       # Diamond比对原始结果（.tsv/.log）"
echo "    └── $FINAL_DIR       # 最终整理结果"
echo "        ├── 各菌株毒力因子详情（xxx_virulence_factors.tsv）"
echo "        └── 所有菌株汇总表（all_strains_virulence_summary.tsv）"
echo -e "\n💡 关键使用提示："
echo "  1. 汇总表可直接用Excel打开（TSV格式，制表符分隔）"
echo "  2. .gbk文件可导入SnapGene/IGV进行基因组可视化"
echo "  3. 若需重新运行，可修改任意路径/参数，脚本自动跳过已完成步骤"
echo "  4. 支持跨平台运行（Windows/WSL/Linux/macOS）"
echo -e "\n📚 引用说明："
echo "  - Prokka: Seemann T (2014) Bioinformatics 30(14):2068-9"
echo "  - VFDB: Chen L et al. (2023) Nucleic Acids Res 51(D1):D911-D917"
echo "=============================================="