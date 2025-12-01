#!/bin/bash
set -euo pipefail

# ==================== 已确认的VFDB路径（不用改）====================
VFDB="/mnt/g/wsl/database/vfdb/VFDB_setA_nt.fas"
BAKTA_DIR="/mnt/d/WSL/disk/projects/VP1/bakta_annotations"
GENOME_DIR="/mnt/d/WSL/disk/projects/VP1/clean_genome"
OUTPUT="virulence_gene_summary.csv"  # 最终汇总表

# 构建VFDB索引（仅需一次，已存在则自动跳过）
echo "🔧 构建VFDB索引..."
makeblastdb -in "$VFDB" -dbtype nucl -out "$VFDB" -quiet 2>/dev/null || echo "✅ 索引已存在"

# 创建结果表头（字段清晰，方便后续筛选）
echo "菌株名,基因ID,基因名,毒力功能,注释来源,BLAST相似度,E值,血清型,危害等级" > "$OUTPUT"

# 遍历所有Bakta注释菌株
for STRAIN_DIR in "$BAKTA_DIR"/*/; do
    STRAIN_NAME=$(basename "$STRAIN_DIR")
    TSV="$STRAIN_DIR/$STRAIN_NAME.tsv"  # Bakta高精度注释表
    FNA="$GENOME_DIR/$STRAIN_NAME.fna"  # 基因组序列文件

    echo -e "\n📌 处理菌株：$STRAIN_NAME"
    # 跳过缺少核心文件的菌株（不报错，仅提示）
    [ ! -f "$TSV" ] && echo "⚠️  缺少注释文件，跳过" && continue
    [ ! -f "$FNA" ] && echo "⚠️  缺少基因组文件，跳过" && continue

    # 第一步：筛选毒力相关候选基因（关键词精准匹配）
    echo "🔍 筛选毒力候选基因..."
    grep -E "virulence|toxin|tdh|trh|hemolysin|adhesin|T3SS|T6SS|invasion" "$TSV" | \
    awk -F "\t" '{
        gene_id = $1;
        gene_name = ($2 == "-") ? "unknown_vir_gene" : $2;  # 无名基因命名
        func = $5;
        print gene_id "\t" gene_name "\t" func
    }' > tmp_cand.tsv
