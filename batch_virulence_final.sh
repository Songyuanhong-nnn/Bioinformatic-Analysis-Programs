#!/bin/bash
set -euo pipefail

# ==================== 固定路径（无需修改，按你的项目目录）====================
ENV_NAME="bakta"  # 已确认的活跃环境
BAKTA_ANNOT_DIR="/mnt/d/WSL/disk/projects/VP1/bakta_annotations"
CLEAN_GENOME_DIR="/mnt/d/WSL/disk/projects/VP1/clean_genome"
OUTPUT_CSV="virulence_extraction_results_final.csv"
# VFDB数据库路径（若你的路径不同，仅修改这一行！）
VFDB_DB="/mnt/d/WSL/disk/projects/VP1/database/VFDB_setA_nt.fas"
# ==============================================================================

# 自动激活bakta环境（防止脚本运行中环境失效）
source activate "$ENV_NAME"

# 核心检查：VFDB数据库是否存在
if [ ! -f "$VFDB_DB" ]; then
    echo "❌ 错误：VFDB数据库文件不存在！"
    echo "   请按以下步骤修复："
    echo "   1. 访问 http://www.mgc.ac.cn/VFs/download.htm"
    echo "   2. 下载 'VFDB_setA_nt.fas'（细菌毒力基因核心库）"
    echo "   3. 解压后放在路径：/mnt/d/WSL/disk/projects/VP1/database/"
    echo "   4. 重新运行本脚本"
    exit 1
fi

# 给VFDB建BLAST索引（仅第一次运行时执行，后续自动跳过）
if [ ! -f "${VFDB_DB}.nhr" ]; then
    echo "🔧 给VFDB数据库建BLAST索引（仅需一次）..."
    makeblastdb -in "$VFDB_DB" -dbtype nucl -out "$VFDB_DB" -quiet
fi

# 创建输出文件（含关键字段，方便后续筛选）
echo "菌株名,基因ID,基因名,毒力功能描述,注释来源,BLAST相似度(%),E值,血清型关联" > "$OUTPUT_CSV"

# 统计待处理菌株数
strain_count=$(ls -l "$BAKTA_ANNOT_DIR" | grep "^d" | wc -l)
echo "===== 开始提取毒力基因（环境：$ENV_NAME，菌株数：$strain_count）====="

# 遍历所有Bakta注释菌株（仅处理完整注释的菌株）
for STRAIN_DIR in "$BAKTA_ANNOT_DIR"/*/; do
    STRAIN_NAME=$(basename "$STRAIN_DIR")
    BAKTA_TSV="$STRAIN_DIR/${STRAIN_NAME}.tsv"  # Bakta高精度注释表
    BAKTA_FNA="$CLEAN_GENOME_DIR/${STRAIN_NAME}.fna"  # 基因组序列

    echo -e "\n=================================================="
    echo "📌 正在处理：$STRAIN_NAME"
    echo "📂 注释文件：$BAKTA_TSV"
    echo "=================================================="

    # 跳过缺少核心文件的菌株（避免报错）
    if [ ! -f "$BAKTA_TSV" ] || [ ! -f "$BAKTA_FNA" ]; then
        echo "⚠️  跳过：缺少Bakta注释文件或基因组文件"
        continue
    fi

    # 1. 精准筛选毒力相关基因（关键词覆盖副溶血性弧菌关键毒力基因）
    echo "🔍 筛选毒力候选基因（tdh/trh/毒素/黏附素等）..."
    grep -E "virulence|toxin|hemolysin|tdh|trh|pil[ABCD]|fimbria|adhesin|invasion|pathogen|T3SS|T6SS|secreted" "$BAKTA_TSV" | \
    awk -F "\t" '{
        gene_id = $1;
        gene_name = ($2 == "-") ? "unknown_vir_gene" : $2;  # 无名基因命名为unknown
        func = $5;
        print gene_id "\t" gene_name "\t" func
    }' > "${STRAIN_NAME}_candidate.tsv"

    # 2. 提取候选基因的核苷酸序列（用于BLAST验证）
    echo "📑 提取候选基因序列..."
    awk -F "\t" 'NR==FNR {gene_ids[$1]=1; next} /^>/ {
        flag=0;
        split(substr($0,2), id_arr, " ");  # 分割序列ID（去掉">"后按空格拆分）
        if (id_arr[1] in gene_ids) flag=1;
        print;
        next
    } flag' "${STRAIN_NAME}_candidate.tsv" "$BAKTA_FNA" > "${STRAIN_NAME}_candidate.fna"

    # 3. BLAST比对VFDB验证（严格阈值：相似度≥70%，E值<1e-10）
    echo "🧬 BLAST比对VFDB数据库（验证毒力基因）..."
    blastn -query "${STRAIN_NAME}_candidate.fna" \
           -db "$VFDB_DB" \
           -outfmt "6 qseqid sseqid pident evalue" \
           -evalue 1e-10 \
           -perc_identity 70 \
           -num_threads 8 \
           -quiet > "${STRAIN_NAME}_blast.tsv"

    # 4. 整合结果（去重，写入汇总表）
    echo "📊 整合并写入结果..."
    join -t $'\t' -1 1 -2 1 \
        <(sort -k1,1 "${STRAIN_NAME}_candidate.tsv") \
        <(sort -k1,1 "${STRAIN_NAME}_blast.tsv") \
    | awk -v strain="$STRAIN_NAME" -F "\t" '{
        # 血清型关联：从菌株名提取（如GCA_XXX_O3K6_XXX → 标注O3K6）
        sero = "";
        if (strain ~ /O3.*K6/ || strain ~ /O4.*K8/ || strain ~ /O1.*K25/) {
            sero = substr(strain, match(strain, /O[0-9]+.*K[0-9]+/), RLENGTH);
        } else {
            sero = "未知血清型";
        }
        # 格式：菌株名,基因ID,基因名,功能,来源,相似度,E值,血清型
        print strain "," $1 "," $2 "," $3 ",Bakta+VFDB," $5 "," $6 "," sero
    }' >> "$OUTPUT_CSV"

    # 清理临时文件（避免占用空间）
    rm -f "${STRAIN_NAME}_candidate.tsv" "${STRAIN_NAME}_candidate.fna" "${STRAIN_NAME}_blast.tsv"
    echo "✅ 完成：$STRAIN_NAME"
done

# 最终统计
result_count=$(wc -l "$OUTPUT_CSV" | awk '{print $1 - 1}')  # 减去表头
echo -e "\n=================================================="
echo "🎉 所有菌株处理完成！"
echo "📄 最终毒力基因汇总表：$OUTPUT_CSV"
echo "📈 共筛选到有效毒力基因：$result_count 个"
echo "📋 汇总表字段：菌株名,基因ID,基因名,毒力功能描述,注释来源,BLAST相似度(%),E值,血清型关联"
echo "=================================================="
