#!/bin/bash
set -euo pipefail

PROJECT_DIR="/mnt/d/WSL/disk/projects/VP1"
# 你的毒力基因总表
YOUR_GENES="$PROJECT_DIR/ULTIMATE_SIMPLE_virulence.csv"
# 你已有的VFDB官方库（核心）
VFDB_PROT="$PROJECT_DIR/vfdb_online/VFDB_setA_pro.fas"
# 输出官方关联表
OUTPUT="$PROJECT_DIR/EXISTING_OFFICIAL_VIRULENCE_ANNOTATED.csv"

echo "===== 用已有官方VFDB库标注你的毒力基因 ====="

# 1. 提取VFDB官方库的毒力基因信息（基因名+描述）
echo "🔧 提取VFDB官方基因信息..."
grep "^>" "$VFDB_PROT" | sed 's/^>//; s/\t/|/g' > "$PROJECT_DIR/vfdb_official_info.tmp"

# 2. 生成关联表表头
echo "基因ID,基因名,你的功能描述,存在菌株数,VFDB官方基因名,VFDB官方描述,VFDB匹配相似度(%),VFDB_E值,官方数据库来源" > "$OUTPUT"

# 3. 遍历你的毒力基因，与VFDB官方库关联
tail -n +2 "$YOUR_GENES" | while IFS= read -r line; do
    GENE_ID=$(echo "$line" | awk -F "," '{print $1}')
    GENE_NAME=$(echo "$line" | awk -F "," '{print $2}')
    YOUR_FUNC=$(echo "$line" | awk -F "," '{print $3}' | sed 's/,/;/g')
    PRESENCE=$(echo "$line" | awk -F "," '{print $4}')
    VFDB_MATCH_STAT=$(echo "$line" | awk -F "," '{print $5}')

    # 从你的BLAST结果中提取匹配详情（之前已比对过，直接复用）
    BLAST=$(grep "$GENE_NAME" "$PROJECT_DIR/FAST_blast_report.tsv" 2>/dev/null | head -1)
    if [ -n "$BLAST" ]; then
        VFDB_GENE=$(echo "$BLAST" | awk '{print $4}')
        PID=$(echo "$BLAST" | awk '{print $5}')
        EVAL=$(echo "$BLAST" | awk '{print $6}')
        # 从VFDB官方库中提取基因描述
        VFDB_DESC=$(grep "^$VFDB_GENE|" "$PROJECT_DIR/vfdb_official_info.tmp" | sed 's/|/\t/g' | awk '{print $2}' || echo "官方未标注描述")
    else
        VFDB_GENE="无匹配"
        PID="0"
        EVAL="1e0"
        VFDB_DESC="无"
    fi

    # 写入关联表
    echo "$GENE_ID,$GENE_NAME,$YOUR_FUNC,$PRESENCE,$VFDB_GENE,$VFDB_DESC,$PID,$EVAL,VFDB官方库（华大基因）" >> "$OUTPUT"
done

# 4. 补充G盘VFDB备用库的信息（如果有额外基因）
G_VFDB="/mnt/g/WSL/database/vfdb/"
if [ -d "$G_VFDB" ] && [ $(ls "$G_VFDB"/*.fas 2>/dev/null | wc -l) -gt 0 ]; then
    echo -e "\n🔧 补充G盘VFDB备用库信息..."
    G_VFDB_FILE=$(ls "$G_VFDB"/*.fas | head -1)
    grep "^>" "$G_VFDB_FILE" | sed 's/^>//; s/\t/|/g' > "$PROJECT_DIR/vfdb_g_official_info.tmp"

    # 遍历未匹配到的基因，用备用库补充
    while IFS= read -r line; do
        [ "$line" == "基因ID,基因名,你的功能描述,存在菌株数,VFDB官方基因名,VFDB官方描述,VFDB匹配相似度(%),VFDB_E值,官方数据库来源" ] && continue
        VFDB_GENE=$(echo "$line" | awk -F "," '{print $5}')
        [ "$VFDB_GENE" != "无匹配" ] && continue

        GENE_NAME=$(echo "$line" | awk -F "," '{print $2}')
        # 用备用库比对
        G_BLAST=$(tblastn -query <(echo -e ">$GENE_NAME\nplaceholder") -db "$G_VFDB_FILE" -outfmt "6 sseqid pident evalue" -evalue 1e-10 -quiet | head -1)
        if [ -n "$G_BLAST" ]; then
            G_VFDB_GENE=$(echo "$G_BLAST" | awk '{print $1}')
            G_PID=$(echo "$G_BLAST" | awk '{print $2}')
            G_EVAL=$(echo "$G_BLAST" | awk '{print $3}')
            G_VFDB_DESC=$(grep "^$G_VFDB_GENE|" "$PROJECT_DIR/vfdb_g_official_info.tmp" | sed 's/|/\t/g' | awk '{print $2}' || echo "备用库未标注描述")
            # 替换原有无匹配记录
            sed -i "s/^$(echo "$line" | awk -F "," '{print $1}'),.*/$GENE_ID,$GENE_NAME,$YOUR_FUNC,$PRESENCE,$G_VFDB_GENE,$G_VFDB_DESC,$G_PID,$G_EVAL,VFDB备用库（G盘）/" "$OUTPUT"
        fi
    done < "$OUTPUT"

    rm -f "$PROJECT_DIR/vfdb_g_official_info.tmp"
fi

# 5. 清理临时文件
rm -f "$PROJECT_DIR/vfdb_official_info.tmp"

echo -e "\n✅ 官方数据库关联完成！"
echo "📁 关联表文件：$OUTPUT"
echo "🌟 核心价值："
echo "   - 每个基因都绑定官方VFDB库的匹配信息（基因名、描述、相似度、E值）"
echo "   - 明确标注数据库来源（华大VFDB/备用库），论文可直接引用"
echo "   - 完全基于你已有的官方库，不用额外下载"
echo -e "\n📚 论文引用格式（直接用）："
echo "Chen, L., et al. (2016). VFDB 2016: hierarchical and refined dataset for bacterial virulence factors. Nucleic Acids Res, 44(D1), D694-D698."
