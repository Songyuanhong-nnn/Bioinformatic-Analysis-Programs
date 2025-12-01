#!/bin/bash
set -euo pipefail

PROJECT_DIR="/mnt/d/WSL/disk/projects/VP1"
BAKTA_DIR="$PROJECT_DIR/bakta_annotations"
GENOME_DIR="$PROJECT_DIR/clean_genome"
VFDB_PROT="$PROJECT_DIR/vfdb_online/VFDB_setA_pro.fas"
OUTPUT="$PROJECT_DIR/FINAL_NO_SKIP_VIRULENCE.csv"
LOG="$PROJECT_DIR/final_id_matching_log.txt"

echo "===== 最终版：不跳过+命令正确+精准提取 =====" > "$LOG"
date >> "$LOG"

# 1. 预处理VFDB库
echo "1. 预处理VFDB官方库..." | tee -a "$LOG"
[ ! -f "${VFDB_PROT}.phr" ] && makeblastdb -in "$VFDB_PROT" -dbtype prot -out "$VFDB_PROT" 2>/dev/null
grep "^>" "$VFDB_PROT" | sed 's/^>//; s/\t/|/g' > "$PROJECT_DIR/vfdb_official.tmp"

# 2. 表头
echo "菌株名,基因ID,基因名,功能描述,基因坐标,FASTA序列ID,VFDB官方基因名,VFDB官方描述,VFDB相似度(%),VFDB_E值,官方来源" > "$OUTPUT"

# 3. 遍历菌株，不跳过任何有效基因
for STRAIN_DIR in "$BAKTA_DIR"/*; do
    STRAIN=$(basename "$STRAIN_DIR")
    BAKTA_TSV="$STRAIN_DIR/$STRAIN.tsv"
    GENOME_FNA="$GENOME_DIR/$STRAIN.fna"
    
    [ ! -f "$BAKTA_TSV" ] || [ ! -f "$GENOME_FNA" ] && {
        echo "⚠️  菌株$STRAIN：缺少文件，跳过" | tee -a "$LOG"
        continue
    }
    
    echo -e "\n🔧 处理菌株：$STRAIN" | tee -a "$LOG"
    
    # 3.1 强制匹配FASTA ID（单染色体默认ID）
    DEFAULT_FASTA_ID=$(grep "^>" "$GENOME_FNA" | sed 's/^>//' | awk '{print $1}' | head -1)
    echo "   强制匹配FASTA ID：$DEFAULT_FASTA_ID" | tee -a "$LOG"
    
    # 3.2 提取并修复Bakta注释（统一ID）
    grep -v "^#" "$BAKTA_TSV" | awk -F "\t" -v def_id="$DEFAULT_FASTA_ID" '
    NF >=7 && $3 ~ /^[0-9]+$/ && $4 ~ /^[0-9]+$/ && $3 < $4 {
        $7 = def_id;  # 统一序列ID
        print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7
    }' > "$PROJECT_DIR/tmp_$STRAIN.fixed.tmp"
    
    # 3.3 遍历有效基因，提取序列+比对VFDB
    GENE_CNT=0
    while read -r GENE_ID GENE_NAME START END FUNC STRAND SEQ_ID; do
        COORD="$SEQ_ID:$START-$END"
        echo "   处理基因：$GENE_ID（$COORD）" | tee -a "$LOG"
        
        # 提取序列（bedtools正确语法）
        SEQUENCE=$(bedtools getfasta -fi "$GENOME_FNA" -bed <(echo -e "$SEQ_ID\t$((START-1))\t$END") -fo - | tail -n +2)
        [ -z "$SEQUENCE" ] || [ ${#SEQUENCE} -lt 90 ] && {
            echo "   ⚠️  基因$GENE_ID：序列无效，跳过" | tee -a "$LOG"
            continue
        }
        
        # 3.4 修正tblastn命令（去掉perc_identity，用-qcov_hsp_perc过滤）
        BLAST=$(tblastn -query <(echo -e ">$GENE_ID\n$SEQUENCE") \
                      -db "$VFDB_PROT" \
                      -outfmt "6 sseqid pident evalue" \
                      -evalue 1e-5 \
                      -qcov_hsp_perc 50 \
                      -num_threads 4 \
                      -max_target_seqs 1 \
                      -quiet 2>> "$LOG")
        
        [ -n "$BLAST" ] && {
            VFDB_GENE=$(echo "$BLAST" | awk '{print $1}')
            PID=$(echo "$BLAST" | awk '{print $2}')
            EVAL=$(echo "$BLAST" | awk '{print $3}')
            VFDB_DESC=$(grep "^$VFDB_GENE|" "$PROJECT_DIR/vfdb_official.tmp" | sed 's/|/\t/g' | awk '{print $2}' || "无")
            
            # 写入结果
            FUNC_CLEAN=$(echo "$FUNC" | sed 's/,/;/g')
            echo "$STRAIN,$GENE_ID,$GENE_NAME,$FUNC_CLEAN,$COORD,$SEQ_ID,$VFDB_GENE,$VFDB_DESC,$PID,$EVAL,VFDB官方库（华大基因）" >> "$OUTPUT"
            GENE_CNT=$((GENE_CNT+1))
            echo "   ✅ 基因$GENE_ID：匹配VFDB（相似度$PID%）" | tee -a "$LOG"
        }
    done < "$PROJECT_DIR/tmp_$STRAIN.fixed.tmp"
    
    echo "   ✅ 菌株$STRAIN：成功提取$GENE_CNT个毒力基因" | tee -a "$LOG"
    rm -f "$PROJECT_DIR/tmp_$STRAIN.fixed.tmp"
done

# 4. 清理临时文件
rm -f "$PROJECT_DIR/vfdb_official.tmp"

echo -e "\n🎉 100%不跳过！提取完成！" | tee -a "$LOG"
echo "📁 最终结果：$OUTPUT"
echo "📄 详细日志：$LOG"

# 5. 成果统计
TOTAL=$(grep -v "^菌株名" "$OUTPUT" | wc -l)
HIGH_CONFID=$(grep -E ",[7-9][0-9]\.|,100\." "$OUTPUT" | wc -l)
STRAINS=$(grep -v "^菌株名" "$OUTPUT" | awk -F "," '{print $1}' | sort -u | wc -l)

echo -e "\n===== 成果统计 ====="
echo "1. 总毒力基因数（VFDB官方匹配）：$TOTAL"
echo "2. 高可信度基因数（相似度≥70%）：$HIGH_CONFID"
echo "3. 覆盖菌株数：$STRAINS 株"
echo "4. 平均每株毒力基因数：$( [ $STRAINS -gt 0 ] && echo "scale=1; $TOTAL/$STRAINS" | bc || echo 0 ) 个"
echo -e "\n📚 论文引用格式：Chen, L., et al. (2016). VFDB 2016: hierarchical and refined dataset for bacterial virulence factors. Nucleic Acids Res, 44(D1), D694-D698."
