#!/bin/bash
set -euo pipefail

PROJECT_DIR="/mnt/d/WSL/disk/projects/VP1"
BAKTA_DIR="$PROJECT_DIR/bakta_annotations"
GENOME_DIR="$PROJECT_DIR/clean_genome"
VFDB_PROT="$PROJECT_DIR/vfdb_online/VFDB_setA_pro.fas"
OUTPUT="$PROJECT_DIR/BATCH_FINAL_VIRULENCE.csv"
LOG="$PROJECT_DIR/batch_process_log.txt"

# 清理之前的残留文件
rm -f "$PROJECT_DIR/tmp_*.tmp" "$PROJECT_DIR/vfdb_official.tmp"

echo "===== 批量快速处理：确保跑完所有菌株 =====" > "$LOG"
date >> "$LOG"

# 1. 预处理VFDB库（只做1次）
echo "1. 预处理VFDB官方库..." | tee -a "$LOG"
[ ! -f "${VFDB_PROT}.phr" ] && makeblastdb -in "$VFDB_PROT" -dbtype prot -out "$VFDB_PROT" 2>/dev/null
grep "^>" "$VFDB_PROT" | sed 's/^>//; s/\t/|/g' > "$PROJECT_DIR/vfdb_official.tmp"

# 2. 只写1次表头（避免重复）
[ ! -f "$OUTPUT" ] && echo "菌株名,基因ID,基因名,功能描述,基因坐标,FASTA序列ID,VFDB官方基因名,VFDB官方描述,VFDB相似度(%),VFDB_E值,官方来源" > "$OUTPUT"

# 3. 批量处理所有菌株，跳过已处理的
for STRAIN_DIR in "$BAKTA_DIR"/*; do
    STRAIN=$(basename "$STRAIN_DIR")
    BAKTA_TSV="$STRAIN_DIR/$STRAIN.tsv"
    GENOME_FNA="$GENOME_DIR/$STRAIN.fna"
    
    # 跳过条件：缺少文件/已处理过该菌株
    if [ ! -f "$BAKTA_TSV" ] || [ ! -f "$GENOME_FNA" ]; then
        echo "⚠️  菌株$STRAIN：缺少文件，跳过" | tee -a "$LOG"
        continue
    fi
    grep -q "^$STRAIN," "$OUTPUT" && {
        echo "⚠️  菌株$STRAIN：已处理，跳过" | tee -a "$LOG"
        continue
    }
    
    echo -e "\n🔧 批量处理菌株：$STRAIN" | tee -a "$LOG"
    
    # 3.1 强制匹配FASTA ID
    DEFAULT_FASTA_ID=$(grep "^>" "$GENOME_FNA" | sed 's/^>//' | awk '{print $1}' | head -1)
    [ -z "$DEFAULT_FASTA_ID" ] && {
        echo "⚠️  菌株$STRAIN：无FASTA ID，跳过" | tee -a "$LOG"
        continue
    }
    echo "   匹配FASTA ID：$DEFAULT_FASTA_ID" | tee -a "$LOG"
    
    # 3.2 提取该菌株所有有效基因（批量生成bed文件）
    BED_FILE="$PROJECT_DIR/tmp_$STRAIN.bed"
    GENE_INFO_FILE="$PROJECT_DIR/tmp_$STRAIN.info"
    grep -v "^#" "$BAKTA_TSV" | awk -F "\t" -v def_id="$DEFAULT_FASTA_ID" '
    NF >=7 && $3 ~ /^[0-9]+$/ && $4 ~ /^[0-9]+$/ && $3 < $4 {
        # 输出bed格式（用于批量提取序列）
        print def_id "\t" $((3-1)) "\t" $4 "\t" $1
        # 输出基因信息（后续匹配用）
        print $1 "\t" $2 "\t" $5
    }' > "$BED_FILE"
    
    # 3.3 批量提取所有基因序列（避免单个提取卡住）
    FASTA_FILE="$PROJECT_DIR/tmp_$STRAIN.seq.fna"
    bedtools getfasta -fi "$GENOME_FNA" -bed "$BED_FILE" -fo "$FASTA_FILE" 2>> "$LOG"
    [ ! -s "$FASTA_FILE" ] && {
        echo "⚠️  菌株$STRAIN：无有效序列，跳过" | tee -a "$LOG"
        rm -f "$BED_FILE" "$FASTA_FILE"
        continue
    }
    
    # 3.4 批量比对VFDB（1次比对该菌株所有基因，更快更稳定）
    BLAST_FILE="$PROJECT_DIR/tmp_$STRAIN.blast"
    tblastn -query "$FASTA_FILE" \
            -db "$VFDB_PROT" \
            -outfmt "6 qseqid sseqid pident evalue" \
            -evalue 1e-3 \
            -qcov_hsp_perc 40 \
            -num_threads 8 \
            -max_target_seqs 1 \
            -quiet 2>> "$LOG" > "$BLAST_FILE"
    
    # 3.5 解析批量比对结果
    GENE_CNT=0
    while read -r BLAST_LINE; do
        QSEQID=$(echo "$BLAST_LINE" | awk '{print $1}')  # 基因ID
        VFDB_GENE=$(echo "$BLAST_LINE" | awk '{print $2}')
        PID=$(echo "$BLAST_LINE" | awk '{print $3}')
        EVAL=$(echo "$BLAST_LINE" | awk '{print $4}')
        
        # 匹配基因名和功能描述
        GENE_INFO=$(grep -w "$QSEQID" "$BED_FILE" | head -1)
        [ -z "$GENE_INFO" ] && continue
        SEQ_ID=$(echo "$GENE_INFO" | awk '{print $1}')
        START=$(echo "$GENE_INFO" | awk '{print $2+1}')  # 还原为原始坐标
        END=$(echo "$GENE_INFO" | awk '{print $3}')
        COORD="$SEQ_ID:$START-$END"
        
        # 从Bakta注释获取基因名和功能
        GENE_NAME=$(grep -w "$QSEQID" "$GENE_INFO_FILE" | awk '{print $2}' || "未知")
        FUNC=$(grep -w "$QSEQID" "$GENE_INFO_FILE" | awk '{print $3}' | sed 's/,/;/g' || "无")
        
        # VFDB官方描述
        VFDB_DESC=$(grep "^$VFDB_GENE|" "$PROJECT_DIR/vfdb_official.tmp" | sed 's/|/\t/g' | awk '{print $2}' || "无")
        
        # 写入结果
        echo "$STRAIN,$QSEQID,$GENE_NAME,$FUNC,$COORD,$SEQ_ID,$VFDB_GENE,$VFDB_DESC,$PID,$EVAL,VFDB官方库（华大基因）" >> "$OUTPUT"
        GENE_CNT=$((GENE_CNT+1))
    done < "$BLAST_FILE"
    
    echo "   ✅ 菌株$STRAIN：成功提取$GENE_CNT个毒力基因" | tee -a "$LOG"
    
    # 清理该菌株临时文件
    rm -f "$BED_FILE" "$FASTA_FILE" "$BLAST_FILE" "$GENE_INFO_FILE"
done

# 4. 清理全局临时文件
rm -f "$PROJECT_DIR/vfdb_official.tmp"

echo -e "\n🎉 批量处理完成！所有菌株已跑完！" | tee -a "$LOG"
echo "📁 最终结果：$OUTPUT"
echo "📄 处理日志：$LOG"

# 5. 成果统计
TOTAL=$(grep -v "^菌株名" "$OUTPUT" | wc -l)
HIGH_CONFID=$(grep -E ",[7-9][0-9]\.|,100\." "$OUTPUT" | wc -l)
STRAINS=$(grep -v "^菌株名" "$OUTPUT" | awk -F "," '{print $1}' | sort -u | wc -l)

echo -e "\n===== 最终成果统计 ====="
echo "1. 总毒力基因数（VFDB官方匹配）：$TOTAL"
echo "2. 高可信度基因数（相似度≥70%）：$HIGH_CONFID"
echo "3. 覆盖菌株数：$STRAINS 株"
echo "4. 平均每株毒力基因数：$( [ $STRAINS -gt 0 ] && echo "scale=1; $TOTAL/$STRAINS" | bc || echo 0 ) 个"
echo -e "\n📚 论文引用格式：Chen, L., et al. (2016). VFDB 2016: hierarchical and refined dataset for bacterial virulence factors. Nucleic Acids Res, 44(D1), D694-D698."
