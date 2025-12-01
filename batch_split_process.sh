#!/bin/bash
set -euo pipefail

PROJECT_DIR="/mnt/d/WSL/disk/projects/VP1"
VFDB_PROT="$PROJECT_DIR/vfdb_online/VFDB_setA_pro.fas"
OUTPUT="$PROJECT_DIR/FINAL_BATCH_VIRULENCE.csv"
LOG="$PROJECT_DIR/batch_split_log.txt"

# 初始化结果文件
[ ! -f "$OUTPUT" ] && echo "菌株名,基因ID,基因名,功能描述,坐标,VFDB基因,VFDB描述,相似度(%),E值,官方来源" > "$OUTPUT"

# 获取所有菌株（排除已处理的）
ALL_STRAINS=($(ls "$PROJECT_DIR/bakta_annotations" | grep -v "^$"))
PROCESSED_STRAINS=($(grep -v "^菌株名" "$OUTPUT" | awk -F "," '{print $1}' | sort -u))
TO_PROCESS=()
for STRAIN in "${ALL_STRAINS[@]}"; do
    if ! [[ " ${PROCESSED_STRAINS[@]} " =~ " $STRAIN " ]]; then
        TO_PROCESS+=("$STRAIN")
    fi
done

echo "===== 分批处理：共${#TO_PROCESS[@]}株待处理 =====" > "$LOG"
date >> "$LOG"

# 每批处理10株
BATCH_SIZE=10
for ((i=0; i<${#TO_PROCESS[@]}; i+=BATCH_SIZE)); do
    BATCH=("${TO_PROCESS[@]:i:BATCH_SIZE}")
    echo -e "\n🔧 处理批次$((i/BATCH_SIZE+1))：${BATCH[*]}" | tee -a "$LOG"
    
    for STRAIN in "${BATCH[@]}"; do
        BAKTA_TSV="$PROJECT_DIR/bakta_annotations/$STRAIN/$STRAIN.tsv"
        GENOME_FNA="$PROJECT_DIR/clean_genome/$STRAIN.fna"
        [ ! -f "$BAKTA_TSV" ] || [ ! -f "$GENOME_FNA" ] && {
            echo "⚠️  $STRAIN：缺少文件，跳过" | tee -a "$LOG"
            continue
        }
        
        # 提取FASTA ID
        DEFAULT_FASTA_ID=$(grep "^>" "$GENOME_FNA" | sed 's/^>//' | awk '{print $1}' | head -1)
        [ -z "$DEFAULT_FASTA_ID" ] && {
            echo "⚠️  $STRAIN：无FASTA ID，跳过" | tee -a "$LOG"
            continue
        }
        
        # 提取有效基因
        TMP_BED="$PROJECT_DIR/tmp_$STRAIN.bed"
        grep -v "^#" "$BAKTA_TSV" | awk -F "\t" -v def_id="$DEFAULT_FASTA_ID" '
        NF >=7 && $3 ~ /^[0-9]+$/ && $4 ~ /^[0-9]+$/ && $3 < $4 {
            print def_id "\t" $((3-1)) "\t" $4 "\t" $1 "\t" $2 "\t" $5
        }' > "$TMP_BED"
        [ ! -s "$TMP_BED" ] && {
            echo "⚠️  $STRAIN：无有效基因，跳过" | tee -a "$LOG"
            rm -f "$TMP_BED"
            continue
        }
        
        # 提取序列+比对
        TMP_SEQ="$PROJECT_DIR/tmp_$STRAIN.seq.fna"
        bedtools getfasta -fi "$GENOME_FNA" -bed <(awk '{print $1 "\t" $2 "\t" $3 "\t" $4}' "$TMP_BED") -fo "$TMP_SEQ" 2>> "$LOG"
        TMP_BLAST="$PROJECT_DIR/tmp_$STRAIN.blast"
        tblastn -query "$TMP_SEQ" -db "$VFDB_PROT" -outfmt "6 qseqid sseqid pident evalue" -evalue 1e-3 -num_threads 8 -max_target_seqs 1 -quiet 2>> "$LOG" > "$TMP_BLAST"
        
        # 解析结果
        GENE_CNT=0
        while read -r BLAST_LINE; do
            QSEQID=$(echo "$BLAST_LINE" | awk '{print $1}')
            VFDB_GENE=$(echo "$BLAST_LINE" | awk '{print $2}')
            PID=$(echo "$BLAST_LINE" | awk '{print $3}')
            EVAL=$(echo "$BLAST_LINE" | awk '{print $4}')
            
            INFO=$(grep -w "$QSEQID" "$TMP_BED")
            SEQ_ID=$(echo "$INFO" | awk '{print $1}')
            START=$(echo "$INFO" | awk '{print $2+1}')
            END=$(echo "$INFO" | awk '{print $3}')
            GENE_NAME=$(echo "$INFO" | awk '{print $5}')
            FUNC=$(echo "$INFO" | awk '{print $6}' | sed 's/,/;/g')
            COORD="$SEQ_ID:$START-$END"
            VFDB_DESC=$(grep "^$VFDB_GENE|" <(grep "^>" "$VFDB_PROT" | sed 's/^>//; s/\t/|/g') | sed 's/|/\t/g' | awk '{print $2}' || "无")
            
            echo "$STRAIN,$QSEQID,$GENE_NAME,$FUNC,$COORD,$VFDB_GENE,$VFDB_DESC,$PID,$EVAL,VFDB官方库" >> "$OUTPUT"
            GENE_CNT=$((GENE_CNT+1))
        done < "$TMP_BLAST"
        
        echo "✅ $STRAIN：找到$GENE_CNT个毒力基因" | tee -a "$LOG"
        rm -f "$TMP_BED" "$TMP_SEQ" "$TMP_BLAST"
    done
done

echo -e "\n🎉 所有批次处理完成！" | tee -a "$LOG"
echo "最终结果：$OUTPUT"
echo "日志文件：$LOG"

# 统计结果
TOTAL=$(grep -v "^菌株名" "$OUTPUT" | wc -l)
STRAINS=$(grep -v "^菌株名" "$OUTPUT" | awk -F "," '{print $1}' | sort -u | wc -l)
echo -e "\n===== 统计 ====="
echo "总毒力基因数：$TOTAL"
echo "覆盖菌株数：$STRAINS"
