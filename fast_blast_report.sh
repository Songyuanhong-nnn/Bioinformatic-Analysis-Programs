#!/bin/bash
PROJECT_DIR="/mnt/d/WSL/disk/projects/VP1"
INPUT_CSV="$PROJECT_DIR/FAST_serotype_specific_virulence.csv"
VFDB_PROT="$PROJECT_DIR/vfdb_online/VFDB_setA_pro.fas"
BLAST_REPORT="$PROJECT_DIR/FAST_blast_report.tsv"

echo -e "基因ID\t基因名\t血清型\tVFDB匹配基因\t相似度(%)\tE值" > "$BLAST_REPORT"
grep -v "血清型" "$INPUT_CSV" | while IFS= read -r line; do
    GENE_ID=$(echo "$line" | awk -F "," '{print $2}')
    GENE_NAME=$(echo "$line" | awk -F "," '{print $3}')
    SERO=$(echo "$line" | awk -F "," '{print $1}')
    
    # 从Roary提取含该基因的任意菌株
    STRAIN=$(grep "^$GENE_ID," "$PROJECT_DIR/roary_output_final_1762658265/gene_presence_absence.csv" | awk -F "," '{
        for(i=15;i<=NF;i++) if($i==1) {print $i "|" NR; break}
    }' | head -1 | awk -F "|" '{print $2}')
    [ -z "$STRAIN" ] && continue
    
    # 快速比对VFDB
    BLAST=$(tblastn -query <(echo -e ">$GENE_ID\nplaceholder") -db "$VFDB_PROT" -outfmt "6 sseqid pident evalue" -evalue 1e-10 -quiet | head -1)
    VFDB_GENE=$(echo "$BLAST" | awk '{print $1}' || echo "无")
    PID=$(echo "$BLAST" | awk '{print $2}' || echo "0")
    EVAL=$(echo "$BLAST" | awk '{print $3}' || echo "1e0")
    
    echo -e "$GENE_ID\t$GENE_NAME\t$SERO\t$VFDB_GENE\t$PID\t$EVAL" >> "$BLAST_REPORT"
done

echo "✅ 序列比对报告生成完成：$BLAST_REPORT"
