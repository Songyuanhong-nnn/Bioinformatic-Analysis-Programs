#!/bin/bash
ROARY_CSV="roary_output/gene_presence_absence.csv"
OUTPUT_DIR="final_markers"
mkdir -p $OUTPUT_DIR

# 筛选高危害特有基因
HEADERS=$(head -n1 $ROARY_CSV | tr ',' '\n')
HIGH=(${HIGH_HAZARD//,/ })
LOW=(${LOW_MEDIUM_HAZARD//,/ })

> $OUTPUT_DIR/high_hazard_unique.txt
tail -n+2 $ROARY_CSV | while IFS=',' read -r gene rest; do
    presence=($(echo $rest | tr ',' '\n'))
    # 高危害全部存在
    all_high=1
    for g in "${HIGH[@]}"; do
        idx=$(echo "$HEADERS" | grep -n "^$g$" | cut -d':' -f1)
        [ -z "$idx" ] && { all_high=0; break; }
        [ "${presence[$((idx-2))]}" -eq 0 ] && { all_high=0; break; }
    done
    # 低/中危害全部缺失
    all_low=1
    for g in "${LOW[@]}"; do
        idx=$(echo "$HEADERS" | grep -n "^$g$" | cut -d':' -f1)
        [ -z "$idx" ] && { all_low=0; break; }
        [ "${presence[$((idx-2))]}" -ne 0 ] && { all_low=0; break; }
    done
    [ $all_high -eq 1 ] && [ $all_low -eq 1 ] && echo $gene >> $OUTPUT_DIR/high_hazard_unique.txt
done

# 提取特有基因序列并关联VFDB毒力功能
python3 - << 'PYSCRIPT'
from Bio import SeqIO
import os

# 读取特有基因
with open('final_markers/high_hazard_unique.txt', 'r') as f:
    unique_genes = [l.strip() for l in f if l.strip()]

# 提取序列
with open('final_markers/unique_genes.faa', 'w') as outf:
    for gid in os.listdir('prokka_annotations'):
        faa_file = f'prokka_annotations/{gid}/{gid}.faa'
        for record in SeqIO.parse(faa_file, 'fasta'):
            if record.id in unique_genes:
                SeqIO.write(record, outf, 'fasta')
                unique_genes.remove(record.id)
                if not unique_genes:
                    break
        if not unique_genes:
            break

# BLAST比对VFDB（复用已有数据库）
os.system('blastp -query final_markers/unique_genes.faa \
                -db vfdb_online/VFDB_setA \
                -out final_markers/vf_blast.blast \
                -outfmt "6 qseqid sseqid pident evalue stitle" \
                -evalue 1e-10 \
                -pident 70 \
                -num_threads 4')

# 整理结果
with open('final_markers/vf_blast.blast', 'r') as f, \
     open('final_markers/marker_virulence_annot.tsv', 'w') as outf:
    outf.write('Marker_Gene\tVFDB_Gene\tIdentity(%)\tE_value\tDescription\n')
    for line in f:
        parts = line.strip().split('\t')
        outf.write(f'{parts[0]}\t{parts[1]}\t{parts[2]}\t{parts[3]}\t{parts[4]}\n')

print("标志物筛选完成！")
PYSCRIPT
