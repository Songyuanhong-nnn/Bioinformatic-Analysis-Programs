#!/bin/bash
# 步骤1：整合血清型到表型文件
awk -F ',' 'NR==FNR{a[$1]=$2;next}{print $0 "," a[$1]}' serotype_final_results.csv phenotype.txt > phenotype_with_serotype.txt

# 步骤2：筛选高毒力专属基因（核心靶标第一步）
export HIGH_RISK_STRAINS=$(awk -F ',' '$2==1{print $1}' phenotype_with_serotype.txt | tr '\n' ',')
export LOW_RISK_STRAINS=$(awk -F ',' '$2==0{print $1}' phenotype_with_serotype.txt | tr '\n' ',')

python3 - << 'PYSCRIPT'
import pandas as pd
# 读取Roary核心矩阵
df = pd.read_csv("roary_output_final/gene_presence_absence.csv", low_memory=False)
# 筛选：所有高毒力菌株都有，所有低毒力菌株都没有的基因
high_cols = [col for col in df.columns if col.strip() in os.environ["HIGH_RISK_STRAINS"].split(',')]
low_cols = [col for col in df.columns if col.strip() in os.environ["LOW_RISK_STRAINS"].split(',')]

# 高毒力全有（值≠0），低毒力全没有（值=0）
df["is_target"] = (df[high_cols].astype(str) != '0').all(axis=1) & (df[low_cols].astype(str) == '0').all(axis=1)
target_genes = df[df["is_target"]]

# 保存候选靶标（包含基因名、注释、序列信息）
target_genes[["Gene", "Annotation", "Sequence"]].to_csv("candidate_virulence_genes.csv", index=False)
print(f"筛选出 {len(target_genes)} 个高毒力专属候选靶标")
PYSCRIPT

# 步骤3：关联VFDB验证（确认是否为已知毒力基因）
blastp -query <(awk -F ',' 'NR>1{print ">"$1"\n"$3}' candidate_virulence_genes.csv) \
       -db vfdb_online/VFDB_setA_protein.faa \
       -evalue 1e-10 -num_threads 8 -outfmt 6 -max_target_seqs 1 \
       -out candidate_vs_vfdb.blastp
echo "VFDB验证完成，可查看 candidate_vs_vfdb.blastp 确认已知毒力基因"
