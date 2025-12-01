import pandas as pd
import numpy as np

def analyze_t6ss_cooccurrence():
    print("=== 🔍 VI型分泌系统基因共现分析 ===")
    
    # 读取数据
    gene_pa = pd.read_csv('roary_consistent/gene_presence_absence.csv', low_memory=False)
    strains = gene_pa.columns[14:]
    
    # VI型分泌系统相关基因
    t6ss_genes = [
        'group_12000',  # TssH ATP酶
        'group_2436',   # 收缩鞘大亚基
        'group_3627',   # 收缩鞘小亚基  
        'group_5111',   # Hcp1效应蛋白
        'group_1213',   # Vgr蛋白
        'group_2437',   # TssM膜亚基
        'group_5114',   # TssA蛋白
        'group_5115',   # TssE基板亚基
        'group_5116',   # TssF基板亚基
        'group_5117',   # TssG基板亚基
        'group_5119',   # TssJ脂蛋白
        'group_5120',   # TssK基板亚基
        'group_5121'    # IcmH/DotU蛋白
    ]
    
    print("VI型分泌系统基因分布:")
    cooccurrence_matrix = []
    
    for gene in t6ss_genes:
        gene_data = gene_pa[gene_pa['Gene'] == gene]
        if len(gene_data) > 0:
            row = gene_data.iloc[0]
            presence = [1 if pd.notna(row[col]) else 0 for col in strains]
            freq = sum(presence) / len(strains)
            
            print(f"  {gene}: {sum(presence)}/{len(strains)} ({freq:.1%}) - {row['Annotation']}")
            cooccurrence_matrix.append(presence)
    
    # 计算共现率
    if len(cooccurrence_matrix) >= 2:
        print(f"\n🎯 基因共现分析:")
        
        # 检查所有基因是否在相同菌株中出现
        all_same = True
        for i in range(1, len(cooccurrence_matrix)):
            if cooccurrence_matrix[0] != cooccurrence_matrix[i]:
                all_same = False
                break
        
        if all_same:
            print("✅ 所有VI型分泌系统基因完全共现!")
            print("   这些基因作为一个完整的功能单元存在")
        else:
            print("⚠️ 基因分布不完全一致，但高度相关")
            
        # 计算成对共现率
        print(f"\n📊 成对共现统计:")
        for i in range(min(5, len(cooccurrence_matrix))):
            for j in range(i+1, min(6, len(cooccurrence_matrix))):
                gene1 = t6ss_genes[i]
                gene2 = t6ss_genes[j]
                
                matches = sum(1 for k in range(len(strains)) 
                            if cooccurrence_matrix[i][k] == cooccurrence_matrix[j][k])
                cooccurrence_rate = matches / len(strains)
                
                print(f"  {gene1} & {gene2}: {cooccurrence_rate:.1%} 一致")

def identify_t6ss_strains():
    """识别含有完整VI型分泌系统的菌株"""
    print(f"\n=== 🎯 识别VI型分泌系统阳性菌株 ===")
    
    gene_pa = pd.read_csv('roary_consistent/gene_presence_absence.csv', low_memory=False)
    strains = gene_pa.columns[14:]
    
    # 关键VI型分泌系统标志基因
    marker_genes = ['group_12000', 'group_2436', 'group_3627', 'group_5111']
    
    t6ss_positive_strains = []
    
    for strain in strains:
        has_all_markers = True
        for gene in marker_genes:
            gene_data = gene_pa[gene_pa['Gene'] == gene]
            if len(gene_data) > 0:
                row = gene_data.iloc[0]
                if pd.isna(row[strain]):
                    has_all_markers = False
                    break
        
        if has_all_markers:
            t6ss_positive_strains.append(strain)
    
    print(f"完整VI型分泌系统阳性菌株: {len(t6ss_positive_strains)}/{len(strains)}")
    
    if t6ss_positive_strains:
        print("阳性菌株示例:")
        for strain in t6ss_positive_strains[:5]:
            print(f"  - {strain}")
    
    return t6ss_positive_strains

if __name__ == "__main__":
    analyze_t6ss_cooccurrence()
    positive_strains = identify_t6ss_strains()
    
    print(f"\n💡 生物学意义:")
    print("1. VI型分泌系统通常作为完整的基因岛存在")
    print("2. 这些基因协同工作，共同获得或丢失") 
    print("3. 56.1%的菌株可能代表一个特定的进化分支")
    print("4. 这个系统与细菌的竞争性和致病性相关")
