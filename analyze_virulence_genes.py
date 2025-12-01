import pandas as pd
import numpy as np

def analyze_virulence_patterns():
    """直接分析毒力基因分布模式"""
    print("=== 毒力基因分布模式分析 ===")
    
    # 读取Roary数据
    gene_pa = pd.read_csv('roary_output_final_1762658265/gene_presence_absence.csv', low_memory=False)
    
    # 寻找注释中包含毒力相关关键词的基因
    virulence_keywords = [
        'toxin', 'hemolysin', 'virulence', 'effector', 'adhesin', 
        'invasion', 'secretion', 'transferase', 'capsule', 'antigen',
        'hemolysin', 'enterotoxin', 'cytotoxin', 'leukocidin'
    ]
    
    print("寻找毒力相关基因...")
    virulence_genes = []
    
    for idx, row in gene_pa.iterrows():
        annotation = str(row['Annotation']).lower()
        if any(keyword in annotation for keyword in virulence_keywords):
            # 计算该基因的分布统计
            strain_columns = gene_pa.columns[14:]
            presence_count = sum(pd.notna(row[col]) for col in strain_columns)
            presence_freq = presence_count / len(strain_columns)
            
            virulence_genes.append({
                'Gene': row['Gene'],
                'Annotation': row['Annotation'],
                'Presence_Frequency': presence_freq,
                'Presence_Count': presence_count,
                'Total_Strains': len(strain_columns)
            })
    
    print(f"找到 {len(virulence_genes)} 个毒力相关基因")
    
    # 按出现频率排序
    virulence_df = pd.DataFrame(virulence_genes)
    if len(virulence_df) > 0:
        virulence_df = virulence_df.sort_values('Presence_Frequency', ascending=False)
        
        print(f"\n=== 前20个毒力相关基因 ===")
        for i, row in virulence_df.head(20).iterrows():
            print(f"{i+1:2d}. {row['Gene']:20} {row['Presence_Frequency']:.1%} ({row['Presence_Count']}/{row['Total_Strains']})")
            print(f"     {row['Annotation']}\n")
        
        # 保存结果
        virulence_df.to_csv('virulence_related_genes.csv', index=False)
        print("✓ 毒力基因列表已保存到: virulence_related_genes.csv")
    
    return virulence_df

def find_genes_with_variable_distribution(gene_pa):
    """寻找分布变异性高的基因（可能代表特异性）"""
    print("\n=== 寻找分布变异性高的基因 ===")
    
    strain_columns = gene_pa.columns[14:]
    variable_genes = []
    
    for idx, row in gene_pa.iterrows():
        if idx % 5000 == 0:
            print(f"处理进度: {idx}/{len(gene_pa)}")
        
        presence_count = sum(pd.notna(row[col]) for col in strain_columns)
        presence_freq = presence_count / len(strain_columns)
        
        # 选择出现频率在20%-80%之间的基因（既不是核心基因也不是罕见基因）
        if 0.2 <= presence_freq <= 0.8:
            variable_genes.append({
                'Gene': row['Gene'],
                'Annotation': row['Annotation'],
                'Presence_Frequency': presence_freq,
                'Presence_Count': presence_count,
                'Total_Strains': len(strain_columns)
            })
    
    variable_df = pd.DataFrame(variable_genes)
    variable_df = variable_df.sort_values('Presence_Frequency', ascending=False)
    
    print(f"找到 {len(variable_df)} 个分布变异性高的基因")
    
    # 显示前20个
    print(f"\n=== 前20个高变异性基因 ===")
    for i, row in variable_df.head(20).iterrows():
        annotation = str(row['Annotation'])[:70] + "..." if len(str(row['Annotation'])) > 70 else row['Annotation']
        print(f"{i+1:2d}. {row['Gene']:20} {row['Presence_Frequency']:.1%}")
        print(f"     {annotation}\n")
    
    variable_df.to_csv('variable_distribution_genes.csv', index=False)
    print("✓ 高变异性基因列表已保存到: variable_distribution_genes.csv")
    
    return variable_df

if __name__ == "__main__":
    # 分析毒力基因
    virulence_df = analyze_virulence_patterns()
    
    # 读取Roary数据用于变异性分析
    gene_pa = pd.read_csv('roary_output_final_1762658265/gene_presence_absence.csv', low_memory=False)
    
    # 寻找高变异性基因
    variable_df = find_genes_with_variable_distribution(gene_pa)
    
    print("\n=== 分析完成 ===")
    print("建议下一步:")
    print("1. 从 virulence_related_genes.csv 中选择候选基因进行验证")
    print("2. 从 variable_distribution_genes.csv 中选择分布特异的基因")
    print("3. 结合文献验证这些基因的功能")
