import pandas as pd
import os

def extract_target_sequences():
    """提取推荐靶标的序列信息"""
    print("=== 提取靶标序列信息 ===")
    
    # 读取Roary数据
    gene_pa = pd.read_csv('roary_output_final_1762658265/gene_presence_absence.csv', low_memory=False)
    
    # 目标基因列表
    target_genes = [
        'group_9360', 'group_3365', 'tdh2', 'vopS', 'epsE_2',
        'tdh1', 'ureC', 'group_6946', 'hlyA', 'group_2253'
    ]
    
    # 查找这些基因的详细信息
    target_info = []
    
    for gene in target_genes:
        gene_data = gene_pa[gene_pa['Gene'] == gene]
        
        if len(gene_data) > 0:
            row = gene_data.iloc[0]
            
            # 计算准确的分布
            strain_columns = gene_pa.columns[14:]
            presence_count = sum(pd.notna(row[col]) for col in strain_columns)
            presence_freq = presence_count / len(strain_columns)
            
            target_info.append({
                'Gene': gene,
                'Annotation': row['Annotation'],
                'Presence_Count': presence_count,
                'Total_Strains': len(strain_columns),
                'Presence_Percentage': f"{presence_freq:.1%}",
                'Non_annotated': row['Non_annotated'],
                'No_isolates': row['No_isolates'],
                'No_sequences': row['No_sequences'],
                'Avg_sequences_per_isolate': row['Avg_sequences_per_isolate']
            })
        else:
            print(f"警告: 未找到基因 {gene}")
    
    # 创建详细报告
    df = pd.DataFrame(target_info)
    df.to_csv('target_genes_detailed_info.csv', index=False)
    
    print(f"\n✓ 成功提取 {len(df)} 个靶标的详细信息")
    print("✓ 保存到: target_genes_detailed_info.csv")
    
    # 显示详细信息
    print(f"\n=== 靶标详细信息 ===")
    for info in target_info:
        print(f"\n🎯 {info['Gene']}")
        print(f"   功能: {info['Annotation']}")
        print(f"   分布: {info['Presence_Count']}/{info['Total_Strains']} ({info['Presence_Percentage']})")
        print(f"   非注释菌株: {info['Non_annotated']}")
        print(f"   独特序列: {info['No_sequences']}")
    
    return df

if __name__ == "__main__":
    extract_target_sequences()
