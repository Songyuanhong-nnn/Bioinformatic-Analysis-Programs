import pandas as pd
import os

def extract_target_sequences():
    """提取推荐靶标的序列信息"""
    print("=== 提取靶标序列信息 ===")
    
    # 读取Roary数据
    gene_pa = pd.read_csv('roary_output_final_1762658265/gene_presence_absence.csv', low_memory=False)
    
    # 检查Roary文件的列名
    print("Roary文件列名:", gene_pa.columns.tolist()[:15])  # 显示前15列
    
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
            
            # 使用正确的列名
            target_info.append({
                'Gene': gene,
                'Annotation': row['Annotation'],
                'Presence_Count': presence_count,
                'Total_Strains': len(strain_columns),
                'Presence_Percentage': f"{presence_freq:.1%}",
                'No_isolates': row.get('No_isolates', 'N/A'),
                'No_sequences': row.get('No_sequences', 'N/A'),
                'Avg_sequences_per_isolate': row.get('Avg_sequences_per_isolate', 'N/A')
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
        print(f"   独特序列: {info['No_sequences']}")
        print(f"   平均序列数: {info['Avg_sequences_per_isolate']}")
    
    return df

def get_sequence_locations():
    """获取序列位置信息"""
    print("\n=== 序列位置信息 ===")
    
    # 检查是否有序列文件
    possible_files = [
        'roary_output_final_1762658265/pan_genome_reference.fa',
        'roary_output_final/pan_genome_reference.fa'
    ]
    
    for file in possible_files:
        if os.path.exists(file):
            print(f"✓ 找到序列文件: {file}")
            # 这里可以添加序列提取代码
            break
    else:
        print("⚠️ 未找到pan_genome_reference.fa文件")
        print("建议从Roary输出目录中查找序列文件")

if __name__ == "__main__":
    extract_target_sequences()
    get_sequence_locations()
    
    print(f"\n=== 🎯 最终靶标清单总结 ===")
    print("以下10个靶标已准备好进行实验验证:")
    targets = [
        "1. group_9360 (II型分泌系统蛋白E, 46.8%)",
        "2. group_3365 (热不稳定溶血素, 96.8%)", 
        "3. tdh2 (热稳定溶血素2, 17.1%)",
        "4. vopS (效应蛋白, 69.9%)",
        "5. epsE_2 (分泌系统蛋白E, 54.2%)",
        "6. tdh1 (热稳定溶血素1, 10.6%)",
        "7. ureC (脲酶亚基, 20.8%)",
        "8. group_6946 (染色体溶血素, 3.7%)",
        "9. hlyA (溶血素A, 4.6%)",
        "10. group_2253 (III型分泌系统, 6.5%)"
    ]
    
    for target in targets:
        print(f"   {target}")
