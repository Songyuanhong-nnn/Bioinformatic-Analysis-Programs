import pandas as pd
import numpy as np
import os

def load_roary_data():
    """
    读取真实的Roary输出数据
    """
    # 尝试不同的可能路径
    possible_paths = [
        'roary_output_final/gene_presence_absence.csv',
        'roary_output_final_1762658265/gene_presence_absence.csv',
        'gene_presence_absence.csv'
    ]
    
    for path in possible_paths:
        if os.path.exists(path):
            print(f"✓ 找到Roary文件: {path}")
            try:
                gene_pa = pd.read_csv(path)
                print(f"  成功读取，共 {len(gene_pa)} 个基因簇")
                return gene_pa, path
            except Exception as e:
                print(f"  读取失败: {e}")
    
    print("❌ 未找到Roary输出文件")
    return None, None

def create_serotype_groups():
    """
    创建血清型分组 - 需要根据你的实际情况修改
    """
    # 这里需要你提供真实的菌株血清型信息
    # 示例：根据你的kaptive结果或已知信息
    
    # 从你的目录结构看，可能有这些信息
    try:
        # 尝试读取kaptive结果
        kaptive_o = pd.read_csv('kaptive_o_serotype_results.tsv', sep='\t')
        kaptive_k = pd.read_csv('kaptive_k_serotype_results.tsv', sep='\t')
        print("✓ 读取Kaptive血清型预测结果")
        
        # 创建血清型分组
        serotype_groups = {}
        for idx, row in kaptive_o.iterrows():
            strain = row['Genome']
            o_type = row.get('Type', 'Unknown')
            confidence = row.get('Match_Confidence', 'Unknown')
            serotype_groups[strain] = {
                'O_serotype': o_type,
                'Confidence': confidence
            }
        
        return serotype_groups
    except Exception as e:
        print(f"⚠️ 无法读取Kaptive结果: {e}")
        return None

def find_high_risk_specific_genes(gene_pa, high_risk_serotypes=None):
    """
    根据真实数据找出高危血清型特有基因
    """
    if high_risk_serotypes is None:
        high_risk_serotypes = ['O3:K6', 'O10:K4']  # 根据文献调整
    
    # 获取菌株列
    strain_columns = gene_pa.columns[14:]
    print(f"分析 {len(strain_columns)} 个菌株")
    
    # 显示前10个菌株名，方便你识别
    print("前10个菌株名称:")
    for i, strain in enumerate(strain_columns[:10]):
        print(f"  {i+1}. {strain}")
    
    # 由于我们没有完整的血清型分组，这里提供两种方法：
    
    print("\n请选择分析方法:")
    print("1. 手动指定高危菌株索引（如前10个为高危组）")
    print("2. 基于菌株名称模式识别")
    
    # 方法1：简单分组（你需要根据实际情况调整）
    high_risk_indices = list(range(min(10, len(strain_columns))))  # 前10个作为高危组示例
    high_risk_strains = [strain_columns[i] for i in high_risk_indices]
    other_strains = [strain for i, strain in enumerate(strain_columns) if i not in high_risk_indices]
    
    print(f"\n分组情况:")
    print(f"高危组: {len(high_risk_strains)} 个菌株")
    print(f"其他组: {len(other_strains)} 个菌株")
    
    # 找出特异性基因
    specific_genes = []
    
    for idx, row in gene_pa.iterrows():
        # 计算在高危组中的出现频率
        high_risk_presence = sum(pd.notna(row[col]) for col in high_risk_strains)
        high_risk_freq = high_risk_presence / len(high_risk_strains) if high_risk_strains else 0
        
        # 计算在其他组中的出现频率
        other_presence = sum(pd.notna(row[col]) for col in other_strains)
        other_freq = other_presence / len(other_strains) if other_strains else 0
        
        # 筛选条件：高危组频率 > 80%，其他组频率 < 20%
        if high_risk_freq >= 0.8 and other_freq <= 0.2:
            gene_info = {
                'Gene': row['Gene'],
                'Annotation': row['Annotation'],
                'HighRisk_Frequency': high_risk_freq,
                'Other_Frequency': other_freq,
                'No_in_HighRisk': high_risk_presence,
                'No_in_Other': other_presence,
                'No_absent_in_HighRisk': len(high_risk_strains) - high_risk_presence,
                'No_absent_in_Other': len(other_strains) - other_presence
            }
            specific_genes.append(gene_info)
    
    return pd.DataFrame(specific_genes)

def main():
    """
    主函数
    """
    print("=== 真实数据分析 ===")
    
    # 1. 读取Roary数据
    gene_pa, data_path = load_roary_data()
    if gene_pa is None:
        print("请确保Roary输出文件存在")
        return
    
    # 2. 尝试读取血清型信息
    serotype_groups = create_serotype_groups()
    
    # 3. 找出特异性基因
    print("\n正在分析特异性基因...")
    specific_genes_df = find_high_risk_specific_genes(gene_pa)
    
    if len(specific_genes_df) == 0:
        print("❌ 未找到符合条件的特异性基因")
        print("建议调整分组策略或筛选阈值")
        return
    
    print(f"\n✓ 找到 {len(specific_genes_df)} 个候选特异性基因")
    
    # 4. 保存结果
    specific_genes_df.to_csv('high_risk_specific_genes_REAL.csv', index=False)
    print("✓ 真实数据结果已保存到: high_risk_specific_genes_REAL.csv")
    
    # 5. 显示统计信息
    print(f"\n=== 统计信息 ===")
    print(f"总候选基因数: {len(specific_genes_df)}")
    
    # 按注释类型统计
    annotation_stats = specific_genes_df['Annotation'].apply(
        lambda x: 'hypothetical' if 'hypothetical' in str(x).lower() else 'annotated'
    ).value_counts()
    
    print("注释情况:")
    for ann_type, count in annotation_stats.items():
        print(f"  {ann_type}: {count}")
    
    # 显示前20个基因
    print(f"\n=== 前20个候选基因 ===")
    for i, row in specific_genes_df.head(20).iterrows():
        annotation_preview = row['Annotation'][:60] + "..." if len(str(row['Annotation'])) > 60 else row['Annotation']
        print(f"{i+1:2d}. {row['Gene']:20} 高危: {row['HighRisk_Frequency']:.1%} 其他: {row['Other_Frequency']:.1%}")
        print(f"     {annotation_preview}")

if __name__ == "__main__":
    main()
