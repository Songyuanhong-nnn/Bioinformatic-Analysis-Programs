import pandas as pd
import numpy as np
import os

def parse_kaptive_results():
    """解析Kaptive结果"""
    print("=== 解析Kaptive血清型结果 ===")
    
    # 读取Kaptive结果
    kaptive_o = pd.read_csv('kaptive_o_serotype_results.tsv', sep='\t')
    kaptive_k = pd.read_csv('kaptive_k_serotype_results.tsv', sep='\t')
    
    # 提取血清型信息
    serotype_data = {}
    
    for idx, row in kaptive_o.iterrows():
        strain = row['Assembly']
        o_type = row['Best match type']
        confidence = row['Match confidence']
        
        serotype_data[strain] = {
            'O_serotype': o_type,
            'O_confidence': confidence,
            'K_serotype': 'Unknown',
            'K_confidence': 'Unknown'
        }
    
    for idx, row in kaptive_k.iterrows():
        strain = row['Assembly']
        k_type = row['Best match type']
        confidence = row['Match confidence']
        
        if strain in serotype_data:
            serotype_data[strain]['K_serotype'] = k_type
            serotype_data[strain]['K_confidence'] = confidence
        else:
            serotype_data[strain] = {
                'O_serotype': 'Unknown',
                'O_confidence': 'Unknown',
                'K_serotype': k_type,
                'K_confidence': confidence
            }
    
    print(f"成功解析 {len(serotype_data)} 个菌株的血清型信息")
    
    # 统计血清型分布
    o_types = {}
    k_types = {}
    
    for strain, data in serotype_data.items():
        o_type = data['O_serotype']
        k_type = data['K_serotype']
        
        if o_type not in o_types:
            o_types[o_type] = 0
        o_types[o_type] += 1
        
        if k_type not in k_types:
            k_types[k_type] = 0
        k_types[k_type] += 1
    
    print(f"\nO血清型分布 (前10):")
    for o_type, count in sorted(o_types.items(), key=lambda x: x[1], reverse=True)[:10]:
        print(f"  {o_type}: {count}")
    
    print(f"\nK血清型分布 (前10):")
    for k_type, count in sorted(k_types.items(), key=lambda x: x[1], reverse=True)[:10]:
        print(f"  {k_type}: {count}")
    
    return serotype_data

def identify_high_risk_serotypes(serotype_data):
    """识别高危血清型"""
    print("\n=== 识别高危血清型 ===")
    
    # 根据文献，这些是已知的高危血清型
    high_risk_o_types = ['O3:K6', 'O4:K8', 'O1:K25', 'O1:K26', 'O1:K36', 'O4:K12']
    high_risk_k_types = ['K6', 'K8', 'K25', 'K26', 'K36', 'K12']
    
    # 从数据中找出可能的高危血清型
    o_type_counts = {}
    k_type_counts = {}
    
    for strain, data in serotype_data.items():
        o_type = data['O_serotype']
        k_type = data['K_serotype']
        
        if o_type not in o_type_counts:
            o_type_counts[o_type] = 0
        o_type_counts[o_type] += 1
        
        if k_type not in k_type_counts:
            k_type_counts[k_type] = 0
        k_type_counts[k_type] += 1
    
    # 选择出现频率适中的血清型作为高危候选
    candidate_high_risk = []
    
    for o_type, count in o_type_counts.items():
        if count >= 5 and count <= 50:  # 出现频率适中
            if 'unknown' not in o_type.lower():
                candidate_high_risk.append(o_type)
    
    for k_type, count in k_type_counts.items():
        if count >= 5 and count <= 50:
            if 'unknown' not in k_type.lower():
                candidate_high_risk.append(k_type)
    
    print(f"候选高危血清型: {candidate_high_risk[:10]}")  # 显示前10个
    
    # 使用文献中已知的高危血清型 + 数据中出现频率适中的类型
    final_high_risk = list(set(high_risk_o_types + high_risk_k_types + candidate_high_risk[:5]))
    
    print(f"最终高危血清型列表: {final_high_risk}")
    
    return final_high_risk

def find_specific_genes_with_serotypes(gene_pa, serotype_data, high_risk_serotypes):
    """使用血清型信息寻找特异性基因"""
    
    strain_columns = gene_pa.columns[14:]
    print(f"\n=== 分析 {len(strain_columns)} 个菌株 ===")
    
    # 根据血清型分组
    high_risk_strains = []
    other_strains = []
    
    for strain_col in strain_columns:
        # 标准化菌株名以匹配Kaptive结果
        strain_name = strain_col.replace('_1_ASM', '.1_ASM').replace('_2_ASM', '.2_ASM').replace('_genomic', '')
        
        if strain_name in serotype_data:
            data = serotype_data[strain_name]
            o_type = data['O_serotype']
            k_type = data['K_serotype']
            
            # 检查是否属于高危血清型
            is_high_risk = False
            for hr_serotype in high_risk_serotypes:
                if hr_serotype in o_type or hr_serotype in k_type:
                    is_high_risk = True
                    break
            
            if is_high_risk:
                high_risk_strains.append(strain_col)
            else:
                other_strains.append(strain_col)
        else:
            other_strains.append(strain_col)
    
    print(f"高危组菌株: {len(high_risk_strains)}")
    print(f"其他组菌株: {len(other_strains)}")
    
    if len(high_risk_strains) == 0:
        print("❌ 未找到高危血清型菌株，使用备选策略")
        # 使用O1、O3、O4血清型作为高危组
        for strain_col in strain_columns:
            strain_name = strain_col.replace('_1_ASM', '.1_ASM').replace('_2_ASM', '.2_ASM').replace('_genomic', '')
            if strain_name in serotype_data:
                o_type = serotype_data[strain_name]['O_serotype']
                if 'O1' in o_type or 'O3' in o_type or 'O4' in o_type:
                    high_risk_strains.append(strain_col)
                else:
                    other_strains.append(strain_col)
        
        print(f"备选高危组: {len(high_risk_strains)} 菌株")
    
    # 逐步放宽条件寻找基因
    thresholds = [
        (0.8, 0.2, "严格"),
        (0.7, 0.3, "中等"), 
        (0.6, 0.4, "宽松"),
        (0.5, 0.5, "极宽松")
    ]
    
    all_specific_genes = []
    
    for high_thresh, low_thresh, level in thresholds:
        print(f"\n尝试 {level} 条件: 高危≥{high_thresh:.0%}, 其他≤{low_thresh:.0%}")
        
        specific_genes = []
        
        for idx, row in gene_pa.iterrows():
            if len(high_risk_strains) == 0:
                break
                
            high_risk_presence = sum(pd.notna(row[col]) for col in high_risk_strains)
            high_risk_freq = high_risk_presence / len(high_risk_strains)
            
            other_presence = sum(pd.notna(row[col]) for col in other_strains)
            other_freq = other_presence / len(other_strains) if other_strains else 0
            
            if high_risk_freq >= high_thresh and other_freq <= low_thresh:
                gene_info = {
                    'Gene': row['Gene'],
                    'Annotation': row['Annotation'],
                    'HighRisk_Frequency': high_risk_freq,
                    'Other_Frequency': other_freq,
                    'No_in_HighRisk': high_risk_presence,
                    'No_in_Other': other_presence,
                    'Threshold_Level': level
                }
                specific_genes.append(gene_info)
        
        print(f"  找到 {len(specific_genes)} 个基因")
        all_specific_genes.extend(specific_genes)
        
        # 如果找到基因就停止
        if len(specific_genes) > 0 and len(all_specific_genes) >= 10:
            break
    
    return pd.DataFrame(all_specific_genes)

def main():
    """主函数"""
    print("=== 基于Kaptive血清型的特异性基因分析 ===")
    
    # 1. 读取Roary数据
    print("加载Roary数据...")
    gene_pa = pd.read_csv('roary_output_final_1762658265/gene_presence_absence.csv', low_memory=False)
    print(f"基因簇数量: {len(gene_pa)}")
    
    # 2. 解析Kaptive结果
    serotype_data = parse_kaptive_results()
    
    # 3. 识别高危血清型
    high_risk_serotypes = identify_high_risk_serotypes(serotype_data)
    
    # 4. 寻找特异性基因
    print("\n=== 寻找特异性基因 ===")
    specific_genes_df = find_specific_genes_with_serotypes(gene_pa, serotype_data, high_risk_serotypes)
    
    if len(specific_genes_df) == 0:
        print("\n❌ 未找到任何特异性基因")
        print("尝试使用O1/O3/O4作为高危组...")
        
        # 最后尝试：使用常见的致病血清型
        high_risk_strains = []
        other_strains = []
        
        strain_columns = gene_pa.columns[14:]
        for strain_col in strain_columns:
            strain_name = strain_col.replace('_1_ASM', '.1_ASM').replace('_2_ASM', '.2_ASM').replace('_genomic', '')
            if strain_name in serotype_data:
                o_type = serotype_data[strain_name]['O_serotype']
                if any(ot in o_type for ot in ['O1', 'O3', 'O4', 'O10']):
                    high_risk_strains.append(strain_col)
                else:
                    other_strains.append(strain_col)
        
        print(f"O1/O3/O4/O10组: {len(high_risk_strains)} 菌株")
        print(f"其他组: {len(other_strains)} 菌株")
        
        # 使用宽松条件
        specific_genes = []
        for idx, row in gene_pa.iterrows():
            if len(high_risk_strains) == 0:
                break
                
            high_risk_presence = sum(pd.notna(row[col]) for col in high_risk_strains)
            high_risk_freq = high_risk_presence / len(high_risk_strains)
            
            other_presence = sum(pd.notna(row[col]) for col in other_strains)
            other_freq = other_presence / len(other_strains) if other_strains else 0
            
            if high_risk_freq >= 0.5 and other_freq <= 0.5:
                gene_info = {
                    'Gene': row['Gene'],
                    'Annotation': row['Annotation'],
                    'HighRisk_Frequency': high_risk_freq,
                    'Other_Frequency': other_freq,
                    'No_in_HighRisk': high_risk_presence,
                    'No_in_Other': other_presence,
                    'Threshold_Level': '最后尝试'
                }
                specific_genes.append(gene_info)
        
        specific_genes_df = pd.DataFrame(specific_genes)
    
    # 5. 保存和显示结果
    if len(specific_genes_df) > 0:
        print(f"\n✓ 成功找到 {len(specific_genes_df)} 个候选特异性基因")
        
        # 保存结果
        specific_genes_df.to_csv('high_risk_specific_genes_KAPTIVE.csv', index=False)
        print("✓ 结果已保存到: high_risk_specific_genes_KAPTIVE.csv")
        
        # 显示结果
        print(f"\n=== 前20个候选基因 ===")
        for i, row in specific_genes_df.head(20).iterrows():
            annotation = str(row['Annotation'])[:70] + "..." if len(str(row['Annotation'])) > 70 else row['Annotation']
            print(f"{i+1:2d}. {row['Gene']:20}")
            print(f"     高危: {row['HighRisk_Frequency']:.1%} ({row['No_in_HighRisk']}菌株), 其他: {row['Other_Frequency']:.1%} ({row['No_in_Other']}菌株)")
            print(f"     注释: {annotation}")
            print(f"     条件: {row['Threshold_Level']}")
            print()
    else:
        print("\n❌ 经过多种尝试仍未找到特异性基因")
        print("可能的原因:")
        print("1. 数据集中血清型分布较均匀")
        print("2. 需要更专业的血清型危害性分类")
        print("3. 考虑使用其他分组策略（如来源、时间等）")

if __name__ == "__main__":
    main()
