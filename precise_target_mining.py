import pandas as pd
import numpy as np
import re

def find_true_virulence_factors():
    """寻找真正的毒力因子"""
    print("=== 精准毒力因子挖掘 ===")
    
    # 读取Roary数据
    gene_pa = pd.read_csv('roary_output_final_1762658265/gene_presence_absence.csv', low_memory=False)
    strain_columns = gene_pa.columns[14:]
    
    # 精准毒力关键词 - 只包括真正的毒力相关功能
    true_virulence_keywords = [
        'hemolysin', 'cytotoxin', 'enterotoxin', 'leukocidin', 
        'adhesin', 'invasin', 'collagenase', 'protease', 'phospholipase',
        'siderophore', 'capsule', 'secretion system', 'type iii', 'type iv',
        'toxin', 'effector', 'virulence factor', 'pathogenicity'
    ]
    
    # 排除的一般代谢基因
    exclude_keywords = [
        'transferase', 'methyltransferase', 'phosphatase', 'kinase',
        'dehydrogenase', 'synthase', 'polymerase', 'reductase',
        'transporter', 'receptor', 'ribosomal', 'rna', 'trna'
    ]
    
    print("筛选真正的毒力因子...")
    true_virulence_genes = []
    
    for idx, row in gene_pa.iterrows():
        annotation = str(row['Annotation']).lower()
        
        # 检查是否包含毒力关键词
        has_virulence = any(keyword in annotation for keyword in true_virulence_keywords)
        
        # 检查是否不包含排除关键词
        not_excluded = not any(exclude in annotation for exclude in exclude_keywords)
        
        if has_virulence and not_excluded:
            # 计算分布
            presence_count = sum(pd.notna(row[col]) for col in strain_columns)
            presence_freq = presence_count / len(strain_columns)
            
            # 只选择分布频率在10%-90%的基因（排除核心和罕见基因）
            if 0.1 <= presence_freq <= 0.9:
                true_virulence_genes.append({
                    'Gene': row['Gene'],
                    'Annotation': row['Annotation'],
                    'Presence_Frequency': presence_freq,
                    'Presence_Count': presence_count,
                    'Total_Strains': len(strain_columns),
                    'Category': 'True_Virulence'
                })
    
    print(f"找到 {len(true_virulence_genes)} 个真正的毒力因子候选")
    return pd.DataFrame(true_virulence_genes)

def analyze_hemolysin_genes():
    """专门分析溶血素基因"""
    print("\n=== 溶血素基因专项分析 ===")
    
    gene_pa = pd.read_csv('roary_output_final_1762658265/gene_presence_absence.csv', low_memory=False)
    strain_columns = gene_pa.columns[14:]
    
    hemolysin_genes = []
    
    for idx, row in gene_pa.iterrows():
        annotation = str(row['Annotation']).lower()
        
        # 精准匹配溶血素相关基因
        if any(keyword in annotation for keyword in ['hemolysin', 'hly', 'thl', 'tly']):
            presence_count = sum(pd.notna(row[col]) for col in strain_columns)
            presence_freq = presence_count / len(strain_columns)
            
            hemolysin_genes.append({
                'Gene': row['Gene'],
                'Annotation': row['Annotation'],
                'Presence_Frequency': presence_freq,
                'Presence_Count': presence_count,
                'Total_Strains': len(strain_columns),
                'Category': 'Hemolysin'
            })
    
    print(f"找到 {len(hemolysin_genes)} 个溶血素相关基因")
    return pd.DataFrame(hemolysin_genes)

def analyze_secretion_system_genes():
    """分析分泌系统基因"""
    print("\n=== 分泌系统基因分析 ===")
    
    gene_pa = pd.read_csv('roary_output_final_1762658265/gene_presence_absence.csv', low_memory=False)
    strain_columns = gene_pa.columns[14:]
    
    secretion_genes = []
    
    for idx, row in gene_pa.iterrows():
        annotation = str(row['Annotation']).lower()
        
        # 分泌系统相关基因
        if any(keyword in annotation for keyword in [
            'type iii', 'type iv', 'type vi', 't3ss', 't4ss', 't6ss',
            'secretion system', 'effector', 'injectisome'
        ]):
            presence_count = sum(pd.notna(row[col]) for col in strain_columns)
            presence_freq = presence_count / len(strain_columns)
            
            secretion_genes.append({
                'Gene': row['Gene'],
                'Annotation': row['Annotation'],
                'Presence_Frequency': presence_freq,
                'Presence_Count': presence_count,
                'Total_Strains': len(strain_columns),
                'Category': 'Secretion_System'
            })
    
    print(f"找到 {len(secretion_genes)} 个分泌系统相关基因")
    return pd.DataFrame(secretion_genes)

def find_high_value_hypotheticals():
    """寻找高价值的假设蛋白"""
    print("\n=== 高价值假设蛋白分析 ===")
    
    gene_pa = pd.read_csv('roary_output_final_1762658265/gene_presence_absence.csv', low_memory=False)
    strain_columns = gene_pa.columns[14:]
    
    # 读取高变异性基因
    try:
        variable_genes = pd.read_csv('variable_distribution_genes.csv')
        hypothetical_variable = variable_genes[
            variable_genes['Annotation'].str.contains('hypothetical', case=False, na=False)
        ]
        print(f"从高变异性基因中找到 {len(hypothetical_variable)} 个假设蛋白")
        return hypothetical_variable
    except:
        print("无法读取高变异性基因文件")
        return pd.DataFrame()

def prioritize_and_rank_targets(*dataframes):
    """对所有候选靶标进行优先排序"""
    print("\n=== 靶标优先排序 ===")
    
    # 合并所有数据框
    all_candidates = pd.concat(dataframes, ignore_index=True)
    
    if len(all_candidates) == 0:
        print("没有候选靶标可排序")
        return pd.DataFrame()
    
    # 评分系统
    scored_candidates = []
    
    for idx, row in all_candidates.iterrows():
        score = 0
        
        # 1. 功能重要性 (40分)
        category_scores = {
            'Hemolysin': 40,
            'Secretion_System': 35,
            'True_Virulence': 30,
            'Hypothetical': 20
        }
        score += category_scores.get(row.get('Category', 'Unknown'), 15)
        
        # 2. 分布特异性 (30分) - 中等频率得分最高
        freq = row['Presence_Frequency']
        if 0.3 <= freq <= 0.7:
            specificity_score = 30 * (1 - abs(freq - 0.5) / 0.2)  # 0.5频率得满分
        else:
            specificity_score = 15
        score += specificity_score
        
        # 3. 注释质量 (20分) - 非假设蛋白得分高
        annotation = str(row['Annotation']).lower()
        if 'hypothetical' not in annotation:
            score += 20
        
        # 4. 基因长度暗示 (10分) - 通过基因名简单判断
        gene_name = str(row['Gene'])
        if any(pattern in gene_name for pattern in ['group_', 'cluster_']):
            score += 5  # 可能是较大的基因簇
        
        scored_candidates.append({
            **row.to_dict(),
            'Score': score
        })
    
    # 创建评分DataFrame并排序
    scored_df = pd.DataFrame(scored_candidates)
    scored_df = scored_df.sort_values('Score', ascending=False)
    scored_df['Rank'] = range(1, len(scored_df) + 1)
    
    return scored_df

def main():
    """主函数"""
    print("=== 精准靶标挖掘分析 ===")
    
    # 1. 各种专项分析
    true_virulence_df = find_true_virulence_factors()
    hemolysin_df = analyze_hemolysin_genes()
    secretion_df = analyze_secretion_system_genes()
    hypothetical_df = find_high_value_hypotheticals()
    
    # 2. 优先排序
    all_candidates = [df for df in [true_virulence_df, hemolysin_df, secretion_df, hypothetical_df] if len(df) > 0]
    
    if all_candidates:
        ranked_targets = prioritize_and_rank_targets(*all_candidates)
        
        if len(ranked_targets) > 0:
            print(f"\n✓ 总共评估 {len(ranked_targets)} 个候选靶标")
            
            # 保存结果
            ranked_targets.to_csv('ranked_candidate_targets.csv', index=False)
            print("✓ 排序结果已保存到: ranked_candidate_targets.csv")
            
            # 显示前30个最佳靶标
            print(f"\n=== 前30个最佳候选靶标 ===")
            for i, row in ranked_targets.head(30).iterrows():
                annotation = str(row['Annotation'])[:80] + "..." if len(str(row['Annotation'])) > 80 else row['Annotation']
                print(f"{row['Rank']:2d}. [{row['Category']}] 分数: {row['Score']:.1f}")
                print(f"    基因: {row['Gene']}")
                print(f"    分布: {row['Presence_Frequency']:.1%} ({row['Presence_Count']}/{row['Total_Strains']})")
                print(f"    注释: {annotation}\n")
            
            # 按类别统计
            print(f"\n=== 按类别统计 ===")
            category_stats = ranked_targets.groupby('Category').size()
            for category, count in category_stats.items():
                top_in_category = ranked_targets[ranked_targets['Category'] == category].head(3)
                print(f"{category}: {count}个基因")
                for _, gene in top_in_category.iterrows():
                    print(f"  - {gene['Gene']} (分数: {gene['Score']:.1f})")
        else:
            print("❌ 没有合格的候选靶标")
    else:
        print("❌ 所有专项分析均未找到候选靶标")
    
    print(f"\n=== 分析建议 ===")
    print("1. 优先验证前10个高分靶标")
    print("2. 重点关注溶血素和分泌系统相关基因")
    print("3. 对高价值假设蛋白进行功能预测")
    print("4. 设计实验验证这些靶标的特异性")

if __name__ == "__main__":
    main()
