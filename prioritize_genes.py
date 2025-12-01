import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
from collections import defaultdict

def prioritize_candidate_genes(specific_genes_df, high_risk_serotypes):
    """
    对候选基因进行多维度评分和优先排序
    """
    scoring_data = []
    
    for idx, gene_row in specific_genes_df.iterrows():
        score = 0
        scores_detail = {}
        
        # 1. 特异性评分 (40%)
        specificity_score = (gene_row['HighRisk_Frequency'] - gene_row['Other_Frequency']) * 40
        score += specificity_score
        scores_detail['specificity'] = specificity_score
        
        # 2. 功能注释评分 (30%)
        annotation = str(gene_row['Annotation']).lower()
        annotation_score = 0
        
        # 高价值功能标签
        high_value_keywords = {
            'membrane': 8, 'surface': 8, 'secreted': 8, 'toxin': 10, 
            'hemolysin': 10, 'adhesin': 9, 'pilin': 7, 'porin': 7,
            'capsule': 9, 'lipopolysaccharide': 9, 'o-antigen': 10,
            'k-antigen': 10, 'virulence': 8, 'effector': 8,
            'antigen': 8, 'export': 6, 'transport': 6
        }
        
        for keyword, points in high_value_keywords.items():
            if keyword in annotation:
                annotation_score += points
        
        annotation_score = min(annotation_score, 30)  # 上限30分
        score += annotation_score
        scores_detail['annotation'] = annotation_score
        
        # 3. 保守性评分 (20%) - 在高危组内保守
        conservation_score = gene_row['HighRisk_Frequency'] * 20
        score += conservation_score
        scores_detail['conservation'] = conservation_score
        
        # 4. 假设蛋白惩罚 (-10分如果是假设蛋白)
        hypothetical_penalty = -10 if 'hypothetical' in annotation else 0
        score += hypothetical_penalty
        scores_detail['hypothetical_penalty'] = hypothetical_penalty
        
        # 5. 基因长度评分 (10%)
        gene_length = estimate_gene_length(gene_row)
        length_score = min(gene_length / 100, 10)  # 每100bp得1分，最高10分
        score += length_score
        scores_detail['length'] = length_score
        
        scoring_data.append({
            'Gene': gene_row['Gene'],
            'Annotation': gene_row['Annotation'],
            'HighRisk_Freq': gene_row['HighRisk_Frequency'],
            'Other_Freq': gene_row['Other_Frequency'],
            'Total_Score': score,
            **scores_detail
        })
    
    # 创建评分DataFrame并排序
    priority_df = pd.DataFrame(scoring_data)
    priority_df = priority_df.sort_values('Total_Score', ascending=False)
    priority_df['Priority_Rank'] = range(1, len(priority_df) + 1)
    
    return priority_df

def estimate_gene_length(gene_row):
    """
    估计基因长度 - 简化版本
    """
    # 可以根据注释信息中的长度信息来估计
    # 这里使用一个简单的估计方法
    annotation = str(gene_row['Annotation']).lower()
    
    if any(word in annotation for word in ['toxin', 'hemolysin']):
        return 800  # 毒素通常较短
    elif any(word in annotation for word in ['membrane', 'surface']):
        return 1200  # 膜蛋白中等长度
    elif any(word in annotation for word in ['transferase', 'synthesis']):
        return 1500  # 合成酶较长
    else:
        return 1000  # 默认长度

def categorize_genes_by_function(priority_df):
    """
    按功能类别对基因进行分组
    """
    categories = {
        '表面抗原相关': ['membrane', 'surface', 'porin', 'pilin', 'adhesin', 'antigen'],
        '毒素与毒力因子': ['toxin', 'hemolysin', 'virulence', 'effector'],
        '荚膜合成相关': ['capsule', 'capsular', 'k-antigen', '多糖', 'polysaccharide'],
        'LPS合成相关': ['lipopolysaccharide', 'lps', 'o-antigen'],
        '分泌系统': ['secreted', 'secretion', 'transporter', 'export'],
        '代谢与适应性': ['metabolism', 'regulation', 'transferase', 'synthesis'],
        '假设蛋白': ['hypothetical', 'unknown', 'putative', 'unnamed']
    }
    
    categorized_genes = {}
    
    for category, keywords in categories.items():
        category_genes = []
        for idx, row in priority_df.iterrows():
            annotation = str(row['Annotation']).lower()
            if any(keyword in annotation for keyword in keywords):
                category_genes.append({
                    'Gene': row['Gene'],
                    'Annotation': row['Annotation'],
                    'Total_Score': row['Total_Score'],
                    'Priority_Rank': row['Priority_Rank']
                })
        categorized_genes[category] = category_genes
    
    return categorized_genes

def main():
    """
    主函数 - 执行完整的基因优先排序流程
    """
    print("=== 候选基因优先排序分析 ===")
    
    # 1. 首先需要读取你的特异性基因数据
    # 如果你已经有 specific_genes_df，可以直接使用
    # 如果没有，我们先创建一个示例数据来测试
    
    try:
        # 尝试读取已有的特异性基因文件
        specific_genes_df = pd.read_csv('high_risk_specific_genes.csv')
        print("✓ 读取 high_risk_specific_genes.csv 成功")
    except:
        print("⚠️ 无法读取 high_risk_specific_genes.csv，创建示例数据用于测试")
        # 创建示例数据
        specific_genes_df = pd.DataFrame({
            'Gene': [f'gene_{i}' for i in range(1, 51)],
            'Annotation': [
                'hypothetical protein',
                'membrane protein', 
                'toxin hemolysin',
                'capsular polysaccharide synthesis protein',
                'O-antigen transferase',
                'virulence factor',
                'surface antigen',
                'porin protein',
                'adhesin protein',
                'secretion system protein'
            ] * 5,
            'HighRisk_Frequency': np.random.uniform(0.8, 1.0, 50),
            'Other_Frequency': np.random.uniform(0.0, 0.2, 50)
        })
    
    print(f"待分析基因数量: {len(specific_genes_df)}")
    
    # 2. 定义高危血清型
    high_risk_serotypes = ['O3K6', 'O10K4', 'O4KUT']
    
    # 3. 执行优先排序
    print("执行多维度评分...")
    priority_genes = prioritize_candidate_genes(specific_genes_df, high_risk_serotypes)
    
    # 4. 功能分类
    print("进行功能分类...")
    gene_categories = categorize_genes_by_function(priority_genes)
    
    # 5. 显示结果
    print(f"\n=== 分析完成 ===")
    print(f"总基因数: {len(priority_genes)}")
    print(f"评分范围: {priority_genes['Total_Score'].min():.1f} - {priority_genes['Total_Score'].max():.1f}")
    
    print(f"\n=== Top 20 候选基因 ===")
    top_20 = priority_genes.head(20)[['Gene', 'Annotation', 'Total_Score', 'Priority_Rank']]
    for idx, row in top_20.iterrows():
        print(f"{row['Priority_Rank']:2d}. {row['Gene']:15} {row['Total_Score']:5.1f}分 - {row['Annotation']}")
    
    print(f"\n=== 功能类别分布 ===")
    for category, genes in gene_categories.items():
        print(f"  {category}: {len(genes)} 个基因")
        if genes and len(genes) <= 5:  # 如果基因数少，显示全部
            for gene in genes[:3]:
                print(f"    - {gene['Gene']} (分数: {gene['Total_Score']:.1f})")
    
    # 6. 保存结果
    priority_genes.to_csv('prioritized_candidate_genes.csv', index=False)
    print(f"\n✓ 结果已保存到: prioritized_candidate_genes.csv")
    
    # 7. 生成推荐列表
    print(f"\n=== 最终推荐靶标 ===")
    # 排除假设蛋白，选择每个类别的前2个
    recommended_genes = []
    for category, genes in gene_categories.items():
        if category not in ['假设蛋白', '代谢与适应性'] and genes:
            top_genes = sorted(genes, key=lambda x: x['Total_Score'], reverse=True)[:2]
            recommended_genes.extend(top_genes)
    
    # 按分数排序并去重
    seen_genes = set()
    final_recommendations = []
    for gene in sorted(recommended_genes, key=lambda x: x['Total_Score'], reverse=True):
        if gene['Gene'] not in seen_genes:
            final_recommendations.append(gene)
            seen_genes.add(gene['Gene'])
    
    print(f"推荐 {len(final_recommendations)} 个高质量候选靶标:")
    for i, gene in enumerate(final_recommendations[:15], 1):
        print(f"{i:2d}. {gene['Gene']:15} {gene['Total_Score']:5.1f}分 - {gene['Annotation']}")

if __name__ == "__main__":
    main()
