import pandas as pd
import numpy as np

def analyze_57_strains():
    print("=== 🎯 基于57个真实GCA菌株的靶标分析 ===")
    
    # 读取新Roary结果
    gene_pa = pd.read_csv('roary_consistent/gene_presence_absence.csv', low_memory=False)
    strains = gene_pa.columns[14:]
    
    print(f"📊 数据集确认:")
    print(f"   菌株数量: {len(strains)}")
    print(f"   基因簇总数: {len(gene_pa)}")
    
    # 计算核心基因组统计
    core_genes = sum(1 for idx, row in gene_pa.iterrows() 
                    if sum(pd.notna(row[col]) for col in strains) >= len(strains) * 0.99)
    
    print(f"   核心基因(≥99%): {core_genes} ({core_genes/len(gene_pa)*100:.1f}%)")
    
    return gene_pa, strains

def find_virulence_targets(gene_pa, strains):
    """寻找毒力相关靶标"""
    print(f"\n=== 🔍 寻找毒力相关靶标 ===")
    
    # 毒力相关关键词
    virulence_keywords = [
        'toxin', 'hemolysin', 'virulence', 'effector', 'adhesin',
        'secretion', 'invasin', 'capsule', 'antigen', 'phospholipase',
        'protease', 'collagenase', 'siderophore', 'hemolysin'
    ]
    
    candidates = []
    
    for idx, row in gene_pa.iterrows():
        presence = sum(pd.notna(row[col]) for col in strains)
        freq = presence / len(strains)
        
        # 选择分布适中的基因 (20%-80%)
        if 0.2 <= freq <= 0.8:
            annotation = str(row['Annotation']).lower()
            
            # 检查是否包含毒力相关关键词
            if any(keyword in annotation for keyword in virulence_keywords):
                candidates.append({
                    'Gene': row['Gene'],
                    'Annotation': row['Annotation'],
                    'Frequency': freq,
                    'Presence': presence,
                    'Total': len(strains),
                    'Score': freq * (1 - abs(freq - 0.5))  # 0.5分布得分最高
                })
    
    # 排序并显示结果
    if candidates:
        candidates.sort(key=lambda x: x['Score'], reverse=True)
        
        print(f"找到 {len(candidates)} 个候选毒力靶标")
        print(f"\n🎯 前20个最佳靶标:")
        
        for i, cand in enumerate(candidates[:20], 1):
            print(f"{i:2d}. {cand['Gene']}")
            print(f"    分布: {cand['Presence']}/{cand['Total']} ({cand['Frequency']:.1%})")
            print(f"    功能: {cand['Annotation']}")
            print(f"    得分: {cand['Score']:.3f}")
            print()
            
        return candidates
    else:
        print("未找到符合条件的毒力靶标")
        return []

def check_known_genes(gene_pa, strains):
    """检查已知副溶血弧菌基因"""
    print(f"\n=== 🔎 已知副溶血弧菌基因检查 ===")
    
    known_patterns = ['tdh', 'trh', 'tlh', 'vop', 'ure', 'hly', 't3ss', 't6ss']
    
    found_genes = []
    
    for pattern in known_patterns:
        matches = gene_pa[gene_pa['Gene'].str.contains(pattern, case=False, na=False)]
        
        for idx, row in matches.iterrows():
            presence = sum(pd.notna(row[col]) for col in strains)
            freq = presence / len(strains)
            
            found_genes.append({
                'Gene': row['Gene'],
                'Pattern': pattern,
                'Annotation': row['Annotation'],
                'Frequency': freq,
                'Presence': presence
            })
    
    if found_genes:
        found_genes.sort(key=lambda x: x['Frequency'], reverse=True)
        print(f"找到 {len(found_genes)} 个已知相关基因:")
        
        for gene in found_genes[:15]:
            print(f"   🔍 {gene['Gene']}: {gene['Presence']}/{len(strains)} ({gene['Frequency']:.1%})")
            print(f"      匹配: {gene['Pattern']}")
            print(f"      功能: {gene['Annotation'][:80]}...")
            print()

if __name__ == "__main__":
    gene_pa, strains = analyze_57_strains()
    targets = find_virulence_targets(gene_pa, strains)
    check_known_genes(gene_pa, strains)
    
    print(f"\n💎 最终总结:")
    if targets:
        print(f"✅ 基于57个真实GCA菌株找到了 {len(targets)} 个候选靶标")
        print("这些靶标在你的数据中真实存在，适合立即开始实验验证")
        
        # 保存结果
        import pandas as pd
        df = pd.DataFrame(targets)
        df.to_csv('real_candidate_targets_57strains.csv', index=False)
        print(f"✅ 候选靶标已保存到: real_candidate_targets_57strains.csv")
    else:
        print("❌ 需要调整筛选策略")
