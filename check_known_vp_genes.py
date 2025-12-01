import pandas as pd

def check_known_vp_virulence_genes():
    """检查已知的副溶血弧菌毒力基因"""
    print("=== 已知副溶血弧菌毒力基因检查 ===")
    
    # 已知的副溶血弧菌关键毒力基因
    known_vp_genes = [
        'tdh',  # 热稳定直接溶血素
        'trh',  # TDH相关溶血素
        'tlh',  # 热不稳定溶血素
        'T3SS1', 'T3SS2',  # III型分泌系统
        'vopC', 'vopQ', 'vopR', 'vopS',  # T3SS效应蛋白
        'vvhA',  # 溶血素
        'rtxA',  # 重复毒素
        'hlyA',  # 溶血素A
        'ureR', 'ureC'  # 脲酶相关
    ]
    
    gene_pa = pd.read_csv('roary_output_final_1762658265/gene_presence_absence.csv', low_memory=False)
    
    found_genes = []
    
    for known_gene in known_vp_genes:
        # 在基因名和注释中搜索
        name_matches = gene_pa[gene_pa['Gene'].str.contains(known_gene, case=False, na=False)]
        annotation_matches = gene_pa[gene_pa['Annotation'].str.contains(known_gene, case=False, na=False)]
        
        all_matches = pd.concat([name_matches, annotation_matches]).drop_duplicates()
        
        if len(all_matches) > 0:
            for idx, row in all_matches.iterrows():
                strain_columns = gene_pa.columns[14:]
                presence_count = sum(pd.notna(row[col]) for col in strain_columns)
                presence_freq = presence_count / len(strain_columns)
                
                found_genes.append({
                    'Known_Gene': known_gene,
                    'Matched_Gene': row['Gene'],
                    'Annotation': row['Annotation'],
                    'Presence_Frequency': presence_freq,
                    'Presence_Count': presence_count
                })
    
    if found_genes:
        found_df = pd.DataFrame(found_genes)
        found_df = found_df.sort_values('Presence_Frequency', ascending=False)
        
        print(f"找到 {len(found_df)} 个已知毒力基因的匹配项")
        print(f"\n=== 已知毒力基因分布 ===")
        for idx, row in found_df.iterrows():
            print(f"{row['Known_Gene']} -> {row['Matched_Gene']}")
            print(f"  注释: {row['Annotation']}")
            print(f"  分布: {row['Presence_Frequency']:.1%} ({row['Presence_Count']}/216)\n")
        
        found_df.to_csv('known_vp_virulence_genes.csv', index=False)
    else:
        print("未找到已知副溶血弧菌毒力基因的精确匹配")

if __name__ == "__main__":
    check_known_vp_virulence_genes()
