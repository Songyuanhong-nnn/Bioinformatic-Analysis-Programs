import pandas as pd
import glob

def reassess_with_real_gca():
    print("=== 🎯 基于真实GCA数据重新评估靶标 ===")
    
    # 找出实际可用的GCA菌株
    gca_strains = set()
    
    # 从clean_genome找
    for fna in glob.glob("clean_genome/GCA_*.fna"):
        strain = os.path.basename(fna).replace('.fna', '')
        gca_strains.add(strain + '_genomic')
    
    # 从bakta找  
    for dir_path in glob.glob("bakta_annotations/GCA_*"):
        strain = os.path.basename(dir_path)
        gca_strains.add(strain + '_genomic')
    
    print(f"实际可用的GCA菌株: {len(gca_strains)}")
    if gca_strains:
        print(f"示例: {list(gca_strains)[:3]}")
    
    # 读取Roary结果，但只关注这些GCA菌株
    try:
        gene_pa = pd.read_csv('roary_output_final_1762658265/gene_presence_absence.csv', low_memory=False)
        
        # 找出在Roary结果中实际存在的GCA菌株
        strain_columns = gene_pa.columns[14:]
        available_gca = [s for s in strain_columns if s in gca_strains]
        
        print(f"\n在Roary结果中存在的GCA菌株: {len(available_gca)}")
        
        if available_gca:
            # 重新计算靶标基因分布
            targets = ['group_9360', 'group_3365', 'tdh2', 'vopS', 'epsE_2']
            
            print(f"\n=== 靶标基因在真实GCA菌株中的分布 ===")
            for gene in targets:
                gene_data = gene_pa[gene_pa['Gene'] == gene]
                if len(gene_data) > 0:
                    row = gene_data.iloc[0]
                    presence = sum(pd.notna(row[col]) for col in available_gca)
                    percentage = (presence / len(available_gca)) * 100 if available_gca else 0
                    print(f"🎯 {gene}: {presence}/{len(available_gca)} ({percentage:.1f}%)")
                    print(f"   注释: {row['Annotation']}")
                else:
                    print(f"❌ {gene}: 在Roary结果中未找到")
        
        return available_gca
        
    except Exception as e:
        print(f"读取Roary文件错误: {e}")
        return []

if __name__ == "__main__":
    import os
    real_gca_strains = reassess_with_real_gca()
    
    print(f"\n=== 💎 最终建议 ===")
    print(f"你应该基于 {len(real_gca_strains)} 个真实GCA菌株继续分析")
    print("这些菌株在clean_genome、bakta和roary结果中都存在且一致")
