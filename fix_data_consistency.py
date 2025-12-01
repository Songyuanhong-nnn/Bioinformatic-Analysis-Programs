import pandas as pd
import os
import glob
import shutil

def analyze_naming_patterns():
    print("=== 🔍 分析命名模式 ===")
    
    # 检查各目录的命名模式
    print("1. clean_genome命名模式:")
    for fna in glob.glob("clean_genome/*.fna")[:3]:
        basename = os.path.basename(fna).replace('.fna', '')
        print(f"   {basename}")
    
    print("\n2. bakta_annotations命名模式:")
    for dir_path in glob.glob("bakta_annotations/*")[:3]:
        basename = os.path.basename(dir_path)
        print(f"   {basename}")
    
    print("\n3. roary_input命名模式:")
    for gff in glob.glob("roary_input/*.gff")[:3]:
        basename = os.path.basename(gff).replace('.gff', '')
        print(f"   {basename}")

def find_matching_strains():
    print("\n=== 🔄 寻找匹配的菌株 ===")
    
    # 收集所有命名变体
    strains = {}
    
    # clean_genome中的菌株 (GCA_000328405.1_ASM32840v1_genomic)
    for fna in glob.glob("clean_genome/GCA_*.fna"):
        name = os.path.basename(fna).replace('.fna', '')
        base_id = name.replace('.', '_')  # 转换为roary格式
        strains[base_id] = strains.get(base_id, {})
        strains[base_id]['clean_genome'] = name
    
    # bakta中的菌株 (GCA_000328405.1_ASM32840v1_genomic)
    for dir_path in glob.glob("bakta_annotations/GCA_*"):
        name = os.path.basename(dir_path)
        base_id = name.replace('.', '_')  # 转换为roary格式
        strains[base_id] = strains.get(base_id, {})
        strains[base_id]['bakta'] = name
    
    # roary_input中的菌株 (GCA_000328405_1_ASM32840v1_genomic)
    for gff in glob.glob("roary_input/GCA_*.gff"):
        name = os.path.basename(gff).replace('.gff', '')
        strains[name] = strains.get(name, {})
        strains[name]['roary_input'] = name
    
    # 统计匹配情况
    matched = 0
    total = 0
    
    print("菌株匹配情况:")
    for base_id, sources in list(strains.items())[:10]:
        has_clean = 'clean_genome' in sources
        has_bakta = 'bakta' in sources  
        has_roary = 'roary_input' in sources
        
        if has_clean and has_bakta and has_roary:
            matched += 1
            status = "✅ 完全匹配"
        else:
            status = "❌ 不匹配"
        
        print(f"   {base_id}: clean={has_clean}, bakta={has_bakta}, roary={has_roary} {status}")
        total += 1
    
    print(f"\n匹配统计: {matched}/{total} 个菌株完全匹配")
    
    return strains

def create_consistent_dataset(strains):
    print("\n=== 🎯 创建一致的数据集 ===")
    
    # 找出完全匹配的菌株
    consistent_strains = []
    
    for base_id, sources in strains.items():
        if 'clean_genome' in sources and 'bakta' in sources and 'roary_input' in sources:
            consistent_strains.append({
                'base_id': base_id,
                'clean_genome': sources['clean_genome'],
                'bakta': sources['bakta'], 
                'roary_input': sources['roary_input']
            })
    
    print(f"找到 {len(consistent_strains)} 个完全匹配的菌株")
    
    if consistent_strains:
        print("\n示例菌株:")
        for strain in consistent_strains[:5]:
            print(f"   {strain['base_id']}")
            print(f"     clean: {strain['clean_genome']}")
            print(f"     bakta: {strain['bakta']}")
            print(f"     roary: {strain['roary_input']}")
    
    return consistent_strains

def reassess_targets_with_consistent_data(consistent_strains):
    print("\n=== 📊 基于一致数据重新评估靶标 ===")
    
    if not consistent_strains:
        print("没有一致的数据可用于分析")
        return
    
    # 获取roary格式的菌株名
    roary_strain_names = [s['roary_input'] + '_genomic' for s in consistent_strains]
    
    print(f"使用 {len(roary_strain_names)} 个一致菌株重新分析")
    print(f"示例: {roary_strain_names[:3]}")
    
    # 读取Roary结果
    try:
        gene_pa = pd.read_csv('roary_output_final_1762658265/gene_presence_absence.csv', low_memory=False)
        
        # 重新计算靶标基因分布
        targets = ['group_9360', 'group_3365', 'tdh2', 'vopS', 'epsE_2']
        
        print(f"\n靶标基因在一致数据集中的分布:")
        for gene in targets:
            gene_data = gene_pa[gene_pa['Gene'] == gene]
            if len(gene_data) > 0:
                row = gene_data.iloc[0]
                
                # 计算在一致菌株中的分布
                presence = 0
                for strain in roary_strain_names:
                    if strain in gene_pa.columns and pd.notna(row[strain]):
                        presence += 1
                
                percentage = (presence / len(roary_strain_names)) * 100
                print(f"🎯 {gene}: {presence}/{len(roary_strain_names)} ({percentage:.1f}%)")
                print(f"   注释: {row['Annotation']}")
            else:
                print(f"❌ {gene}: 在Roary结果中未找到")
                
    except Exception as e:
        print(f"读取Roary文件错误: {e}")

if __name__ == "__main__":
    analyze_naming_patterns()
    strains = find_matching_strains()
    consistent_strains = create_consistent_dataset(strains)
    reassess_targets_with_consistent_data(consistent_strains)
    
    print(f"\n=== 💡 最终建议 ===")
    if consistent_strains:
        print(f"基于 {len(consistent_strains)} 个完全匹配的菌株继续分析")
        print("这些菌株在clean_genome、bakta和roary_input中都存在且标识一致")
    else:
        print("需要修复数据标识不一致的问题")
