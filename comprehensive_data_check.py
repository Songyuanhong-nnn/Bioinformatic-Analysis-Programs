import pandas as pd
import os
import glob

def check_all_data_sources():
    print("=== 🔍 全面数据状态检测 ===")
    
    # 1. 检测原始基因组文件
    print("\n1. 📁 原始基因组文件:")
    gca_fna = glob.glob("clean_genome/GCA_*.fna")
    gcf_fna = glob.glob("clean_genome/GCF_*.fna") 
    all_fna = glob.glob("clean_genome/*.fna")
    
    print(f"   GCA文件: {len(gca_fna)}")
    print(f"   GCF文件: {len(gcf_fna)}")
    print(f"   总文件: {len(all_fna)}")
    
    # 2. 检测Bakta注释结果
    print("\n2. 🔬 Bakta注释结果:")
    bakta_gca = glob.glob("bakta_annotations/GCA_*")
    bakta_gcf = glob.glob("bakta_annotations/GCF_*")
    bakta_all = glob.glob("bakta_annotations/*/")
    
    print(f"   GCA注释: {len(bakta_gca)}")
    print(f"   GCF注释: {len(bakta_gcf)}")
    print(f"   总注释: {len(bakta_all)}")
    
    # 3. 检测Prokka注释结果
    print("\n3. 🧬 Prokka注释结果:")
    prokka_gca = glob.glob("prokka_annotations/GCA_*")
    prokka_gcf = glob.glob("prokka_annotations/GCF_*") 
    prokka_all = glob.glob("prokka_annotations/*/")
    
    print(f"   GCA注释: {len(prokka_gca)}")
    print(f"   GCF注释: {len(prokka_gcf)}")
    print(f"   总注释: {len(prokka_all)}")
    
    # 4. 检测Roary输入GFF文件
    print("\n4. 📊 Roary输入GFF文件:")
    roary_gca = glob.glob("roary_input/GCA_*.gff")
    roary_gcf = glob.glob("roary_input/GCF_*.gff")
    roary_all = glob.glob("roary_input/*.gff")
    
    print(f"   GCA GFF: {len(roary_gca)}")
    print(f"   GCF GFF: {len(roary_gcf)}")
    print(f"   总GFF: {len(roary_all)}")
    
    # 5. 检测Roary输出中的菌株
    print("\n5. 🎯 Roary输出菌株分析:")
    try:
        gene_pa = pd.read_csv('roary_output_final_1762658265/gene_presence_absence.csv', low_memory=False)
        strain_columns = gene_pa.columns[14:]
        
        gca_strains = [s for s in strain_columns if 'GCA_' in s]
        gcf_strains = [s for s in strain_columns if 'GCF_' in s]
        other_strains = [s for s in strain_columns if 'GCA_' not in s and 'GCF_' not in s]
        
        print(f"   GCA菌株: {len(gca_strains)}")
        print(f"   GCF菌株: {len(gcf_strains)}")
        print(f"   其他菌株: {len(other_strains)}")
        print(f"   总菌株: {len(strain_columns)}")
        
        # 显示前5个菌株示例
        print(f"   菌株示例: {strain_columns[:3].tolist()}")
        
    except Exception as e:
        print(f"   读取Roary文件错误: {e}")
    
    return {
        'clean_genome': {'gca': len(gca_fna), 'gcf': len(gcf_fna), 'total': len(all_fna)},
        'bakta': {'gca': len(bakta_gca), 'gcf': len(bakta_gcf), 'total': len(bakta_all)},
        'prokka': {'gca': len(prokka_gca), 'gcf': len(prokka_gcf), 'total': len(prokka_all)},
        'roary_input': {'gca': len(roary_gca), 'gcf': len(roary_gcf), 'total': len(roary_all)},
        'roary_output': {'gca': len(gca_strains), 'gcf': len(gcf_strains), 'total': len(strain_columns)}
    }

def check_data_consistency(results):
    print("\n=== 🔄 数据一致性分析 ===")
    
    # 检查各阶段数据量是否匹配
    stages = ['clean_genome', 'bakta', 'prokka', 'roary_input', 'roary_output']
    
    print("各阶段GCA数据量:")
    for stage in stages:
        if stage in results:
            print(f"   {stage}: {results[stage]['gca']}")
    
    print(f"\n一致性状态:")
    gca_counts = [results[stage]['gca'] for stage in stages if stage in results]
    
    if len(set(gca_counts)) == 1:
        print("   ✅ 所有阶段GCA数据量一致")
    else:
        print("   ⚠️  GCA数据量不一致")
        for stage in stages:
            if stage in results:
                print(f"      {stage}: {results[stage]['gca']}")

def find_actual_working_set():
    print("\n=== 🎯 确定实际工作数据集 ===")
    
    # 找出真正可用的GCA数据集
    sources = {
        'clean_genome': set([os.path.basename(f).replace('.fna', '') for f in glob.glob("clean_genome/GCA_*.fna")]),
        'bakta': set([os.path.basename(d) for d in glob.glob("bakta_annotations/GCA_*")]),
        'prokka': set([os.path.basename(d) for d in glob.glob("prokka_annotations/GCA_*")]),
        'roary_input': set([os.path.basename(f).replace('.gff', '') for f in glob.glob("roary_input/GCA_*.gff")])
    }
    
    print("各源GCA标识数量:")
    for source, ids in sources.items():
        print(f"   {source}: {len(ids)}")
    
    # 找出所有源中都存在的GCA
    common_gca = set.intersection(*[ids for ids in sources.values() if ids])
    
    print(f"\n共同存在的GCA数量: {len(common_gca)}")
    if common_gca:
        print(f"示例: {list(common_gca)[:3]}")
    
    return common_gca

if __name__ == "__main__":
    results = check_all_data_sources()
    check_data_consistency(results)
    common_gca = find_actual_working_set()
    
    print(f"\n=== 💡 建议 ===")
    print(f"基于检测结果，你的实际工作数据集包含 {len(common_gca)} 个GCA基因组")
    print("下一步应该基于这个一致的数据集重新分析")
