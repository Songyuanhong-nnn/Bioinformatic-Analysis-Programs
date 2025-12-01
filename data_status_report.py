import pandas as pd
import glob
import os

def generate_final_report():
    print("=== 📊 最终数据状态报告 ===")
    
    # 收集所有数据源信息
    sources = {
        'clean_genome': set([os.path.basename(f).replace('.fna', '') for f in glob.glob("clean_genome/GCA_*.fna")]),
        'bakta_annotations': set([os.path.basename(d) for d in glob.glob("bakta_annotations/GCA_*")]),
        'prokka_annotations': set([os.path.basename(d) for d in glob.glob("prokka_annotations/GCA_*")]),
        'roary_input_gff': set([os.path.basename(f).replace('.gff', '') for f in glob.glob("roary_input/GCA_*.gff")])
    }
    
    print("各数据源GCA数量:")
    for source, ids in sources.items():
        print(f"   {source}: {len(ids)}")
    
    # 找出真正一致的数据集
    non_empty_sources = {k: v for k, v in sources.items() if v}
    if non_empty_sources:
        consistent_dataset = set.intersection(*non_empty_sources.values())
        print(f"\n✅ 一致的数据集大小: {len(consistent_dataset)} 个GCA基因组")
        
        if consistent_dataset:
            print(f"一致数据集示例:")
            for strain in list(consistent_dataset)[:5]:
                print(f"   - {strain}")
    else:
        print("❌ 没有找到可用的GCA数据")
    
    print(f"\n🎯 推荐行动:")
    print("1. 基于一致的数据集继续分析")
    print("2. 如果数量合适(如50-100个)，直接使用现有Roary结果但只关注这些菌株")
    print("3. 如果追求完美一致性，用一致数据集重新运行Roary")

if __name__ == "__main__":
    generate_final_report()
