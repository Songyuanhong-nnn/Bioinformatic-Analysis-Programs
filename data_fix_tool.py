import os
import glob
import shutil

def create_consistent_roary_input():
    """创建一致的Roary输入数据集"""
    print("=== 🛠️ 创建一致的Roary输入数据集 ===")
    
    # 创建干净的目录
    os.makedirs('roary_input_consistent', exist_ok=True)
    
    # 基于clean_genome中的57个GCA文件
    clean_genomes = glob.glob("clean_genome/GCA_*.fna")
    print(f"基于clean_genome中的 {len(clean_genomes)} 个GCA文件")
    
    matched_count = 0
    
    for fna in clean_genomes:
        clean_name = os.path.basename(fna).replace('.fna', '')  # GCA_000328405.1_ASM32840v1_genomic
        roary_name = clean_name.replace('.', '_')  # GCA_000328405_1_ASM32840v1_genomic
        
        # 查找对应的GFF文件
        source_gff1 = f"roary_input/{roary_name}.gff"
        source_gff2 = f"bakta_annotations/{clean_name}/{clean_name}.gff3"
        
        target_gff = f"roary_input_consistent/{roary_name}.gff"
        
        # 优先使用bakta的GFF，其次使用roary_input中的
        if os.path.exists(source_gff2):
            shutil.copy2(source_gff2, target_gff)
            matched_count += 1
            print(f"✅ 使用bakta: {clean_name}")
        elif os.path.exists(source_gff1):
            shutil.copy2(source_gff1, target_gff) 
            matched_count += 1
            print(f"✅ 使用roary_input: {roary_name}")
        else:
            print(f"❌ 未找到GFF: {clean_name}")
    
    print(f"\n创建了 {matched_count} 个一致的GFF文件")
    return matched_count

def recommend_next_steps(matched_count):
    print(f"\n=== 🚀 下一步建议 ===")
    
    if matched_count >= 50:
        print(f"✅ 成功匹配 {matched_count} 个菌株，数量充足")
        print("建议:")
        print("1. 使用 roary_input_consistent/ 重新运行Roary")
        print("2. 命令: roary -p 20 -f roary_consistent -e -n -v roary_input_consistent/*.gff")
        print("3. 基于新结果继续靶标分析")
    elif matched_count >= 30:
        print(f"⚠️ 匹配 {matched_count} 个菌株，数量尚可")
        print("可以继续分析，但建议补充更多数据")
    else:
        print(f"❌ 只匹配 {matched_count} 个菌株，数量不足")
        print("需要检查数据完整性")

if __name__ == "__main__":
    matched = create_consistent_roary_input()
    recommend_next_steps(matched)
