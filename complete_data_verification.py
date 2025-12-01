import pandas as pd
import os
import glob

def complete_verification():
    print("=== 🔍 完整数据验证 ===")
    
    # 检查所有Roary输出目录
    roary_dirs = glob.glob("roary_*")
    print("找到的Roary输出目录:")
    for dir_path in roary_dirs:
        if os.path.isdir(dir_path):
            size = sum(os.path.getsize(os.path.join(dir_path, f)) for f in os.listdir(dir_path) if os.path.isfile(os.path.join(dir_path, f)))
            print(f"  📁 {dir_path}: {size} 字节")
    
    # 检查每个目录的gene_presence_absence.csv
    print(f"\n📊 各数据集详细信息:")
    
    datasets = {}
    
    for dir_path in roary_dirs:
        gene_pa_file = os.path.join(dir_path, "gene_presence_absence.csv")
        if os.path.exists(gene_pa_file):
            try:
                gene_pa = pd.read_csv(gene_pa_file, low_memory=False)
                strains = len(gene_pa.columns[14:])
                genes = len(gene_pa)
                datasets[dir_path] = {
                    'strains': strains,
                    'genes': genes,
                    'file': gene_pa_file
                }
                print(f"  ✅ {dir_path}: {strains}菌株, {genes}基因簇")
            except Exception as e:
                print(f"  ❌ {dir_path}: 读取错误 - {e}")
    
    # 推荐最佳数据集
    print(f"\n🎯 推荐使用的最佳数据集:")
    if datasets:
        # 优先选择57菌株的数据集
        best_dataset = None
        for dir_path, info in datasets.items():
            if info['strains'] == 57:
                best_dataset = (dir_path, info)
                break
        
        if not best_dataset:
            # 如果没有57菌株的，选择菌株数最接近57的
            best_dataset = min(datasets.items(), key=lambda x: abs(x[1]['strains'] - 57))
        
        dir_path, info = best_dataset
        print(f"  推荐: {dir_path}")
        print(f"    菌株: {info['strains']}")
        print(f"    基因簇: {info['genes']}")
        print(f"    文件: {info['file']}")
        
        return dir_path, info
    else:
        print("  ❌ 未找到可用的数据集")
        return None, None

def check_pan_genome_files(roary_dir):
    """检查pan_genome_reference.fa等关键文件"""
    print(f"\n📁 检查 {roary_dir} 的关键文件:")
    
    key_files = {
        "pan_genome_reference.fa": "泛基因组序列",
        "summary_statistics.txt": "统计摘要", 
        "clustered_proteins": "蛋白聚类"
    }
    
    all_exist = True
    for file, description in key_files.items():
        path = os.path.join(roary_dir, file)
        if os.path.exists(path):
            size = os.path.getsize(path) if os.path.isfile(path) else "目录"
            print(f"  ✅ {file}: {description} - 存在 ({size})")
        else:
            print(f"  ❌ {file}: {description} - 不存在")
            all_exist = False
    
    return all_exist

def extract_from_correct_dataset(roary_dir):
    """从正确数据集中提取序列"""
    print(f"\n🧬 从 {roary_dir} 提取序列:")
    
    pan_genome_file = os.path.join(roary_dir, "pan_genome_reference.fa")
    gene_pa_file = os.path.join(roary_dir, "gene_presence_absence.csv")
    
    if not os.path.exists(pan_genome_file):
        print(f"  ❌ 未找到 {pan_genome_file}")
        return
    
    if not os.path.exists('real_candidate_targets_57strains.csv'):
        print("  ❌ 请先运行靶标分析")
        return
    
    targets = pd.read_csv('real_candidate_targets_57strains.csv')
    top_targets = targets.head(5)['Gene'].tolist()
    
    print(f"  提取前5个靶标: {top_targets}")
    
    # 提取序列
    sequences = {}
    current_gene = None
    current_seq = []
    
    with open(pan_genome_file, 'r') as f:
        for line in f:
            if line.startswith('>'):
                # 保存前一个基因
                if current_gene and current_gene in top_targets:
                    sequences[current_gene] = ''.join(current_seq)
                    print(f"  ✅ 提取: {current_gene}")
                
                # 开始新基因
                current_gene = line[1:].strip().split()[0]
                current_seq = []
            else:
                current_seq.append(line.strip())
        
        # 保存最后一个基因
        if current_gene and current_gene in top_targets:
            sequences[current_gene] = ''.join(current_seq)
            print(f"  ✅ 提取: {current_gene}")
    
    if sequences:
        # 保存序列
        with open('target_sequences_correct.fasta', 'w') as f:
            for gene, seq in sequences.items():
                f.write(f'>{gene}\n{seq}\n')
        
        print(f"  💾 序列保存到: target_sequences_correct.fasta")
        
        # 显示序列信息
        print(f"\n  📊 提取的序列:")
        for gene, seq in sequences.items():
            print(f"    {gene}: {len(seq)} bp")
    else:
        print("  ❌ 未找到靶标序列")

if __name__ == "__main__":
    roary_dir, info = complete_verification()
    
    if roary_dir and info:
        files_ok = check_pan_genome_files(roary_dir)
        
        if files_ok:
            print(f"\n🎉 数据集 {roary_dir} 完整可用！")
            extract_from_correct_dataset(roary_dir)
        else:
            print(f"\n⚠️ 数据集 {roary_dir} 缺少关键文件")
    else:
        print("\n❌ 没有可用的完整数据集")
