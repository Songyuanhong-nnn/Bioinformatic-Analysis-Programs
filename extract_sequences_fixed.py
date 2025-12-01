import pandas as pd
import os

def extract_target_sequences():
    print("=== 🧬 提取靶标基因真实序列 ===")
    
    # 读取pan_genome_reference.fa
    pan_genome_file = "roary_consistent/pan_genome_reference.fa"
    
    if not os.path.exists(pan_genome_file):
        print("❌ 未找到pan_genome_reference.fa")
        print("请检查文件是否存在:")
        print("ls -la roary_consistent/pan_genome_reference.fa")
        return
    
    # 读取候选靶标
    target_file = "real_candidate_targets_57strains.csv"
    if not os.path.exists(target_file):
        print("❌ 请先运行靶标分析")
        return
    
    targets = pd.read_csv(target_file)
    top_targets = targets.head(10)['Gene'].tolist()  # 取前10个
    
    print(f"提取前10个靶标的序列: {top_targets}")
    
    # 从pan_genome_reference.fa中提取序列
    sequences = {}
    current_gene = None
    current_seq = []
    
    with open(pan_genome_file, 'r') as f:
        for line in f:
            if line.startswith('>'):
                # 保存前一个基因
                if current_gene and current_gene in top_targets:
                    sequences[current_gene] = ''.join(current_seq)
                    print(f"✅ 提取: {current_gene}")
                
                # 开始新基因
                current_gene = line[1:].strip().split()[0]  # 取基因名
                current_seq = []
            else:
                current_seq.append(line.strip())
        
        # 保存最后一个基因
        if current_gene and current_gene in top_targets:
            sequences[current_gene] = ''.join(current_seq)
            print(f"✅ 提取: {current_gene}")
    
    # 输出结果
    print(f"\n🎯 成功提取 {len(sequences)} 个靶标序列")
    
    if sequences:
        # 保存到文件
        with open('target_sequences.fasta', 'w') as f:
            for gene, seq in sequences.items():
                f.write(f'>{gene}\n{seq}\n')
        
        print(f"💾 序列已保存到: target_sequences.fasta")
        
        # 显示序列信息
        print(f"\n📊 序列信息:")
        for gene, seq in sequences.items():
            print(f"   {gene}: {len(seq)} bp")
            
    else:
        print("❌ 未找到任何靶标序列")
        print("可能的原因:")
        print("1. 基因名不匹配")
        print("2. pan_genome_reference.fa不包含这些基因")
        print("3. 文件格式问题")

if __name__ == "__main__":
    extract_target_sequences()
