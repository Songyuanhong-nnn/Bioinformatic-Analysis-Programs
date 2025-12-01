import pandas as pd
import os

def extract_target_sequences():
    print("=== 🧬 提取靶标基因真实序列 ===")
    
    # 读取pan_genome_reference.fa
    pan_genome_file = "roary_consistent/pan_genome_reference.fa"
    
    if not os.path.exists(pan_genome_file):
        print("❌ 未找到pan_genome_reference.fa")
        return
    
    # 读取候选靶标
    target_file = "real_candidate_targets_57strains.csv"
    if not os.path.exists(target_file):
        print("❌ 请先运行靶标分析")
        return
    
    targets = pd.read_csv(target_file)
    top_targets = targets.head(5)['Gene'].tolist()
    
    print(f"提取前5个靶标的序列: {top_targets}")
    
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
                
                # 开始新基因
                current_gene = line[1:].strip().split()[0]  # 取基因名
                current_seq = []
            else:
                current_seq.append(line.strip())
        
        # 保存最后一个基因
        if current_gene and current_gene in top_targets:
            sequences[current_gene] = ''.join(current_seq)
    
    # 输出结果
    print(f"\n✅ 成功提取 {len(sequences)} 个靶标序列")
    for gene, seq in sequences.items():
        print(f"\n🎯 {gene}:")
        print(f"   序列长度: {len(seq)} bp")
        print(f"   前50bp: {seq[:50]}...")
        print(f"   后50bp: ...{seq[-50:]}" if len(seq) > 100 else "")
    
    # 保存到文件
    with open('target_sequences_57strains.fasta', 'w') as f:
        for gene, seq in sequences.items():
            f.write(f'>{gene}\n{seq}\n')
    
    print(f"\n💾 序列已保存到: target_sequences_57strains.fasta")

if __name__ == "__main__":
    extract_target_sequences()
