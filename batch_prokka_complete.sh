#!/bin/bash
# 遍历所有整理后的基因组（raw_genomes目录下的FASTA文件）
for genome_fasta in raw_genomes/*.fna; do
    # 提取基因组ID（去掉路径和后缀，统一格式）
    genome_id=$(basename "$genome_fasta" .fna | sed 's/\./_/g')
    
    # 跳过已完成注释的基因组（避免重复）
    if [ -d "prokka_annotations/$genome_id" ]; then
        echo "✅ 已注释：$genome_id，跳过"
        continue
    fi
    
    echo "🚀 正在注释：$genome_id"
    # Prokka核心命令（修复3个问题：去掉多余反斜杠、添加输入文件、参数顺序正确）
    prokka --outdir prokka_annotations/"$genome_id" \
           --prefix "$genome_id" \
           --kingdom Bacteria \
           --genus Vibrio \
           --species parahaemolyticus \
           --gcode 11 \
           --cpus 8 \
           --force \
           "$genome_fasta"  # 关键：指定要注释的FASTA文件路径
done

# 收集所有GFF文件（用于后续Roary泛基因组分析）
mkdir -p roary_input
cp prokka_annotations/*/*.gff roary_input/ 2>/dev/null
echo "📥 已收集所有注释后的GFF文件到 roary_input/ 目录"
