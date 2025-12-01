#!/bin/bash
# 副溶血性弧菌Prokka注释重试脚本（完整版本）
for genome_fasta in raw_genomes/*.fna; do
    # 提取基因组ID并处理特殊字符
    genome_id=$(basename "$genome_fasta" .fna | sed 's/\./_/g')
    annot_dir="prokka_annotations/$genome_id"
    
    # 严格判定注释成功：存在GFF和FAA文件，且FAA序列数≥1000
    faa_file=$(find "$annot_dir" -name "*.faa" 2>/dev/null | head -n1)
    gff_file=$(find "$annot_dir" -name "*.gff" 2>/dev/null | head -n1)
    
    if [ -n "$gff_file" ] && [ -n "$faa_file" ]; then
        faa_count=$(grep -c "^>" "$faa_file" 2>/dev/null)
        if [ "$faa_count" -ge 1000 ]; then
            echo "✅ 已成功：$genome_id（CDS：$faa_count），跳过"
            continue
        fi
    fi
    
    # 注释失败：删除旧目录并重新运行Prokka
    echo "🔄 重新注释：$genome_id"
    rm -rf "$annot_dir"
    prokka --outdir "$annot_dir" \
           --prefix "$genome_id" \
           --kingdom Bacteria \
           --genus Vibrio \
           --species parahaemolyticus \
           --gcode 11 \
           --cpus 12 \
           --force \
           "$genome_fasta"
done

# 收集所有成功的GFF文件
mkdir -p roary_input
rm -rf roary_input/*.gff
cp prokka_annotations/*/*.gff roary_input/ 2>/dev/null
final_gff=$(ls roary_input/*.gff | wc -l)
total_genome=$(ls raw_genomes/*.fna | wc -l)
echo -e "\n📊 注释任务结束！"
echo "成功注释：$final_gff 个 | 总基因组：$total_genome 个"
echo "GFF文件已保存到：roary_input/"
