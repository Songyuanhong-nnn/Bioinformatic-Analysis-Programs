#!/bin/bash
# 精准提取数据，兼容已存在的注释目录
echo -e "菌株名,毒力基因数,O血清型,K血清型" > final_results.csv

# 检查Kaptive结果文件（若不存在，提示先运行血清型分析）
if [ ! -f "kaptive_o_serotype_results.tsv" ]; then
  echo "错误：未找到O抗原分型结果！请先运行Kaptive O抗原分析命令"
  exit 1
fi
if [ ! -f "kaptive_k_serotype_results.tsv" ]; then
  echo "错误：未找到K抗原分型结果！请先运行Kaptive K抗原分析命令"
  exit 1
fi

# 处理Kaptive结果，建立菌株-血清型映射（兼容特殊字符）
awk -F'\t' 'NR>1 {gsub(/\./,"_",$1); print $1","$2}' kaptive_o_serotype_results.tsv > o_map.tmp
awk -F'\t' 'NR>1 {gsub(/\./,"_",$1); print $1","$2}' kaptive_k_serotype_results.tsv > k_map.tmp

# 动态搜索所有Bakta注释的annotation.gff3（绝对路径，不遗漏）
gff_files=$(find /mnt/d/WSL/disk/projects/VP1/bakta_annotations -name "annotation.gff3")

if [ -z "$gff_files" ]; then
  echo "警告：未找到annotation.gff3文件！请确认注释目录是否为bakta_annotations"
  exit 1
fi

# 遍历每个注释文件，提取数据
for gff in $gff_files; do
  # 提取菌株名（从目录名截取，替换特殊字符适配Kaptive结果）
  strain_raw=$(basename $(dirname "$gff"))
  strain=$(echo "$strain_raw" | sed 's/\./_/g')
  
  # 统计毒力基因数（筛选virulence/toxin相关条目，避免漏检）
  vir_count=$(grep -c -E "virulence|toxin|pathogen" "$gff")
  
  # 匹配O血清型（无则填Unknown）
  o_sero=$(awk -v s="$strain" -F',' '$1==s {print $2}' o_map.tmp)
  o_sero=${o_sero:-"Unknown"}
  
  # 匹配K血清型（无则填Unknown）
  k_sero=$(awk -v s="$strain" -F',' '$1==s {print $2}' k_map.tmp)
  k_sero=${k_sero:-"Unknown"}
  
  # 写入结果（还原原始菌株名，便于核对）
  echo -e "\"$strain_raw\",\"$vir_count\",\"$o_sero\",\"$k_sero\"" >> final_results.csv
done

# 清理临时文件
rm -f o_map.tmp k_map.tmp

echo "数据提取完成！"
echo "结果文件：final_results.csv"
echo "包含字段：菌株名、毒力基因数、O血清型、K血清型"
