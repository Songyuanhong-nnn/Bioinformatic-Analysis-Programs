#!/bin/bash
# 全局搜索提取，不依赖固定目录结构
echo -e "菌株名,毒力基因数,O血清型,K血清型" > final_results.csv

# 检查Kaptive结果文件
if [ ! -f "kaptive_o_serotype_results.tsv" ] || [ ! -f "kaptive_k_serotype_results.tsv" ]; then
  echo "错误：Kaptive血清型结果文件缺失！请确认已执行O/K抗原分型"
  exit 1
fi

# 处理Kaptive菌株名（统一替换特殊字符，便于匹配）
awk -F'\t' 'NR>1 {
  gsub(/\./,"_",$1); 
  gsub(/\-/,"_",$1); 
  print $1","$2
}' kaptive_o_serotype_results.tsv > o_map.tmp

awk -F'\t' 'NR>1 {
  gsub(/\./,"_",$1); 
  gsub(/\-/,"_",$1); 
  print $1","$2
}' kaptive_k_serotype_results.tsv > k_map.tmp

# 全局搜索所有annotation.gff3（项目目录下所有层级）
gff_files=$(find /mnt/d/WSL/disk/projects/VP1 -name "annotation.gff3")

if [ -z "$gff_files" ]; then
  echo "严重错误：未在项目目录下找到任何annotation.gff3文件！"
  echo "请先执行以下命令确认Bakta是否真的注释成功："
  echo "find /mnt/d/WSL/disk/projects/VP1 -name \"*.gff3\""
  exit 1
fi

# 遍历所有找到的注释文件
for gff in $gff_files; do
  # 提取菌株名（从文件路径中截取最可能的菌株目录名）
  # 假设路径格式如：xxx/菌株名/annotation.gff3，取倒数第二级目录名
  strain_raw=$(basename $(dirname "$gff"))
  # 统一菌株名格式（匹配Kaptive处理后的格式）
  strain=$(echo "$strain_raw" | sed 's/\./_/g; s/\-/_/g')
  
  # 统计毒力基因数（精准筛选功能描述含毒力相关的条目）
  vir_count=$(grep -c -E "virulence|toxin|pathogenic|pathogen|virulent" "$gff")
  
  # 匹配O血清型
  o_sero=$(awk -v s="$strain" -F',' '$1==s {print $2}' o_map.tmp)
  o_sero=${o_sero:-"Unknown"}
  
  # 匹配K血清型
  k_sero=$(awk -v s="$strain" -F',' '$1==s {print $2}' k_map.tmp)
  k_sero=${k_sero:-"Unknown"}
  
  # 写入结果（保留原始菌株名，便于后续核对）
  echo -e "\"$strain_raw\",\"$vir_count\",\"$o_sero\",\"$k_sero\"" >> final_results.csv
done

# 清理临时文件
rm -f o_map.tmp k_map.tmp

echo "数据提取完成！"
echo "找到的注释文件数量：$(echo "$gff_files" | wc -w)"
echo "结果文件：final_results.csv"
echo "提示：若部分菌株血清型为Unknown，可能是Kaptive未分型或菌株名匹配失败"
