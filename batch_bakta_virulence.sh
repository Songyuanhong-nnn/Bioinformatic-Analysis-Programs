#!/bin/bash
# 批量用 bakta 注释 VP1 目录下的无后缀菌株文件（CPxxxx.1 格式），提取毒力因子

# 菌株文件目录（直接在 VP1 目录下，无需子文件夹）
INPUT_DIR="./"
# 输出结果目录（自动创建，避免覆盖原有文件）
OUTPUT_DIR="./bakta_results_virulence"
mkdir -p $OUTPUT_DIR

# 批量处理所有无后缀的 CPxxxx.1 格式文件（精准匹配你的菌株文件）
for strain_file in $INPUT_DIR/CP*.1; do
  # 跳过不存在的文件（避免通配符不匹配时出错）
  [ -f "$strain_file" ] || continue
  
  # 获取菌株名称（比如从 "CP003972.1" 提取出 "CP003972.1" 作为名称）
  strain_name=$(basename "$strain_file")
  echo "正在处理菌株：$strain_name"
  
  # 用 bakta 注释（--virulence 提取毒力因子，--threads 8 用8线程加速，可根据CPU调整）
  bakta --virulence --threads 8 --output "$OUTPUT_DIR/$strain_name" "$strain_file"
  
  # 单独提取毒力因子结果到统一文件夹，方便后续汇总
  mkdir -p ./virulence_results_final
  if [ -f "$OUTPUT_DIR/$strain_name/virulence.tsv" ]; then
    cp "$OUTPUT_DIR/$strain_name/virulence.tsv" "./virulence_results_final/${strain_name}_virulence.tsv"
    echo "毒力因子结果已保存：./virulence_results_final/${strain_name}_virulence.tsv"
  else
    echo "警告：$strain_name 未检测到毒力因子，跳过结果复制"
  fi
  
  echo "$strain_name 处理完成！"
  echo "----------------------------------------"
done

echo "所有菌株处理完成！结果汇总："
echo "1. 完整 bakta 注释结果：$OUTPUT_DIR"
echo "2. 毒力因子提取结果：./virulence_results_final"
