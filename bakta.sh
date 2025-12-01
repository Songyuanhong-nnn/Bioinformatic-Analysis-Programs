#!/bin/bash
conda activate bakta
# ===================== 固定参数（直接用，无需修改）=====================
DB_PATH="/mnt/g/WSL/database/bakta/db"
OUTPUT_ROOT="./bakta_results"
THREADS=12
GENUS="Vibrio"
SPECIES="vulnificus"
GRAM="-"
FNA_DIR="./"

# ===================== 前置检查（确保环境正常）=====================
# 检查输入目录
if [ ! -d "$FNA_DIR" ]; then
    echo "❌ 错误：目录 $FNA_DIR 不存在！"
    exit 1
fi
# 检查 .fna 文件数量
fna_count=$(ls -1 "$FNA_DIR"/*.fna 2>/dev/null | wc -l)
if [ $fna_count -eq 0 ]; then
    echo "❌ 错误：$FNA_DIR 下没找到 .fna 文件！"
    exit 1
fi
# 创建输出目录
mkdir -p "$OUTPUT_ROOT"
# 删除旧的失败记录文件
rm -f "$OUTPUT_ROOT/failed_samples.txt"

# ===================== 核心逻辑：完整性检测+选择性重跑=====================
count=0
success_count=0
skip_count=0
fail_count=0

for fna_file in "$FNA_DIR"/*.fna; do
    [ -f "$fna_file" ] || continue
    count=$((count + 1))
    sample_name=$(basename "$fna_file" | cut -d'_' -f1-2)
    sample_output="$OUTPUT_ROOT/$sample_name"
    gbff_file="$sample_output/${sample_name}.gbff"
    tsv_file="$sample_output/${sample_name}.tsv"
    faa_file="$sample_output/${sample_name}.faa"

    echo -e "\n========================================"
    echo "正在处理样本（$count/$fna_count）：$sample_name"
    echo "输入文件：$fna_file"
    echo "输出目录：$sample_output"
    echo "========================================"

    # 🔍 精准完整性检测（3重验证，保留已成功样本）
    is_complete=0
    if [ -f "$gbff_file" ] && [ -f "$tsv_file" ] && [ -f "$faa_file" ]; then
        # 验证文件大小（过滤空文件/不完整文件）
        gbff_size=$(wc -c < "$gbff_file")
        tsv_size=$(wc -c < "$tsv_file")
        faa_size=$(wc -c < "$faa_file")
        if [ $gbff_size -gt 500000 ] && [ $tsv_size -gt 20000 ] && [ $faa_size -gt 50000 ]; then
            # 验证 gbff 结尾是否有 "//"（Bakta 完整注释标志）
            tail -5 "$gbff_file" | grep -q "//"
            if [ $? -eq 0 ]; then
                is_complete=1
            fi
        fi
    fi

    # ✅ 注释完整，直接跳过
    if [ $is_complete -eq 1 ]; then
        echo "✅ 样本 $sample_name 注释完整，跳过处理！"
        skip_count=$((skip_count + 1))
        continue
    fi

    # ⚠️ 注释不完整/未运行，删除旧结果并重跑
    if [ -d "$sample_output" ]; then
        echo "⚠️ 样本 $sample_name 注释不完整，删除旧结果后重跑..."
        rm -rf "$sample_output"
    fi

    # 直接调用 bakta（无需路径，依赖当前激活的 Conda 环境）
    bakta \
        --db "$DB_PATH" \
        --output "$sample_output" \
        --prefix "$sample_name" \
        --genus "$GENUS" \
        --species "$SPECIES" \
        --gram "$GRAM" \
        --threads "$THREADS" \
        --skip-plot \
        --complete \
        --force \
        "$fna_file"

    # 🔍 重跑后再次验证完整性
    if [ $? -eq 0 ] && [ -f "$gbff_file" ] && [ -f "$tsv_file" ] && [ -f "$faa_file" ]; then
        tail -5 "$gbff_file" | grep -q "//"
        if [ $? -eq 0 ]; then
            echo "✅ 样本 $sample_name 处理成功！"
            success_count=$((success_count + 1))
        else
            echo "❌ 样本 $sample_name 处理后仍不完整！"
            echo "$sample_name" >> "$OUTPUT_ROOT/failed_samples.txt"
            fail_count=$((fail_count + 1))
        fi
    else
        echo "❌ 样本 $sample_name 处理失败！"
        echo "$sample_name" >> "$OUTPUT_ROOT/failed_samples.txt"
        fail_count=$((fail_count + 1))
    fi
done

# ===================== 处理总结（清晰统计结果）=====================
echo -e "\n========================================"
echo "📊 处理结果统计："
echo "总样本数：$count"
echo "跳过（已完整）：$skip_count 个"
echo "处理成功：$success_count 个"
echo "处理失败：$fail_count 个"
echo "========================================"
echo "结果保存目录：$OUTPUT_ROOT"
if [ -f "$OUTPUT_ROOT/failed_samples.txt" ]; then
    echo "❌ 失败样本列表：$OUTPUT_ROOT/failed_samples.txt"
else
    echo "✅ 所有样本处理成功！"
fi
echo "========================================"
