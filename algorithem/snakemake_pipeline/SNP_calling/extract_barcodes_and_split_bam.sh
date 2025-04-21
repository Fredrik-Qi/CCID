#!/bin/bash

# 检查参数
if [ "$#" -ne 2 ]; then
    echo "用法: $0 barcode.txt 输出目录"
    exit 1
fi

# 读取参数
barcode_file="$1"
output_dir="$2"

# 创建输出目录（如果不存在）
mkdir -p "$output_dir"

# 自动检测 Sort_Population 列号
sort_col=$(awk -F'\t' 'NR==1 {for(i=1;i<=NF;i++) if($i=="Sort_Population") print i}' "$barcode_file")
echo "Sort_Population 列号: $sort_col"

# 确保找到了 Sort_Population 列
if [ -z "$sort_col" ]; then
    echo "错误：未找到 Sort_Population 列"
    exit 1
fi

# 获取所有唯一的 Sort Population 值
sort_values=$(awk -F'\t' -v col="$sort_col" 'NR>1 {print $col}' "$barcode_file" | sort -u)

# 获取 XC_DNBPE 所在的列号
col_num=$(awk -F'\t' 'NR==1 {for(i=1;i<=NF;i++) if($i=="XC_DNBPE") print i}' "$barcode_file")

# 确保找到 XC_DNBPE 列
if [ -z "$col_num" ]; then
    echo "错误：未找到 XC_DNBPE 列"
    exit 1
fi

# 遍历每个 Sort Population 值，进行数据处理
for value in $sort_values; do
    # 生成文件名
    subset_file="$output_dir/${value}_subset_stats_all_0.txt"
    output_file="$output_dir/${value}_subset_stats_all_0_cellBClist.txt"

    # 提取对应的行存入 subset_file
    awk -F'\t' -v val="$value" -v col="$sort_col" '$col == val' "$barcode_file" > "$subset_file"

    # 提取 XC_DNBPE 列，并存入 output_file
    awk -F'\t' -v col="$col_num" '{print $col}' "$subset_file" > "$output_file"

    echo "已生成: $output_file"
done
