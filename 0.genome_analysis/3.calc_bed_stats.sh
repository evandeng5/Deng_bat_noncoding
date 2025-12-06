#!/bin/bash

# 1. 定义输出文件夹名称
OUT_DIR="filtered_bed_le1000"

# 2. 创建输出目录 (如果不存在)
if [ ! -d "$OUT_DIR" ]; then
    mkdir -p "$OUT_DIR"
    echo "Created output directory: $OUT_DIR"
else
    echo "Output directory already exists: $OUT_DIR"
fi

echo "Starting filtering process (Removing regions > 1000bp)..."
echo "-----------------------------------------------------"
printf "%-30s %-15s %-15s\n" "File Name" "Original Count" "Filtered Count"

# 3. 遍历当前目录下所有的 .bed 文件
for file in *.bed; do
    # 排除可能不存在的情况
    [ -e "$file" ] || continue

    # 定义输出文件路径
    out_file="$OUT_DIR/$file"

    # --- 核心步骤：使用 awk 过滤 ---
    # 逻辑：($3 - $2) 计算长度
    # 条件：<= 1000 (即保留小于等于1000的，去除大于1000的)
    # 动作：默认打印符合条件的行
    awk '($3 - $2) <= 1000' "$file" > "$out_file"

    # --- (可选) 统计并打印过滤前后的行数，方便您确认效果 ---
    orig_count=$(wc -l < "$file")
    filt_count=$(wc -l < "$out_file")
    
    printf "%-30s %-15s %-15s\n" "$file" "$orig_count" "$filt_count"
done

echo "-----------------------------------------------------"
echo "All done! Filtered files are in the '$OUT_DIR' folder."