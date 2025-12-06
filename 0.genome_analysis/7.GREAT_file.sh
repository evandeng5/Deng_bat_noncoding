for file in ./*significant.bed; do
    filename="${file%_significant.bed}"
    # 使用 awk 提取第一、四、五列，并保持列之间使用 TAB 分隔
    awk '{print $1 "\t" $4 "\t" $5}' "$file" > "${filename}_Great.bed"

    echo "Processed: ${file}"
done
