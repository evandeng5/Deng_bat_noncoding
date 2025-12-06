#!/bin/bash

# 遍历当前目录下所有的.gff文件
for file in ./*.gff; do
    # 获取文件名并去掉扩展名
    filename="${file%.gff}"

    # 提取第六列的 log10(p值) 并转回原始p值
    awk '{print $6}' "$file" | awk '{print 10^(-$1)}' > "${filename}_pvalues.tmp"

    # 使用 R 计算校正后的p值，并将结果输出到中间文件
    Rscript -e "
        pvalues <- read.table('${filename}_pvalues.tmp')\$V1
        padj <- p.adjust(pvalues, method='fdr')
        write.table(padj, file='${filename}_padj.tmp', row.names=FALSE, col.names=FALSE, quote=FALSE)
    "

    # 将校正后的p值追加到原始文件的最后一列，并筛选显著的行 (p.adjust <= 0.05)
    paste "$file" "${filename}_padj.tmp" | awk '$NF <= 0.05' > "${filename}_significant.bed"

    # 重命名染色体名称：将NC_0+替换为chr，去掉 .11
    sed -i 's/NC_0\+\([0-9]\+\)\.[0-9]\+/chr\1/' "${filename}_significant.bed"

    # 清理中间文件
    rm "${filename}_pvalues.tmp" "${filename}_padj.tmp"

    echo "Processed (adjusted p-values): ${file}"
done