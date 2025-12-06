#!/bin/bash

# --- 配置部分 ---
# Conda 环境名称
CONDA_ENV_NAME="genome"

# GFF 文件路径
# 注意：这里必须是您之前拆分 GFF 所在的目录
export GFF_DIR="/home/dxy/shared_data/species_runned/Homo_sapiens/gff_split_by_chr"

# 输出文件前缀
export OUTPUT_PREFIX="4d_codon_6sp_"

# 并行任务数量 (设置为服务器 CPU 核心数，这里设为 20)
NUM_JOBS=20

# --- 环境检查 ---

# 激活 Conda 环境 (确保脚本可以找到 msa_view)
# 注意：在脚本中 source activate 比 conda activate 更稳定
source ~/tools/miniconda3/bin/activate "$CONDA_ENV_NAME" 2>/dev/null || conda activate "$CONDA_ENV_NAME"

# 检查 parallel 是否安装
if ! command -v parallel &> /dev/null
then
    echo "GNU parallel is not installed in the $CONDA_ENV_NAME environment."
    echo "Please run: conda install -c conda-forge parallel"
    exit 1
fi

echo "Starting parallel processing with $NUM_JOBS jobs..."

# --- 核心执行部分 (修正版) ---

# 解释：
# 1. 我们 export 了 GFF_DIR 和 OUTPUT_PREFIX，所以 parallel 内部的子 shell 可以直接读取它们，不需要拼接。
# 2. 整个命令块用单引号 '...' 包裹，防止外部 Shell 提前解析变量。
# 3. {} 和 {.} 是 parallel 的特殊符号，会被自动替换。

ls NC_*.maf | parallel -j "$NUM_JOBS" '
    # 1. 获取输入文件名
    INPUT_NAME={}
    
    # 2. 提取染色体 ID (去除 .subset.maf 后缀)
    # 例如：NC_000001.11.subset.maf -> NC_000001.11
    CHR_ID=$(echo "$INPUT_NAME" | sed "s/\.subset\.maf//")
    
    # 3. 构建该染色体对应的 GFF 文件路径
    # 使用环境变量 $GFF_DIR
    CHR_GFF_FILE="${GFF_DIR}/${CHR_ID}.gff"
    
    # 4. 检查 GFF 文件是否存在
    if [ ! -f "$CHR_GFF_FILE" ]; then
        echo "Error: GFF file not found for $CHR_ID at $CHR_GFF_FILE"
        exit 1
    fi

    # 5. 构建输出文件名
    # {.} 表示去除最后一个后缀(.maf)，变成 .subset
    # 结果示例: 4d_codon_6sp_NC_000001.11.subset.ss
    OUTPUT_FILE="${OUTPUT_PREFIX}{.}.ss"

    echo "Processing $INPUT_NAME with GFF $CHR_GFF_FILE ..."

    # 6. 执行 msa_view
    msa_view "$INPUT_NAME" --4d --features "$CHR_GFF_FILE" > "$OUTPUT_FILE"
'

echo "Batch processing complete."