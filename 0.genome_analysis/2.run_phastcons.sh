#!/bin/bash

# --- 1. 设置路径变量 ---
# 模型文件路径
MOD_FILE="/home/dxy/2024_bat_comparative_genome/1.analysis/1.accelerated_detect/nonconserved-all-4d_6sp.mod"

# 输入 MAF 文件夹
MAF_DIR="/home/dxy/2024_bat_comparative_genome/1.analysis/1.accelerated_detect/1.maf/6_sp_maf"

# 输出目录 (相对路径，会在当前目录下创建)
WIG_DIR="wig"
BED_DIR="bed"

# 创建输出目录
mkdir -p "$WIG_DIR"
mkdir -p "$BED_DIR"

# 设置并行任务数
JOBS=30

# --- 2. 执行并行 PhastCons ---

echo "Starting phastCons with $JOBS jobs..."

# 关键修正：使用 ls 获取绝对路径，但在 parallel 内部处理文件名
ls "$MAF_DIR"/NC_*.subset.maf | parallel -j "$JOBS" '
    # 1. 获取输入的绝对路径
    INPUT_FILE={}
    
    # 2. 【关键修正】提取纯文件名 (去除路径)
    # 例如: /home/.../NC_000001.11.subset.maf -> NC_000001.11.subset.maf
    FILENAME=$(basename "$INPUT_FILE")
    
    # 3. 提取序列 ID (去除 .subset.maf 后缀)
    # 例如: NC_000001.11.subset.maf -> NC_000001.11
    SEQ_ID="${FILENAME%.subset.maf}"
    
    # 4. 定义输出文件路径 (使用相对路径)
    # 注意：变量引用需要正确嵌套引号
    OUT_WIG="'"$WIG_DIR"'/${SEQ_ID}.wig"
    OUT_BED="'"$BED_DIR"'/${SEQ_ID}.bed"
    
    echo "Processing ${SEQ_ID} ..."
    
    # 5. 运行 phastCons
    phastCons "$INPUT_FILE" \
        "'"$MOD_FILE"'" \
        --target-coverage 0.3 \
        --expected-length 45 \
        --rho 0.3 \
        --msa-format MAF \
        --seqname "$SEQ_ID" \
        --most-conserved "$OUT_BED" \
        > "$OUT_WIG"
'

echo "All phastCons jobs finished."