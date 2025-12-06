#!/usr/bin/env bash
set -euo pipefail

# 用法：bash batch_maf_subset.sh [IN_DIR] [OUT_DIR] [SPECIES_LIST] [THREADS]
IN_DIR="${1:-.}"
OUT_DIR="${2:-subset_maf}"
SPECIES_LIST="${3:-species.lst}"
THREADS="${4:-1}"

# 1) 基本检查
command -v mafSpeciesSubset >/dev/null || { echo "[ERR] mafSpeciesSubset 未找到"; exit 1; }
[[ -s "$SPECIES_LIST" ]] || { echo "[ERR] 物种列表不存在或为空：$SPECIES_LIST"; exit 1; }

mkdir -p "$OUT_DIR"

# 2) 收集待处理文件（支持 .maf 与 .maf.gz）
# 修改：查找符号链接 (-type l)
mapfile -t FILES < <(find "$IN_DIR" -maxdepth 1 -type l \( -name "*.maf" -o -name "*.maf.gz" \) | sort)
[[ ${#FILES[@]} -gt 0 ]] || { echo "[WARN] 未找到 .maf/.maf.gz 文件"; exit 0; }

# 3) 定义单文件处理函数
run_one() {
  local f="$1"
  local fn base out
  fn="${f##*/}"
  base="$fn"
  base="${base%.maf.gz}"   # 先去掉 .maf.gz
  base="${base%.maf}"      # 再尝试去掉 .maf
  out="${OUT_DIR}/${base}.subset.maf"

  echo "[INFO] ${fn} -> ${out}"
  mafSpeciesSubset "$f" "$SPECIES_LIST" "$out"
}

export -f run_one
export OUT_DIR SPECIES_LIST

# 4) 串行/并行执行
if (( THREADS > 1 )) && command -v parallel >/dev/null; then
  printf "%s\0" "${FILES[@]}" | parallel -0 -j "$THREADS" run_one {}
else
  for f in "${FILES[@]}"; do
    run_one "$f"
  done
fi

echo "[DONE] 输出目录：$OUT_DIR"

