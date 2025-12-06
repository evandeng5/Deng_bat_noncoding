#!/bin/bash
set -e

THREADS=10
export THREADS

# 让通配在空集合时不报错
shopt -s nullglob

find . -maxdepth 1 -type d -not -name "." | parallel -j "$THREADS" --no-run-if-empty --progress '
    set -euo pipefail
    shopt -s expand_aliases

    dir="{}"
    accessionID=$(basename "$dir")
    SRA_FILE="$dir/${accessionID}.sra"

    if [ -f "$SRA_FILE" ]; then
        echo "[INFO] SRA file for $accessionID already exists, skip prefetch."
    else
        echo "[INFO] Fetching $accessionID..."
        prefetch -O "$dir" "$accessionID" || { echo "[ERR] Prefetch failed for $accessionID"; exit 1; }
    fi

    if [ -f "$SRA_FILE" ]; then
        echo "[INFO] Converting $accessionID to FASTQ..."
        fasterq-dump --split-files -O "$dir" "$SRA_FILE" || { echo "[ERR] fasterq-dump failed for $accessionID"; exit 1; }

        # 用 find 枚举，避免 *.fastq.gz 未匹配时的字面量残留
        mapfile -t fastq_files < <(find "$dir" -maxdepth 1 -type f \( -name "*.fastq" -o -name "*.fastq.gz" \) | sort)
        num_fastq_files=${#fastq_files[@]}

        # 尝试匹配配对文件
        R1=""
        R2=""

        # 优先找 *_1 / *_2（含 .gz）
        mapfile -t r1s < <(find "$dir" -maxdepth 1 -type f \( -name "*_1.fastq" -o -name "*_1.fastq.gz" \) | sort)
        mapfile -t r2s < <(find "$dir" -maxdepth 1 -type f \( -name "*_2.fastq" -o -name "*_2.fastq.gz" \) | sort)

        if [[ ${#r1s[@]} -ge 1 && ${#r2s[@]} -ge 1 ]]; then
            R1="${r1s[0]}"
            R2="${r2s[0]}"
        fi

        if [[ -n "$R1" && -n "$R2" ]]; then
            echo "[INFO] Paired-end detected for $accessionID → ($R1 , $R2)"
            fastp -i "$R1" -I "$R2" \
                  -o "$dir/${accessionID}_1_clean.fastq" -O "$dir/${accessionID}_2_clean.fastq" \
                  -l 60 -w 10 \
                  --html "$dir/${accessionID}_fastp_report.html" \
                  --json "$dir/${accessionID}_fastp_report.json" \
                  || { echo "[ERR] fastp (PE) failed for $accessionID"; exit 1; }
            echo "[DONE] fastp paired-end for $accessionID"
        else
            # 单端：允许 fasterq-dump 只产出 *_1.fastq 的情况
            if [[ $num_fastq_files -eq 1 ]]; then
                SE_IN="${fastq_files[0]}"
                echo "[INFO] Single-end detected for $accessionID → $SE_IN"
                fastp -i "$SE_IN" -o "$dir/${accessionID}_clean.fastq" \
                      -l 15 \
                      -h "$dir/${accessionID}_fastp_report.html" \
                      -j "$dir/${accessionID}_fastp_report.json" \
                      || { echo "[ERR] fastp (SE) failed for $accessionID"; exit 1; }
                echo "[DONE] fastp single-end for $accessionID"
            elif [[ $num_fastq_files -eq 2 ]]; then
                # 兼容某些单端流程把一个文件压缩成 .gz 的场景（解压后仍是 1 个）
                # 但更稳妥的做法是只在明确只有 1 个时当 SE，这里仅提示
                echo "[WARN] Found 2 FASTQs for $accessionID:"
                printf "       %s\n" "${fastq_files[@]}"
                echo "[WARN] If this is single-end with one .gz backup, remove the extra file and rerun."
            else
                echo "[WARN] Unexpected FASTQ count ($num_fastq_files) for $accessionID:"
                printf "       %s\n" "${fastq_files[@]}"
                echo "[WARN] Skip fastp for $accessionID due to ambiguous inputs."
            fi
        fi
    else
        echo "[ERR] Prefetch failed for $accessionID, skip fasterq-dump and fastp."
    fi
'

echo "All tasks completed successfully!"
