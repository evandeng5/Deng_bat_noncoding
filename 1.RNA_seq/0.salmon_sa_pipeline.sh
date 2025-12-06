#!/usr/bin/env bash
set -euo pipefail

# -------------------------------
# Salmon SA decoy-aware pipeline
# Author: you
# -------------------------------
# Requirements:
#   - salmon >= 1.10
#   - gffread (optional; only if --build-from-gtf is used)
#
# Sample manifest (TSV), header required:
#   sample_id  r1                          r2
#   SRR13417529 /path/SRR13417529_1.fq.gz  /path/SRR13417529_2.fq.gz
#   SRR13417530 /path/SRR13417530_1.fq.gz  /path/SRR13417530_2.fq.gz
# For SE reads: put r2 as "-" (dash)
#
# Example:
#   ./salmon_sa_pipeline.sh \
#     --transcripts transcripts.fa \        # 或者 --build-from-gtf 配合 --gtf --genome
#     --genome genome.fa \
#     --gtf annotation.gtf \                # 若提供，将用于 -g 生成基因层面结果
#     --samples samples.tsv \
#     --outdir quants \
#     --index sa_index \
#     --threads 16 \
#     --kmer 31
#
# 仅有 GFF/GTF + 基因组、没有转录本FASTA 时：
  # ./salmon_sa_pipeline.sh \
  #   --build-from-gtf \
  #   --gtf annotation.gff \                # gtf 或 gff 均可（gffread 可自动识别）
  #   --genome genome.fa \
  #   --samples samples.tsv \
  #   --outdir quants \
  #   --index sa_index

show_help() {
  cat <<'EOF'
Usage: salmon_sa_pipeline.sh [options]

Required (Choose one transcript source):
  --transcripts  FILE   Transcripts FASTA (cDNA优先; CDS也可)
  --build-from-gtf       Use gffread to build transcripts from --gtf + --genome

Required genome/annotation:
  --genome      FILE   Genome FASTA (as decoys, and for gffread if enabled)
  --gtf         FILE   GTF/GFF for gene mapping (-g) and/or gffread extraction

Samples:
  --samples     FILE   Tab-separated manifest with header: sample_id  r1  r2
                       For SE data, set r2 to "-"

General:
  --outdir      DIR    Output directory (default: quants)
  --index       DIR    Salmon index dir (default: sa_index)
  --kmer        INT    k-mer size (default: 31)
  --threads     INT    Threads (default: 8)
  --libtype     STR    Salmon -l (default: A; auto-detect)
  --name        STR    Optional run name to tag logs

Flags:
  --skip-index          Skip index if already exists
  --help                Show this help

Notes:
  1) 推荐 cDNA 作转录本（更易捕获UTR reads）；只有 GFF/GTF+genome 时加 --build-from-gtf 自动用 gffread 提取。
  2) 若提供 --gtf，定量时会自动加 -g，生成 gene-level quant（quant.genes.sf）。
  3) decoys.txt 将从 --genome 中自动构建；gentrome.fa = transcripts.fa + genome.fa
EOF
}

# --------------- defaults ---------------
TRANSCRIPTS=""
BUILD_FROM_GTF="false"
GENOME=""
GTF=""
SAMPLES=""
OUTDIR="quants"
INDEX="sa_index"
KMER=31
THREADS=8
LIBTYPE="A"
SKIP_INDEX="false"
RUNNAME=""

# --------------- parse args ---------------
while [[ $# -gt 0 ]]; do
  case "$1" in
    --transcripts) TRANSCRIPTS="$2"; shift 2;;
    --build-from-gtf) BUILD_FROM_GTF="true"; shift 1;;
    --genome) GENOME="$2"; shift 2;;
    --gtf) GTF="$2"; shift 2;;
    --samples) SAMPLES="$2"; shift 2;;
    --outdir) OUTDIR="$2"; shift 2;;
    --index) INDEX="$2"; shift 2;;
    --kmer) KMER="$2"; shift 2;;
    --threads) THREADS="$2"; shift 2;;
    --libtype) LIBTYPE="$2"; shift 2;;
    --skip-index) SKIP_INDEX="true"; shift 1;;
    --name) RUNNAME="$2"; shift 2;;
    --help|-h) show_help; exit 0;;
    *) echo "[ERR] Unknown option: $1"; show_help; exit 1;;
  esac
done

# --------------- sanity checks ---------------
if [[ "$BUILD_FROM_GTF" != "true" && -z "$TRANSCRIPTS" ]]; then
  echo "[ERR] Provide --transcripts OR --build-from-gtf."
  exit 1
fi

if [[ -z "$GENOME" ]]; then
  echo "[ERR] --genome is required."
  exit 1
fi

if [[ -z "$SAMPLES" ]]; then
  echo "[ERR] --samples manifest is required."
  exit 1
fi

if [[ ! -f "$GENOME" ]]; then
  echo "[ERR] genome not found: $GENOME"
  exit 1
fi

if [[ ! -f "$SAMPLES" ]]; then
  echo "[ERR] samples manifest not found: $SAMPLES"
  exit 1
fi

mkdir -p "$OUTDIR" logs

# --------------- build transcripts if needed ---------------
TRANSC_FA=""

if [[ "$BUILD_FROM_GTF" == "true" ]]; then
  if ! command -v gffread >/dev/null 2>&1; then
    echo "[ERR] gffread not found but --build-from-gtf requested."
    exit 1
  fi
  if [[ -z "$GTF" || ! -f "$GTF" ]]; then
    echo "[ERR] --gtf is required and must exist when --build-from-gtf is used."
    exit 1
  fi
  TRANSC_FA="transcripts.from_gtf.fa"
  if [[ -s "$TRANSC_FA" ]]; then
    echo "[INFO] Found existing $TRANSC_FA, skipping gffread."
  else
    echo "[INFO] Running gffread to extract transcripts..."
    # gffread 自动识别 gtf/gff；-w 写出转录本cDNA序列
    gffread "$GTF" -g "$GENOME" -w "$TRANSC_FA"
  fi
else
  if [[ ! -f "$TRANSCRIPTS" ]]; then
    echo "[ERR] transcripts FASTA not found: $TRANSCRIPTS"
    exit 1
  fi
  TRANSC_FA="$TRANSCRIPTS"
fi

# --------------- decoys & gentrome ---------------
DECOYS="decoys.txt"
GENTROME="gentrome.fa"

if [[ -s "$DECOYS" ]]; then
  echo "[INFO] $DECOYS exists. Using it."
else
  echo "[INFO] Building $DECOYS from $GENOME"
  grep "^>" "$GENOME" | cut -d " " -f1 | sed 's/^>//' > "$DECOYS"
fi

if [[ -s "$GENTROME" ]]; then
  echo "[INFO] $GENTROME exists. Using it."
else
  echo "[INFO] Building $GENTROME (transcripts + genome; transcripts first)"
  cat "$TRANSC_FA" "$GENOME" > "$GENTROME"
fi

# --------------- salmon index ---------------
if [[ -d "$INDEX" && "$SKIP_INDEX" == "true" ]]; then
  echo "[INFO] Index $INDEX exists and --skip-index set. Skipping indexing."
else
  echo "[INFO] Building Salmon SA index at $INDEX (k=$KMER, threads=$THREADS)"
  salmon index -t "$GENTROME" -d "$DECOYS" -i "$INDEX" -k "$KMER" -p "$THREADS"
fi

# --------------- detect header & columns ---------------
header=$(head -n1 "$SAMPLES" | tr -d '\r')
IFS=$'\t' read -r h1 h2 h3 <<< "$header"
if [[ "$h1" != "sample_id" || "$h2" != "r1" || "$h3" != "r2" ]]; then
  echo "[ERR] samples manifest must have header: sample_id<tab>r1<tab>r2"
  exit 1
fi

# --------------- quant loop ---------------
echo "[INFO] Starting quantification..."
tail -n +2 "$SAMPLES" | while IFS=$'\t' read -r sid r1 r2; do
  [[ -z "$sid" ]] && continue
  out="$OUTDIR/$sid"
  log="logs/${sid}${RUNNAME:+.$RUNNAME}.log"

  if [[ -s "$out/quant.sf" ]]; then
    echo "[INFO] $sid exists -> $out/quant.sf, skip."
    continue
  fi

  echo "[INFO] Quantifying $sid ..."
  mkdir -p "$out"

  # Build command
  cmd=(salmon quant -i "$INDEX" -l "$LIBTYPE" -p "$THREADS" -o "$out" --validateMappings)

  # gene-level if GTF present
  if [[ -n "$GTF" && -f "$GTF" ]]; then
    cmd+=(-g "$GTF")  # Salmon 会输出 quant.genes.sf
  fi

  # PE vs SE
  if [[ "$r2" == "-" ]]; then
    # single-end: Salmon 需 --readFiles; 推荐强制给出 --fldMean/SD? 这里留默认
    cmd+=(-r "$r1")
  else
    cmd+=(-1 "$r1" -2 "$r2")
  fi

  # Run & log
  {
    echo "[CMD] ${cmd[*]}"
    "${cmd[@]}"
    echo "[DONE] $sid"
  } &> "$log"

done

echo "[ALL DONE] Outputs in $OUTDIR. Index in $INDEX."
