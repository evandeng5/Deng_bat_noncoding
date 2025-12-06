#!/usr/bin/env bash
set -euo pipefail

# ====== 可配置参数 ======
MAF_PATH="/home/dxy/2024_bat_comparative_genome/2.data/WGS/evolver21species_9_9.maf"
BED_GLOB="./*.bed"
OUT_MAF_DIR="./maf_output"
OUT_FASTA_DIR="./fasta_output"
LOG_TSV="./logs/qc_summary.tsv"

# 覆盖率与计数阈值
COV_LOW=0.8
COV_HIGH=1.2
OUTGROUP_MIN=5
YIN_MIN=2
YANG_MIN=2
BATS_TOTAL_MIN=10

# === 物种 ===
OUTGROUP=(Mus_musculus_house_mouse Human Sorex_araneus_European_shrew Bos_taurus_cattle Canis_lupus_familiaris_dog Equus_caballus_horse)
YIN=(Cynopterus_sphinx_Indian_short_nosed_fruit_bat Rousettus_aegyptiacus_Egyptian_rousette Aselliscus_stoliczkanus_Stoliczkas_tridentbat Rhinolophus_ferrumequinum_greater_horseshoe_bat)
YANG=(Pteronotus_mesoamericanus_bats Desmodus_rotundus_common_vampire_bat Phyllostomus_discolor_pale_spear_nosed_bat Sturnira_hondurensis_bats Artibeus_lituratus_great_fruit_eatingbat Rhynchonycteris_naso_proboscis_bat Tadarida_brasiliensis_Brazilian_free_tailed_bat Antrozous_pallidus_pallid_bat Eptesicus_fuscus_big_brown_bat Myotis_myotis_bats Myotis_yumanensis_bats)

# 目录先创建
mkdir -p "$OUT_MAF_DIR" "$OUT_FASTA_DIR" \
         ./fasta_qualified ./fasta_outgroup_lowcov \
         ./fasta_noyin_lowquality ./fasta_noyang_lowquality \
         ./fasta_bats_total_low ./logs

# 将数组转成字符串导出（平行子进程可见）
export OUTGROUP_STR="${OUTGROUP[*]}"
export YIN_STR="${YIN[*]}"
export YANG_STR="${YANG[*]}"
printf 'OUTGROUP_STR=%s\nYIN_STR=%s\nYANG_STR=%s\n' "$OUTGROUP_STR" "$YIN_STR" "$YANG_STR" > ./logs/species_env.snapshot

# 写日志表头（和 Python 的列完全一致）
if [[ ! -s "$LOG_TSV" ]]; then
  echo -e "chrom\tstart\tend\tref_len\tok_outgroup\tok_yin\tok_yang\tok_bats_total\tcnt_outgroup_ok\tcnt_yin_ok\tcnt_yang_ok\tcnt_bats_ok\ttotal_species_in_fasta\tverdict\treason\tfasta_path\tcoverage_mode\tref_non_gap" > "$LOG_TSV"
fi

process_bed_line() {
  # 还原为本地数组
  set -f
  IFS=' ' read -r -a OUTGROUP_ARR <<< "${OUTGROUP_STR:-}"
  IFS=' ' read -r -a YIN_ARR       <<< "${YIN_STR:-}"
  IFS=' ' read -r -a YANG_ARR      <<< "${YANG_STR:-}"
  set +f

  local chrom="$1" start="$2" end="$3"
  local ref_len=$((end - start))
  local region="Human.${chrom}:${start}-${end}"
  local maf_out="${OUT_MAF_DIR}/${region}.maf"
  local fasta_out="${OUT_FASTA_DIR}/${region}.fasta"

  taffy view -i "$MAF_PATH" -r "$region" -m -o "$maf_out" 2>> ./logs/taffy.err

  # 验证 maf
  if [[ ! -s "$maf_out" ]] || ! grep -q '^##maf' "$maf_out"; then
    echo -e "${chrom}\t${start}\t${end}\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\tfail_maf_empty_or_noheader\tNA\t$maf_out\tNA\tNA" >> "$LOG_TSV"
    grep -m1 '^s Human' "$maf_out" >> "./logs/qc_stderr.${chrom}_${start}_${end}.log" 2>/dev/null || true
    return 14
  fi

  if ! msa_view "$maf_out" --missing-as-indels > "$fasta_out" 2>> "./logs/qc_stderr.${chrom}_${start}_${end}.log"; then
    echo -e "${chrom}\t${start}\t${end}\t${ref_len}\t0\t0\t0\t0\t0\t0\t0\t0\t0\tfail_msa_view\tNA\t$maf_out\tNA\tNA" >> "$LOG_TSV"
    return 15
  fi

  python3 qc_fasta.py \
    --fasta "$fasta_out" \
    --ref-len "$ref_len" \
    --cov-low "$COV_LOW" \
    --cov-high "$COV_HIGH" \
    --min-outgroup "$OUTGROUP_MIN" \
    --min-yin "$YIN_MIN" \
    --min-yang "$YANG_MIN" \
    --min-bats-total "$BATS_TOTAL_MIN" \
    --outgroup "${OUTGROUP_ARR[@]}" \
    --yin "${YIN_ARR[@]}" \
    --yang "${YANG_ARR[@]}" \
    --log-tsv "$LOG_TSV" \
    --coverage-mode species_insert_ok  \
    # --debug \
    1>> ./logs/qc_stdout.log \
    2>> "./logs/qc_stderr.${chrom}_${start}_${end}.log"

  local rc=$?
  case "$rc" in
    0)  mv -f "$fasta_out" ./fasta_qualified/ ;;
    10) mv -f "$fasta_out" ./fasta_outgroup_lowcov/ ;;
    11) mv -f "$fasta_out" ./fasta_noyin_lowquality/ ;;
    12) mv -f "$fasta_out" ./fasta_noyang_lowquality/ ;;
    13) mv -f "$fasta_out" ./fasta_bats_total_low/ ;;
    *)  mv -f "$fasta_out" ./fasta_bats_total_low/ ;;
  esac
}

# 仅在定义之后导出函数 & 变量
export -f process_bed_line
export MAF_PATH OUT_MAF_DIR OUT_FASTA_DIR LOG_TSV
export COV_LOW COV_HIGH OUTGROUP_MIN YIN_MIN YANG_MIN BATS_TOTAL_MIN
export OUTGROUP_STR YIN_STR YANG_STR

# 无匹配时不把字面量 *.bed 当成文件
shopt -s nullglob

# 并行跑；--env 确保变量在子进程可见
for bed in $BED_GLOB; do
  awk 'BEGIN{OFS="\t"} $1!~/^#/ && NF>=3 {print $1,$2,$3}' "$bed" | \
  parallel --env OUTGROUP_STR,YIN_STR,YANG_STR,MAF_PATH,OUT_MAF_DIR,OUT_FASTA_DIR,LOG_TSV,COV_LOW,COV_HIGH,OUTGROUP_MIN,YIN_MIN,YANG_MIN,BATS_TOTAL_MIN \
           -j 500 --colsep '\t' process_bed_line {1} {2} {3}
done
