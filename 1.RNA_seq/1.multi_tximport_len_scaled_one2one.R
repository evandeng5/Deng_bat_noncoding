#!/usr/bin/env Rscript
# Rscript multi_tximport_len_scaled_one2one.R --root_dir . --species_info species_info.tsv --orthogroups /home/xydeng/1.unfinished/1.unfinished/12.gene2counts/Orthogroups.tsv --keep_genes RYR1,RYR2,RYR3,CAMK2D,CAMK2B --out_prefix salmon_lenScaled_o2o
suppressPackageStartupMessages({
  library(optparse)
  library(tximport)
  library(readr)
  library(dplyr)
  library(stringr)
})

## ======================= 命令行参数 =======================

option_list <- list(
  make_option(c("-r", "--root_dir"), type="character",
              help="包含各物种目录的根目录（每个物种目录下有 quants/*/quant.sf）", metavar="DIR"),
  make_option(c("-s", "--species_info"), type="character",
              help="species_info.tsv，列至少包含: species, gff；可选: group"),
  make_option(c("-o", "--orthogroups"), type="character",
              help="Orthogroups.tsv，用于 one2one 判定"),
  make_option(c("-k", "--keep_genes"), type="character",
              default="RYR1,RYR2,RYR3",
              help="额外强制保留基因名，逗号分隔 [default: %default]"),
  make_option(c("-p", "--out_prefix"), type="character",
              default="one2one_lenScaled",
              help="输出前缀 [default: %default]"),
  make_option(c("--countsFromAbundance"), type="character",
              default="lengthScaledTPM",
              help="传给 tximport 的 countsFromAbundance 参数 [default: %default]")
)

opt <- parse_args(OptionParser(option_list=option_list))

if (is.null(opt$root_dir) || is.null(opt$species_info) || is.null(opt$orthogroups)) {
  cat("ERROR: 缺少必要参数\n",
      "示例:\n",
      "  Rscript multi_tximport_len_scaled_one2one.R \\\n",
      "    --root_dir 13.11_8_salmon_8_species \\\n",
      "    --species_info species_info.tsv \\\n",
      "    --orthogroups Orthogroups.tsv \\\n",
      "    --keep_genes RYR1,RYR2,RYR3,CAMK2D,CAMK2B \\\n",
      "    --out_prefix salmon_lenScaled_o2o\n")
  quit(status=1)
}

## ======================= 工具函数 =======================

# 统一基因名：去掉 gene- / gene_ 前缀，去掉首尾空格
normalize_gene <- function(x) {
  x <- trimws(x)
  x <- sub("^gene[-_]", "", x)
  x
}

# 从 attributes 列中按优先顺序提取属性（向量化）
get_attr <- function(x, keys) {
  res <- rep(NA_character_, length(x))
  for (k in keys) {
    # 支持 GFF (key=val) 或 GTF (key "val") 等格式
    m <- stringr::str_match(x, paste0(k, "[ =\"]([^;\" ]+)"))
    hit <- !is.na(m[, 2])
    res[hit & is.na(res)] <- m[hit, 2]
  }
  res
}

# 不依赖 GenomicFeatures，直接从 GFF/GTF 解析 tx2gene
make_tx2gene <- function(gff_path) {
  if (!file.exists(gff_path)) {
    stop("[ERR] GFF/GTF 文件不存在: ", gff_path)
  }
  message("[INFO] 从注释文件直接解析 tx2gene: ", gff_path)

  df <- readr::read_tsv(
    gff_path,
    comment = "#",
    col_names = FALSE,
    col_types = readr::cols(.default = "c")
  )

  if (ncol(df) < 9) {
    stop("[ERR] 注释文件列数不足(<9)，请检查是否为标准 GFF/GTF: ", gff_path)
  }

  # 优先用 transcript/mRNA 行
  tx_rows <- df[df$X3 %in% c("transcript", "mRNA"), ]

  # 某些注释没有 transcript 行，只在 exon 里有 transcript_id/gene_id
  if (nrow(tx_rows) == 0) {
    tx_rows <- df[df$X3 %in% c("exon", "CDS"), ]
  }

  if (nrow(tx_rows) == 0) {
    stop("[ERR] 未找到 transcript/mRNA/exon/CDS 行，无法构建 tx2gene: ", gff_path)
  }

  tx_id   <- get_attr(tx_rows$X9, c("ID", "transcript_id", "transcript", "Name"))
  gene_id <- get_attr(tx_rows$X9, c("gene_id", "Parent", "gene", "gene_name"))

  tx2gene <- data.frame(
    tx   = tx_id,
    gene = gene_id,
    stringsAsFactors = FALSE
  )

  tx2gene <- tx2gene[!is.na(tx2gene$tx) & !is.na(tx2gene$gene), ]
  tx2gene <- unique(tx2gene)
  message("[INFO] tx2gene 行数: ", nrow(tx2gene))

  # 去重 + 排除空白
  tx2gene <- tx2gene[tx2gene$tx != "", ]
  tx2gene <- tx2gene[!duplicated(tx2gene$tx), ]

  message("[INFO] tx2gene 行数 (去重后): ", nrow(tx2gene))

  if (nrow(tx2gene) == 0) {
    stop("[ERR] 未能从注释构建任何 tx2gene 映射，请检查 ", gff_path)
  }

  
  tx2gene
}

## ======================= 读取 species_info =======================

sp_info <- readr::read_tsv(opt$species_info, col_types = cols())

required_cols <- c("species", "gff")
if (!all(required_cols %in% colnames(sp_info))) {
  stop("species_info.tsv 必须包含列: species, gff；可选列: group")
}

if (!"group" %in% colnames(sp_info)) {
  sp_info$group <- sp_info$species
}

current_species <- unique(sp_info$species)
message("[INFO] Species (", length(current_species), "): ",
        paste(current_species, collapse = ", "))

## ======================= 读取 Orthogroups，判定 one2one =======================

message("[INFO] 读取 Orthogroups 用于 one2one 判定: ", opt$orthogroups)

og_con <- file(opt$orthogroups, "r")
on.exit(close(og_con), add = TRUE)

header <- strsplit(readLines(og_con, n = 1), "\t")[[1]]
if (length(header) < 2) {
  stop("[ERR] Orthogroups.tsv 表头格式异常")
}

# 与 species_info 重叠的物种列
sp_cols_idx <- which(header %in% current_species)
sp_cols_names <- header[sp_cols_idx]

if (length(sp_cols_idx) == 0) {
  stop("[ERR] Orthogroups.tsv 中找不到任何与 species_info.tsv 一致的物种列")
}

missing_in_og <- setdiff(current_species, sp_cols_names)
if (length(missing_in_og) > 0) {
  message("[WARN] 下列物种在 Orthogroups.tsv 中无对应列，将不用于 one2one 判定: ",
          paste(missing_in_og, collapse = ", "))
}

og2genes <- list()  # og -> list(species -> genes[标准名])
gene2og  <- list()  # 标准基因名 -> og 集合

while (length(line <- readLines(og_con, n = 1)) > 0) {
  fields <- strsplit(line, "\t")[[1]]
  if (length(fields) < 1) next

  og <- strsplit(fields[1], ":")[[1]][1]
  mp <- list()

  ok <- TRUE
  for (j in seq_along(sp_cols_idx)) {
    idx <- sp_cols_idx[j]
    sp  <- sp_cols_names[j]
    if (idx > length(fields)) {
      ok <- FALSE
      break
    }
    cell <- trimws(fields[idx])
    if (nzchar(cell)) {
      genes <- strsplit(cell, ",")[[1]]
      genes <- trimws(genes[nzchar(genes)])
      if (length(genes) > 0) {
        genes <- normalize_gene(genes)  # 关键：去掉 gene- 等前缀
        mp[[sp]] <- genes
      }
    }
  }

  if (!ok || length(mp) == 0) next

  og2genes[[og]] <- mp

  # 建立 gene(标准名) -> og 映射
  for (gs in unlist(mp, use.names = FALSE)) {
    if (!nzchar(gs)) next
    if (is.null(gene2og[[gs]])) {
      gene2og[[gs]] <- og
    } else {
      gene2og[[gs]] <- unique(c(gene2og[[gs]], og))
    }
  }
}

# 判定 one2one：
# 对于 Orthogroups 中包含的这些 sp_cols_names 物种，每个 OG 要求每个物种正好 1 个基因
one2one_ogs <- character(0)
for (og in names(og2genes)) {
  mp <- og2genes[[og]]
  if (!all(sp_cols_names %in% names(mp))) next
  if (all(vapply(mp[sp_cols_names], length, integer(1)) == 1L)) {
    one2one_ogs <- c(one2one_ogs, og)
  }
}
one2one_ogs <- unique(one2one_ogs)

one2one_genes <- unique(unlist(
  lapply(one2one_ogs, function(og) unlist(og2genes[[og]], use.names = FALSE)),
  use.names = FALSE
))

keep_extra <- unique(trimws(unlist(strsplit(opt$keep_genes, ","))))
keep_extra <- keep_extra[nzchar(keep_extra)]
keep_extra_norm <- normalize_gene(keep_extra)

message("[INFO] one2one OGs = ", length(one2one_ogs))
message("[INFO] one2one genes (normalized) = ", length(one2one_genes))
message("[INFO] Extra keep genes (normalized) = ",
        paste(keep_extra_norm, collapse = ", "))

## ======================= tximport：各物种汇总 =======================

all_counts <- NULL
coldata <- data.frame()

for (i in seq_len(nrow(sp_info))) {
  sp   <- sp_info$species[i]
  grp  <- sp_info$group[i]
  gff  <- sp_info$gff[i]

  sp_quants_dir <- file.path(opt$root_dir, sp, "quants")
  if (!dir.exists(sp_quants_dir)) {
    stop("[ERR] 目录不存在: ", sp_quants_dir)
  }

  sample_dirs <- list.dirs(sp_quants_dir, full.names = TRUE, recursive = FALSE)
  files <- file.path(sample_dirs, "quant.sf")
  exists_mask <- file.exists(files)

  if (!any(exists_mask)) {
    stop("[ERR] 在 ", sp_quants_dir, " 未找到任何 quant.sf")
  }

  files <- files[exists_mask]
  sample_names_raw <- basename(dirname(files))
  sample_names <- paste(sp, sample_names_raw, sep = "_")  # 防止不同物种同名样本
  names(files) <- sample_names

  message("[INFO] 物种 ", sp, ": 样本数 = ", length(files))

  tx2gene <- make_tx2gene(gff)

  ## -------- 对齐 quant.sf 与 tx2gene 中的转录本 ID --------
  # 读取一个 quant.sf，用于检查 ID 是否匹配
  example_quant <- files[1]
  qhead <- readr::read_tsv(example_quant, n_max = 200, col_types = cols())
  salmon_tx <- qhead$Name

  # 1) 直接匹配
  direct_overlap <- length(intersect(salmon_tx, tx2gene$tx))

  if (direct_overlap == 0) {
    # 2) 尝试只去掉版本号匹配（两边都去）
    salmon_core <- sub("\\.[0-9]+$", "", salmon_tx)
    tx_core     <- sub("\\.[0-9]+$", "", tx2gene$tx)
    core_overlap <- length(intersect(salmon_core, tx_core))

    if (core_overlap > 0) {
      message("[INFO] 使用去版本号的转录本 ID 对齐: ", sp)
      tx2gene$tx <- tx_core
    } else {
      stop("[ERR] ", sp, ": quant.sf 中的 Name 与 tx2gene 完全不匹配。\n",
           "  请检查该物种 Salmon index 使用的转录本 fasta 是否和 gff 同源。")
    }
  }

  # 最终 sanity check
  if (length(intersect(salmon_tx, tx2gene$tx)) == 0) {
    stop("[ERR] ", sp, ": 对齐后仍然没有任何转录本 ID 重叠，终止。")
  }

  ## -------- 安全 tximport，失败时做详细诊断 --------
  txi <- tryCatch({
    tximport(
      files,
      type = "salmon",
      tx2gene = tx2gene,
      countsFromAbundance = opt$countsFromAbundance,
      dropInfReps = TRUE
    )
  }, error = function(e) {
    message("\n[ERR] tximport 在物种 ", sp, " 失败，开始诊断...\n")

    # 1) 检查每个 quant.sf 的 Name 列长度 & 是否有重复
    for (f in files) {
      q <- try(readr::read_tsv(f, show_col_types = FALSE), silent = TRUE)
      if (inherits(q, "try-error")) {
        message("[ERR] 无法读取 quant.sf: ", f)
      } else {
        message("[DBG] 文件: ", f,
                " | 转录本数: ", nrow(q),
                " | 唯一转录本数: ", length(unique(q$Name)),
                " | 是否有重复: ", ifelse(anyDuplicated(q$Name) > 0, "YES", "NO"))
      }
    }

    # 2) 和第一个样本的 Name 做对比，看看谁不一致
    q0 <- readr::read_tsv(files[1], show_col_types = FALSE)
    base_ids <- q0$Name
    for (f in files[-1]) {
      qf <- readr::read_tsv(f, show_col_types = FALSE)
      if (!identical(base_ids, qf$Name)) {
        message("[DBG] 样本(ID顺序/集合不一致): ", f)
        message("      base_ids 长度 = ", length(base_ids),
                ", 该样本长度 = ", length(qf$Name))
        message("      仅在 base_ids 中的前5个: ",
                paste(head(setdiff(base_ids, qf$Name), 5), collapse = ", "))
        message("      仅在该样本中的前5个: ",
                paste(head(setdiff(qf$Name, base_ids), 5), collapse = ", "))
      }
    }

    # 3) 检查和 tx2gene 的交集情况
    q_all <- unique(unlist(lapply(files, function(f) {
      readr::read_tsv(f, show_col_types = FALSE)$Name
    })))
    inter_n <- length(intersect(q_all, tx2gene$tx))

    message("[DBG] tx2gene$tx 数量: ", nrow(tx2gene))
    message("[DBG] 所有 quant.sf 唯一转录本 ID 数量: ", length(q_all))
    message("[DBG] 二者交集数量: ", inter_n)

    if (inter_n > 0) {
      message("[DBG] 交集示例: ",
              paste(head(intersect(q_all, tx2gene$tx), 10), collapse = ", "))
    } else {
      message("[DBG] 完全无交集，说明 GFF 与 Salmon index 可能不匹配，或 tx2gene 构建规则有问题。")
    }

    # 最后真正抛出原始错误，防止脚本悄悄继续
    stop(e)
  })



  sp_counts <- as.data.frame(txi$counts)
  sp_counts$gene_id <- rownames(sp_counts)
  rownames(sp_counts) <- NULL
  sp_counts <- sp_counts[, c("gene_id", setdiff(colnames(sp_counts), "gene_id"))]

  # 保存 TPM 矩阵（txi$abundance）
  sp_tpm <- as.data.frame(txi$abundance)
  sp_tpm$gene_id <- rownames(sp_tpm)
  rownames(sp_tpm) <- NULL
  sp_tpm <- sp_tpm[, c("gene_id", setdiff(colnames(sp_tpm), "gene_id"))]

  # 合并各物种 TPM 矩阵
  if (exists("all_tpm")) {
    all_tpm <- full_join(all_tpm, sp_tpm, by = "gene_id")
  } else {
    all_tpm <- sp_tpm
  }



  if (is.null(all_counts)) {
    all_counts <- sp_counts
  } else {
    all_counts <- full_join(all_counts, sp_counts, by = "gene_id")
  }

  coldata <- rbind(
    coldata,
    data.frame(
      sample  = sample_names,
      species = sp,
      group   = grp,
      stringsAsFactors = FALSE
    )
  )
}

# NA 补 0
all_counts[is.na(all_counts)] <- 0
count_cols <- setdiff(colnames(all_counts), "gene_id")

## ======================= gene_meta & one2one 子矩阵 =======================

# 从 gene_id 中提取用于匹配的部分，再 normalize
# 如果 gene_id 中有 '|', 取最后一段，否则就是自身
symbol_raw  <- ifelse(grepl("\\|", all_counts$gene_id),
                      sub(".*\\|", "", all_counts$gene_id),
                      all_counts$gene_id)
symbol_norm <- normalize_gene(symbol_raw)

pick_one_og <- function(gs) {
  ogs <- gene2og[[gs]]
  if (is.null(ogs) || length(ogs) == 0) return(NA_character_)
  sort(ogs)[1]
}

gene_meta <- data.frame(
  gene_id    = all_counts$gene_id,
  symbol     = symbol_norm,
  is_one2one = symbol_norm %in% one2one_genes,
  orthogroup = vapply(symbol_norm, pick_one_og, character(1)),
  stringsAsFactors = FALSE
)

# one2one + 额外保留基因
keep_set <- unique(c(one2one_genes, keep_extra_norm))
sub_idx <- symbol_norm %in% keep_set
counts_one2one <- all_counts[sub_idx, , drop = FALSE]

## ======================= TPM 过滤低表达基因 =======================

# all_tpm: gene_id + 所有样本的 TPM
# counts_one2one: 已经筛过 one2one(+keep) 的 counts

# 对齐 one2one 基因的 TPM 行
tpm_mat <- all_tpm[match(counts_one2one$gene_id, all_tpm$gene_id), , drop = FALSE]

# 防御：去掉 gene_id 列，只留数值矩阵
tpm_mat_vals <- as.matrix(tpm_mat[, -1, drop = FALSE])
rownames(tpm_mat_vals) <- counts_one2one$gene_id

# 阈值设置
tpm_threshold <- 1      # TPM > 1
min_samples   <- 20      # 至少 3 个 SRR 样本满足条件（你可以改成 10 等）

# 每个基因在多少个样本中 TPM > 阈值
expr_count <- rowSums(tpm_mat_vals > tpm_threshold, na.rm = TRUE)

# 保留满足条件的基因
keep_filter <- expr_count >= min_samples

counts_filtered <- counts_one2one[keep_filter, , drop = FALSE]

# 输出过滤后矩阵
filtered_out <- paste0(opt$out_prefix, "_counts_one2one_TPM1.tsv")
write_tsv(counts_filtered, filtered_out)
message("[OK] 低表达过滤后矩阵: ", filtered_out,
        " (genes=", nrow(counts_filtered), 
        ", TPM>", tpm_threshold, " in ≥", min_samples, " samples)")


# 计算每个物种的平均 TPM
tpm_mat_df <- as.data.frame(tpm_mat)
coldata_df <- as.data.frame(coldata)

tpm_species_mean <- sapply(unique(coldata_df$species), function(sp) {
  sp_samples <- coldata_df$sample[coldata_df$species == sp]
  valid_samples <- intersect(sp_samples, colnames(tpm_mat_df))
  if (length(valid_samples) > 0) {
    rowMeans(tpm_mat_df[, valid_samples, drop = FALSE], na.rm = TRUE)
  } else {
    rep(NA_real_, nrow(tpm_mat_df))
  }
})
tpm_species_mean <- as.data.frame(tpm_species_mean)
tpm_species_mean$gene_id <- counts_one2one$gene_id
write_tsv(tpm_species_mean, paste0(opt$out_prefix, "_TPM_species_mean.tsv"))
message("[OK] 每物种平均 TPM 输出: ", opt$out_prefix, "_TPM_species_mean.tsv")

## ======================= 输出 =======================

counts_all_out <- paste0(opt$out_prefix, "_counts_all.tsv")
coldata_out    <- paste0(opt$out_prefix, "_coldata.tsv")
gene_meta_out  <- paste0(opt$out_prefix, "_gene_meta.tsv")
counts_o2o_out <- paste0(opt$out_prefix, "_counts_one2one.tsv")

write_tsv(all_counts, counts_all_out)
write_tsv(coldata,    coldata_out)
write_tsv(gene_meta,  gene_meta_out)
write_tsv(counts_one2one, counts_o2o_out)

# 输出全物种 TPM 矩阵
tpm_all_out <- paste0(opt$out_prefix, "_tpm_all.tsv")
write_tsv(all_tpm, tpm_all_out)
message("[OK] 全集 TPM 矩阵: ", tpm_all_out,
        " (genes=", nrow(all_tpm), ", samples=", length(colnames(all_tpm)) - 1, ")")

message("[OK] 全集矩阵: ", counts_all_out,
        " (genes=", nrow(all_counts), ", samples=", length(count_cols), ")")
message("[OK] coldata: ", coldata_out,
        " (rows=", nrow(coldata), ")")
message("[OK] gene_meta: ", gene_meta_out)
message("[OK] one2one(+keep) 子矩阵: ", counts_o2o_out,
        " (genes=", nrow(counts_one2one), ")")
