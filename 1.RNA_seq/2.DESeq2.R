#!/usr/bin/env Rscript
#Rscript DESeq2.R --counts 11_8_ChangeDog_counts_one2one.tsv  --coldata 11_8_ChangeDog_coldata.tsv  --design "~ group" --contrast "group:gliding,non_gliding" --outdir DESeq2_RVR_ctrl_vs_case
suppressPackageStartupMessages({
  library(DESeq2)
  library(readr)
  library(dplyr)
  library(ggplot2)
})

# ---------- tiny arg parser ----------
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL, required = FALSE) {
  if (flag %in% args) {
    i <- which(args == flag)[1]
    if (i == length(args)) stop("Missing value for ", flag)
    return(args[i + 1])
  } else {
    if (required) stop("Missing required flag: ", flag)
    return(default)
  }
}

counts_path <- get_arg("--counts", required = TRUE)
coldata_path <- get_arg("--coldata", required = TRUE)
design_str   <- get_arg("--design", default = "~ group")
contrast_str <- get_arg("--contrast", default = NA)  # e.g. "group:ctrl,case"
outdir       <- get_arg("--outdir", default = "DESeq2_out")

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

message("== Inputs ==")
message("counts  : ", counts_path)
message("coldata : ", coldata_path)
message("design  : ", design_str)
message("contrast: ", ifelse(is.na(contrast_str), "(auto if 2 levels)", contrast_str))
message("outdir  : ", outdir)

# ---------- read inputs ----------
counts_df <- read_tsv(counts_path, show_col_types = FALSE)
id_col <- colnames(counts_df)[1]  # 支持 gene_id 或 orthogroup

# make rownames
row_ids <- make.unique(counts_df[[id_col]])
cts <- counts_df[,-1, drop = FALSE]
cts <- as.matrix(cts)
mode(cts) <- "numeric"
cts[is.na(cts)] <- 0
cts <- round(cts)
rownames(cts) <- row_ids    # ← 关键

if (!all(sapply(counts_df[-1], is.numeric))) {
  message("Coercing counts to numeric (round to integer).")
}


coldata <- read_tsv(coldata_path, show_col_types = FALSE)
stopifnot("sample" %in% colnames(coldata))
rownames(coldata) <- coldata$sample

# 检查列名一致
if (!all(colnames(cts) %in% rownames(coldata))) {
  missing <- setdiff(colnames(cts), rownames(coldata))
  stop("coldata 里缺少以下样本：", paste(missing, collapse = ", "))
}
# coldata 重新按 counts 列顺序排列
coldata <- coldata[colnames(cts), , drop = FALSE]

# ---------- build DESeq2 object ----------
# 设计公式
design_formula <- as.formula(design_str)

dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData   = coldata,
                              design    = design_formula)

# 推荐 size factor 类型（避免全 0 的基因影响）
dds <- estimateSizeFactors(dds, type = "poscounts")

# 预过滤（可选：极低表达基因）
keep <- rowSums(counts(dds) >= 5) >= 2
dds <- dds[keep, ]

# 拟合
dds <- DESeq(dds)

# ---------- results / contrasts ----------
# 解析 contrast，如 "group:ctrl,case" 表示 case vs ctrl
do_results <- TRUE
if (is.na(contrast_str)) {
  # 若未指定，对使用设计中的第一个因子列尝试自动二元对比
  # 从 design 提取变量名
  vars <- all.vars(design_formula)
  vars <- setdiff(vars, "(Intercept)")
  if (length(vars) == 0) {
    message("No variable in design; skip results.")
    do_results <- FALSE
  } else {
    fac <- vars[1]
    if (!is.factor(coldata[[fac]])) {
      coldata[[fac]] <- factor(coldata[[fac]])
      colData(dds)[[fac]] <- coldata[[fac]]
    }
    levs <- levels(coldata[[fac]])
    if (length(levs) == 2) {
      contrast <- c(fac, levs[1], levs[2]) # 默认第二个对第一个：lev2 vs lev1？→ DESeq2 需要 c(factor, ref, target) 表示 target/ref 的 log2FC，文档是 c(factor, levelN, levelD)
      # 注意：DESeq2 的 contrast 定义为 c(factor, numeratorLevel, denominatorLevel)
      # 我们按常见书写：contrast = c(fac, levs[2], levs[1]) => log2(lev2/lev1)
      contrast <- c(fac, levs[2], levs[1])
      message("Auto contrast = ", paste(contrast, collapse = " "))
    } else {
      message("Design factor '", fac, "' has ", length(levs), " levels; please pass --contrast to specify.")
      do_results <- FALSE
    }
  }
} else {
  parts <- strsplit(contrast_str, ":", fixed = TRUE)[[1]]
  stopifnot(length(parts) == 2)
  fac <- parts[1]
  levs2 <- strsplit(parts[2], ",", fixed = TRUE)[[1]]
  stopifnot(length(levs2) == 2)
  # DESeq2 expects c(factor, numeratorLevel, denominatorLevel) = log2(numerator/denominator)
  contrast <- c(fac, levs2[2], levs2[1])  # 让 "A,B" 表示 B vs A（与上一步脚本保持直觉：后者对前者）
  if (!fac %in% colnames(coldata)) stop("contrast 因子不在 coldata 中：", fac)
  if (!all(levs2 %in% levels(factor(coldata[[fac]]))))
    stop("contrast 的水平不在 coldata$", fac, " 中：", paste(levs2, collapse = ", "))
  message("User contrast = ", paste(contrast, collapse = " "))
}

if (do_results) {
  res <- results(dds, contrast = contrast)
  # 如果有 apeglm，做 LFC 收缩（对二元对比有效）
  shrunk <- tryCatch({
    suppressPackageStartupMessages(library(apeglm))
    lfcShrink(dds, contrast = contrast, type = "apeglm")
  }, error = function(e) {
    message("[warn] apeglm shrinkage failed: ", e$message, " → use unshrunken results")
    res
  })
  res_tbl <- as.data.frame(shrunk) %>%
    tibble::rownames_to_column(id_col) %>%
    arrange(padj)
  out_res <- file.path(outdir, "DESeq2_results_shrunk.csv")
  write.csv(res_tbl, out_res, row.names = FALSE)
  message("✔ results written: ", out_res)
} else {
  message("Skip results (no valid contrast).")
}

# ---------- PCA & exports ----------
vsd <- vst(dds, blind = FALSE)
p <- plotPCA(vsd, intgroup = intersect(c("group","species"), colnames(coldata)))
ggsave(file.path(outdir, "PCA.png"), plot = p, width = 6, height = 5, dpi = 200)

# 归一化计数与 size factor
sf <- sizeFactors(dds)
write_tsv(tibble(sample = names(sf), sizeFactor = sf),
          file.path(outdir, "size_factors.tsv"))

nc <- counts(dds, normalized = TRUE)
nc_out <- data.frame(id = rownames(nc), nc, check.names = FALSE)
names(nc_out)[1] <- id_col
write_tsv(nc_out, file.path(outdir, "normalized_counts.tsv"))

# VST 矩阵
vst_mat <- assay(vsd)
vst_out <- data.frame(id = rownames(vst_mat), vst_mat, check.names = FALSE)
names(vst_out)[1] <- id_col
write_tsv(vst_out, file.path(outdir, "vst_matrix.tsv"))

message("Done. Outputs in: ", outdir)
