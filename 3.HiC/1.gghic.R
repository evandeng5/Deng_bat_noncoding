# ==============================================================================
# 1. 加载依赖
# ==============================================================================
suppressPackageStartupMessages({
  library(gghic)
  library(ggplot2)
  library(dplyr)
  library(tibble)
  library(HiCExperiment)
  library(GenomicRanges)
  library(scales)
  library(rtracklayer)
  library(GenomeInfoDb)
})

# ==============================================================================
# 2. 输入文件
# ==============================================================================
cool_file_path <- "ERR13363410.allValidPairs.10kb.cool"
gtf_plot_file  <- "Myotis_nattereri_GCA_964212035.2.gtf"
tad_file_path  <- "ERR13363410_TADs_domains.bed"

target_gene_symbol <- "CAMK2D"   # 当前用 gene_id / gene_name 都可以，只要和下面的列一致
flank_size <- 300000

# ==============================================================================
# 3. 自动定位基因位置
# ==============================================================================
message("正在从 GTF 文件中查找目标基因位置...")

gtf_data <- import(gtf_plot_file, format = "gtf")
# ---- 新增：确保有 gene_symbol 这一列 ----
gtf_mcols <- mcols(gtf_data)

if (!"gene_symbol" %in% names(gtf_mcols)) {
  if ("gene_name" %in% names(gtf_mcols)) {
    gtf_mcols$gene_symbol <- gtf_mcols$gene_name
  } else if ("gene_id" %in% names(gtf_mcols)) {
    gtf_mcols$gene_symbol <- gtf_mcols$gene_id
  } else {
    stop(
      "GTF 中既没有 gene_symbol、gene_name 也没有 gene_id；",
      "无法为 gghic 生成基因符号列，请检查 GTF。"
    )
  }
  mcols(gtf_data) <- gtf_mcols
  
  # 写出一个带 gene_symbol 的临时 GTF 供 geom_annotation 使用
  gtf_plot_file <- sub("\\.gtf$", "_with_symbol.gtf", gtf_plot_file)
  export(gtf_data, gtf_plot_file, format = "gtf")
} else {
  gtf_plot_file <- gtf_plot_file
}
# 按需要改这里：gene_id / gene_name / symbol
target_gene_gr <- subset(
  gtf_data,
  type == "gene" & gene_id == target_gene_symbol
)

if (length(target_gene_gr) == 0) {
  stop(paste(
    "错误：在 GTF 文件中找不到基因", target_gene_symbol,
    "。请检查基因名或选择正确的列（gene_id / gene_name / symbol）。"
  ))
}

target_chrom <- as.character(seqnames(target_gene_gr))[1]
gene_start   <- start(target_gene_gr)[1]
gene_end     <- end(target_gene_gr)[1]

plot_start <- max(1, gene_start - flank_size)
plot_end   <- gene_end + flank_size

message(paste0("找到基因: ", target_gene_symbol))
message(paste0("位置: ", target_chrom, ":", gene_start, "-", gene_end))
message(paste0("绘图区域: ", target_chrom, ":", plot_start, "-", plot_end))

# ==============================================================================
# 4. 检查 Cool 文件染色体名称是否匹配
# ==============================================================================
cf <- HiCExperiment::CoolFile(cool_file_path)
avail_info  <- HiCExperiment::availableChromosomes(cf)
chrom_names <- as.character(GenomeInfoDb::seqnames(avail_info))

if (!target_chrom %in% chrom_names) {
  stop(paste0(
    "染色体名称不匹配！\n",
    "GTF 中的名称: ", target_chrom, "\n",
    "Cool 文件中的名称示例: ", paste(head(chrom_names), collapse = ", "), "\n",
    "请手动统一名称（例如在 GTF 或 cool 中加/去 'chr' 前缀）。"
  ))
}

# ==============================================================================
# 5. 读取 TAD
# ==============================================================================
message("正在读取 TAD 数据...")
df_tad <- read.table(tad_file_path, header = FALSE, stringsAsFactors = FALSE)
if (ncol(df_tad) >= 3) colnames(df_tad)[1:3] <- c("seqnames", "start", "end")

df_tad_filtered <- df_tad |>
  dplyr::filter(
    seqnames == target_chrom,
    end  >= plot_start,
    start <= plot_end
  )

# ==============================================================================
# 6. 构建 ChromatinContacts 并绘图（Hi-C + TAD + 基因轨道）
# ==============================================================================
message("正在构建 ChromatinContacts 对象并导入 Hi-C 数据...")

cc <- ChromatinContacts(
  cool_file_path,
  focus = paste0(target_chrom, ":", plot_start, "-", plot_end)
)
cc <- import(cc)


message("开始绘图（包含 RYR2 基因模型）...")

p <- gghic(cc) +
  # ---- TAD 轨迹 ----
geom_tad2(
  data = df_tad_filtered,
  aes(seqnames = seqnames, start = start, end = end),
  colour   = "#2c3e50",
  linetype = "dashed",
  size     = 1.0,
  alpha    = 0.8
) +
# ---- 基因注释轨道：只画目标基因 ----



geom_annotation(
  gtf_path     = gtf_plot_file,
  chrom_prefix = FALSE,    # 染色体名里没有 chr
  style        = "arrow",
  fill         = "#4682B4",
  colour       = "#4682B4",
  fontsize     = 8
) +
  ggtitle(paste0("Bat Genome (", target_chrom, ") - ", target_gene_symbol))

print(p)
