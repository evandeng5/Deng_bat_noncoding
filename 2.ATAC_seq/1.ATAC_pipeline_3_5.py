import os
import argparse
import subprocess
from concurrent.futures import ThreadPoolExecutor
import logging
import time


# 配置日志格式
logging.basicConfig(
    format="[%(asctime)s] %(message)s",  # 时间 + 日志内容
    datefmt="%Y-%m-%d %H:%M",  # 时间格式
    level=logging.INFO  # 设定日志级别
)

def log_message(message):
    logging.info(message)



def run_command(command):
    """运行 shell 命令"""
    log_message(f"Running command: {command}")
    try:
        subprocess.run(command, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        log_message(f"Error: {e}")

def process_sample(accession_id, qc_tool, reference, prefix, outputdir,hal,species,haltarget):
    # 创建样本目录，避免线程争抢
    sample_dir = os.path.join(".", accession_id)
    os.makedirs(sample_dir, exist_ok=True)
    
    # 下载 SRA 文件
    log_message(f"1.Fetching {accession_id}...")
    run_command(f"prefetch  {accession_id} -p  --max-size 50G")
    sra_file = os.path.join(sample_dir, f"{accession_id}.sra")
    fastq_r1 = os.path.join(sample_dir, f"{accession_id}_val_1.fq")
    fastq_r2 = os.path.join(sample_dir, f"{accession_id}_val_2.fq")

    # 检查 SRA 文件是否成功下载
    if not os.path.isfile(sra_file):
        log_message(f"Error: Prefetch failed for {accession_id}, skipping conversion.")
        return
    
    # 转换 SRA 为 FASTQ
    log_message(f"2.Converting {accession_id} to FASTQ...")
    run_command(f"fasterq-dump --split-files -O {sample_dir} {sra_file}")
    
    # 质控处理
    if os.path.isfile(fastq_r1) and os.path.isfile(fastq_r2):
        if qc_tool == "fastp":
            log_message(f"Running double-end fastp for {accession_id}...")
            run_command(
                f"fastp -i {fastq_r1} -I {fastq_r2} "
                f"-o {sample_dir}/{prefix}_val_1.fq -O {sample_dir}/{prefix}_val_2.fq "
                f"-l 60 -w 10 --html {sample_dir}/{accession_id}_fastp_report.html "
                f"--json {sample_dir}/{accession_id}_fastp_report.json"
            )
        elif qc_tool == "trim_galore":
            log_message(f"Running trim_galore for {accession_id}...")
            run_command(f"trim_galore --fastqc -j 4 -q 25 --phred33 --basename {prefix} --length 35 --paired -o {sample_dir} {fastq_r1} {fastq_r2}")
        log_message(f"{qc_tool} processing completed for {accession_id}")
    elif os.path.isfile(fastq_r1):
        run_command(f"trim_galore --fastqc --length 25 --basename {prefix} --quality 25 {fastq_r1} -j 4 -o {sample_dir} {fastq_r1}")
        log_message(f"Running single-end fastp for {accession_id}...")
        
        
    #bwa比对    
    log_message("3.bwa aligning...")
    # log_message(f"Indexing reference genome: {reference}")
    # run_command(f"bwa index {reference}")
    run_command(f"bwa mem -t 8  {reference} {sample_dir}/{prefix}_val_1.fq {sample_dir}/{prefix}_val_2.fq   >  {sample_dir}/{prefix}_ATAC_seq.sam ")
    #samtools转格式
    log_message("4.Running samtools...")
    run_command(f"samtools view -bS {sample_dir}/{prefix}_ATAC_seq.sam > {sample_dir}/{prefix}_ATAC_seq.bam")
    run_command(f"samtools sort {sample_dir}/{prefix}_ATAC_seq.bam -o {sample_dir}/{prefix}_ATAC_seq_sorted.bam -@ 8")
    
    #picard去重复
    log_message("5.Running picard remove duplicates... ")
    test_RG = subprocess.run(f"samtools view -H {sample_dir}/{prefix}_ATAC_seq_sorted.bam | grep '@RG'", 
                            shell=True, capture_output=True, text=True)
    if test_RG.stdout.strip():  # 如果有输出，则表示 @RG 存在
        log_message(f"Read group found in {sample_dir}/{prefix}_ATAC_seq_sorted.bam")
        run_command(f"""java -jar /home/xydeng/software/biosoftware/picard/build/libs/picard.jar MarkDuplicates \
                        I={sample_dir}/{prefix}_ATAC_seq_sorted.bam \
                        O={sample_dir}/{prefix}_ATAC_seq_sorted_RG_dedup.bam \
                        M={sample_dir}/{prefix}_ATAC_seq_sorted_RG_dedup.txt \
                        REMOVE_DUPLICATES=true """)
    else:
        log_message(f"No read group found, adding read group first...")
        run_command(f"""java -jar /home/xydeng/software/biosoftware/picard/build/libs/picard.jar AddOrReplaceReadGroups \
                        I={sample_dir}/{prefix}_ATAC_seq_sorted.bam \
                        O={sample_dir}/{prefix}_ATAC_seq_sorted_RG.bam \
                        RGID=2 RGLB=lib2 RGPL=illumina RGPU=unit2 RGSM=sample2""")
            
        run_command(f"""java -jar /home/xydeng/software/biosoftware/picard/build/libs/picard.jar MarkDuplicates \
                        I={sample_dir}/{prefix}_ATAC_seq_sorted_RG.bam \
                        O={sample_dir}/{prefix}_ATAC_seq_sorted_RG_dedup.bam \
                        M={sample_dir}/{prefix}_ATAC_seq_sorted_RG_dedup.txt \
                        REMOVE_DUPLICATES=true """)

    
    #ATAC偏移校正
    log_message(f"6.Correcting ATAC shift")    
    run_command(f"samtools index {sample_dir}/{prefix}_ATAC_seq_sorted_RG_dedup.bam")
    
    run_command(f"""alignmentSieve --numberOfProcessors 15 \
                --ATACshift  \
                --bam {sample_dir}/{prefix}_ATAC_seq_sorted_RG_dedup.bam \
                -o {sample_dir}/{prefix}_ATAC_seq_sorted_RG_dedup_shifted.bam """)
    
    
    # macs3 peak calling
    log_message(f"7.Running peak calling...")
    run_command(f"faCount {reference} -summary > genome.faCount.report")
    cmd = "awk 'NR==2 {print $3+$4+$5+$6}' genome.faCount.report"
    genome_size = subprocess.check_output(cmd, shell=True, text=True).strip()
    log_message(f"Successfully obtained genome size: {genome_size}")

    
    run_command(f"""macs3 callpeak \
                -t {sample_dir}/{prefix}_ATAC_seq_sorted_RG_dedup_shifted.bam \
                -f BAM \
                -g {genome_size} \
                --nomodel \
                --shift \
                -100 \
                --extsize 200 \
                --keep-dup all\
                -n {prefix} \
                --outdir {sample_dir}/{outputdir} \
                -B""")
    
    if hal and species:
        log_message(f"8.Running halLiftover...")
        run_command(f"""halLiftover --bedType 5 \
                    {hal} \
                    {species}   \
                    {sample_dir}/{outputdir}/{prefix}_peaks.narrowPeak \
                    {haltarget} {sample_dir}/{outputdir}/Human_show_{prefix}_peaks.narrowPeak""")
        run_command(f"""halLiftover --bedType 4 \
                    {hal} \
                    {species}   \
                    {sample_dir}/{outputdir}/{prefix}_treat_pileup.bdg \
                    {haltarget} {sample_dir}/{outputdir}/Human_show_{prefix}_treat_pileup.bdg""")

    else:
        log_message("hal file not provided, skipping halLiftover step.")
    
    log_message(f"Processing completed for {accession_id}")
    

    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="ATAC-seq pipeline, from SRA to peak calling" )
    parser.add_argument("-c", "--qc_tool", choices=["fastp", "trim_galore"], default="fastp", help="Select quality control tools: fastp 或 trim_galore")
    parser.add_argument("-t", "--threads", type=int, default=10, help="Number of threads for parallel processing")
    parser.add_argument("-i", "--input", nargs="+", required=True, help="One or more accession IDs input, separated by spaces")
    parser.add_argument("-r", "--reference", required=True, help="Reference genome file path")
    parser.add_argument("-p","--prefix",required=True,help="Output file prefix")
    parser.add_argument("-o","--outputdir",default="output",help="Peakfile Output folder")
    parser.add_argument("-a","--alignment",help="Whole genome alignment file")
    parser.add_argument("-s","--species",help="Species coordinates to be converted")
    parser.add_argument("-ht","--haltarget",default="Human",help="Reference species name")
    args = parser.parse_args()
    # 检查依赖关系
    if args.alignment and not args.species:
        parser.error("The -a/--alignment option has been filled in, but -s/--species and --haltarget reference cannot be empty, please provide species conversion information。")
    
    start_time = time.time()    
    log_message("🚀🚀🚀 Processing ATAC-seq pipeline ...")
    run_command(f"bwa index {args.reference} ")
    
    with ThreadPoolExecutor(max_workers=args.threads) as executor:
        for accession_id in args.input:
            executor.submit(process_sample, accession_id, args.qc_tool, args.reference, args.prefix,args.outputdir,args.alignment,args.species,args.haltarget)
    
    end_time = time.time()
    elapsed_time = end_time - start_time
    log_message(f"All tasks completed successfully! ⏳ Total runtime: {elapsed_time:.2f} seconds.")
