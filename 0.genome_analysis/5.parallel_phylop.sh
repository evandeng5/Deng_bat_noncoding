find . -maxdepth 1 -name "*.bed" | parallel '
    # {} 是完整的输入文件名，例如 ./file1.bed
    # {.} 是不带后缀的文件名，例如 ./file1
    # {/.} 是不带路径和后缀的文件名，例如 file1

    INPUT_BED={}; 
    BASENAME={/.}; # 移除路径和.bed后缀
    
    phyloP --method LRT --mode ACC \
           --features $INPUT_BED \
           --branch Antrozous_pallidus_pallid_bat-Aselliscus_stoliczkanus_Stoliczkas_tridentbat \
           --msa-format MAF \
           --gff-scores /home/dxy/2024_bat_comparative_genome/1.analysis/1.accelerated_detect/nonconserved-all-4d_21sp.mod \
           /home/dxy/2024_bat_comparative_genome/1.analysis/1.accelerated_detect/1.maf/21_sp_maf/${BASENAME}.maf \
           > ${BASENAME}.gff
'