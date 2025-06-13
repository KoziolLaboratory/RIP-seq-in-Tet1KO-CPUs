library(exomePeak2)

setwd('...')

# Input 样本
input_bam <- c("./01.Data_AllBAM/BAM/TetCon_Input1.bam",
               "./01.Data_AllBAM/BAM/TetCon_Input2.bam",
               "./01.Data_AllBAM/BAM/TetCon_Input3.bam",
               "./01.Data_AllBAM/BAM/TetKO_Input1.bam",
               "./01.Data_AllBAM/BAM/TetKO_Input2.bam",
               "./01.Data_AllBAM/BAM/TetKO_Input3.bam")

# IP 样本
ip_bam <- c("./01.Data_AllBAM/BAM/TetCon_IP1.bam",
            "./01.Data_AllBAM/BAM/TetCon_IP2.bam",
            "./01.Data_AllBAM/BAM/TetCon_IP3.bam",
            "./01.Data_AllBAM/BAM/TetKO_IP1.bam",
            "./01.Data_AllBAM/BAM/TetKO_IP2.bam",
            "./01.Data_AllBAM/BAM/TetKO_IP3.bam")

# 变量名
name <- c("Con-1","Con-2","Con-3","KO-1","KO-2","KO-3")



# Build Guitar Coordinates
library(Guitar)
txdb <- makeTxDbFromGFF(file      = ".../Mus_musculus.GRCm38.102.gtf",
                        format    = "gtf",
                        dataSource= "Ensembl",
                        organism  = "Mus musculus")
# 注释文件
gtf = ".../Mus_musculus.GRCm38.102.gtf"

# 两组比较
result <- exomePeak2(   gff = gtf,
                       txdb = txdb,
            bam_ip        = ip_bam[1:3],
            bam_input     = input_bam[1:3],
            bam_ip_treated    = ip_bam[4:6],
            bam_input_treated = input_bam[4:6],
            save_dir          = "...exomepeak_Compare/",
            experiment_name   = 'Compare',
            strandness        = "1st_strand",
            parallel          = 1,
            p_cutoff          = 1e-10,
            fragment_length   = 150)