#!/bin/bash
#SBATCH -J Motif
#SBATCH -p q_cn
#SBATCH -o JobOut.Motif.%j
#SBATCH -n 1
#SBATCH -c 6
#SBATCH --mem=40g
scratch60/0.CodeBase/Tools_Homer/20250515_Liang_Tet1_RIP-seq/Run_Homer.sh
# conda activate Tools_Homer # 先在Terminal中激活环境，然后再提交任务

.../Tools_Homer/bin/perl                               \
.../Tools_Homer/bin/findMotifsGenome.pl                \
.../Tools_Homer/20250515_Liang_Tet1_RIP-seq/Common_Con_KO/formatted_Common_Con_KO.bed          \
.../Tools_Homer/data/genomes/mm10   ./Common_Con_KO/   -size 75   -len 8,10,12

.../Tools_Homer/bin/perl                               \
.../Tools_Homer/bin/findMotifsGenome.pl                \
.../Tools_Homer/20250515_Liang_Tet1_RIP-seq/Unique_Con/formatted_Unique_Con.bed          \
.../Tools_Homer/data/genomes/mm10   ./Unique_Con/   -size 75   -len 8,10,12

.../Tools_Homer/bin/perl                               \
.../Tools_Homer/bin/findMotifsGenome.pl                \
.../Tools_Homer/20250515_Liang_Tet1_RIP-seq/Unique_KO/formatted_Unique_KO.bed          \
.../Tools_Homer/data/genomes/mm10   ./Unique_KO/   -size 75   -len 8,10,12