#!/bin/bash
#SBATCH -J Genome
#SBATCH -p q_fat_2
#SBATCH -o JobOut.Genome.%j
#SBATCH -n 1
#SBATCH -c 6
#SBATCH --mem=256g

##############
### Genome ###
##############
for i in "${!read1_array[@]}"
do
hisat2 -p 6 --dta -x ./Mouse/GRCm38/release-102/Mus_musculus.GRCm38.dna -1 ${read1_array[$i]} -2 ${read2_array[$i]} -S "${read1_array[$i]}.sam"
done

###############
### SAM2BAM ###
###############
# Convert SAM file to BAM file
# Sort the BAM file
# Index the sorted BAM file
for i in "${!read_array[@]}"
do
  samtools view -S -b "${read_array[$i]}" > "${read_array[$i]%.sam}.bam"
  samtools sort "${read_array[$i]%.sam}.bam" -o "${read_array[$i]%.sam}_sorted.bam"
  samtools index "${read_array[$i]%.sam}_sorted.bam"
done