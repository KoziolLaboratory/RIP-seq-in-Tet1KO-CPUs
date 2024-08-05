#!/bin/bash
#SBATCH -J Preprocess
#SBATCH -p q_cn
#SBATCH -o JobOut.Preprocess.%j
#SBATCH -n 1
#SBATCH --mem=128g

##########
### QC ###
##########
module load sratoolkit
module load fastqc
module load fastp
fastqc -t 12 ./*.fq.gz
module load multiqc
multiqc .

###############
### Adapter ###
###############
# List the files and save them to Read1 and Read2
# Read the file contents into an array
# Loop through each pair of files
module load cutadapt
ls ./*_1.fq.gz > Read1
ls ./*_2.fq.gz > Read2
mapfile -t read1_array < Read1
mapfile -t read2_array < Read2

for i in "${!read1_array[@]}"
do
  echo "Processing pair ${i}:"
  echo "  Read1: ${read1_array[$i]}"

  cutadapt \
    -a file:./Adapter.fasta \
    -o "RemovedAdapter_$(basename ${read1_array[$i]})" \
    "${read1_array[$i]}"
done

##################
### Remove10bp ###
##################
module load cutadapt
  cutadapt \
    -u 10  \
    -o "Removed10bp_$(basename ${read1_array[$i]})" \
    "${read1_array[$i]}"

##################
### BowtierRNA ###
##################
bowtie2 \
--very-sensitive-local \
--no-unal \
-I 1 -X 1000 -p 12 \
-x $index \
-1 ${read1_array[$i]} \
-2 ${read2_array[$i]} \
--un-conc-gz "RemovedrRNA_${read1_array[$i]}_%.fq.gz" \
2> "sample_Map2rRNAStat_${read1_array[$i]}.xls"

###################
### trim_galore ###
###################
  trim_galore \
    -q 25 \
    --phred33 \
    --length 36 \
    --stringency 3 \
    --paired "${read1_array[$i]}" "${read2_array[$i]}" \
    --gzip -o ./Trim_Galore






