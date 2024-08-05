#!/bin/bash
#SBATCH -J PeakCalling
#SBATCH -p q_cn
#SBATCH -o JobOut.PeakCalling.%j
#SBATCH -n 1
#SBATCH --mem=64g

inputdir="/path/to/input"
outputdir="/path/to/output"

###################
### PeakCalling ###
###################

macs2 callpeak -t $inputdir/TetCon_IP1.bam  -c $inputdir/TetCon_Input1.bam  --outdir $outputdir -n TetCon_1 -g mm --shift 0 --extsize 150 -q 0.01
macs2 callpeak -t $inputdir/TetCon_IP2.bam  -c $inputdir/TetCon_Input2.bam  --outdir $outputdir -n TetCon_2 -g mm --shift 0 --extsize 150 -q 0.01
macs2 callpeak -t $inputdir/TetCon_IP3.bam  -c $inputdir/TetCon_Input3.bam  --outdir $outputdir -n TetCon_3 -g mm --shift 0 --extsize 150 -q 0.01
macs2 callpeak -t $inputdir/TetKO_IP1.bam   -c $inputdir/TetKO_Input1.bam   --outdir $outputdir -n  TetKO_1 -g mm --shift 0 --extsize 150 -q 0.01
macs2 callpeak -t $inputdir/TetKO_IP2.bam   -c $inputdir/TetKO_Input2.bam   --outdir $outputdir -n  TetKO_2 -g mm --shift 0 --extsize 150 -q 0.01
macs2 callpeak -t $inputdir/TetKO_IP3.bam   -c $inputdir/TetKO_Input3.bam   --outdir $outputdir -n  TetKO_3 -g mm --shift 0 --extsize 150 -q 0.01