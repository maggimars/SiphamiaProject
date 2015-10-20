#!/bin/bash

#SBATCH --job-name=flex_bact
#SBATCH --partition=compute
#SBATCH --time=70:30:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=12
#SBATCH --nodes=1
#SBATCH --mail-user=margaret.marsbrisbin@oist.jp
#SBATCH --mail-type=BEGIN,FAIL,END
#SBATCH --input=none
#SBATCH --output=trimB_%j.out
#SBATCH --error=trimB_%j.err
#SBATCH --array 29-32
. $HOME/.bashrc

DATA=/home/m/maggi-brisbin/SiphComp/reads

flexbar -t $HOME/SiphComp/trimmed/B${SLURM_ARRAY_TASK_ID}_trim_ -r $DATA/B${SLURM_ARRAY_TASK_ID}_1.fastq \
-p $DATA/B${SLURM_ARRAY_TASK_ID}_2.fastq -f fastq -a $HOME/SiphComp/adapterfiles/B${SLURM_ARRAY_TASK_ID}adapters.fasta -ao 1
