#!/bin/bash

#SBATCH --job-name=bt_lo
#SBATCH --partition=compute
#SBATCH --time=70:30:00
#SBATCH --mem=64G
#SBATCH --nodes=1
#SBATCH --mail-user=margaret.marsbrisbin@oist.jp
#SBATCH --mail-type=BEGIN,FAIL,END
#SBATCH --input=none
#SBATCH --output=bt_lo_%j.out
#SBATCH --error=bt_lo_%j.err
#SBATCH --array 2-5
. $HOME/.bashrc

DATA=/home/m/maggi-brisbin/SiphComp/

bowtie2 -t -x $DATA/siphgenome/indexgen/siphgenome -1 $DATA/trimmed/LO${SLURM_ARRAY_TASK_ID}_trim__1.fastq \
-2 $DATA/trimmed/LO${SLURM_ARRAY_TASK_ID}_trim__2.fastq -S $DATA/aligned/LO${SLURM_ARRAY_TASK_ID}.sam
