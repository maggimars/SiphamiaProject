#!/bin/bash

#SBATCH --job-name=unmapped
#SBATCH --partition=compute
#SBATCH --time=70:30:00
#SBATCH --mem=64G
#SBATCH --nodes=1
#SBATCH --mail-user=margaret.marsbrisbin@oist.jp
#SBATCH --mail-type=BEGIN,FAIL,END
#SBATCH --input=none
#SBATCH --output=unmapped_%j.out
#SBATCH --error=unmapped_%j.err
#SBATCH --array 2-5
. $HOME/.bashrc

DATA=/home/m/maggi-brisbin/SiphComp/aligned

samtools view -b -f 13 $DATA/LO${SLURM_ARRAY_TASK_ID}_sorted.bam >$DATA/LO${SLURM_ARRAY_TASK_ID}_unmapped.bam
