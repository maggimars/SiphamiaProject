#!/bin/bash

#SBATCH --job-name=bamtofast
#SBATCH --partition=compute
#SBATCH --time=150:30:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=12
#SBATCH --nodes=1
#SBATCH --mail-user=margaret.marsbrisbin@oist.jp
#SBATCH --mail-type=BEGIN,FAIL,END
#SBATCH --input=none
#SBATCH --output=bamtofast_%j.out
#SBATCH --error=bamtofast_%j.err
#SBATCH --array 2-5
. $HOME/.bashrc

DATA=/home/m/maggi-brisbin/SiphComp

bedtools bamtofastq -i $DATA/aligned/LO${SLURM_ARRAY_TASK_ID}.qsort.bam \
-fq $DATA/fishfiles/LO${SLURM_ARRAY_TASK_ID}.end1.fastq \
-fq2 $DATA/fishfiles/LO${SLURM_ARRAY_TASK_ID}.end2.fastq
