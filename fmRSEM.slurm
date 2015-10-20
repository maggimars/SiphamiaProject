#!/bin/bash

#SBATCH --job-name=fmRSEM
#SBATCH --partition=compute
#SBATCH --time=150:30:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=12
#SBATCH --nodes=1
#SBATCH --mail-user=margaret.marsbrisbin@oist.jp
#SBATCH --mail-type=BEGIN,FAIL,END
#SBATCH --input=none
#SBATCH --output=fmRSEM_%j.out
#SBATCH --error=fmRSEM_%j.err
#SBATCH --array 2-5
. $HOME/.bashrc

DATA=/home/m/maggi-brisbin/SiphComp/trimmed

rsem-calculate-expression --bowtie2 --bowtie2-sensitivity-level very_sensitive --seed-length 19 --num-threads 12 \
--paired-end $DATA/M${SLURM_ARRAY_TASK_ID}_trim__1.fastq $DATA/M${SLURM_ARRAY_TASK_ID}_trim__2.fastq \
$HOME/SiphComp/fishfiles/ref_transcript/ts_bt2_ref $HOME/SiphComp/fishRSEM/M${SLURM_ARRAY_TASK_ID}
