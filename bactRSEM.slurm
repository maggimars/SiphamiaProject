#!/bin/bash

#SBATCH --job-name=bactRSEM
#SBATCH --partition=compute
#SBATCH --time=150:30:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=12
#SBATCH --nodes=1
#SBATCH --mail-user=margaret.marsbrisbin@oist.jp
#SBATCH --mail-type=BEGIN,FAIL,END
#SBATCH --input=none
#SBATCH --output=bactRSEM_%j.out
#SBATCH --error=bactRSEM_%j.err
#SBATCH --array 29-32
. $HOME/.bashrc


DATA=/home/m/maggi-brisbin/SiphComp

rsem-calculate-expression --bowtie2 --bowtie2-sensitivity-level very_sensitive --seed-length 19 --num-threads 12 \
--paired-end $DATA/trimmed/B${SLURM_ARRAY_TASK_ID}_trim__1.fastq $DATA/trimmed/B${SLURM_ARRAY_TASK_ID}_trim__2.fastq \
$DATA/siphgenome/bt2ref/out/rsem_bt2spike_ref $DATA/RSEM/B${SLURM_ARRAY_TASK_ID}
