#!/bin/bash

#SBATCH --job-name=floRSEM
#SBATCH --partition=compute
#SBATCH --time=150:30:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=12
#SBATCH --nodes=1
#SBATCH --mail-user=margaret.marsbrisbin@oist.jp
#SBATCH --mail-type=BEGIN,FAIL,END
#SBATCH --input=none
#SBATCH --output=floRSEM_%j.out
#SBATCH --error=floRSEM_%j.err
#SBATCH --array 2-5
. $HOME/.bashrc

DATA=/home/m/maggi-brisbin/SiphComp/fishfiles

rsem-calculate-expression --bowtie2 --bowtie2-sensitivity-level very_sensitive --seed-length 19 --num-threads 12 \
--paired-end $DATA/LO${SLURM_ARRAY_TASK_ID}.end1.fastq $DATA/LO${SLURM_ARRAY_TASK_ID}.end2.fastq \
$DATA/ref_transcript/ts_bt2_ref $HOME/SiphComp/fishRSEM/LO${SLURM_ARRAY_TASK_ID}
