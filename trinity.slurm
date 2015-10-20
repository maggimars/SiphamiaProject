#!/bin/bash

#SBATCH --job-name=trinity
#SBATCH --partition=compute
#SBATCH --time=150:30:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=12
#SBATCH --nodes=1
#SBATCH --mail-user=margaret.marsbrisbin@oist.jp
#SBATCH --mail-type=BEGIN,FAIL,END
#SBATCH --input=none
#SBATCH --output=trinity_%j.out
#SBATCH --error=trinity_%j.err
. $HOME/.bashrc

LDATA=/home/m/maggi-brisbin/SiphComp/fishfiles
MDATA=/home/m/maggi-brisbin/SiphComp/trimmed

Trinity --seqType fq --max_memory 64G  \
     --left $LDATA/LO2.end1.fastq,$LDATA/LO3.end1.fastq,$LDATA/LO4.end1.fastq,$LDATA/LO5.end1.fastq,$MDATA/M2_trim__1.fastq,$MDATA/M3_trim__1.fastq,$MDATA/M4_trim__1.fastq,$MDATA/M5_trim__1.fastq \
     --right $LDATA/LO2.end2.fastq,$LDATA/LO3.end2.fastq,$LDATA/LO4.end2.fastq,$LDATA/LO5.end2.fastq,$MDATA/M2_trim__2.fastq,$MDATA/M3_trim__2.fastq,$MDATA/M4_trim__2.fastq,$MDATA/M5_trim__2.fastq  \
     --CPU 12
