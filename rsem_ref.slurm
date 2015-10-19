#!/bin/bash

#SBATCH --job-name=rsem_ref
#SBATCH --partition=compute
#SBATCH --time=25:30:00
#SBATCH --mem-per-cpu=5G
#SBATCH --ntasks=1
#SBATCH --mail-user=margaret.marsbrisbin@oist.jp
#SBATCH --mail-type=BEGIN,FAIL,END
#SBATCH --input=none
#SBATCH --output=rsemref_%j.out
#SBATCH --error=rsemref_%j.err
. $HOME/.bashrc

rsem-prepare-reference --bowtie2 /home/m/maggi-brisbin/SiphComp/siphgenome/spike/refwspike.fa \
/home/m/maggi-brisbin/SiphComp/siphgenome/bt2ref/out/resem_bt2spike_ref
