#!/bin/bash

#SBATCH --job-name=trans_genemap
#SBATCH --partition=compute
#SBATCH --time=150:30:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=12
#SBATCH --nodes=1
#SBATCH --mail-user=margaret.marsbrisbin@oist.jp
#SBATCH --mail-type=BEGIN,FAIL,END
#SBATCH --input=none
#SBATCH --output=tgm_%j.out
#SBATCH --error=tgm_%j.err

. $HOME/.bashrc

DATA=/home/m/maggi-brisbin/SiphComp/fishfiles/trinity_out_dir
REF=/home/m/maggi-brisbin/SiphComp/fishfiles/ref_transcript

extract-transcript-to-gene-map-from-trinity $DATA/Trinity.fasta $REF/ref_transcript_
