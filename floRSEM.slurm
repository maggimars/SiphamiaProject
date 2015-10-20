#!/bin/bash

#SBATCH --job-name=ref_trans
#SBATCH --partition=compute
#SBATCH --time=150:30:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=12
#SBATCH --nodes=1
#SBATCH --mail-user=margaret.marsbrisbin@oist.jp
#SBATCH --mail-type=BEGIN,FAIL,END
#SBATCH --input=none
#SBATCH --output=ref_trans_%j.out
#SBATCH --error=ref_trans_%j.err
. $HOME/.bashrc

DATA=/home/m/maggi-brisbin/SiphComp/fishfiles
    
. $HOME/trinityrnaseq-2.1.0/util/align_and_estimate_abundance.pl --transcripts $DATA/trinity_out_dir/Trinity.fasta \
--est_method RSEM --aln_method bowtie --trinity_mode --prep_reference --output_prefix $DATA/ref_transcript/ref_transcript

