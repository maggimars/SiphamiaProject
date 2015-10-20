# SiphamiaProject
RNA-seq Differential Expression project with Siphamia tubifer and Photobacterium mandapamensis

# Workflow
## Bacteria
1. Add Spike-in sequences to P. mandapamensis genome
2. Create reference genome index with RSEM (rsem_ref.slurm)
3. Trim barcodes from sample files with Flexbar (trim_bact.slurm, trim_lo.slurm)
4. Calculate differential expression with RSEM (bactRSEM.slurm, loRSEM.slurm)
5. Analyse differential expression with edgeR (bacteria.genes.R)
6. GO enrichment

## Fish
2. Trim barecodes from fish muscle sequences with Flexbar (flexM.slurm)
1. Create bowtie2 reference index for P. mandapamensis genome
2. Align mixed fish and bacteria (Light Organ) samples to bacteria genome with bowtie2 (align.slurm)
4. Create file of unaligned sequences - these are the fish LO sequences - Samtools (sambam.slurm, bamsort.slurm, unmapped.slurm)
5. Convert LO files to fastq files - bedtools (sort_unmapped.slurm, bamtofast.slurm)
6. Construct a reference transcriptome with Trinity (trinity.slurm)
7. Add spike-in sequences to Trinity file
8. Creat reference index with RSEM (trans_genemap.slurm, reftrans.slurm)
9. Calculate differential expression with RSEM (floRSEM.slurm, fmRSEM.slurm)
10. Analyse differential expression with edgeR
11. Filter Trinity.fasta file to remove genes with less than 1 FPKM averaged across samples and keep only the longest isoform for each gene - Biopython (https://gist.github.com/maggimars/dde254a81a453bd30b54, (https://gist.github.com/maggimars/e746b0fed46a5c6222d1)
12. GO enrichment
