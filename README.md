# SiphamiaProject
RNA-seq Differential Expression project with Siphamia tubifer and Photobacterium mandapamensis

# Workflow
## Bacteria
1. Add Spike-in sequences to P. mandapamensis genome
2. Create reference genome index with RSEM (rsem_ref.slurm)
3. Trim barcodes from sample files with Flexbar (trim_bact.slurm, trim_lo.slurm)
4. Calculate differential expression with RSEM 
5. Analyse differential expression with edgeR
6. GO enrichment

## Fish
1. Create bowtie2 reference index for P. mandapamensis genome
2. Trim barecodes from fish muscle sequences with Flexbar
3. Align mixed fish and bacteria (Light Organ) samples to bacteria genome with bowtie2
4. Create file of unaligned sequences - these are the fish LO sequences - Samtools
5. Convert LO files to fastq files - bedtools
6. Construct a reference transcriptome with Trinity
7. Add spike-in sequences to Trinity file
8. Creat reference index with RSEM 
9. Calculate differential expression with RSEM 
10. Analyse differential expression with edgeR
11. Filter Trinity.fasta file to remove genes with less than 1 FPKM averaged across samples and keep only the longest isoform for each gene - Biopython
12. GO enrichment
