#!/bin/bash
#SBATCH --job-name=featurecounts
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=bem6982@nyu.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=12:00:00
#SBATCH --output=/gpfs/scratch/bem6982/bioinformatics/assn_03_gene_exp/logs/gene_exp_%j.log
#SBATCH -p cpu_medium

module purge

## Load required modules
module load subread/1.6.3

## Run featureCounts on bam files
	## Unstranded
featureCounts -M -s 0 -B -a /gpfs/share/apps/iGenomes/Homo_sapiens/UCSC/hg38/Annotation/Genes/genes.gtf -o ./outputs/counts.txt ./bams/SRR7049609.mapping.sorted.bam ./bams/SRR7049610.mapping.sorted.bam ./bams/SRR7049611.mapping.sorted.bam ./bams/SRR7049612.mapping.sorted.bam ./bams/SRR7049615.mapping.sorted.bam ./bams/SRR7049616.mapping.sorted.bam

