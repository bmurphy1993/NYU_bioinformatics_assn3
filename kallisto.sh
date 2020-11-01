#!/bin/bash
#SBATCH --job-name=kallisto
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=bem6982@nyu.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=12:00:00
#SBATCH --output=/gpfs/scratch/bem6982/bioinformatics/assn_03_gene_exp/logs/kallisto_%j.log
#SBATCH -p cpu_medium

module purge

## Load required modules
module load kallisto/0.44.0

## Create Homo Sapiens Index file
kallisto index -i HomoSapiens /gpfs/data/courses/bmscga2604/gencode.v34.pc_transcripts.fa

##Generate pseudoalignment 
kallisto quant -i HomoSapiens -o /gpfs/scratch/bem6982/bioinformatics/assn_03_gene_exp/outputs/$1 -b 100 --bias -t 4  /gpfs/scratch/bem6982/bioinformatics/assn_03_gene_exp/fastqs/$1.sra_1.fastq /gpfs/scratch/bem6982/bioinformatics/assn_03_gene_exp/fastqs/$1.sra_2.fastq  
##Flags:
##i - re-index 
##-b bootstrap 100
##--bias learns parameters for a model of sequences specific bias and corrects the abundances accordlingly.
##-t, --threads=INT 

