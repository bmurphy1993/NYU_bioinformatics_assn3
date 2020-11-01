# NYU_bioinformatics_assn3

# Investigating the cellular response to exogenous DNA

1. Run 'featurecounts_all.sh on high performance cluster (HPC)
2. Run 'kallisto.sh' on HPC
3. Secure copy ("scp") outputs from 'featurecounts_all.sh' and 'kallisto.sh' from HPC to local machine
  - 'kallisto.sh' outputs 6 folders, one for each sample, each with three files: 'abundance.h5', 'abundance.tsv', 'run_info.json'.
  - copy entire folders for each of the kallisto outputs to local machine
4. Run 'transcript-to-gene-import_Group_7.R' on RStudio (change file paths to local directory)
5. Run 'DeSeq2Analysis_Group_7.R' on RStudio (change file paths to local directory)
6. Run 'DeSeq2Analysis_kallisto_Group_7.R' on RStudio (change file paths to local directory)
