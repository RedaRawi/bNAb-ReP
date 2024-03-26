# bNAb-ReP
An Automated Pipeline for HIV-1 Resistance Prediction to 33 Neutralizing Antibodies

## Installation

bNAb-ReP has been successfully tested on Linux systems


### Environment
- Download and install conda environment for 64-bit linux or Mac (https://conda.io/miniconda.html) using Python 3.7
  - Install in command line: ./Miniconda3-latest-Linux-x86_64.sh (Linux)
  - Open a new terminal window or source your bash with: source ~/.bashrc
- Create bNAb-ReP environment in command line: conda create --name bNAb-ReP
- Activate bNAb-ReP environment in command line: source activate bNAb-ReP
- Install require packages by running the following in command line:
  - conda install -c r r
  - conda install -c bioconda r-bio3d
  - conda install -c r r-rcurl
  - conda install -c r r-jsonlite
  - conda install -c cidermole jdk8 (Linux) or conda install -c cyclus java-jdk (Mac)
  - conda install -c bioconda mafft
  - conda install -c r r-data.table
  - conda install -c r r-foreach
  - conda install -c r r-doparallel 
  - conda install -c r r-rocr
  - conda install -c conda-forge r-geosphere
  - conda install -c conda-forge readline
- Install R library h2o (version 3.16.0.2) (https://cran.r-project.org/web/packages/h2o/index.html) (see "Old sources") manually by:
  #### The following two commands remove any previously installed H2O packages for R.
  if ("package:h2o" %in% search()) { detach("package:h2o", unload=TRUE) }

  if ("h2o" %in% rownames(installed.packages())) { remove.packages("h2o") }
  #### Next, we download packages that H2O depends on.
  pkgs <- c("RCurl","jsonlite")

  for (pkg in pkgs) {
  if (! (pkg %in% rownames(installed.packages()))) { install.packages(pkg) }
  }

  #### Now we download, install and initialize the H2O package for R.
  install.packages("h2o", type="source", repos="https://h2o-release.s3.amazonaws.com/h2o/rel-wheeler/2/R")

  #### Finally, let's load H2O and start up an H2O cluster
  library(h2o)

  h2o.init()

## Run 
bNAb-ReP can be run in the command line

### 4 input arguments are necessary to run the bNAb-ReP
  1.  bNAb, choose one of 33 bNAbs (10-1074, 2F5, 2G12, 35O22, 3BNC117, 4E10, 8ANC195, CH01, DH270.1, DH270.5, DH270.6, HJ16, NIH45-46, PG16, PG9, PGDM1400, PGT121, PGT128, PGT135,   PGT145, PGT151, VRC-CH31, VRC-PG04, VRC01, VRC03, VRC07, VRC13, VRC26.08, VRC26.25, VRC29.03, VRC34.01, VRC38.01, b12)
  2.  HIV-1 Env sequence(s) in FASTA format (https://en.wikipedia.org/wiki/FASTA_format).
      (Please use ".fasta" for file name extension)
  3.  Path to installed MAFFT software (e.g. "/usr/local/bin/mafft")
  4.  Output prefix (e.g. "OUTPUT")

### Dependencies
  1. Please download directories alignments and models
  2. Copy directories alignments and models in the current path, where you are intending to execute the script

### Execute in the command line
R --vanilla < run_bNAb-ReP_v1.1-3.R VRC01 testing.fasta mafft OUTPUT

### Result
Results will be saved in OUTPUT_probabilities.txt. For each sequence in the input alignment file a probability value between 0-1 is calculated, with higher values corresponding to neutralization sensitivity.



### Example/Test

#### HXB2 and BG505 testing sequences
R --vanilla < run_bNAb-ReP_v1.1-3.R VRC01 testing.fasta mafft OUTPUT

# Contact
Reda Rawi: reda.rawi@nih.gov
