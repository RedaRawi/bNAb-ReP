# bNAb-ReP
An Automated Pipeline for HIV-1 Resistance Prediction to 33 Neutralizing Antibodies

## Installation

bNAb-ReP was tested on Linux and Mac

### Environment
- Download and install conda environment for 64-bit linux or Mac (https://conda.io/miniconda.html) using Python 3.6
  - Install in command line: ./Miniconda3-latest-Linux-x86_64.sh (Linux) or Miniconda3-latest-MacOSX-x86_64 (Mac)
- Create bNAb-ReP environment in command line: conda create --name bNAb-ReP
- Activate bNAb-ReP environment in command line: source activate bNAb-ReP
- Install require packages by running the following in command line:
  - conda install -c r r r=3.4.1
  - conda install -c conda-forge readline
  - conda install -c bioconda r-bio3d
  - conda install -c r r-rcurl
  - conda install -c r r-jsonlite
  - conda install -c cidermole jdk8
  - conda install -c bioconda mafft
- Install R library h2o (version 3.16.0.2) (https://cran.r-project.org/web/packages/h2o/index.html) (see "Old sources") manually by:
  - Download h2o package (sepecific version)
  - Open R in command line by typing: R
  - Run the following command in R with the path to the downloaded directory in the first argument: install.packages( "path/h2o_3.16.0.2.tar.gz", type = "source", repos = NULL )
  

## Run 
bNAb-ReP can be run in the command line

### 4 input arguments are necessary to run the bNAb-ReP
  1.  bNAb, choose one of 33 bNAbs (10-1074, 2F5, 2G12, 35O22, 3BNC117, 4E10, 8ANC195, CH01, DH270.1, DH270.5, DH270.6, HJ16, NIH45-46, PG16, PG9, PGDM1400, PGT121, PGT128, PGT135, PGT145, PGT151, VRC-CH31, VRC-PG04, VRC01, VRC03, VRC07, VRC13, VRC26.08, VRC26.25, VRC29.03, VRC34.01, VRC38.01, b12)
  2.  HIV-1 Env sequence(s) in FASTA format (https://en.wikipedia.org/wiki/FASTA_format).
      Please use ".fasta" for file name extension 
  3.  Path to installed MAFFT software (e.g. "/usr/local/bin/mafft")
  4.  Output prefix (e.g. "OUTPUT")

### Dependencies
  1. Please download directories alignments and models
  2. Copy directories alignments and models in the current path, where you are intending to execute the script

### Execute in the command line
R --vanilla < run_bNAb-ReP_v1.1-3.R VRC01 testing.fasta /usr/local/bin/mafft OUTPUT

### Result
Results will be saved in OUTPUT_probabilities.txt. For each sequence in the input alignment file a probability value between 0-1 is calculated, with higher values corresponding to neutralization sensitivity.



### Example/Test

#### HXB2 and BG505 testing sequences
R --vanilla < run_bNAb-ReP_v1.1-3.R VRC01 testing.fasta /usr/bin/mafft OUTPUT

# Contact
Reda Rawi: reda.rawi@nih.gov
