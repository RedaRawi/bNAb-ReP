# bNAb-ReP
An Automated Pipeline for HIV-1 Env Resistance Prediction to 33 Neutralizing Antibodies

## Motivation
Apply bNAb-ReP to predict neutralization susceptibility of HIV-1 Env sequences

## Installation

### Requirements
- R (https://www.r-project.org)
  - R libraries
    - bio3d
    - h2o (version 3.16.0.2) (https://cran.r-project.org/web/packages/h2o/index.html) (see "Old sources")
- MAFFT (version 7) (https://mafft.cbrc.jp/alignment/software/)

## Run 
bNAb-ReP can be run in the command line

4 input arguments are necessary to run the CRISPro.
  1.  bNAb, choose one of 33 bNAbs (10-1074, 2F5, 2G12, 35O22, 3BNC117, 4E10, 8ANC195, CH01, DH270.1, DH270.5, DH270.6, HJ16, NIH45-46, PG16, PG9, PGDM1400, PGT121, PGT128, PGT135, PGT145, PGT151, VRC-CH31, VRC-PG04, VRC01, VRC03, VRC07, VRC13, VRC26.08, VRC26.25, VRC29.03, VRC34.01, VRC38.01, b12)
  2.  HIV-1 Env sequence(s) in FASTA format (https://en.wikipedia.org/wiki/FASTA_format).
      Please use ".fasta" for file name extension 
  3.  Path to installed MAFFT software (e.g. "/usr/local/bin/mafft")
  4.  Output prefix (e.g. "OUTPUT")
  
### Execute in the command line
R --vanilla < CRISPro_v.1.1-2.R 1.5 10 0.0005 8 rama8000-transpro.data mutate.py OUTPUT 1 target.pdb A alternative.pdb A


### Result
- Results will be saved in OUTPUT_DF.csv
- Explanation of Results
  - Option 1 (with target and alternative conformation input) columns include:
    - Residue number
    - Chain in target conformation
    - Amino acid in target conformation
    - Phi angle in target conformation
    - Psi angle in target conformation
    - Compatibility of trans proline angles in target conformation
    - Helical secondary structure in target conformation
    - DSSP secondary structure prediction in target conformation
    - Number of clashes in target conformation
    - Chain in alternative conformation
    - Amino acid in alternative conformation
    - Phi angle in alternative conformation
    - Psi angle in alternative conformation
    - Trans proline angle in alternative conformation
    - Helical secondary structure in alternative conformation
    - DSSP secondary structure prediction in alternative conformation
    - Number of clashes in alternative conformation
    - Compatible with target conformation
    - Compatible with alternative conformation
    - Compatible with target and not with alternative conformation
  - Option 2 (with targe conformation input only) columns include:
    - Chain in target conformation
    - Residue number in target conformation
    - Amino acid in target conformation
    - Phi angle in target conformation
    - Psi angle in target conformation
    - Trans proline angle in target conformation
    - Helical secondary structure in target conformation
    - DSSP secondary structure prediction in target conformation
    - Number of clashes in target conformation
    - Compatible with target conformation
  - Option 3 (with alternative conformation input only) columns include: 
    - Chain in alternative conformation
    - Residue number in alternative conformation
    - Amino acid in alternative conformation
    - Phi angle in alternative conformation
    - Psi angle in alternative conformation
    - Trans proline angle in alternative conformation
    - Helical secondary structure in alternative conformation
    - DSSP secondary structure prediction in alternative conformation
    - Number of clashes in alternative conformation
    - Compatible with alternative conformation


### Example/Test

#### RSV F example
Copy all required files into a directory (Files: CRISPro_v.1.1-2.R, rama8000-transpro.data, mutate.py, 4jhw_trimer.pdb, and 3rrr_trimer.pdb)

- Option 1: R --vanilla < CRISPro_v.1.1-2.R 1.5 10 0.0005 8 rama8000-transpro.data mutate.py RSV_F_option-1 1 4jhw_trimer.pdb A 3rrr_trimer.pdb A
- Option 2: R --vanilla < CRISPro_v.1.1-2.R 1.5 10 0.0005 8 rama8000-transpro.data mutate.py RSV_F_option-2 2 4jhw_trimer.pdb A 
- Option 3: R --vanilla < CRISPro_v.1.1-2.R 1.5 10 0.0005 8 rama8000-transpro.data mutate.py RSV_F_option-3 3 3rrr_trimer.pdb A

# Contact
Reda Rawi: reda.rawi@nih.gov
