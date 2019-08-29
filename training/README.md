## Generate Training data
See example files for VRC01 and VRC34.01 bNAb in example directory

### Please install all dependencies (see main GitHub page)

### 3 Input arguments required
1. Training sequences in FASTA alignment file format.
2. Neutralization file: Each line contains either a 0 (resistant) or 1 (sensitive) for the corresponding sequences in the input sequences. 
3. Path to installed MAFFT software (e.g. "/usr/local/bin/mafft")

### Execute Preprocessing script in the command line
R --vanilla < bNAb-ReP_preprocess_v.1.1-1.R VRC34.01_IC50_50_alignment.fasta VRC34.01_IC50_50_neutralization.txt /usr/local/bin/mafft


## Perform GBM model training

### 4 Input arguments required
1. Number of cores (CPUs); Required for multiprocessing
2. Memory (in GB); Please provide number of GB in memory you want the training to allocate/use
3. (Full) Path to the current directory where the scripts are executed (output will be saved there) (e.g. "/home/rawir/bNAb-ReP")
4. (Full) Path to the training data, which was preprocessed in the previous step (e.g. "/home/rawir/bNAb-ReP/Training.txt")

### Execute GBM model training script in the command line
R --vanilla < bNAb-ReP_GBM_training.R 8 32 /home/rawir/bNAb-ReP /home/rawir/bNAb-ReP/Training_VRC34.01_IC50_50.txt


## Perform GBM model testing

### 6 Input arguments required
1. Reference alignment (used in training)
2. Path to full model trained (file name should be: "full")
3. File with adjusted cutoff
4. New test alignment
5. Path to installed MAFFT software (e.g. "/usr/local/bin/mafft")
6. Output prefix

### Execute GBM model training script in the command line
R --vanilla < bNAb-ReP_GBM_training.R VRC34.01_IC50_50_alignment.fasta /home/rawir/bNAb-ReP/full final_cutoff_cutoff_balanced.txt NewTestSequences.fasta /usr/local/bin/mafft Output
